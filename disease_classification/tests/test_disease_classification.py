"""Check fold isolation, genus filters, RF alignment and saved figure inputs."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import biom
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import disease_common as common
import disease_data as data
import disease_model as model


class DiseaseClassificationTests(unittest.TestCase):
    def test_metadata_rejects_duplicate_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "meta.tsv"
            pd.DataFrame({"sample_id": ["s1", "s1"]}).to_csv(path, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                common.load_metadata(path)

    def test_labels_exclude_unknown_status(self):
        meta = pd.DataFrame({"sample_id": ["a", "b", "c", "d"],
                             "status_cleaned": ["HC", "Disease", "Other", None]})
        with patch.object(common, "load_metadata", return_value=meta):
            self.assertEqual(common.load_labels().to_dict(), {"a": 0, "b": 1})

    def test_auc_handles_missing_and_single_class(self):
        self.assertEqual(model.safe_auc([0, 1, np.nan], [0.1, 0.9, 0.5]), 1)
        self.assertTrue(np.isnan(model.safe_auc([1, 1], [0.1, 0.9])))

    def test_mean_and_prevalence_filters_differ(self):
        table = biom.Table(np.array([[0.995] * 4, [0.005] * 4, [0] * 4]),
                           ["common", "rare", "absent"], ["a", "b", "c", "d"])
        meta = pd.DataFrame({"sample_id": ["d", "b", "a", "c"],
                             "bioproject_id": ["P2", "P1", "P1", "P2"]})
        self.assertEqual(data.select_genera(table, meta, "mean").tolist(), ["common"])
        self.assertEqual(data.select_genera(table, meta, "prevalence").tolist(), ["common", "rare"])
        self.assertEqual(data.parse_genus("d__B; g__", True), "g__")
        self.assertIsNone(data.parse_genus("d__B; g__", False))

    def test_fit_aligns_test_features_and_never_trains_on_test_labels(self):
        train = biom.Table(np.array([[9, 8, 1, 2], [1, 2, 9, 8]]), ["a", "b"], ["s1", "s2", "s3", "s4"])
        test = biom.Table(np.array([[8, 2], [2, 8]]), ["b", "a"], ["s5", "s6"])
        labels = pd.Series([0, 0, 1, 1], index=train.ids(axis="sample"))
        ids, scores, importance = model.fit_fold(train, test, labels, jobs=1)
        reference = RandomForestClassifier(**common.RF_OPTIONS, n_jobs=1)
        reference.fit(np.array([[0.9, 0.1], [0.8, 0.2], [0.1, 0.9], [0.2, 0.8]]), labels)
        np.testing.assert_array_equal(ids, ["s5", "s6"])
        np.testing.assert_allclose(scores, reference.predict_proba([[0.2, 0.8], [0.8, 0.2]])[:, 1])
        self.assertAlmostEqual(importance.importance.sum(), 1)

    def test_prepare_writes_disjoint_project_splits(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            records = [{"sample_id": f"{p}_{g}_{i}", "bioproject_id": p,
                        "status_cleaned": g, "site_cleaned": "saliva"}
                       for p in ("P1", "qiita_123") for g in ("HC", "Disease Name") for i in range(16)]
            meta = pd.DataFrame(records)
            low = dict(sample_id="low", bioproject_id="P1", status_cleaned="HC", site_cleaned="saliva")
            source_meta = pd.concat([meta, pd.DataFrame([low])], ignore_index=True)
            source_meta.to_csv(root / "metadata_oral_qc_dedup.tsv", sep="\t", index=False)
            table = biom.Table(np.array([[1194] * 64 + [6], [6] * 64 + [0]]),
                               ["f1", "f2"], source_meta.sample_id.tolist())
            with biom.util.biom_open(str(root / "biom_oral_qc_dedup.biom"), "w") as handle:
                table.to_hdf5(handle, "test")
            taxonomy = root / "taxonomy.tsv"
            pd.DataFrame({"Feature ID": ["f1", "f2"], "Taxon": ["g__Common", "g__Rare"]}).to_csv(taxonomy, sep="\t", index=False)
            split_dir = root / "splits"
            with patch.object(data, "INPUT_DIR", root), patch.object(common, "INPUT_DIR", root), \
                 patch.object(data, "TAXONOMY_PATH", taxonomy), patch.object(data, "DATA_DIR", split_dir), \
                 patch.object(common, "DATA_DIR", split_dir):
                for scheme, n_features in (("mean", 1), ("prevalence", 2)):
                    data.prepare_splits(scheme)
                    manifest = common.read_manifest(scheme)
                    self.assertEqual(len(manifest), 2)
                    for row in manifest.itertuples(index=False):
                        name = common.fold_name(row)
                        train = biom.load_table(str(common.split_path(scheme, "train", name)))
                        test = biom.load_table(str(common.split_path(scheme, "test", name)))
                        self.assertFalse(set(train.ids(axis="sample")) & set(test.ids(axis="sample")))
                        held = meta.set_index("sample_id").loc[test.ids(axis="sample"), "bioproject_id"]
                        self.assertEqual(set(held), {row.study})
                        self.assertEqual(train.shape, (n_features, 32))
                        self.assertEqual(test.shape, (n_features, 32))
                        self.assertNotIn("low", train.ids(axis="sample"))
                        np.testing.assert_allclose(test.sum(axis="sample"), 1)

    def test_auc_summary_uses_prevalence_predictions_and_folds(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = pd.DataFrame([dict(site="saliva", disease="Disease", study="P1")])
            manifest.to_csv(root / "run_info_prevalence.csv", index=False)
            # The alternative filter has different predictions and an extra fold.
            pd.concat([manifest, manifest.assign(study="P2")]).to_csv(
                root / "run_info_mean.csv", index=False)
            predictions = pd.DataFrame({"sampel_id": ["a", "b"],
                                        "group": ["saliva_Disease_P1"] * 2,
                                        "pred_prob": [0.1, 0.9]})
            predictions.to_csv(root / "RF_results_lsoP.csv", index=False)
            predictions.assign(pred_prob=[0.9, 0.1]).to_csv(root / "RF_results_lso.csv", index=False)
            with patch.object(common, "DATA_DIR", root), patch.object(common, "RESULT_DIR", root), \
                 patch.object(model, "RESULT_DIR", root), \
                 patch.object(model, "load_labels", return_value=pd.Series([0, 1], index=["a", "b"])):
                result = model.auc_tables()
            self.assertEqual(len(result), 1)
            self.assertEqual(result.iloc[0].test_study, "P1")
            self.assertEqual(result.iloc[0].auc, 1.0)

    def test_saved_auc_tables_match_predictions_and_metadata(self):
        saved = []
        with patch.object(pd.DataFrame, "to_csv", lambda frame, path, **kw: saved.append((path, frame.copy()))):
            model.auc_tables()
        for path, calculated in saved:
            pd.testing.assert_frame_equal(calculated.reset_index(drop=True), pd.read_csv(path), atol=1e-14, rtol=0)


if __name__ == "__main__":
    unittest.main()
