"""Scientific invariants for age features, grouped CV, and HC summaries."""

from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from age_common import read_predictions, site_predictions
from age_model import ModelConfig, build_xy_pool, fit_fold, load_genus_table, split_pool
from age_stats import disease_age_gap_stats
from run_age_prediction import main


class AgePredictionTests(unittest.TestCase):
    def test_count_loader_preserves_genus_index_and_first_sample(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "counts.tsv"
            path.write_text("\ts1\ts2\ts3\nA\t2\t3\t4\nB\t5\t7\t9\n")
            features, samples, matrix = load_genus_table({"s1", "s3"}, path)
            self.assertEqual(features.tolist(), ["A", "B"])
            self.assertEqual(samples.tolist(), ["s1", "s3"])
            np.testing.assert_array_equal(matrix, [[2, 5], [4, 9]])
            with self.assertRaisesRegex(ValueError, "No sample IDs overlap"):
                load_genus_table({"missing"}, path)

    def test_features_use_hc_training_vocabulary_and_test_order(self):
        metadata = pd.DataFrame({"sample_id": ["h1", "h2", "d1", "t1", "t2"],
                                 "status_cleaned": ["HC", "HC", "Disease", "HC", "Disease"],
                                 "age_cleaned": [10, 20, 30, 40, 50]})
        counts = np.array([[2, 2, 0], [1, 3, 0], [0, 0, 9], [1, 1, 8], [3, 1, 0]], dtype=float)
        result = build_xy_pool(counts, metadata.sample_id.to_numpy(), metadata.iloc[:3],
                               metadata.iloc[[4, 3]], metadata.set_index("sample_id").age_cleaned)
        x_train, y_train, x_test, y_test = result
        np.testing.assert_allclose(x_train, np.log10([[.5 + 1e-6, .5 + 1e-6], [.25 + 1e-6, .75 + 1e-6]]))
        np.testing.assert_allclose(x_test, np.log10([[.75 + 1e-6, .25 + 1e-6], [.1 + 1e-6, .1 + 1e-6]]))
        np.testing.assert_array_equal(y_train, [10, 20])
        np.testing.assert_array_equal(y_test, [50, 40])

    def test_grouped_folds_and_insufficient_hc(self):
        pool = pd.DataFrame({"sample_id": [f"s{i}" for i in range(24)],
                             "bioproject_id": np.repeat(["p1", "p2", "p3"], 8),
                             "status_cleaned": ["Disease"] * 24})
        config = ModelConfig(folds=5, trees=8, jobs=1)
        held_out = []
        for train, test in split_pool(pool, "project", config):
            self.assertFalse(set(train.bioproject_id) & set(test.bioproject_id))
            held_out.extend(test.sample_id)
            self.assertIsNone(fit_fold(None, None, train, test, None, config))
        self.assertCountEqual(held_out, pool.sample_id)
        self.assertEqual(len(list(split_pool(pool.iloc[:8], "project", config))), 0)

    def test_prediction_site_alignment_and_unique_test_samples(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            pd.DataFrame({"sample_id": ["001", "002", "003"], "true_age": [20, 30, 40],
                          "pred_age": [22, 25, 80]}).to_csv(
                root / "age_predict_test_1_oral_mucosa_genus_range_1_80.csv", index=False)
            metadata = pd.DataFrame({"sample_id": ["001", "002", "003"],
                                     "status_cleaned": ["HC", "HC", "Disease"],
                                     "bioproject_id": ["p"] * 3, "site_cleaned": ["oral mucosa"] * 3})
            loaded = site_predictions("oral mucosa", metadata, root)
            self.assertEqual(len(loaded), 3)
            self.assertEqual(loaded.sample_id.tolist(), ["001", "002", "003"])
            np.testing.assert_array_equal(loaded.delta_age, [2, -5, 40])
            self.assertEqual(len(site_predictions("tongue", metadata, root)), 0)
            pd.DataFrame({"sample_id": ["001"], "true_age": [20], "pred_age": [22]}).to_csv(
                root / "age_predict_test_2_oral mucosa_genus_range_1_80.csv", index=False)
            with self.assertRaisesRegex(ValueError, "Repeated out-of-fold"):
                read_predictions(root, "genus", "oral mucosa")

    def test_disease_statistics_hc_anchor(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            ids = [f"s{i}" for i in range(24)]
            meta = pd.DataFrame({"sample_id": ids, "site_cleaned": ["saliva"] * 24,
                                 "status_cleaned": ["HC"] * 12 + ["Disease"] * 12,
                                 "bioproject_id": ["p1"] * 24})
            pd.DataFrame({"sample_id": ids, "true_age": [30] * 24,
                          "pred_age": [32] * 12 + [37] * 12}).to_csv(
                root / "age_predict_test_1_saliva_genus_range_1_80.csv", index=False)
            stats = disease_age_gap_stats(meta, root, n_boot=30)
            self.assertEqual(stats.hc_n.iloc[0], 12)
            self.assertEqual(stats.hc_med.iloc[0], 2)
            self.assertEqual(stats.delta.iloc[0], 5)
            self.assertEqual(stats.ci_lo.iloc[0], 5)
            self.assertEqual(stats.ci_hi.iloc[0], 5)

    def test_training_cli_for_both_schemes(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            ids = [f"s{i}" for i in range(36)]
            metadata = pd.DataFrame({"sample_id": ids, "age_cleaned": np.arange(36) + 10,
                                     "site_cleaned": ["saliva"] * 36,
                                     "bioproject_id": np.repeat(["p1", "p2", "p3"], 12),
                                     "status_cleaned": (["HC"] * 9 + ["Disease"] * 3) * 3})
            metadata.to_csv(root / "meta.tsv", sep="\t", index=False)
            pd.DataFrame(np.random.RandomState(42).randint(1, 20, (4, 36)),
                         index=["A", "B", "C", "D"], columns=ids).to_csv(root / "counts.tsv", sep="\t")
            for scheme, token in [("project", "genus"), ("random", "random")]:
                output = root / scheme
                main(["--scheme", scheme, "--sites", "saliva", "--metadata", str(root / "meta.tsv"),
                      "--genus-table", str(root / "counts.tsv"), "--output-dir", str(output),
                      "--trees", "8", "--folds", "3", "--jobs", "1"])
                predictions = read_predictions(output, token, "saliva")
                self.assertCountEqual(predictions.sample_id, ids)
                self.assertEqual(predictions.fold.nunique(), 3)
                main(["--scheme", scheme, "--sites", "saliva", "--metadata", str(root / "meta.tsv"),
                      "--genus-table", str(root / "counts.tsv"), "--output-dir", str(output),
                      "--trees", "8", "--folds", "2", "--jobs", "1"])
                self.assertEqual(read_predictions(output, token, "saliva").fold.nunique(), 2)


if __name__ == "__main__":
    unittest.main(buffer=True)
