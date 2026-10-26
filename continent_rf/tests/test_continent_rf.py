"""Small scientific invariants; run with unittest, no cluster or real data needed."""

from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from biom import Table
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from rf_cv import CVConfig, run_cv, stratified_project_fold_split
from rf_data import aggregate_genus, load_dataset, parse_genus
from run_continent_rf import main


class ContinentRFTests(unittest.TestCase):
    def test_whole_project_never_leaks_across_folds(self):
        positive = ["shared", "p1", "p2"]
        negative = ["shared", "n1", "n2"]
        tests = []
        for train, test in stratified_project_fold_split(
            positive, negative, 3, np.random.RandomState(43)
        ):
            self.assertFalse(set(train) & set(test))
            self.assertEqual(set(train) | set(test), set(positive + negative))
            tests.extend(test)
        self.assertCountEqual(tests, set(positive + negative))

    def test_genus_sum_normalization_and_zero_rows(self):
        table = Table(np.array([[2, 0, 0], [3, 1, 0], [5, 3, 0], [10, 0, 7]]),
                      ["a", "b", "c", "unassigned"], ["s1", "s2", "s3"])
        result = aggregate_genus(table, {"a": "A", "b": "A", "c": "B"})
        np.testing.assert_array_equal(result.toarray(), [[.5, .5], [.25, .75], [0, 0]])
        self.assertTrue(sparse.isspmatrix_csr(result))
        self.assertEqual(parse_genus("k__Bacteria; g__ A; s__x"), "A")
        self.assertIsNone(parse_genus("k__Bacteria; g__"))
        with self.assertRaisesRegex(ValueError, "No BIOM features"):
            aggregate_genus(table, {})

    def test_ovo_ova_labels_and_held_out_predictions(self):
        regions = ["North America", "Asia", "Europe"]
        rows = [dict(sample_id=f"{r}_{p}_{s}", region=r, project=f"{r}_{p}",
                     unit=f"{r}_{p}", site="Saliva")
                for r in regions for p in range(3) for s in range(2)]
        meta = pd.DataFrame(rows)
        meta["_matrix_index"] = np.arange(len(meta))
        matrix = sparse.csr_matrix(np.repeat(np.eye(3), 6, axis=0))
        config = CVConfig(repeats=2, trees=8, jobs=1)
        with tempfile.TemporaryDirectory() as temp:
            for task in ["ovo", "ova"]:
                path = Path(temp) / f"{task}.tsv"
                folds, summaries = run_cv(meta, matrix, "oral", task, path, config)
                self.assertEqual(len(summaries), 3)
                self.assertEqual(len(folds), 18)
                pred = pd.read_csv(path, sep="\t")
                identity = ["target_region"] if task == "ova" else ["region_1", "region_2"]
                self.assertTrue((pred.groupby(identity + ["repeat", "sample_id"]).size() == 1).all())
                self.assertTrue((pred.groupby(identity + ["repeat", "project"])["fold"].nunique() == 1).all())
                target = pred["target_region"] if task == "ova" else pred["region_1"]
                np.testing.assert_array_equal(pred["true_label"], (pred["region"] == target).astype(int))
                self.assertTrue(all(row["mean_auc"] == 1 for row in summaries))

    def test_runner_alignment_and_empty_outputs(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            table = Table(np.array([[0, 4, 1, 3], [2, 0, 1, 1]]), ["a", "b"], ["s4", "s2", "s1", "s3"])
            biom_path = root / "counts.biom"
            biom_path.write_text(table.to_json("test"))
            metadata = root / "meta.tsv"
            pd.DataFrame({"sample_id": ["s1", "s2", "s3", "s4"],
                          "region": ["Asia", "Asia", "Europe", "Europe"],
                          "project": ["p1", "p1", "p2", "p2"]}).to_csv(metadata, sep="\t", index=False)
            meta, matrix = load_dataset("gut", metadata, biom_path, min_study_samples=1)
            self.assertEqual(meta["sample_id"].tolist(), ["s4", "s2", "s1", "s3"])
            np.testing.assert_array_equal(matrix.toarray(), [[0, 1], [1, 0], [.5, .5], [.75, .25]])
            args = ["--datasets", "gut", "--level", "feature", "--tasks", "ova", "ovo",
                    "--gut-metadata", str(metadata), "--gut-biom", str(biom_path),
                    "--output-dir", str(root / "output"), "--min-study-samples", "1"]
            # Only one project per region: all tasks are skipped, but TSVs retain headers.
            main(args)
            for task in ["ova", "ovo"]:
                for kind in ["summary", "fold_metrics", "predictions"]:
                    frame = pd.read_csv(root / "output/gut" / f"region_{task}_{kind}.tsv", sep="\t")
                    self.assertTrue(frame.empty)
                    self.assertIn("dataset", frame.columns)
            config = json.loads((root / "output/gut/run_config.json").read_text())
            self.assertEqual(config["datasets"], ["gut"])
            self.assertEqual(config["split_mode"], "project")

    def test_genus_cli_trains_both_datasets_and_tasks(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            regions = ["North America", "Asia", "Europe"]
            meta = pd.DataFrame([
                {"sample_id": f"s{i}_{j}", "region": region, "project": f"p{i}_{j}"}
                for i, region in enumerate(regions) for j in range(2)
            ])
            table = Table(np.array([[1, 2, 3, 4, 5, 6], [3, 2, 1, 6, 5, 4], [1] * 6]),
                          ["a", "b", "unassigned"], meta["sample_id"].tolist())
            (root / "counts.biom").write_text(table.to_json("test"))
            meta.to_csv(root / "gut.tsv", sep="\t", index=False)
            oral = meta.rename(columns={"region": "continent_group", "project": "bioproject_id"})
            oral["site_cleaned"] = "saliva"
            oral.to_csv(root / "oral.tsv", sep="\t", index=False)
            pd.DataFrame({"Feature ID": ["a", "b", "unassigned"],
                          "Taxon": ["g__A", "g__B", "g__"]}).to_csv(
                root / "taxonomy.tsv", sep="\t", index=False)
            runner = Path(__file__).resolve().parents[1] / "scripts/run_continent_rf.py"
            command = [sys.executable, str(runner), "--tasks", "ova", "ovo",
                       "--split-mode", "unit", "--repeats", "1", "--trees", "8",
                       "--jobs", "1", "--min-study-samples", "1",
                       "--taxonomy", str(root / "taxonomy.tsv"), "--output-dir", str(root / "out")]
            for dataset in ["oral", "gut"]:
                command += [f"--{dataset}-metadata", str(root / f"{dataset}.tsv"),
                            f"--{dataset}-biom", str(root / "counts.biom")]
            subprocess.run(command, check=True, capture_output=True, text=True)
            for dataset in ["oral", "gut"]:
                for task in ["ova", "ovo"]:
                    summary = pd.read_csv(root / f"out/{dataset}/region_{task}_summary.tsv", sep="\t")
                    self.assertEqual(len(summary), 3)
                    self.assertEqual(set(summary["dataset"]), {dataset})
                    self.assertTrue((summary["n_evaluations"] == 2).all())
                config = json.loads((root / f"out/{dataset}/run_config.json").read_text())
                self.assertEqual(config["level"], "genus")
                self.assertEqual(config["split_unit"], "region||project")


if __name__ == "__main__":
    with redirect_stdout(io.StringIO()):
        unittest.main()
