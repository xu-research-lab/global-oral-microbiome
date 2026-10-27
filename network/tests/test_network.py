"""Validate paired cohorts, pooling, module roles and final table contents."""
from contextlib import redirect_stdout
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import networkx as nx
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import network_common as common
import prepare_network as prepare
import pool_network as pool
from module_members import member_table


class NetworkTests(unittest.TestCase):
    def test_project_pairing_uses_all_sites(self):
        metadata = pd.DataFrame({
            "bioproject_id": ["P1", "P1", "P2", "P3", "P3"],
            "status_cleaned": ["Healthy Control", "Dental Caries", "Dental Caries", "Healthy Control", "OSCC"],
            "site_cleaned": ["Saliva", "Tooth", "Tooth", "Saliva", None],
        })
        self.assertEqual(prepare.select_projects(metadata).index.tolist(), ["P1"])

    def test_unit_filter_requires_detection_in_five_samples(self):
        counts = pd.DataFrame([[1] * 5 + [0], [2] * 5 + [0], [1, 1, 0, 0, 0, 0]],
                              index=["A", "B", "rare"])
        filtered = prepare.filter_unit(counts)
        self.assertEqual(filtered.shape, (2, 5))
        self.assertEqual(filtered.index.tolist(), ["A", "B"])
        self.assertIsNone(prepare.filter_unit(counts.iloc[:, :4]))

    def test_abundance_uses_all_group_samples(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp) / "biom_paired/HC_genus"
            folder.mkdir(parents=True)
            pd.DataFrame([[1], [1]], index=["A", "B"], columns=["s1"]).to_csv(folder / "P1.tsv", sep="\t")
            pd.DataFrame([[2]], index=["B"], columns=["s2"]).to_csv(folder / "P2.tsv", sep="\t")
            self.assertEqual(common.load_mean_abundance("HC", tmp), {"A": 0.25, "B": 0.75})

    def test_cartography_marks_connectors_and_combined_roles(self):
        graph = nx.Graph([("a", "b"), ("b", "c"), ("c", "a"), ("b", "d"), ("d", "e"), ("d", "f")])
        nx.set_node_attributes(graph, {n: 1 if n in "abc" else 2 for n in graph}, "rank")
        self.assertEqual(common.cartography_role_nodes(graph), {"b": common.ROLE_BLUE, "d": common.ROLE_PURPLE})

    def test_cache_fingerprint_follows_contents_and_sample_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            first, second = Path(tmp) / "first.csv", Path(tmp) / "second.csv"
            first.write_text("matrix A")
            second.write_text(first.read_text())
            unit = dict(study="P1", n=10, path=str(first))
            expected = pool.input_fingerprint([unit])
            self.assertEqual(expected, pool.input_fingerprint([dict(unit, path=str(second))]))
            self.assertNotEqual(expected, pool.input_fingerprint([dict(unit, n=11)]))
            first.write_text("matrix B")
            self.assertNotEqual(expected, pool.input_fingerprint([unit]))

    def test_prepare_invalidates_only_changed_cohort_correlations(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            records = [dict(sample_id=f"{p}_{g}_{i}", bioproject_id=p,
                            status_cleaned=g, site_cleaned="Saliva")
                       for p in ("P1", "P2") for g in ("Healthy Control", "Dental Caries") for i in range(5)]
            metadata = root / "metadata.tsv"
            counts_path = root / "counts.tsv"
            pd.DataFrame(records).to_csv(metadata, sep="\t", index=False)
            counts = pd.DataFrame([np.ones(20), np.full(20, 2)], index=["A", "B"],
                                  columns=[r["sample_id"] for r in records])
            counts.to_csv(counts_path, sep="\t")
            output = root / "result"
            with patch.object(prepare, "RESULT_DIR", output), redirect_stdout(io.StringIO()):
                prepare.prepare_units(metadata, counts_path)
                units = pd.read_csv(output / "units.csv")
                self.assertEqual(len(units), 4)
                self.assertTrue(units.kept.all())
                self.assertTrue(all(not Path(path).is_absolute() for path in units.path))
                correlation = output / "per_site_sparcc/HC_genus/P1.csv"
                correlation.parent.mkdir(parents=True)
                correlation.write_text("cached")
                prepare.prepare_units(metadata, counts_path)
                self.assertTrue(correlation.exists())
                counts.iloc[0, 0] = 3
                counts.to_csv(counts_path, sep="\t")
                prepare.prepare_units(metadata, counts_path)
                self.assertFalse(correlation.exists())

    def test_saved_backbone_edges_match_cohort_statistics(self):
        units = pd.read_csv(common.RESULT_DIR / "units.csv", dtype={"study": str})
        with redirect_stdout(io.StringIO()):
            for group in common.GROUPS:
                selected = pool.discover_units(group, units)
                features = pool.collect_fid(selected)
                accumulated, cohorts = pool.build_stats(selected, features)
                consensus = pd.read_csv(common.RESULT_DIR / "backbone_netmoss" / f"{group}_genus" /
                                        "pooled_union.csv", index_col=0).reindex(index=features, columns=features).to_numpy()
                missing = ~np.isfinite(consensus)
                consensus[missing] = pool.consensus_pool_py(accumulated)[missing]
                signs = np.sign(consensus)
                np.fill_diagonal(signs, 1)
                actual, _, _, _ = pool.build_tables(accumulated, cohorts, len(features), features,
                                                     0.20, 5, 0.60, 0.20, 0.025, signs)
                saved = pd.read_csv(common.RESULT_DIR / "backbone_v2" / f"{group}_genus_uw" / "backbone_edges.csv")
                pd.testing.assert_frame_equal(actual.loc[actual.final, saved.columns].reset_index(drop=True),
                                              saved, atol=1e-12, rtol=1e-10)

    def test_saved_wide_table_matches_networks_and_abundance(self):
        actual = member_table(common.prepare_graphs())
        saved = pd.read_csv(common.RESULT_DIR / "module_networks/all_module_members_table_wide.csv").fillna("")
        pd.testing.assert_frame_equal(actual, saved, check_exact=True)


if __name__ == "__main__":
    unittest.main()
