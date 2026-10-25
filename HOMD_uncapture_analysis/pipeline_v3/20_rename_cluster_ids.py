#!/usr/bin/env python
"""Rename the treecluster cluster IDs in every table.

Two names are distinguished:

  n_asvs >= 2   a cluster proper          "21907"   -> "c21907"
  n_asvs == 1   a single ASV              "-28239"  -> "s28239"

The cluster IDs come from treecluster.  A cluster of two or more ASVs is
numbered with a positive integer; an ASV that clusters with nothing else is
numbered with a negative one, and that sign is what made the ID read as a
negative number in R and in a table's first column.  Prefixing "c" and turning
the sign into "s" gives both kinds a name that is a plain identifier, and says
which kind it is.

The 5,143 single-ASV clusters that also contain reference sequences carry a
positive ID and no "-" to replace; by the rule they keep their number
unchanged.  They are listed in the map with action "unchanged".

The rewrite is textual: only the ID field of each line is touched, every other
field is copied byte for byte, so no numeric formatting changes.  Running the
script twice is a no-op -- an ID that already starts with "c" or "s" followed
by digits is left alone.

The upstream scripts (01_treecluster_cluster.py, 02_cluster_aggregate.py) still
write numeric IDs and compare them as numbers, so the rename is a
post-processing step: run this script again whenever the cluster tables are
regenerated, before anything that writes a figure or a model table.

Output: cluster_id_rename_map.tsv (old, new, n_asvs, action)
"""
import csv
import os
import re
import sys

PIPE = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
MAP_OUT = os.path.join(PIPE, "cluster_id_rename_map.tsv")

# cluster_id -> number of ASVs in the cluster
n_asvs = {}


def take(path, id_col, n_col):
    with open(os.path.join(PIPE, path), newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if id_col not in rd.fieldnames or n_col not in rd.fieldnames:
            sys.exit("missing column in %s: %s / %s" % (path, id_col, n_col))
        for row in rd:
            cid, n = row[id_col], row[n_col]
            if n in ("", "NA"):
                continue
            n = int(float(n))
            if cid in n_asvs and n_asvs[cid] != n:
                sys.exit("conflicting n_asvs for %s: %d vs %d (%s)"
                         % (cid, n_asvs[cid], n, path))
            n_asvs[cid] = n


# The cluster table is the authority for every cluster that holds at least one
# ASV.  asv_clusters.tsv is not used: its cluster_size counts the ASVs before
# the chimera screen and disagrees with n_asvs wherever an ASV was dropped
# (113 vs 112 for 21907), which is a difference that could flip a single-ASV
# cluster into a multi-ASV one.  clusters.tsv adds the clusters that hold no ASV
# at all (they keep their number); where the two overlap they agree.
take("02_clusters_aggregated.tsv", "cluster_id", "n_asvs")
take("clusters.tsv", "cluster_id", "n_asv")


def renamed(cid):
    """The new ID, or the old one when the rule leaves it alone."""
    if re.fullmatch(r"[cs][0-9]+", cid):
        return cid                                   # already renamed
    n = n_asvs.get(cid)
    if n is None:
        return cid
    if n >= 2:
        return "c" + cid
    if n == 1 and cid.startswith("-"):
        return "s" + cid[1:]
    return cid


# ---- the tables -------------------------------------------------------------
# (file, column carrying the ID).  `cluster_id` is also the name of the first
# column of the two matrices, so the row names are covered by the same entry.
TABLES = [
    ("02_clusters_aggregated.tsv", "cluster_id"),
    ("05_candidate_novel_with_phylum.tsv", "cluster_id"),
    ("07_master_asv_table.tsv", "cluster_id"),
    ("08_supported_clusters.tsv", "cluster_id"),
    ("12d_anchors_supported_v3.tsv", "cluster_id"),
    ("12d_annotation_isolation_sources.tsv", "cluster_id"),
    ("12d_headline_profile_v3.tsv", "cluster_id"),
    ("12f_headline_clusters_v3.tsv", "cluster_id"),
    ("14_cluster_counts_dedup.tsv", "cluster_id"),
    ("14_cluster_relabund_dedup.tsv", "cluster_id"),
    ("14_supported_cluster_annotation.tsv", "cluster_id"),
    ("15_cluster_cohort_detection.tsv", "cluster_id"),
    ("16_cluster_differential_wilcox.tsv", "cluster_id"),
    ("16_cluster_differential_wilcox_sig.tsv", "cluster_id"),
    ("17_maaslin3_all_results.tsv", "feature"),
    ("17_maaslin3_significant.tsv", "feature"),
    ("17f_ctrl_only_results.tsv", "feature"),
    ("asv_clusters.tsv", "cluster_id"),
    ("clusters.tsv", "cluster_id"),
    (os.path.join("analysis_report", "derep_to_cluster_allcand.tsv"), "cluster_id"),
    (os.path.join("analysis_report", "figs",
                  "local_tree_cluster_annotation_check.tsv"), "cluster_id"),
]

changed = {}
for rel, want in TABLES:
    path = os.path.join(PIPE, rel)
    if not os.path.exists(path):
        print("[00b] SKIP (absent): %s" % rel, file=sys.stderr)
        continue
    tmp = path + ".renaming"
    n_rows = 0
    hits = 0
    with open(path, newline="") as fin, open(tmp, "w", newline="") as fout:
        header = fin.readline()
        cols = header.rstrip("\n").split("\t")
        if want not in cols:
            os.remove(tmp)
            sys.exit("%s has no %s column" % (rel, want))
        idx = cols.index(want)
        fout.write(header)
        for line in fin:
            n_rows += 1
            fields = line.rstrip("\n").split("\t")
            if idx < len(fields):
                old = fields[idx]
                new = renamed(old)
                if new != old:
                    fields[idx] = new
                    hits += 1
                    changed[old] = new
                elif old in n_asvs and old not in changed and not re.fullmatch(
                        r"[cs][0-9]+", old):
                    # kept deliberately (a positive single-ASV cluster)
                    changed[old] = old
            fout.write("\t".join(fields) + "\n")
    if hits == 0 and rel.endswith(".tsv"):
        # nothing to do: another run of an already-renamed table
        os.replace(tmp, path)
        print("[00b] %-55s already renamed (0 rows touched)" % rel)
        continue
    os.replace(tmp, path)
    print("[00b] %-55s %6d rows, %6d IDs renamed" % (rel, n_rows, hits))

# ---- map --------------------------------------------------------------------
# written only on the run that did the renaming: a second run has nothing to
# record and would otherwise empty the file
if not changed:
    print("[00b] nothing to rename; map left as it is")
    sys.exit(0)
with open(MAP_OUT, "w", newline="") as fh:
    fh.write("old_id\tnew_id\tn_asvs\taction\n")
    for old in sorted(changed, key=lambda s: (len(s), s)):
        new = changed[old]
        n = n_asvs.get(old, "")
        act = ("unchanged" if new == old else
               ("cluster" if new.startswith("c") else "singleton"))
        fh.write("%s\t%s\t%s\t%s\n" % (old, new, n, act))
print("[00b] map written: %s (%d IDs)" % (MAP_OUT, len(changed)))
