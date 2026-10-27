"""Assign modules using the best of 50 unweighted Louvain restarts."""
import numpy as np
import pandas as pd
import networkx as nx
from networkx.algorithms.community import louvain_communities, modularity
from network_common import RESULT_DIR, GROUPS

def louvain_stats(G, weight=None, seeds=50):
    """Mean/sd of Q and median k over `seeds` Louvain restarts."""
    qs, ks, parts = [], [], []
    for s in range(seeds):
        comms = louvain_communities(G, weight=weight, seed=s)
        qs.append(modularity(G, comms, weight=weight))
        ks.append(len(comms))
        parts.append(comms)
    best = parts[int(np.argmax(qs))]
    return float(np.mean(qs)), float(np.std(qs)), int(np.median(ks)), best


def infer_modules(seeds=50):
    if seeds < 1:
        raise ValueError("At least one Louvain restart is required")
    output = RESULT_DIR / "modularity"
    output.mkdir(parents=True, exist_ok=True)
    for group in GROUPS:
        edges = pd.read_csv(RESULT_DIR / "backbone_v2" / f"{group}_genus_uw" / "backbone_edges.csv")
        graph = nx.Graph()
        graph.add_edges_from(zip(edges.taxon_i.astype(str), edges.taxon_j.astype(str)))
        if not graph.number_of_edges():
            raise ValueError(f"Empty backbone for {group}")
        q, _, _, partition = louvain_stats(graph, seeds=seeds)
        assignments = {node: module for module, nodes in enumerate(
            sorted(partition, key=len, reverse=True), 1) for node in nodes}
        nodes = sorted(graph, key=lambda node: (assignments[node], -graph.degree(node), node))
        pd.DataFrame([dict(taxon_id=node, module=assignments[node], degree=graph.degree(node))
                      for node in nodes]).to_csv(output / f"modules_{group}.csv", index=False)
        print(f"[modules] {group}: {len(partition)} modules; Q={q:.4f}")
