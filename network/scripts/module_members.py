"""Build the wide module-member table from the same graphs as the figure."""
from pathlib import Path

import pandas as pd

from network_common import RESULT_DIR, ROLE_BLUE, ROLE_PURPLE

TABLE_ORDER = ("HC", "PD", "OSCC", "DC")


def member_table(graphs):
    rows = []
    for group in TABLE_ORDER:
        graph = graphs[group]
        hubs, connectors = {}, {}
        for node, color in graph.graph["role_nodes"].items():
            module = graph.nodes[node]["rank"]
            if color != ROLE_BLUE:
                hubs[module] = node
            if color in (ROLE_BLUE, ROLE_PURPLE):
                connectors[module] = node
        members = pd.DataFrame([
            dict(member=str(node).removeprefix("g__"), module=data["rank"],
                 degree=graph.degree(node), abundance=data["abu"] * 100)
            for node, data in graph.nodes(data=True)
        ]).sort_values(["module", "degree", "abundance"], ascending=[True, False, False])
        for module, sub in members.groupby("module", sort=False):
            rows.append({
                "module": f"{group}-M{module}",
                "health_status": "healthy" if group == "HC" else group,
                "module hub": hubs.get(module, ""),
                "main connector": connectors.get(module, ""),
                "member (degree, mean rel. abund.)": "; ".join(
                    f"{r.member}({r.degree},{r.abundance:.6f}%)" for r in sub.itertuples(index=False)),
            })
    return pd.DataFrame(rows)


def write_members(graphs, output_dir=None):
    output_dir = Path(output_dir or RESULT_DIR / "module_networks")
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "all_module_members_table_wide.csv"
    member_table(graphs).to_csv(path, index=False)
    print(f"[members] -> {path}")
