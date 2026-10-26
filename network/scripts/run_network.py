"""Prepare paired cohorts, infer networks, or regenerate the final figure/table."""
import argparse
from pathlib import Path

from network_common import RESULT_DIR


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("plot", "table", "report", "prepare", "sparcc", "pool", "modules", "all"))
    parser.add_argument("--output-dir", type=Path, help="Figure/table output directory")
    parser.add_argument("--counts", type=Path, help="Source genus count TSV for prepare")
    parser.add_argument("--metadata", type=Path, help="Source metadata TSV for prepare")
    parser.add_argument("--unit", type=int, help="One-based SparCC unit index; omit to run all units")
    args = parser.parse_args()
    if args.stage in ("prepare", "all"):
        from prepare_network import prepare_units
        prepare_units(args.metadata, args.counts)
    if args.stage in ("sparcc", "all"):
        from prepare_network import run_sparcc
        run_sparcc(args.unit)
    if args.stage in ("pool", "all"):
        from pool_network import run_pooling
        run_pooling()
    if args.stage in ("modules", "all"):
        from infer_modules import infer_modules
        infer_modules()
    if args.stage in ("plot", "table", "report", "all"):
        from network_common import prepare_graphs
        graphs = prepare_graphs()
        output_dir = args.output_dir or RESULT_DIR / "module_networks"
        if args.stage in ("table", "report", "all"):
            from module_members import write_members
            write_members(graphs, output_dir)
        if args.stage in ("plot", "report", "all"):
            from plot_module_networks import make_combined
            make_combined(graphs, str(output_dir))


if __name__ == "__main__":
    main()
