"""Entry point for prevalence-filtered genus-level LOSO figures."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("plot", "auc", "heatmap", "summarize", "prepare", "fit", "all"))
    parser.add_argument("--scheme", choices=("mean", "prevalence", "both"),
                        help="Filter for prepare/fit (default: prevalence); all uses prevalence")
    parser.add_argument("--jobs", type=int, default=16)
    parser.add_argument("--output-dir", type=Path, help="Figure output directory; inputs remain in the result directory")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if args.scheme is not None and args.stage not in ("prepare", "fit"):
        parser.error("--scheme applies to prepare and fit")
    from disease_common import RESULT_DIR
    from disease_model import auc_tables, fit_models
    scheme = args.scheme or "prevalence"
    schemes = ("mean", "prevalence") if scheme == "both" else (scheme,)
    if args.stage in ("prepare", "all"):
        from disease_data import prepare_splits
        for scheme in schemes:
            prepare_splits(scheme)
    if args.stage in ("fit", "all"):
        for scheme in schemes:
            fit_models(scheme, args.jobs)
    if args.stage in ("summarize", "all"):
        auc_tables()
    output_dir = (args.output_dir or RESULT_DIR).resolve()
    if args.stage in ("plot", "auc", "all"):
        from plot_loso_auc import plot_auc
        plot_auc(output_dir)
    if args.stage in ("plot", "heatmap", "all"):
        rscript = os.environ.get("RSCRIPT", shutil.which("Rscript") or "Rscript")
        subprocess.run([rscript, str(Path(__file__).with_name("plot_feature_heatmap.R")),
                        str(output_dir)], check=True)


if __name__ == "__main__":
    main()
