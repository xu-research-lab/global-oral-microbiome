#!/usr/bin/env python3
"""
Extract PCoA coordinates and eigenvalues from QIIME2 PCoAResults qza file.
Saves as RData-compatible TSV files for downstream db-RDA analysis.

Usage: python extract_pcoa_from_qza.py <qza_file> <output_dir> [max_axes]
  qza_file:    Path to bray_curtis_pcoa_results.qza
  output_dir:  Directory to save output files
  max_axes:    Number of top PCoA axes to extract (default: 1000)
"""

import sys
import os
import zipfile
import numpy as np

def read_non_empty(f):
    """Read next non-empty line from file."""
    while True:
        line = f.readline()
        if not line:
            return None
        decoded = line.decode().strip()
        if decoded:
            return decoded

def extract_pcoa(qza_path, output_dir, max_axes=1000):
    os.makedirs(output_dir, exist_ok=True)

    with zipfile.ZipFile(qza_path) as z:
        names = [n for n in z.namelist() if n.endswith('ordination.txt')]
        if not names:
            print("ERROR: No ordination.txt found in qza file")
            sys.exit(1)
        ord_path = names[0]
        print(f"Found ordination file: {ord_path}")

        with z.open(ord_path) as f:
            # --- Eigenvalues ---
            header = read_non_empty(f)
            print(f"Eigenvals header: {header}")
            n_axes_total = int(header.split('\t')[1])

            eig_line = read_non_empty(f)
            eigenvalues = np.array([float(x) for x in eig_line.split('\t')])
            print(f"Eigenvalues loaded: {len(eigenvalues)} values")

            # --- Proportion explained ---
            prop_header = read_non_empty(f)
            print(f"Proportion explained header: {prop_header}")

            prop_line = read_non_empty(f)
            prop_explained = np.array([float(x) for x in prop_line.split('\t')])
            print(f"Proportion explained loaded: {len(prop_explained)} values")

            # --- Species ---
            species_header = read_non_empty(f)
            print(f"Species header: {species_header}")
            species_parts = species_header.split('\t')
            n_species = int(species_parts[1]) if len(species_parts) > 1 else 0

            if n_species > 0:
                for _ in range(n_species):
                    read_non_empty(f)

            # --- Site ---
            site_header = read_non_empty(f)
            print(f"Site header: {site_header}")
            site_parts = site_header.split('\t')
            n_samples = int(site_parts[1])
            n_site_axes = int(site_parts[2])
            print(f"Site data: {n_samples} samples x {n_site_axes} axes")

            actual_axes = min(max_axes, n_site_axes)
            print(f"Extracting top {actual_axes} axes for {n_samples} samples...")

            # Read site coordinates
            sample_ids = []
            coordinates = []

            line_count = 0
            while True:
                line = f.readline()
                if not line:
                    break
                decoded = line.decode().strip()
                if not decoded:
                    continue
                parts = decoded.split('\t')
                sample_ids.append(parts[0])
                coords = [float(x) for x in parts[1:actual_axes+1]]
                coordinates.append(coords)
                line_count += 1
                if line_count % 10000 == 0:
                    print(f"  Read {line_count} samples...")

            print(f"Total samples read: {len(sample_ids)}")

    # Save eigenvalues
    eig_file = os.path.join(output_dir, 'pcoa_eigenvalues.tsv')
    with open(eig_file, 'w') as fout:
        fout.write('Axis\tEigenvalue\tVarianceExplained\tCumulativeVariance\n')
        total_eig = sum(eigenvalues[eigenvalues > 0])
        cum_var = 0
        for i in range(actual_axes):
            eig = eigenvalues[i]
            prop = eig / total_eig * 100 if total_eig > 0 else 0
            cum_var += prop
            fout.write(f'{i+1}\t{eig}\t{prop}\t{cum_var}\n')
    print(f"Saved eigenvalues: {eig_file}")

    # Save PCoA coordinates
    coord_file = os.path.join(output_dir, 'pcoa_coordinates.tsv')
    coord_header = 'SampleID\t' + '\t'.join([f'PCoA{i+1}' for i in range(actual_axes)])
    with open(coord_file, 'w') as fout:
        fout.write(coord_header + '\n')
        for sid, coords in zip(sample_ids, coordinates):
            coords_str = '\t'.join([f'{c}' for c in coords])
            fout.write(f'{sid}\t{coords_str}\n')
    print(f"Saved coordinates: {coord_file}")

    # Summary
    print(f"\n=== Summary ===")
    print(f"Samples: {len(sample_ids)}")
    print(f"Axes extracted: {actual_axes}")
    print(f"Top 3 eigenvalues: {eigenvalues[:3]}")
    cum_var_10 = sum(eigenvalues[:10]) / total_eig * 100
    cum_var_50 = sum(eigenvalues[:50]) / total_eig * 100
    cum_var_100 = sum(eigenvalues[:100]) / total_eig * 100
    cum_var_all = sum(eigenvalues[:actual_axes]) / total_eig * 100
    print(f"Cumulative variance - top 10: {cum_var_10:.2f}%")
    print(f"Cumulative variance - top 50: {cum_var_50:.2f}%")
    print(f"Cumulative variance - top 100: {cum_var_100:.2f}%")
    print(f"Cumulative variance - top {actual_axes}: {cum_var_all:.2f}%")
    print(f"Negative eigenvalues: {sum(eigenvalues < 0)}")

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python extract_pcoa_from_qza.py <qza_file> <output_dir> [max_axes]")
        sys.exit(1)

    qza = sys.argv[1]
    outdir = sys.argv[2]
    max_a = int(sys.argv[3]) if len(sys.argv) > 3 else 1000

    extract_pcoa(qza, outdir, max_a)
