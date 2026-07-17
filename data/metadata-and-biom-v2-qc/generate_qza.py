#!/usr/bin/env python3
"""
Generate QZA files from the QC-passed BIOM tables.

Requires: qiime2-2022.2 environment (conda activate qiime2-2022.2)

Outputs:
  - biom_all_qc.qza               (FeatureTable[Frequency])
  - genus_table_all_qc.qza         (FeatureTable[Frequency] — requires taxonomy)
  - genus_table_relabund_all_qc.qza (FeatureTable[RelativeFrequency])
"""
import os
import subprocess
from pathlib import Path

OUT_DIR = '/beegfs/huangxiaochang/global_oral/global-oral-microbiome/data/metadata-and-biom-v2-qc'
BIOM_ALL_QC = f'{OUT_DIR}/biom_all_qc.biom'
QZA_ALL_QC = f'{OUT_DIR}/biom_all_qc.qza'
QZA_GENUS = f'{OUT_DIR}/genus_table_all_qc.qza'
QZA_GENUS_RELABUND = f'{OUT_DIR}/genus_table_relabund_all_qc.qza'

# GG2 taxonomy files
GG2_TAXONOMY_TSV = '/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv'


def run_qiime(cmd_args, description):
    print(f'  {description}...')
    cmd = ['qiime'] + cmd_args
    print(f'    Running: {" ".join(cmd)}')
    subprocess.run(cmd, check=True)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # Step 1: Import BIOM → FeatureTable[Frequency] QZA
    print('=== Step 1: Import BIOM → FeatureTable[Frequency] ===')
    if os.path.exists(QZA_ALL_QC):
        os.remove(QZA_ALL_QC)
    run_qiime(
        ['tools', 'import', '--input-path', BIOM_ALL_QC,
         '--output-path', QZA_ALL_QC, '--type', 'FeatureTable[Frequency]'],
        'Import BIOM to FeatureTable[Frequency]'
    )

    # Step 2: Import taxonomy TSV as FeatureData[Taxonomy]
    print('\n=== Step 2: Import taxonomy and collapse to genus ===')
    tax_qza = f'{OUT_DIR}/gg2_taxonomy.qza'
    if os.path.exists(tax_qza):
        os.remove(tax_qza)
    run_qiime(
        ['tools', 'import', '--input-path', GG2_TAXONOMY_TSV,
         '--output-path', tax_qza, '--type', 'FeatureData[Taxonomy]'],
        'Import taxonomy TSV to QZA'
    )

    # Collapse by taxonomy level 6 (genus)
    if os.path.exists(QZA_GENUS):
        os.remove(QZA_GENUS)
    run_qiime(
        ['taxa', 'collapse', '--i-table', QZA_ALL_QC,
         '--i-taxonomy', tax_qza, '--p-level', '6',
         '--o-collapsed-table', QZA_GENUS],
        'Collapse to genus level'
    )
    print(f'    Wrote {QZA_GENUS}')

    # Step 3: Relative abundance transform
    print('\n=== Step 3: Relative abundance transform ===')
    if os.path.exists(QZA_GENUS_RELABUND):
        os.remove(QZA_GENUS_RELABUND)
    run_qiime(
        ['feature-table', 'relative-frequency',
         '--i-table', QZA_GENUS,
         '--o-relative-frequency-table', QZA_GENUS_RELABUND],
        'Convert to relative abundance'
    )
    print(f'    Wrote {QZA_GENUS_RELABUND}')

    print('\n=== Done ===')
    print(f'  {QZA_ALL_QC}')
    print(f'  {QZA_GENUS}')
    print(f'  {QZA_GENUS_RELABUND}')


if __name__ == '__main__':
    main()
