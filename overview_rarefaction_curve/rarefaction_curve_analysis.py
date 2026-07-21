import pandas as pd
import numpy as np
import biom
import random
import os
import multiprocessing

# === Path Configuration ===
INPUT_DIR = "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/01-overview-sampling-distribution/02-patch-pipeline-map/input/"
OUTPUT_DIR = "/home/yanglun/YSD/global_oral/global_oral_v2/oral-v2-02/01-overview-sampling-distribution/02-patch-pipeline-map/output/"
BIOM_FILE = os.path.join(INPUT_DIR, "biom_oral_qc.biom")
META_FILE = os.path.join(INPUT_DIR, "metadata_oral_qc.tsv")
OUT_FILE = os.path.join(OUTPUT_DIR, "rarefaction_pcofg.tsv")
TAX_FILE = os.path.join(INPUT_DIR, "taxonomy.tsv")

# Global variables for multiprocessing. 
# Using globals allows Linux 'fork' to share memory without expensive pickling.
MAT = None
FEATURE_TO_TAX = None
KEEP_SAMPLE_INDICES = None

def get_tax_name(tax_list, prefix):
    """Extract the name of a specific taxonomic level (e.g., 'p__' for Phylum) from a taxonomy list."""
    for t in tax_list:
        if t.startswith(prefix):
            val = t[len(prefix):].strip()
            # If the level is incompletely annotated, skip it so it isn't counted as a "known unique" taxon
            if val != "":
                return t 
    return None

def process_step(n_samples):
    """Calculate unique taxa for a specific number of samples (designed for parallel execution)."""
    # Extract the first 'n_samples' from the randomly shuffled sequence
    s_indices = KEEP_SAMPLE_INDICES[:n_samples]

    # Slice the sparse matrix for these specific samples
    sub_mat = MAT[:, s_indices]

    # Find features (rows) that have an abundance > 0 in at least one selected sample.
    # getnnz(axis=1) is highly optimized and counts non-zero values across columns.
    present_features = sub_mat.getnnz(axis=1).nonzero()[0]

    seen_phylum = set()
    seen_class = set()
    seen_order = set()
    seen_family = set()
    seen_genus = set()

    # Update the sets of discovered taxa
    for f_idx in present_features:
        tax = FEATURE_TO_TAX[f_idx]
        if tax['Phylum']: seen_phylum.add(tax['Phylum'])
        if tax['Class']: seen_class.add(tax['Class'])
        if tax['Order']: seen_order.add(tax['Order'])
        if tax['Family']: seen_family.add(tax['Family'])
        if tax['Genus']: seen_genus.add(tax['Genus'])

    return {
        'Samples': n_samples,
        'Phylum': len(seen_phylum),
        'Class': len(seen_class),
        'Order': len(seen_order),
        'Family': len(seen_family),
        'Genus': len(seen_genus)
    }

def main():
    print("1. Loading metadata...")
    # Use native Python to parse line by line, skipping extra columns that cause errors
    valid_samples = set()
    with open(META_FILE, 'r', encoding='utf-8', errors='ignore') as f:
        # Read the header and locate the column index for sample_id
        header = f.readline().strip('\n').split('\t')
        if 'sample_id' not in header:
            raise ValueError("Column 'sample_id' not found in metadata file")
        run_idx = header.index('sample_id')
        
        for line in f:
            parts = line.strip('\n').split('\t')
            # Extract sample_id as long as the row length is sufficient
            if len(parts) > run_idx:
                run_id = parts[run_idx].strip()
                if run_id:
                    valid_samples.add(run_id)

    print("2. Loading BIOM table (this may take a minute)...")
    table = biom.load_table(BIOM_FILE)

    # Filter samples that exist in both metadata and the biom table
    table_samples = set(table.ids(axis='sample'))
    keep_samples = list(valid_samples.intersection(table_samples))
    print(f"   Found {len(keep_samples)} valid matched samples.")

    print("3. Loading external Taxonomy file...")
    # Read the independent taxonomy.tsv and store it in a dictionary for fast querying
    tax_dict_map = {}
    with open(TAX_FILE, 'r', encoding='utf-8', errors='ignore') as f:
        _ = f.readline() # Skip the first header row (Feature ID, Taxon, Confidence, etc.)
        for line in f:
            parts = line.strip('\n').split('\t')
            if len(parts) >= 2:
                feat_id = parts[0].strip()
                taxon_str = parts[1].strip()
                tax_dict_map[feat_id] = taxon_str

    # Randomly shuffle sample order to simulate the random sampling accumulation process
    random.seed(42)
    random.shuffle(keep_samples)

    print("4. Parsing taxonomy for each feature...")
    feature_ids = table.ids(axis='observation')
    feature_to_tax = []

    # Pre-parse the Phylum, Class, Order, Family, and Genus for each OTU/ASV
    for fid in feature_ids:
        tax_dict = {'Phylum': None, 'Class': None, 'Order': None, 'Family': None, 'Genus': None}

        # Get the taxonomy string directly from the newly loaded dictionary
        if fid in tax_dict_map:
            tax_str = tax_dict_map[fid]
            tax = [t.strip() for t in tax_str.split(';')]

            tax_dict['Phylum'] = get_tax_name(tax, 'p__')
            tax_dict['Class'] = get_tax_name(tax, 'c__')
            tax_dict['Order'] = get_tax_name(tax, 'o__')
            tax_dict['Family'] = get_tax_name(tax, 'f__')
            tax_dict['Genus'] = get_tax_name(tax, 'g__')

        feature_to_tax.append(tax_dict)

    print("5. Preparing data for multiprocessing...")
    # Declare globals so worker processes can inherit them
    global MAT, FEATURE_TO_TAX, KEEP_SAMPLE_INDICES
    
    # Convert to sparse matrix to significantly improve slicing computation speed
    MAT = table.matrix_data.tocsc()
    FEATURE_TO_TAX = feature_to_tax
    
    sample_id_to_idx = {sid: idx for idx, sid in enumerate(table.ids(axis='sample'))}
    KEEP_SAMPLE_INDICES = [sample_id_to_idx[sid] for sid in keep_samples]

    total_samples = len(keep_samples)
    
    # Define the sample checkpoints to record the results.
    steps = sorted(list(set(
        [1, 5, 10, 50, 100, 200, 500, 800] +
        list(range(1000, 10000, 1000)) +
        list(range(10000, total_samples, 5000)) +
        [total_samples]
    )))

    # Fetch CPU cores assigned by Slurm, default to available system cores if not in Slurm
    num_cores = int(os.environ.get("SLURM_CPUS_PER_TASK", multiprocessing.cpu_count()))
    print(f"6. Calculating accumulation curve data using {num_cores} parallel processes...")
    
    # Launch multiprocessing pool
    with multiprocessing.Pool(processes=num_cores) as pool:
        results = pool.map(process_step, steps)

    print("7. Saving results...")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df_results = pd.DataFrame(results)
    df_results = df_results.sort_values('Samples')
    df_results.to_csv(OUT_FILE, sep='\t', index=False)
    print(f"   Done! Data successfully saved to {OUT_FILE}")

if __name__ == "__main__":
    main()