import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from fafbseg import flywire
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr
import os

# Define MN IDs (from find_disinhibition_paths.py)
mn12v = [720575940661338497,720575940632281927]
mn12d = [720575940630233404,720575940618539810]
mn11v = [720575940645527918,720575940615197747]
mn11d = [720575940630868793,720575940618165019]
mn10 = [720575940635360924]
mn_all = mn12v + mn12d + mn11v + mn11d + mn10

def get_upstream_connections(target_ids, min_weight=1):
    """
    Query upstream connections for target_ids.
    Returns a DataFrame with columns: pre, post, weight, pre_nt
    """
    print(f"Querying upstream connections for {len(target_ids)} targets...")
    
    # 1. Get connectivity
    try:
        # Using get_connectivity with upstream=True
        # Note: min_score is cleft score, usually 30 or 50. 
        # The user mentioned "synaptic weight", which is usually synapse count.
        # We can filter by weight later or use filtered=True if available.
        conn = flywire.get_connectivity(
            target_ids,
            upstream=True,
            downstream=False,
            filtered=True, # Use filtered if available to remove low confidence
            min_score=50
        )
    except Exception as e:
        print(f"Error fetching connectivity: {e}")
        return pd.DataFrame()

    if conn.empty:
        print("No connections found.")
        return pd.DataFrame()

    # Filter by min_weight (synapse count)
    if min_weight > 0:
        conn = conn[conn['weight'] >= min_weight]
    
    if conn.empty:
        print(f"No connections with weight >= {min_weight}.")
        return pd.DataFrame()

    # Save raw connectivity
    conn.to_csv(f"raw_connectivity_{len(target_ids)}_targets.csv", index=False)

    # 2. Get annotations for source neurons (pre)
    pre_ids = conn['pre'].unique().tolist()
    print(f"Fetching annotations for {len(pre_ids)} source neurons...")
    
    try:
        # Assuming search_annotations is available as per user's existing code
        # If list is too long, we might need to chunk it, but let's try direct first
        annotations = flywire.search_annotations(pre_ids)
    except Exception as e:
        print(f"Error fetching annotations: {e}")
        return pd.DataFrame()

    if annotations.empty:
        print("No annotations found.")
        return pd.DataFrame()
    
    # Save annotations
    annotations.to_csv(f"annotations_{len(pre_ids)}_neurons.csv", index=False)

    # Merge annotations to get NT
    # Expecting 'root_id' and 'top_nt' in annotations
    if 'top_nt' not in annotations.columns:
        print("Annotations missing 'top_nt' column.")
        return pd.DataFrame()

    conn_annotated = conn.merge(
        annotations[['root_id', 'top_nt']],
        left_on='pre',
        right_on='root_id',
        how='left'
    )
    
    return conn_annotated

def apply_signed_weight(row):
    nt = str(row['top_nt']).lower()
    weight = row['weight']
    
    if nt in ['ach', 'acetylcholine']:
        return weight # Positive
    elif nt in ['gaba', 'glutamate']:
        return -weight # Negative
    else:
        return 0 # Others

def save_connectivity_matrix(df, title_suffix="", min_weight_filter=0, reindex_cols=None):
    """
    df: DataFrame with 'pre', 'post', 'signed_weight'
    reindex_cols: list of column IDs to ensure are present (fills with 0)
    """
    if df.empty:
        print("No data to save.")
        return None

    # Pivot to matrix: rows=pre, cols=post
    # If there are multiple connections (shouldn't be for unique pre-post pairs in simple connectivity), sum them
    matrix = df.pivot_table(index='pre', columns='post', values='signed_weight', aggfunc='sum', fill_value=0)
    
    # Reindex columns if provided (to ensure all targets are present)
    if reindex_cols is not None:
        # Ensure reindex_cols are same type as columns if possible, usually both are ints or strings
        # If matrix columns are ints and reindex_cols are ints, it works.
        matrix = matrix.reindex(columns=reindex_cols, fill_value=0)

    # Filter matrix by absolute weight if needed
    if min_weight_filter > 0:
        # Filter rows (sources) that have at least one connection with abs(weight) >= min_weight_filter
        mask = matrix.abs().max(axis=1) >= min_weight_filter
        matrix = matrix.loc[mask]
        
        if matrix.empty:
            print(f"No connections left after filtering with min_weight={min_weight_filter}")
            return None

    print(f"Saving matrix for {matrix.shape[0]} sources and {matrix.shape[1]} targets...")
    
    try:
        filename = f"matrix_{title_suffix.replace(' ', '_')}_min{min_weight_filter}.csv"
        matrix.to_csv(filename)
        print(f"Saved matrix to {filename}")
    except Exception as e:
        print(f"Error saving matrix: {e}")
        
    return matrix

def compare_input_similarity(mat_direct, mat_effective):
    """
    Compare the similarity structure of inputs to MNs between Direct (Premotor->MN) 
    and Effective (Upstream->MN) pathways.
    """
    print("\n--- Comparing Input Similarity (RSA) ---")
    
    if mat_direct is None or mat_effective is None:
        print("One of the matrices is None.")
        return

    # 1. Align columns (MNs)
    common_mns = mat_direct.columns.intersection(mat_effective.columns)
    mat_dir = mat_direct[common_mns].fillna(0)
    mat_eff = mat_effective[common_mns].fillna(0)
    
    # Filter columns with non-zero variance for correlation metric
    valid_cols = (mat_dir.std() > 0) & (mat_eff.std() > 0)
    mat_dir_rsa = mat_dir.loc[:, valid_cols]
    mat_eff_rsa = mat_eff.loc[:, valid_cols]
    
    if mat_dir_rsa.shape[1] < 3:
        print(f"Not enough valid MNs (with non-zero variance) for RSA. Valid: {mat_dir_rsa.shape[1]}")
    else:
        print(f"Using {mat_dir_rsa.shape[1]} MNs for RSA comparison.")
        # Compute distance matrices (MN x MN)
        # Transpose: Rows=MNs
        d_dir = pdist(mat_dir_rsa.T, metric='correlation')
        d_eff = pdist(mat_eff_rsa.T, metric='correlation')
        
        # Correlate
        r, p = pearsonr(d_dir, d_eff)
        print(f"RSA Correlation (Input Pattern Similarity): r={r:.3f}, p={p:.3e}")
    
    # 2. Total Drive Correlation (using all common MNs)
    drive_dir = mat_dir.sum(axis=0)
    drive_eff = mat_eff.sum(axis=0)
    
    if len(drive_dir) > 1:
        r_drive, p_drive = pearsonr(drive_dir, drive_eff)
        print(f"Total Drive Correlation (n={len(drive_dir)}): r={r_drive:.3f}, p={p_drive:.3e}")
    else:
        print("Not enough MNs for drive correlation.")

def calculate_effective_connectivity(df_upstream, df_downstream):
    """
    Calculate effective connectivity A -> B -> C.
    df_upstream: A -> B (cols: pre, post, weight/signed_weight)
    df_downstream: B -> C (cols: pre, post, signed_weight)
    
    Returns: DataFrame of effective weights (A -> C)
    """
    print("Calculating effective connectivity...")
    
    # 1. Create Matrix A -> B
    # We assume df_upstream has 'weight' and we want to treat it as signed based on NT if available, 
    # or just use 'weight' if we know they are all Ach (positive).
    # The user said "Only find the ach upstream", so they are Ach -> Positive.
    # So we can just use 'weight' from df_upstream as the value.
    
    # Pivot upstream: index=pre(A), columns=post(B)
    # Sum weights if multiple edges
    mat_ab = df_upstream.pivot_table(index='pre', columns='post', values='weight', aggfunc='sum', fill_value=0)
    
    # 2. Create Matrix B -> C
    # df_downstream has 'signed_weight'
    mat_bc = df_downstream.pivot_table(index='pre', columns='post', values='signed_weight', aggfunc='sum', fill_value=0)
    
    # 3. Align matrices on B (columns of AB, rows of BC)
    # Find common B neurons
    common_b = mat_ab.columns.intersection(mat_bc.index)
    print(f"Number of intermediate neurons (B) matching: {len(common_b)}")
    
    if len(common_b) == 0:
        print("No matching intermediate neurons.")
        return pd.DataFrame()
        
    mat_ab_aligned = mat_ab[common_b]
    mat_bc_aligned = mat_bc.loc[common_b]
    
    # 4. Matrix Multiplication
    # (A x B) @ (B x C) -> (A x C)
    mat_ac = mat_ab_aligned.dot(mat_bc_aligned)
    
    return mat_ac

def main():
    # --- Part 1: Premotor -> MN ---
    print("\n--- Part 1: Analyzing Premotor -> MN ---")
    
    min_weight_query = 5 
    file_mn = "premotor_to_mn_processed.csv"
    
    if os.path.exists(file_mn):
        print(f"Loading existing file: {file_mn}")
        df_mn = pd.read_csv(file_mn)
    else:
        # B -> C
        df_mn = get_upstream_connections(mn_all, min_weight=min_weight_query)
        
        if df_mn.empty:
            print("No premotor neurons found.")
            return

        # Set signed weight for B -> C
        df_mn['signed_weight'] = df_mn.apply(apply_signed_weight, axis=1)
        
        # Save processed data
        df_mn.to_csv(file_mn, index=False)
        print(f"Saved {file_mn}")
    
    # Filter out MNs from sources (pre) in df_mn
    df_mn = df_mn[~df_mn['pre'].isin(mn_all)]
    print(f"Filtered df_mn (excluding MN sources): {len(df_mn)} connections.")

    # Save B->C matrix
    # Ensure all MNs are columns
    mat_direct = save_connectivity_matrix(df_mn, title_suffix="Premotor to MN", min_weight_filter=min_weight_query, reindex_cols=mn_all)
    
    if mat_direct is None:
        print("No direct connectivity matrix generated.")
        return

    # --- Part 2: Ach Upstream -> Premotor -> MN ---
    print("\n--- Part 2: Analyzing Ach Upstream -> Premotor -> MN ---")
    
    # 1. Identify Premotor Neurons (B)
    # Use the rows of mat_direct to ensure consistency
    premotor_ids = mat_direct.index.tolist()
    print(f"Found {len(premotor_ids)} premotor neurons (from direct matrix).")
    
    file_upstream = "upstream_to_premotor_all.csv"
    
    if os.path.exists(file_upstream):
        print(f"Loading existing file: {file_upstream}")
        df_upstream = pd.read_csv(file_upstream)
    else:
        # 2. Query Upstream of Premotor (A -> B)
        # This might be large, so we might want to batch it if it fails, but let's try direct.
        df_upstream = get_upstream_connections(premotor_ids, min_weight=min_weight_query)
        
        if not df_upstream.empty:
            df_upstream.to_csv(file_upstream, index=False)
            print(f"Saved {file_upstream}")
    
    if df_upstream.empty:
        print("No upstream connections to premotor found.")
        return
        
    # 3. Filter for Ach Upstream (A is Ach)
    # Check if 'top_nt' exists
    if 'top_nt' in df_upstream.columns:
        ach_mask = df_upstream['top_nt'].isin(['ach', 'acetylcholine'])
        df_upstream_ach = df_upstream[ach_mask].copy()
        print(f"Filtered for Ach upstream: {len(df_upstream_ach)} connections from {df_upstream_ach['pre'].nunique()} Ach neurons.")
    else:
        print("No NT info for upstream, cannot filter for Ach.")
        return

    if df_upstream_ach.empty:
        print("No Ach upstream neurons found.")
        return

    # Filter out MNs from sources (pre) in df_upstream_ach
    df_upstream_ach = df_upstream_ach[~df_upstream_ach['pre'].isin(mn_all)]
    print(f"Filtered Ach upstream (excluding MNs): {len(df_upstream_ach)} connections.")

    if df_upstream_ach.empty:
        print("No Ach upstream neurons found after excluding MNs.")
        return

    # Save Matrix for Ach Upstream -> Premotor
    # Since we filtered for Ach, signed_weight is just weight
    df_upstream_ach['signed_weight'] = df_upstream_ach['weight']
    # Ensure columns match the premotor neurons in mat_direct
    save_connectivity_matrix(df_upstream_ach, title_suffix="Ach Upstream to Premotor", min_weight_filter=min_weight_query, reindex_cols=premotor_ids)

    # Save Matrix for Ach Upstream -> Ach Premotor
    ach_premotor_ids = df_mn[df_mn['top_nt'].isin(['ach', 'acetylcholine'])]['pre'].unique().tolist()
    ach_premotor_ids = [x for x in ach_premotor_ids if x in premotor_ids]
    print(f"Found {len(ach_premotor_ids)} Ach Premotor neurons for sub-matrix.")
    
    df_upstream_ach_ach = df_upstream_ach[df_upstream_ach['post'].isin(ach_premotor_ids)]
    save_connectivity_matrix(df_upstream_ach_ach, title_suffix="Ach Upstream to Ach Premotor", min_weight_filter=min_weight_query, reindex_cols=ach_premotor_ids)

    # 4. Calculate Effective Connectivity (A -> C)
    # A (Ach) -> B (Any) -> C (MN)
    # Weight A->B is positive (Ach). Weight B->C is signed.
    mat_effective = calculate_effective_connectivity(df_upstream_ach, df_mn)
    
    if not mat_effective.empty:
        # Save Matrix
        filename = f"matrix_effective_AchUpstream_to_MN_min{min_weight_query}.csv"
        mat_effective.to_csv(filename)
        print(f"Saved effective connectivity matrix to {filename}")
        
        # Compare Similarity
        compare_input_similarity(mat_direct, mat_effective)
        
        # Also save a filtered version
        # Filter rows with low max effective weight
        min_eff_weight = 20 # Arbitrary threshold for "effective"
        mask = mat_effective.abs().max(axis=1) >= min_eff_weight
        mat_effective_filtered = mat_effective.loc[mask]
        
        if not mat_effective_filtered.empty:
            filename_filt = f"matrix_effective_AchUpstream_to_MN_min{min_weight_query}_filtered{min_eff_weight}.csv"
            mat_effective_filtered.to_csv(filename_filt)
            print(f"Saved filtered effective matrix to {filename_filt}")

if __name__ == "__main__":
    main()
