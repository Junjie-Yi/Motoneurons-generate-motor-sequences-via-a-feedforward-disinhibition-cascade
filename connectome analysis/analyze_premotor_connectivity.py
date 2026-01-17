import pandas as pd
from fafbseg import flywire
import numpy as np
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
        df_mn = get_upstream_connections(mn_all, min_weight=min_weight_query)
        
        if df_mn.empty:
            print("No premotor neurons found.")
            return

        df_mn['signed_weight'] = df_mn.apply(apply_signed_weight, axis=1)
        
        # Save processed data
        df_mn.to_csv(file_mn, index=False)
        print(f"Saved {file_mn}")
    
    df_mn = df_mn[~df_mn['pre'].isin(mn_all)]
    print(f"Filtered df_mn (excluding MN sources): {len(df_mn)} connections.")

    mat_direct = save_connectivity_matrix(df_mn, title_suffix="Premotor to MN", min_weight_filter=min_weight_query, reindex_cols=mn_all)
    
    if mat_direct is None:
        print("No direct connectivity matrix generated.")
        return

if __name__ == "__main__":
    main()
