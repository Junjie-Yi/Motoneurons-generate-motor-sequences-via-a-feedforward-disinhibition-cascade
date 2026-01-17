import pandas as pd
from fafbseg import flywire
import os

# Neuron IDs copied from find_disinhibition_paths.py
post12v = [720575940643867296,720575940633772257,720575940639769715,720575940629804158,720575940627359109,720575940625002936,720575940632028216,
              720575940630232380,720575940626163658,720575940611548878,720575940621047632,720575940615410911]
post12d = [720575940643867296,720575940627359109,720575940623149581,720575940632028216,720575940611548878,720575940621047632,720575940633772257,
           720575940630683385,720575940619016059,720575940624289789]

mn12v = [720575940661338497,720575940632281927]
mn12d = [720575940630233404,720575940618539810]
mn11v = [720575940645527918,720575940615197747]
mn11d = [720575940630868793,720575940618165019]
mn_all = mn12v+mn12d+mn11v+mn11d

def find_paths_logic(source_ids, target_ids, min_synapse_num, nt_type='ACH'):
    """
    Finds pathways from source_ids -> intermediate (filtered by NT) -> target_ids.
    Returns a DataFrame with the path details.
    """
    sources = source_ids
    targets = target_ids
    
    if not sources or not targets:
        return pd.DataFrame()

    print(f"... Fetching downstream connections from {len(sources)} sources")
    try:
        conn_A_B = flywire.get_connectivity(sources, upstream=False, downstream=True, filtered=True, min_score=50)
        conn_A_B = conn_A_B[conn_A_B['weight'] >= min_synapse_num]
        conn_A_B = conn_A_B[conn_A_B['pre'].isin(sources)]
    except Exception as e:
        print(f"Error fetching downstream connectivity: {e}")
        return pd.DataFrame()
        
    if conn_A_B.empty:
        return pd.DataFrame()
        
    potential_intermediates = conn_A_B['post'].unique().tolist()
    
    print(f"... Checking annotations for {len(potential_intermediates)} potential intermediates")
    try:
        potential_intermediates = [int(x) for x in potential_intermediates]
        annotations = flywire.search_annotations(potential_intermediates)
    except Exception as e:
        print(f"Error fetching annotations: {e}")
        return pd.DataFrame()
        
    if annotations.empty or 'top_nt' not in annotations.columns:
        return pd.DataFrame()
    
    if nt_type.upper() == 'ACH':
        target_nt = 'acetylcholine'
    elif nt_type.upper() == 'GABA':
        target_nt = 'gaba'
    else:
        target_nt = nt_type.lower()
        
    # Get intermediates that match NT and are NOT in the MN list
    filtered_intermediates = annotations[annotations['top_nt'] == target_nt]['root_id'].tolist()
    filtered_intermediates = [x for x in filtered_intermediates if x not in mn_all]
    
    if not filtered_intermediates:
        return pd.DataFrame()
        
    conn_A_B_filtered = conn_A_B[conn_A_B['post'].isin(filtered_intermediates)]
    
    print(f"... Fetching upstream connections to {len(targets)} targets")
    try:
        conn_B_C = flywire.get_connectivity(targets, upstream=True, downstream=False, filtered=True, min_score=50)
        conn_B_C = conn_B_C[conn_B_C['weight'] >= min_synapse_num]
    except Exception as e:
        print(f"Error fetching upstream connectivity: {e}")
        return pd.DataFrame()
        
    if conn_B_C.empty:
        return pd.DataFrame()
        
    conn_B_C_filtered = conn_B_C[conn_B_C['pre'].isin(filtered_intermediates)]
    
    if conn_B_C_filtered.empty:
        return pd.DataFrame()

    # Merge to get full paths: Source -> Intermediate -> Target
    merged = pd.merge(
        conn_A_B_filtered, 
        conn_B_C_filtered, 
        left_on='post', 
        right_on='pre', 
        suffixes=('_A_B', '_B_C')
    )
    
    paths = merged[[
        'pre_A_B', 'post_A_B', 'post_B_C', 'weight_A_B', 'weight_B_C'
    ]].copy()
    
    paths.columns = ['source_id', 'intermediate_id', 'target_id', 'weight_source_to_inter', 'weight_inter_to_target']
    return paths

def run_query(source_name, source_ids, target_name, target_ids, min_synapse=5):
    print(f"\nQuerying: {source_name} -> ACH -> {target_name} (min_synapse_num={min_synapse})")
    df = find_paths_logic(source_ids, target_ids, min_synapse, 'ACH')
    
    if df.empty:
        print("No paths found.")
        print("-" * 50)
        return
        
    # 1. Number of intermediate neurons
    # Count unique intermediate IDs in the valid paths
    intermediate_neurons = df['intermediate_id'].nunique()
    
    # 2. Source -> Intermediate stats
    # We need unique edges (Source, Intermediate)
    edges_source_inter = df[['source_id', 'intermediate_id', 'weight_source_to_inter']].drop_duplicates()
    total_conn_source_inter = len(edges_source_inter)
    total_weight_source_inter = edges_source_inter['weight_source_to_inter'].sum()
    
    # 3. Intermediate -> Target stats
    # We need unique edges (Intermediate, Target) to ensure weight is counted only once per pair
    edges_inter_target = df[['intermediate_id', 'target_id', 'weight_inter_to_target']].drop_duplicates()
    total_conn_inter_target = len(edges_inter_target)
    total_weight_inter_target = edges_inter_target['weight_inter_to_target'].sum()
    
    print("-" * 50)
    print(f"Statistics for {source_name} -> ACH -> {target_name}:")
    print(f"1. Intermediate Neurons (Count):     {intermediate_neurons}")
    print(f"2. {source_name:10} -> Intermediate: {total_conn_source_inter} connections, {total_weight_source_inter} total weight")
    print(f"3. Intermediate -> {target_name:10}: {total_conn_inter_target} connections, {total_weight_inter_target} total weight")
    print("-" * 50)

if __name__ == "__main__":
    # You can adjust this threshold
    min_syn = 5 
    
    # Query 1-3: post12v -> ACH -> mn12d
    run_query("post12v", post12v, "mn12d", mn12d, min_syn)
    
    # Query 4: post12d -> ACH -> mn11v
    run_query("post12d", post12d, "mn11v", mn11v, min_syn)
