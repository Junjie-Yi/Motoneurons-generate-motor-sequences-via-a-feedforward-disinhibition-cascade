import pandas as pd
import matplotlib.pyplot as plt
from neuprint import Client, fetch_adjacencies
from fafbseg import flywire
import os

# --- Configuration ---
# IDs
fafb_id = 720575940618539810
mcns_id = 10619
banc_id = 720575941501585899

# Output files
output_plot = 'upstream_accumulation_comparison.svg'
output_csv_prefix = 'upstream_partners_'

# Thresholds
mcns_threshold = 5
banc_threshold = 3
fafb_threshold = 5

# --- Helper function to get data ---
def get_data(dataset_name, file_suffix, query_func):
    filename = f'{output_csv_prefix}{file_suffix}.csv'
    if os.path.exists(filename):
        print(f"Loading {dataset_name} from local file: {filename}")
        return pd.read_csv(filename)
    else:
        print(f"Querying {dataset_name}...")
        df = query_func()
        # Save raw(ish) results before returning
        df.to_csv(filename, index=False)
        print(f"Saved {dataset_name} data to {filename}")
        return df

# --- Query Functions ---

def query_mcns():
    # Token from previous context
    c = Client('neuprint.janelia.org', dataset='male-cns:v0.9', token='eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6IjQ3OTkyMTEwM0BxcS5jb20iLCJsZXZlbCI6Im5vYXV0aCIsImltYWdlLXVybCI6Imh0dHBzOi8vbGgzLmdvb2dsZXVzZXJjb250ZW50LmNvbS9hL0FDZzhvY0xfWWVGbjNjX3hqWm9uRGFOVUdSN0tzVkJ6NjZfTDM5QnRodW15ZjB2SEc2V2NnUT1zOTYtYz9zej01MD9zej01MCIsImV4cCI6MTk0NjM4NDg5NX0.L5odyYSB76TzaDUl99UMYqkATJqOrLfQkRmmb5BNCgk')
    neuron_df, conn_df = fetch_adjacencies(sources=None, targets=[mcns_id])
    if not conn_df.empty:
        data = conn_df.groupby('bodyId_pre')['weight'].sum().reset_index()
        data.columns = ['partner_id', 'weight']
        return data
    return pd.DataFrame(columns=['partner_id', 'weight'])

def query_banc():
    banc_csv = 'connections_princeton.csv'
    if os.path.exists(banc_csv):
        df_raw = pd.read_csv(banc_csv)
        df_raw.columns = df_raw.columns.str.strip()
        partners = df_raw[df_raw['post_root_id'] == banc_id]
        data = partners.groupby('pre_root_id')['syn_count'].sum().reset_index()
        data.columns = ['partner_id', 'weight']
        return data
    else:
        print(f"Warning: {banc_csv} not found.")
        return pd.DataFrame(columns=['partner_id', 'weight'])

def query_fafb():
    try:
        conn = flywire.get_connectivity([fafb_id], filtered=False, min_score=50, upstream=True)
        if not conn.empty:
            data = conn.groupby('pre')['weight'].sum().reset_index()
            data.columns = ['partner_id', 'weight']
            return data
    except Exception as e:
        print(f"Error querying FAFB: {e}")
    return pd.DataFrame(columns=['partner_id', 'weight'])

# --- Main Execution ---

# 1. Get Data (Local or Query)
mcns_data = get_data('MCNS', 'MCNS', query_mcns)
banc_data = get_data('BANC', 'BANC', query_banc)
fafb_data = get_data('FAFB', 'FAFB', query_fafb)

# 2. Apply Thresholds
print("\nApplying thresholds...")
mcns_data = mcns_data[mcns_data['weight'] >= mcns_threshold]
print(f"MCNS (Threshold {mcns_threshold}): {len(mcns_data)} partners remaining.")

banc_data = banc_data[banc_data['weight'] >= banc_threshold]
print(f"BANC (Threshold {banc_threshold}): {len(banc_data)} partners remaining.")

fafb_data = fafb_data[fafb_data['weight'] >= fafb_threshold]
print(f"FAFB (Threshold {fafb_threshold}): {len(fafb_data)} partners remaining.")

# 3. Visualization
print("\nGenerating visualization...")
plt.figure(figsize=(5, 1.8))

# Helper to convert 0-255 RGB to 0-1
def rgb(r, g, b):
    return (r/255.0, g/255.0, b/255.0)

datasets = [
    ('MCNS', mcns_data, rgb(157, 215, 157)),    # Blue-ish
    ('BANC', banc_data, rgb(242, 225, 4)),   # User example (Orange-ish)
    ('FAFB', fafb_data, rgb(253, 200, 151))      # Green-ish
]

for name, data, color in datasets:
    if data.empty:
        continue
        
    # Sort weights descending
    sorted_weights = data['weight'].sort_values(ascending=False).reset_index(drop=True)
    
    # Calculate accumulated weight
    accumulated_weight = sorted_weights.cumsum()
    
    # Plot
    plt.plot(accumulated_weight.index + 1, accumulated_weight.values, label=name, color=color, linewidth=2)

#plt.xlabel('Rank of Upstream Partner (Sorted by Weight)')
#plt.ylabel('Accumulated Connection Weight')
#plt.title('Accumulated Upstream Connection Weight Comparison')
#plt.legend()
plt.grid(False)
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xscale('linear') # Explicitly set to linear as requested

plt.savefig(output_plot)
print(f"Plot saved to {output_plot}")
# plt.show() # Cannot show in this environment
