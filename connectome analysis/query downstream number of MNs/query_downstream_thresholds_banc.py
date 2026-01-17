import pandas as pd

# Load the CSV file
csv_file = 'connections_princeton.csv'
print(f"Loading {csv_file}...")
df = pd.read_csv(csv_file)

# Clean column names (remove whitespace)
df.columns = df.columns.str.strip()

# Define the neuron IDs to query (pre_root_id)
# Using the ID found in the first few rows of the CSV as an example
neuron_ids = [
    720575941474569521,720575941501585899,720575941573843464,720575941581804383,720575941652694769,720575941617728383,720575941504654274
]

# Define thresholds to test
thresholds = range(2, 11)  # Testing thresholds from 2 to 10

results = []

print(f"Analyzing downstream connections for {len(neuron_ids)} neurons...")

# Filter the dataframe for the specific source neurons first to optimize
# We group by pre_root_id and post_root_id to sum syn_count across different neuropils if needed
# However, the CSV seems to have one row per connection per neuropil. 
# Usually, for a total connection weight, we sum syn_count for the same (pre, post) pair.
filtered_df = df[df['pre_root_id'].isin(neuron_ids)]

# Group by pre and post to get total synapse count per connection (summing over neuropils)
connections = filtered_df.groupby(['pre_root_id', 'post_root_id'])['syn_count'].sum().reset_index()

for threshold in thresholds:
    # Filter connections that meet the threshold
    valid_connections = connections[connections['syn_count'] >= threshold]
    
    # Count unique downstream partners (post_root_id)
    num_downstream = valid_connections['post_root_id'].nunique()
    
    results.append({
        'threshold': threshold,
        'downstream_count': num_downstream
    })
    
    # Optional: Print progress for specific steps
    if threshold % 5 == 0:
        print(f"Threshold {threshold}: {num_downstream} partners")

# Create DataFrame for results
results_df = pd.DataFrame(results)

# Print the results
print("\nResults:")
print(results_df)

# Optionally save to CSV
output_file = 'downstream_counts_from_csv.csv'
results_df.to_csv(output_file, index=False)
print(f"\nResults saved to {output_file}")
