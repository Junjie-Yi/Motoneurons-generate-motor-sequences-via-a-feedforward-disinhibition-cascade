import pandas as pd
from neuprint import Client, fetch_adjacencies

# Configure the client
# NOTE: Please replace 'YOUR_TOKEN_HERE' with your actual neuPrint auth token
# and ensure the dataset name is correct for the male-CNS dataset you are using.
c = Client('neuprint.janelia.org', dataset='male-cns:v0.9', token='eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6IjQ3OTkyMTEwM0BxcS5jb20iLCJsZXZlbCI6Im5vYXV0aCIsImltYWdlLXVybCI6Imh0dHBzOi8vbGgzLmdvb2dsZXVzZXJjb250ZW50LmNvbS9hL0FDZzhvY0xfWWVGbjNjX3hqWm9uRGFOVUdSN0tzVkJ6NjZfTDM5QnRodW15ZjB2SEc2V2NnUT1zOTYtYz9zej01MD9zej01MCIsImV4cCI6MTk0NjM4NDg5NX0.L5odyYSB76TzaDUl99UMYqkATJqOrLfQkRmmb5BNCgk')

# IDs to query (using the ones found in your workspace context)
neuron_ids = [
    533967,
    10619,
    14237,
    440899,
    11393,
    11269
]

thresholds = range(2, 11) # 2 to 10
results = []

print(f"Querying downstream for {len(neuron_ids)} neurons...")

for threshold in thresholds:
    print(f"Processing threshold: {threshold}")
    
    # Fetch adjacencies where these neurons are the source (upstream)
    # We want to find their downstream partners
    # min_total_weight applies the synapse threshold
    neurons_df, conn_df = fetch_adjacencies(
        sources=neuron_ids,
        targets=None, # Find all downstream
        min_total_weight=threshold
    )
    
    if not conn_df.empty:
        # Count unique downstream bodies (bodyId_post)
        num_downstream = conn_df['bodyId_post'].nunique()
    else:
        num_downstream = 0
        
    results.append({
        'threshold': threshold,
        'downstream_count': num_downstream
    })

# Create DataFrame and export
results_df = pd.DataFrame(results)
output_file = 'downstream_counts_by_threshold.csv'
results_df.to_csv(output_file, index=False)

print(f"Results exported to {output_file}")
print(results_df)
