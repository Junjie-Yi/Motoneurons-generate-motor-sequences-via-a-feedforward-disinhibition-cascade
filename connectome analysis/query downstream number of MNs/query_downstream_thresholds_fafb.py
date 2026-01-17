import pandas as pd
import numpy as np
from fafbseg import flywire

np1363 = [720575940630233404,720575940618539810,720575940645527918,720575940615197747,720575940630868793,720575940618165019]
min_cleftScore = 50

# thresholds to test
thresholds = [2,3,4,5,6,7,8,9,10]
# store the length of upstream for each threshold
upstream_lengths = []

motor_conn = flywire.get_connectivity(np1363, filtered=False, min_score=min_cleftScore, upstream=False)
for min_synapseNum in thresholds:
    motor_downstream = motor_conn[motor_conn['weight'] >= min_synapseNum]
    downstream = list(set(list(motor_downstream['post'])))
    upstream_lengths.append(len(downstream))

# Save results to a CSV and print the array for quick view
results_df = pd.DataFrame({
    'min_synapseNum': thresholds,
    'upstream_length': upstream_lengths
})
results_df.to_csv('upstream_lengths_cleftScore50.csv', index=False)
print(upstream_lengths)
