import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

# Use LaTeX fonts for publication-quality appearance
matplotlib.rcParams.update({
    'font.family': 'serif',
    'text.usetex': False,  # Set to True if you have LaTeX installed
    'font.size': 11,
    'axes.labelsize': 11,
    'axes.titlesize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
})

# Load data
df = pd.read_csv("../result/libox/out_libox_all_0516.csv")  # Adjust if your file has a different name
filtered_df = df.copy()

# Map index type for cleaner labels
index_type_map = {
    'libox': 'LiBox',  # Your algorithm
    'alexol': 'ALEX+',
    'lippol': 'LIPP+',
    'xindex': 'XIndex',
    'artolc': 'ART',
    'btreeolc': 'B+tree',
}
filtered_df['index_type'] = filtered_df['index_type'].map(index_type_map).fillna(filtered_df['index_type'])

# Workload labels
def get_workload_type(row):
    if row['read_ratio'] == 1.0:
        return 'RO'
    elif row['read_ratio'] == 0.5 and row['insert_ratio'] == 0.5:
        return 'BAL'
    elif row['insert_ratio'] == 1.0:
        return 'WO'
    else:
        return 'other'

# Trace label mapping with alphabetical labels
trace_map = {
    'covid': '(b) Covid', 
    'fb': '(c) FB', 
    'genome': '(d) Genome', 
    'osm': '(e) OSM', 
    'longitudes-200M': '(a) Longitudes'
}


filtered_df['trace'] = filtered_df['key_path'].apply(lambda x: next((t for t in trace_map if t in x), 'unknown'))
filtered_df['workload'] = filtered_df.apply(get_workload_type, axis=1)

# Publication-quality color scheme for scientific papers
# Making LiBox (your algorithm) stand out with a distinctive color
color_mapping = {
    'LiBox': '#e41a1c',    # Bright red for your algorithm to stand out
    'ALEX+': '#377eb8',    # Blue
    'LIPP+': '#4daf4a',    # Green
    'XIndex': '#984ea3',   # Purple
    'ART': '#ff7f00',  # Orange
    'B+tree': '#000000',   # Black
}

# Get actual indices in the dataset
actual_indices = filtered_df['index_type'].unique()
# Make sure LiBox is first, then the others in a logical order
desired_order = ['LiBox', 'ALEX+', 'LIPP+', 'ART', 'XIndex', 'B+tree']
ordered_index_types = [idx for idx in desired_order if idx in actual_indices]

# Distinctive markers
markers = ['o', 's', '^', 'D', 'v', '*']

# Define traces with alphabetical order
traces = ['longitudes-200M', 'covid', 'fb', 'genome', 'osm']
# Filter to only traces that are in your data
traces = [t for t in traces if t in trace_map]
trace_titles = [trace_map[t] for t in traces]
workloads = ['RO', 'BAL', 'WO']

# Improved figure size and layout for publication
fig, axes = plt.subplots(nrows=3, ncols=len(traces), figsize=(16, 8), sharex=True)

# Collect legend handles
legend_handles = {}

for i, workload in enumerate(workloads):
    for j, trace in enumerate(traces):
        ax = axes[i, j] if len(traces) > 1 else axes[i]
        subset = filtered_df[(filtered_df['workload'] == workload) & (filtered_df['trace'] == trace)]
        
        for idx, index_name in enumerate(ordered_index_types):
            data = subset[subset['index_type'] == index_name]
            if not data.empty:
                color = color_mapping.get(index_name, 'gray')
                marker = markers[idx % len(markers)]
                # Store the line object
                line, = ax.plot(data['thread_num'], data['throughput'] / 1e6,
                          color=color, marker=marker, 
                          markersize=6, linewidth=1.8)
                
                # Store handle for legend
                legend_handles[index_name] = line

        # Add workload label to leftmost plots
        if j == 0:
            ax.set_ylabel(f"{workload}\nThroughput (Mop/s)", fontsize=11)
        
        # Add trace name to bottom row
        if i == 2:
            ax.set_xlabel(f"Number of threads", fontsize=11)
            
            # Set xticks for your 12 threads
            ax.set_xticks([1, 2, 4, 6, 8, 10, 12])
            ax.set_xlim(0, 13)  # Give a bit of margin
            
            # Add trace titles to the bottom of each column
            ax.annotate(trace_titles[j], xy=(0.5, -0.35), xycoords='axes fraction', 
                      ha='center', va='center', fontsize=11)
        
        # Clean, scientific paper style
        ax.grid(True, linestyle='--', alpha=0.4, color='#CCCCCC')
        ax.set_ylim(bottom=0)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

# Create legend
if legend_handles:
    final_handles = []
    final_labels = []
    for index_name in ordered_index_types:
        if index_name in legend_handles:
            final_handles.append(legend_handles[index_name])
            final_labels.append(index_name)
    
    if final_handles:
        fig.legend(final_handles, final_labels, 
                  loc='upper center', 
                  bbox_to_anchor=(0.5, 1.0), 
                  ncol=len(final_handles),
                  fontsize=11, 
                  frameon=False,  # No frame for cleaner look
                  handletextpad=0.5, 
                  columnspacing=1.0, 
                  handlelength=1.5)

# Layout adjustments
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.subplots_adjust(bottom=0.15, wspace=0.3, hspace=0.3)

# Save the figure with publication-quality settings
plt.savefig("gre_libox.pdf", bbox_inches='tight', dpi=300)
plt.savefig("gre_libox.png", bbox_inches='tight', dpi=300)
plt.close()