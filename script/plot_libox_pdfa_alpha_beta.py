import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import os
from datetime import datetime
import subprocess
import tempfile

class ThroughputPlotter:
    """
    A comprehensive plotting class for index structure performance analysis with PDF/A output.
    """
    
    def __init__(self, dpi=600, pdfa_compliance=True):
        self.dpi = dpi
        self.pdfa_compliance = pdfa_compliance
        self.setup_matplotlib()
        
        # Index type mapping
        self.index_type_map = {
            'libox': 'LiBox',
            'alexol': 'ALEX+',
            'lippol': 'LIPP+',
            'xindex': 'XIndex',
            'artolc': 'ART+',
            'btreeolc': 'B+tree',
            'loft': 'LOFT',
            'alex': 'ALEX+',
            'lipp': 'LIPP+',
            'artunsync': 'ART+'
        }
        
        # Color scheme
        self.color_mapping = {
            'LiBox': '#e41a1c',    # Bright red
            'ALEX+': '#377eb8',    # Blue
            'LIPP+': '#4daf4a',    # Green
            'XIndex': '#984ea3',   # Purple
            'ART+': '#ff7f00',      # Orange
            'B+tree': '#333333',   # Dark gray
            'LOFT': '#8c564b',     # Rich burgundy
            'ALEX': '#377eb8',     # Different blue #1f77b4
            'LIPP': '#4daf4a',     # Different green #2ca02c
            'ART': '#ff7f00',     # Different red: #d62728
        }
        
        # Marker mapping
        self.marker_mapping = {
            'LiBox': 'o',      # Circle
            'ALEX+': 's',      # Square
            'LIPP+': '^',      # Triangle up
            'XIndex': 'D',     # Diamond
            'ART+': 'v',        # Triangle down
            'B+tree': '*',     # Star
            'LOFT': 'X',       # X
            'ALEX': 's',       # Pentagon  'p'
            'LIPP': '^',       # Hexagon 'h'
            'ART': 'v',       # Triangle left '<'
        }
        
        # Trace mapping
        self.trace_map = {
            'msr_web': '(b) Msr_web', 
            'w048': '(a) W048', 
            'longitudes-200M': '(c) Longitudes',
            'longitudes': '(c1) Longitudes',
            'fb': '(d) Facebook', 
            'genome': '(e) Genome', 
            'osm': '(f) OSM',
            'planet': 'Planet',
            'fb_100' : 'FB_100',
            'fb_80' : 'FB_80',
            'fb_60' : 'FB_60',
            'fb_40' : 'FB_40',
            'fb_20' : 'FB_20',
        }
        
    def setup_matplotlib(self):
        """Configure matplotlib for publication-quality PDF/A figures."""
        matplotlib.rcParams.update({
            'font.family': 'serif',
            'text.usetex': False,           
            'font.size': 13,
            'axes.labelsize': 13,
            'axes.titlesize': 13,
            'xtick.labelsize': 9,
            'ytick.labelsize': 9,
            'legend.fontsize': 13,
            'axes.linewidth': 0.8,
            'lines.linewidth': 1.8,
            'lines.markersize': 6,
            'lines.markeredgewidth': 0.8,
            'figure.dpi': self.dpi,
            'savefig.dpi': self.dpi,
            'savefig.format': 'pdf',
            'savefig.bbox': 'tight',
            # PDF/A specific settings
            'pdf.fonttype': 42,  # Embed TrueType fonts (required for PDF/A)
            'ps.fonttype': 42,   # Also for EPS output
            'pdf.use14corefonts': False,  # Don't use Type 1 fonts
            'font.serif': ['DejaVu Serif', 'Times New Roman', 'serif'],  # Use embeddable fonts
        })
    
    def load_and_preprocess_data(self, input_file):
        """Load and preprocess the CSV data."""
        print(f"Loading data from: {input_file}")
        df = pd.read_csv(input_file)
        
        # Map index types
        df['index_type'] = df['index_type'].map(self.index_type_map).fillna(df['index_type'])
        
        # Add trace information
        df['trace'] = df['key_path'].apply(lambda x: next((t for t in self.trace_map if t in x), 'unknown'))
        # print(f"test 0: {df['trace']}")
        
        # Add workload type
        df['workload'] = df.apply(self.get_workload_type, axis=1)
        
        return df
    
    def get_workload_type(self, row):
        """Determine workload type from read/insert ratios."""
        read_ratio = row.get('read_ratio', 0)
        insert_ratio = row.get('insert_ratio', 0)
        scan_ratio = row.get('scan_ratio', 0)
        
        if read_ratio == 1.0 and insert_ratio == 0:
            return 'RO'
        elif read_ratio == 0.8 and insert_ratio == 0.2:
            return 'R80'
        elif read_ratio == 0.6 and insert_ratio == 0.4:
            return 'R60'
        elif read_ratio == 0.5 and insert_ratio == 0.5:
            return 'BAL'
        elif read_ratio == 0.4 and insert_ratio == 0.6:
            return 'R40'
        elif read_ratio == 0.2 and insert_ratio == 0.8:
            return 'R20'
        elif insert_ratio == 1.0:
            return 'WO'
        elif scan_ratio == 1.0:
            return 'SO' # scanonly
        elif scan_ratio == 0.95 and insert_ratio == 0.05:
            return 'YCSBE' # YCSB-E
        else:
            return f"R:{read_ratio:.1f}/W:{insert_ratio:.1f}"
    
    def get_ordered_index_types(self, df):
        """Get ordered list of index types present in data."""
        actual_indices = df['index_type'].unique()
        print(f"actual indices: {actual_indices}")
        # desired_order = ['LiBox', 'ALEX+', 'ART+', 'LIPP+'] # for memory vs traces
        desired_order = ['LiBox', 'ALEX+', 'LIPP+', 'ART+', 'XIndex', 'B+tree', 'LOFT'] # for basic type and read ratio
        # desired_order = ['LiBox', 'ALEX+', 'XIndex', 'B+tree', 'LOFT'] # for range search type
        return [idx for idx in desired_order if idx in actual_indices]
    
    def get_traces(self, df):
        """Get ordered list of traces present in data."""
        traces = ['w048', 'msr_web', 'longitudes-200M', 'fb', 'genome', 'osm'] # basic and range test
        #memory_test_traces = ["fb_100","fb_60", "genome_20","genome_80","longitudes-200M_40","msr_web_100","msr_web_60","osm_20","osm_80","w048_20","w048_80","fb_20", "fb_80", "genome_40","longitudes-200M_100","longitudes-200M_60","msr_web_20", "msr_web_80","osm_40","w048_40","fb_40","genome_100","genome_60","longitudes-200M_20", "longitudes-200M_80","msr_web_40", "osm_100" "osm_60","w048_100","w048_60"]
        #traces = ['w048_10_10', 'w048_20_10', 'w048_30_10', 'w048_40_10', 'w048_50_10', 'w048_50_20', 'w048_50_30', 'w048_50_40', 'w048_50_50', 'w048_60_10', 'w048_70_10', 'w048_80_10']
        # traces = ['fb_10_10', 'fb_20_10', 'fb_30_10', 'fb_40_10', 'fb_50_10', 'fb_50_20', 'fb_50_30', 'fb_50_40', 'fb_50_50', 'fb_60_10', 'fb_70_10', 'fb_80_10']
        # traces = ['genome_10_10', 'genome_20_10', 'genome_30_10', 'genome_40_10', 'genome_50_10', 'genome_50_20', 'genome_50_30', 'genome_50_40', 'genome_50_50', 'genome_60_10', 'genome_70_10', 'genome_80_10']
        # traces = ['osm_10_10', 'osm_20_10', 'osm_30_10', 'osm_40_10', 'osm_50_10', 'osm_50_20', 'osm_50_30', 'osm_50_40', 'osm_50_50', 'osm_60_10', 'osm_70_10', 'osm_80_10']
        return [t for t in traces if t in df['trace'].values]
    
    def style_axis(self, ax):
        """Apply consistent styling to axis."""
        ax.grid(True, linestyle='--', alpha=0.3, color='#CCCCCC')
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        for spine in ['left', 'bottom']:
            ax.spines[spine].set_linewidth(0.8)
            ax.spines[spine].set_color('#444444')
    
    def add_legend(self, fig, legend_handles, ordered_index_types):
        """Add legend to figure."""
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
                          frameon=True,
                          fancybox=False,
                          edgecolor='#DDDDDD',
                          handletextpad=0.5, 
                          columnspacing=1.2, 
                          handlelength=1.5)

    def plot_basic_throughput(self, df, output_file):
        """Generate basic throughput plots (original function)."""
        traces = self.get_traces(df)
        trace_titles = [self.trace_map.get(t, t) for t in traces]
        workloads = ['RO', 'BAL', 'WO']
        ordered_index_types = self.get_ordered_index_types(df)
        
        fig, axes = plt.subplots(nrows=3, ncols=len(traces), figsize=(16, 8), sharex=True)
        if len(traces) == 1:
            axes = axes.reshape(len(workloads), 1)
        
        legend_handles = {}
        
        for i, workload in enumerate(workloads):
            for j, trace in enumerate(traces):
                ax = axes[i, j]
                subset = df[(df['workload'] == workload) & (df['trace'] == trace)]
                
                for index_name in ordered_index_types:
                    data = subset[subset['index_type'] == index_name]
                    if not data.empty:
                        color = self.color_mapping.get(index_name, 'gray')
                        marker = self.marker_mapping.get(index_name, 'o')
                        
                        line, = ax.plot(data['thread_num'], data['throughput'] / 1e6,
                                      color=color, marker=marker, linestyle='-',
                                      markersize=6, linewidth=1.8, markeredgecolor='black',
                                      markeredgewidth=0.1, label=index_name)
                        
                        legend_handles[index_name] = line
                
                if j == 0:
                    ax.set_ylabel(f"{workload}\nThroughput (Mop/s)", fontsize=11)
                
                if i == 2:
                    ax.set_xlabel("Number of threads", fontsize=11)
                    ax.set_xticks([1, 8, 16, 24, 32, 40, 48, 56, 64, 72, 84])
                    ax.set_xlim(0, 84)
                    ax.annotate(trace_titles[j], xy=(0.5, -0.35), xycoords='axes fraction', 
                              ha='center', va='center', fontsize=11)
                
                self.style_axis(ax)
                ax.set_ylim(bottom=0)
        
        self.add_legend(fig, legend_handles, ordered_index_types)
        plt.tight_layout(rect=[0, 0.02, 1, 0.92])
        plt.subplots_adjust(bottom=0.15, wspace=0.25, hspace=0.35)
        
        self.save_figure(fig, output_file)

    def plot_read_ratio_throughput(self, df, output_file):
        """Generate read ratio vs throughput plots."""
        traces = self.get_traces(df)
        trace_titles = [self.trace_map.get(t, t) for t in traces]
        ordered_index_types = self.get_ordered_index_types(df)
        
        # Read ratio mapping
        read_ratios = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
        read_ratio_labels = ['0%', '20%', '40%', '60%', '80%', '100%']
        workload_map = {0: 'WO', 0.2: 'R20', 0.4: 'R40', 0.6: 'R60', 0.8: 'R80', 1.0: 'RO'}
        
        fig, axes = plt.subplots(nrows=1, ncols=6, figsize=(16, 3.5), sharex=True)
        axes = axes.flatten()
        
        legend_handles = {}
        
        for j, trace in enumerate(traces):
            if j >= len(axes):
                break
                
            ax = axes[j]
            # print(f"workloads: {df['workload'].values}")
            
            for index_name in ordered_index_types:
                # print(f"index_name: {index_name}")
                throughputs = []
                for read_ratio in read_ratios:
                    workload = workload_map[read_ratio]
                    data = df[(df['trace'] == trace) & 
                             (df['workload'] == workload) & 
                             (df['index_type'] == index_name)]
                    # print(f"read ratio: {read_ratio} workload: {workload}")
                    
                    if not data.empty:
                        # print(f"workload: {df['workload'].values}; current workload: {workload}; read_ratio: {read_ratio}")
                        # Use thread_num=40 or average if multiple threads exist
                        if 40 in data['thread_num'].values:
                            throughput = data[data['thread_num'] == 40]['throughput'].iloc[0] / 1e6
                        else:
                            throughput = data['throughput'].mean() / 1e6
                        throughputs.append(throughput)
                    else:
                        throughputs.append(None)
                
                # Filter out None values for plotting
                valid_indices = [i for i, t in enumerate(throughputs) if t is not None]
                valid_ratios = [read_ratios[i] * 100 for i in valid_indices]
                valid_throughputs = [throughputs[i] for i in valid_indices]
                
                if valid_throughputs:
                    color = self.color_mapping.get(index_name, 'gray')
                    marker = self.marker_mapping.get(index_name, 'o')
                    
                    line, = ax.plot(valid_ratios, valid_throughputs,
                                  color=color, marker=marker, linestyle='-',
                                  markersize=6, linewidth=1.8, markeredgecolor='black',
                                  markeredgewidth=0.1, label=index_name)
                    
                    legend_handles[index_name] = line
            
            ax.set_xlabel("Read Ratio (%)", fontsize=11)
            if j % 3 == 0:  # Left column
                ax.set_ylabel("Throughput (Mop/s)", fontsize=11)
            # Add trace titles to the bottom of each subplot
            ax.annotate(trace_titles[j], xy=(0.5, -0.25), xycoords='axes fraction', 
                      ha='center', va='center', fontsize=11)
            
            ax.set_xticks([0, 20, 40, 60, 80, 100])
            ax.set_xlim(-5, 105)
            self.style_axis(ax)
            ax.set_ylim(bottom=0)
        
        # Hide unused subplots
        for k in range(len(traces), len(axes)):
            axes[k].set_visible(False)
        
        self.add_legend(fig, legend_handles, ordered_index_types)
        plt.tight_layout(rect=[0, 0.02, 1, 0.92])
        plt.subplots_adjust(bottom=0.15, wspace=0.25, hspace=0.35)
        
        self.save_figure(fig, output_file)

    def plot_range_search_throughput(self, df, output_file):
        """Generate range search throughput plots for thread_num=40."""
        traces = self.get_traces(df)
        trace_titles = [self.trace_map.get(t, t) for t in traces]
        workloads = ['YCSBE'] # 'SO' for scanonly, 'YCSBE' for YCSB-E
        # workloads = ['SO', 'YCSBE'] # 'SO' for scanonly, 'YCSBE' for YCSB-E
        ordered_index_types = self.get_ordered_index_types(df)
        print(f"index: {ordered_index_types}")
        
        fig, axes = plt.subplots(nrows=1, ncols=len(traces), figsize=(16, 3.5), sharex=True)
        if len(traces) == 1:
            axes = axes.reshape(len(workloads), 1)
        
        legend_handles = {}
        
        for i, workload in enumerate(workloads):
            for j, trace in enumerate(traces):
                # ax = axes[i, j] # for scanonly and ycsb-e
                ax = axes[j] # for scanonly
                subset = df[(df['workload'] == workload) & (df['trace'] == trace)]
                
                for index_name in ordered_index_types:
                    data = subset[subset['index_type'] == index_name]
                    if not data.empty:
                        color = self.color_mapping.get(index_name, 'gray')
                        marker = self.marker_mapping.get(index_name, 'o')
                        
                        line, = ax.plot(data['thread_num'], data['throughput'] / 1e6,
                                      color=color, marker=marker, linestyle='-',
                                      markersize=6, linewidth=1.8, markeredgecolor='black',
                                      markeredgewidth=0.1, label=index_name)
                        
                        legend_handles[index_name] = line
                
                if j == 0:
                    ax.set_ylabel(f"{workload}\nThroughput (Mop/s)", fontsize=11)
                
                ax.set_xlabel("Number of threads", fontsize=11)
                ax.set_xticks([1, 8, 16, 24, 32, 40, 48, 56, 64, 72, 84])
                ax.set_xlim(0, 84)
                ax.annotate(trace_titles[j], xy=(0.5, -0.25), xycoords='axes fraction', 
                            ha='center', va='center', fontsize=11)
                
                self.style_axis(ax)
                ax.set_ylim(bottom=0)
        
        self.add_legend(fig, legend_handles, ordered_index_types)
        plt.tight_layout(rect=[0, 0.02, 1, 0.92])
        plt.subplots_adjust(bottom=0.15, wspace=0.25, hspace=0.35)
        
        self.save_figure(fig, output_file)

    def plot_memory_vs_traces(self, df, output_file):
        """Generate memory consumption vs percentages of inserted keys plots."""
        
        # Extract percentage from key_path
        import re
        
        def extract_trace_info(key_path):
            """Extract trace name and percentage from path like /path/to/fb_60.csv -> ('fb', 60)"""
            match = re.search(r'/([^/]+)_(\d+)\.csv$', key_path)
            if match:
                trace_name = match.group(1)
                percentage = int(match.group(2))
                return trace_name, percentage
            return None, None
        
        # Create a copy of the dataframe to avoid modifying the original
        df_processed = df.copy()
        
        # Extract trace info
        trace_info = df_processed['key_path'].apply(extract_trace_info)
        df_processed['trace_name'] = trace_info.apply(lambda x: x[0] if x[0] is not None else 'unknown')
        df_processed['percentage'] = trace_info.apply(lambda x: x[1] if x[1] is not None else 0)
        
        # Filter out rows where trace extraction failed
        df_processed = df_processed[df_processed['trace_name'] != 'unknown']
        
        # Get unique traces from the processed data
        available_traces = sorted(df_processed['trace_name'].unique())
        print(f"Available traces: {available_traces}")
        
        # Define the 6 main traces we want to plot
        target_traces = ['w048', 'msr_web', 'longitudes-200M', 'fb', 'genome', 'osm']
        
        # Filter to only include traces that exist in our data
        traces_to_plot = [trace for trace in target_traces if trace in available_traces]
        print(f"Traces to plot: {traces_to_plot}")
        
        # Trace titles mapping
        trace_titles = {
            'w048': '(a) W048', 
            'msr_web': '(b) Msr_web', 
            'longitudes-200M': '(c) Longitudes',
            'fb': '(d) Facebook', 
            'genome': '(e) Genome', 
            'osm': '(f) OSM'
        }
        
        ordered_index_types = self.get_ordered_index_types(df_processed)
        print(f"Index types: {ordered_index_types}")
        
        percentages = [20, 40, 60, 80, 100]
        
        # Create subplots
        fig, axes = plt.subplots(nrows=1, ncols=6, figsize=(16, 3.5), sharex=True)
        axes = axes.flatten()
        
        legend_handles = {}
        
        for j, trace in enumerate(traces_to_plot):
            if j >= len(axes):
                break
                
            ax = axes[j]
            
            for index_name in ordered_index_types:
                memories = []
                valid_percentages = []
                
                for percentage in percentages:
                    # Filter data based on trace name, percentage, and index_type
                    data = df_processed[
                        (df_processed['trace_name'] == trace) & 
                        (df_processed['percentage'] == percentage) & 
                        (df_processed['index_type'] == index_name)
                    ]
                    
                    if not data.empty:
                        memory = data['memory_consumption'].iloc[0] / 1e9  # Convert to GB
                        memories.append(memory)
                        valid_percentages.append(percentage)
                        # print(f"Trace: {trace}, Percentage: {percentage}, Index: {index_name}, Memory: {memory:.2f} GB")
                
                if memories and valid_percentages:
                    color = self.color_mapping.get(index_name, 'gray')
                    marker = self.marker_mapping.get(index_name, 'o')
                    
                    line, = ax.plot(valid_percentages, memories,
                                color=color, marker=marker, linestyle='-',
                                markersize=6, linewidth=1.8, markeredgecolor='black',
                                markeredgewidth=0.1, label=index_name)
                    
                    legend_handles[index_name] = line
            
            # Set labels and formatting
            ax.set_xlabel("Inserted Keys (%)", fontsize=11)
            if j % 3 == 0:  # Left column
                ax.set_ylabel("Memory Consumption (GB)", fontsize=11)
                
            # Add trace title
            trace_title = trace_titles.get(trace, trace)
            # ax.set_title(trace_title, fontsize=12, pad=10)
            ax.annotate(trace_title, xy=(0.5, -0.25), xycoords='axes fraction', 
                      ha='center', va='center', fontsize=11)
            
            ax.set_xticks([20, 40, 60, 80, 100])
            ax.set_xlim(15, 105)
            self.style_axis(ax)
            ax.set_ylim(bottom=0)
        
        # Hide unused subplots
        for k in range(len(traces_to_plot), len(axes)):
            axes[k].set_visible(False)
        
        # ordered_index_types = ['LiBox', 'ALEX+', 'ART+', 'LIPP+'] 
        print(f"ordered_index_types: {ordered_index_types}")
        # Add legend at the top
        self.add_legend(fig, legend_handles, ordered_index_types)
        
        plt.tight_layout(rect=[0, 0.02, 1, 0.92])
        plt.subplots_adjust(bottom=0.1, wspace=0.25, hspace=0.35)
        
        self.save_figure(fig, output_file)

    def plot_alpha_beta_analysis(self, df, output_file, trace_filter='genome'):
        """Generate Alpha and Beta analysis plots for LiBox only."""
        # Filter for specified trace pattern, thread_num=40, and LiBox only
        trace_pattern = trace_filter  # e.g., 'osm' will match 'osm_10_10', 'osm_20_10', etc.
        # print(f"test1: {df['trace']}")
        # print(f"test2: {df[df['trace'].str.contains(trace_pattern, na=False)]['trace']}")

        df_filtered = df[(df['trace'].str.contains(trace_pattern, na=False)) & 
                        (df['thread_num'] == 40) & 
                        (df['index_type'] == 'LiBox')]
        
        print(f"Filtering for traces containing '{trace_pattern}'")
        print(f"Found traces: {sorted(df_filtered['trace'].unique())}")
        
        # Extract alpha and beta from key_path
        def extract_alpha_beta(key_path):
            """Extract alpha and beta values from key_path like w048_30_10.csv"""
            try:
                # Get filename without extension
                filename = key_path.split('/')[-1].replace('.csv', '')
                # print(f"filename: {filename}")
                # Split by underscore and get the last two parts (beta_alpha)
                parts = filename.split('_')
                if len(parts) >= 3:
                    beta = int(parts[-2])
                    alpha = int(parts[-1])
                    return alpha, beta
            except:
                pass
            return None, None
        
        # Add alpha and beta columns
        df_filtered = df_filtered.copy()
        # print(f"key_path: {df_filtered['key_path'].values}")
        df_filtered[['alpha', 'beta']] = df_filtered['key_path'].apply(
            lambda x: pd.Series(extract_alpha_beta(x))
        )
        
        # Remove rows where alpha or beta couldn't be extracted
        df_filtered = df_filtered.dropna(subset=['alpha', 'beta'])
        
        # Debug: Print available alpha/beta values
        print("Available alpha values:", sorted(df_filtered['alpha'].unique()))
        print("Available beta values:", sorted(df_filtered['beta'].unique()))
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4)) # (9, 4))
        
        # Alpha analysis (subplot 1) - vary alpha, keep beta constant
        # Find the most common beta value or use a specific one
        if not df_filtered.empty:
            common_beta = df_filtered['beta'].mode().iloc[0] if not df_filtered['beta'].mode().empty else 30
            alpha_data = df_filtered[df_filtered['beta'] == common_beta]
            
            if not alpha_data.empty:
                alpha_values = sorted(alpha_data['alpha'].unique())
                
                ro_throughputs = []
                wo_throughputs = []
                memories = []
                valid_alphas = []
                
                for alpha in alpha_values:
                    alpha_subset = alpha_data[alpha_data['alpha'] == alpha]
                    
                    ro_data = alpha_subset[alpha_subset['workload'] == 'RO']
                    wo_data = alpha_subset[alpha_subset['workload'] == 'WO']
                    
                    if not ro_data.empty:
                        ro_throughputs.append(ro_data['throughput'].iloc[0] / 1e6)
                        memories.append(ro_data['memory_consumption'].iloc[0] / 1e9)
                        valid_alphas.append(alpha)
                    elif not wo_data.empty:
                        # If no RO data, still add the alpha but with None for RO
                        ro_throughputs.append(None)
                        memories.append(wo_data['memory_consumption'].iloc[0] / 1e9)
                        valid_alphas.append(alpha)
                        
                    if not wo_data.empty:
                        if len(wo_throughputs) < len(valid_alphas):
                            wo_throughputs.append(wo_data['throughput'].iloc[0] / 1e6)
                        else:
                            wo_throughputs[-1] = wo_data['throughput'].iloc[0] / 1e6
                    else:
                        if len(wo_throughputs) < len(valid_alphas):
                            wo_throughputs.append(None)
                
                # Plot throughput on left y-axis, memory on right y-axis
                ax1_twin = ax1.twinx()
                
                if valid_alphas:
                    # Plot RO and WO throughput
                    valid_ro = [(alpha, tput) for alpha, tput in zip(valid_alphas, ro_throughputs) if tput is not None]
                    valid_wo = [(alpha, tput) for alpha, tput in zip(valid_alphas, wo_throughputs) if tput is not None]
                    valid_mem = [(alpha, mem) for alpha, mem in zip(valid_alphas, memories) if mem is not None]
                    
                    if valid_ro:
                        alphas_ro, tputs_ro = zip(*valid_ro)
                        # ax1.plot(alphas_ro, tputs_ro, 'o-', color='blue', label='RO Throughput', linewidth=2, markersize=8)
                        ax1.plot(alphas_ro, tputs_ro, 'o-', color='#e41a1c', label='RO Throughput', 
                                linewidth=1.8, markersize=6, markeredgecolor='black', markeredgewidth=0.1)
                    
                    if valid_wo:
                        alphas_wo, tputs_wo = zip(*valid_wo)
                        # ax1.plot(alphas_wo, tputs_wo, 's-', color='red', label='WO Throughput', linewidth=2, markersize=8)
                        ax1.plot(alphas_wo, tputs_wo, 'o-', color='#377eb8', label='WO Throughput', 
                                linewidth=1.8, markersize=6, markeredgecolor='black', markeredgewidth=0.1)
                    
                    if valid_mem:
                        alphas_mem, mems = zip(*valid_mem)
                        # ax1_twin.plot(alphas_mem, mems, '^-', color='green', label='Index size', linewidth=2, markersize=8)
                        ax1_twin.plot(alphas_mem, mems, '^--', color='#000000', label='Index Size', 
                                     linewidth=1.8, markersize=6, markeredgecolor='black', markeredgewidth=0.1)
                
                ax1.set_xlabel(r'Overflow ratio threshold $\alpha$ (%)', fontsize=11)
                ax1.set_ylabel('Throughput (Mop/s)', color='black', fontsize=11)
                ax1_twin.set_ylabel('Index size (GB)', color='black', fontsize=11)
                ax1.set_ylim(bottom=0, top=240)
                ax1_twin.set_ylim(bottom=0, top=7)
                ax1.set_xticks([5, 10, 20, 30, 40, 50])
                
                # ax1.annotate('Overflow ratio threshold 𝛼(%) vs Performance & Memory', xy=(0.5, -0.2), xycoords='axes fraction', 
                #            ha='center', va='center', fontsize=11)
                
                # Add legends
                # lines1, labels1 = ax1.get_legend_handles_labels()
                # lines2, labels2 = ax1_twin.get_legend_handles_labels()
                # ax1.legend(lines1 + lines2, labels1 + labels2, loc='best')
        
        # Beta analysis (subplot 2) - vary beta, keep alpha constant
        if not df_filtered.empty:
            common_alpha = df_filtered['alpha'].mode().iloc[0] if not df_filtered['alpha'].mode().empty else 10
            beta_data = df_filtered[df_filtered['alpha'] == common_alpha]
            
            if not beta_data.empty:
                beta_values = sorted(beta_data['beta'].unique())
                
                ro_throughputs = []
                wo_throughputs = []
                memories = []
                valid_betas = []
                
                for beta in beta_values:
                    beta_subset = beta_data[beta_data['beta'] == beta]
                    
                    ro_data = beta_subset[beta_subset['workload'] == 'RO']
                    wo_data = beta_subset[beta_subset['workload'] == 'WO']
                    
                    if not ro_data.empty:
                        ro_throughputs.append(ro_data['throughput'].iloc[0] / 1e6)
                        memories.append(ro_data['memory_consumption'].iloc[0] / 1e9)
                        valid_betas.append(beta)
                    elif not wo_data.empty:
                        ro_throughputs.append(None)
                        memories.append(wo_data['memory_consumption'].iloc[0] / 1e9)
                        valid_betas.append(beta)
                        
                    if not wo_data.empty:
                        if len(wo_throughputs) < len(valid_betas):
                            wo_throughputs.append(wo_data['throughput'].iloc[0] / 1e6)
                        else:
                            wo_throughputs[-1] = wo_data['throughput'].iloc[0] / 1e6
                    else:
                        if len(wo_throughputs) < len(valid_betas):
                            wo_throughputs.append(None)
                
                # Plot throughput on left y-axis, memory on right y-axis
                ax2_twin = ax2.twinx()
                
                if valid_betas:
                    # Plot RO and WO throughput
                    valid_ro = [(beta, tput) for beta, tput in zip(valid_betas, ro_throughputs) if tput is not None]
                    valid_wo = [(beta, tput) for beta, tput in zip(valid_betas, wo_throughputs) if tput is not None]
                    valid_mem = [(beta, mem) for beta, mem in zip(valid_betas, memories) if mem is not None]
                    
                    if valid_ro:
                        betas_ro, tputs_ro = zip(*valid_ro)
                        # ax2.plot(betas_ro, tputs_ro, 'o-', color='blue', label='RO Throughput', linewidth=2, markersize=8)
                        ax2.plot(betas_ro, tputs_ro, 'o-', color='#e41a1c', label='RO Throughput', 
                                linewidth=1.8, markersize=6, markeredgecolor='black', markeredgewidth=0.1)
                    
                    if valid_wo:
                        betas_wo, tputs_wo = zip(*valid_wo)
                        # ax2.plot(betas_wo, tputs_wo, 's-', color='red', label='WO Throughput', linewidth=2, markersize=8)
                        ax2.plot(betas_wo, tputs_wo, 'o-', color='#377eb8', label='WO Throughput', 
                                linewidth=1.8, markersize=6, markeredgecolor='black', markeredgewidth=0.1)
                    
                    if valid_mem:
                        betas_mem, mems = zip(*valid_mem)
                        # ax2_twin.plot(betas_mem, mems, '^-', color='green', label='Index size', linewidth=2, markersize=8)
                        ax2_twin.plot(betas_mem, mems, '^--', color='#000000', label='Index Size', 
                                     linewidth=1.8, markersize=6, markeredgecolor='black', markeredgewidth=0.1)
                
                ax2.set_xlabel(r'Underflow ratio threshold $\beta$ (%)', fontsize=11)
                ax2.set_ylabel('Throughput (Mop/s)', color='black', fontsize=11)
                ax2_twin.set_ylabel('Index size (GB)', color='black', fontsize=11)
                ax2.set_ylim(bottom=0, top=240)
                ax2_twin.set_ylim(bottom=0, top=7)
                
                # Add title to the bottom of subplot
                # ax2.annotate('Beta vs Performance & Memory', xy=(0.5, -0.2), xycoords='axes fraction', 
                #            ha='center', va='center', fontsize=11)
                
                # Add legends
                # lines1, labels1 = ax2.get_legend_handles_labels()
                # lines2, labels2 = ax2_twin.get_legend_handles_labels()
                # ax2.legend(lines1 + lines2, labels1 + labels2, loc='best')
        
        self.style_axis(ax1)
        self.style_axis(ax2)
        
        # Create a single legend at the top for both subplots
        # Collect handles from both subplots
        all_handles = []
        all_labels = []
        
        # Get handles from ax1
        lines1, labels1 = ax1.get_legend_handles_labels()
        all_handles.extend(lines1)
        all_labels.extend(labels1)
        
        # Get handles from ax1_twin if it exists
        if 'ax1_twin' in locals():
            lines2, labels2 = ax1_twin.get_legend_handles_labels()
            all_handles.extend(lines2)
            all_labels.extend(labels2)
        
        # Remove duplicates while preserving order
        unique_handles = []
        unique_labels = []
        for handle, label in zip(all_handles, all_labels):
            if label not in unique_labels:
                unique_handles.append(handle)
                unique_labels.append(label)
        
        if unique_handles:
            fig.legend(unique_handles, unique_labels, 
                      loc='upper center', 
                      bbox_to_anchor=(0.5, 1.0), 
                      ncol=len(unique_handles),
                      fontsize=11, 
                      frameon=True,
                      fancybox=False,
                      edgecolor='#DDDDDD',
                      handletextpad=0.5, 
                      columnspacing=1.2, 
                      handlelength=1.5)
        
        plt.tight_layout(rect=[0, 0.02, 1, 0.92])
        plt.subplots_adjust(bottom=0.15, wspace=0.25, hspace=0.35)
        self.save_figure(fig, output_file)

    def convert_to_pdfa_with_ghostscript(self, input_pdf, output_pdf):
        """Convert PDF to PDF/A using Ghostscript."""
        try:
            cmd = [
                'gs', 
                '-dPDFA=1',  # PDF/A-1b compliance
                '-dBATCH', 
                '-dNOPAUSE',
                '-sColorConversionStrategy=UseDeviceIndependentColor',
                '-sDEVICE=pdfwrite',
                '-dPDFACompatibilityPolicy=1',  # Convert non-compliant elements
                f'-sOutputFile={output_pdf}',
                input_pdf
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            print(f"Successfully converted to PDF/A: {output_pdf}")
            return True
        except subprocess.CalledProcessError as e:
            print(f"Warning: Ghostscript conversion failed: {e}")
            print(f"Error output: {e.stderr}")
            return False
        except FileNotFoundError:
            print("Warning: Ghostscript not found. PDF/A conversion skipped.")
            print("Install Ghostscript for full PDF/A compliance.")
            return False

    def save_figure(self, fig, output_file):
        """Save figure in PDF/A format with proper metadata."""
        # Create output directory if it doesn't exist
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        print(f"Saving figure to: {output_file}")
        
        if self.pdfa_compliance:
            # Create PDF/A compliant metadata
            pdf_metadata = {
                'Title': 'Index Performance Analysis',
                'Author': 'Performance Analysis System',
                'Subject': 'Database Index Throughput and Memory Analysis',
                'Creator': 'matplotlib',
                'Producer': 'matplotlib PDF/A backend',
                'CreationDate': datetime.now(),
                'Keywords': 'database, index, performance, throughput'
            }
            
            # Use PdfPages for better PDF/A compliance
            with PdfPages(output_file, metadata=pdf_metadata) as pdf:
                pdf.savefig(fig, bbox_inches='tight', dpi=self.dpi)
                
                # Add PDF/A specific information
                d = pdf.infodict()
                d['Title'] = pdf_metadata['Title']
                d['Author'] = pdf_metadata['Author'] 
                d['Subject'] = pdf_metadata['Subject']
                d['Keywords'] = pdf_metadata['Keywords']
                d['Creator'] = pdf_metadata['Creator']
                d['Producer'] = pdf_metadata['Producer']
        else:
            # Standard PDF save
            plt.savefig(output_file, bbox_inches='tight', dpi=self.dpi)
        
        # Try to convert to PDF/A using Ghostscript for maximum compliance
        if self.pdfa_compliance:
            temp_file = output_file.replace('.pdf', '_temp.pdf')
            os.rename(output_file, temp_file)
            
            if self.convert_to_pdfa_with_ghostscript(temp_file, output_file):
                os.remove(temp_file)
                print(f"PDF/A conversion completed: {output_file}")
            else:
                # Fallback to original file if conversion fails
                os.rename(temp_file, output_file)
                print(f"Using matplotlib PDF (PDF/A features applied): {output_file}")
        
        # Also save as PNG
        png_output = output_file.replace('.pdf', '.png')
        plt.savefig(png_output, bbox_inches='tight', dpi=300, metadata=None)
        print(f"Also saved as: {png_output}")
        
        plt.close()
        print("Figure generation complete.")

    def plot_single_trace_range_search(self, df, output_file, trace_name='osm'):
        """Generate range search throughput plot for a single trace with YCSBE workload."""
        workload = 'YCSBE'
        ordered_index_types = self.get_ordered_index_types(df)
        
        # Filter data for the specified trace and workload
        subset = df[(df['trace'] == trace_name) & (df['workload'] == workload)]
        
        if subset.empty:
            print(f"No data found for trace '{trace_name}' with workload '{workload}'")
            return
        
        # Create single plot
        fig, ax = plt.subplots(figsize=(7, 4))
        
        legend_handles = {}
        
        for index_name in ordered_index_types:
            data = subset[subset['index_type'] == index_name]
            if not data.empty:
                color = self.color_mapping.get(index_name, 'gray')
                marker = self.marker_mapping.get(index_name, 'o')
                
                line, = ax.plot(data['thread_num'], data['throughput'] / 1e6,
                            color=color, marker=marker, linestyle='-',
                            markersize=6, linewidth=1.8, markeredgecolor='black',
                            markeredgewidth=0.1, label=index_name)
                
                legend_handles[index_name] = line
        
        # Set labels and formatting
        ax.set_xlabel("Number of threads", fontsize=13)
        ax.set_ylabel("Throughput (Mop/s)", fontsize=13)
        
        # Get trace title from mapping, fallback to trace_name if not found
        trace_title = self.trace_map.get(trace_name, trace_name.title())
        # ax.set_title(f"{trace_title} - {workload} Workload", fontsize=14, pad=15)
        ax.annotate(trace_title, xy=(0.5, -0.25), xycoords='axes fraction', 
                      ha='center', va='center', fontsize=11)
        
        ax.set_xticks([1, 8, 16, 24, 32, 40, 48, 56, 64, 72, 84])
        ax.set_xlim(0, 84)
        
        self.style_axis(ax)
        ax.set_ylim(bottom=0)
        # Add main title at the top of the figure
        fig.suptitle(f"YCSB-E Performance", fontsize=16, y=0.95)
        
        # Add legend
        if legend_handles:
            final_handles = []
            final_labels = []
            for index_name in ordered_index_types:
                if index_name in legend_handles:
                    final_handles.append(legend_handles[index_name])
                    final_labels.append(index_name)
            
            if final_handles:
                ax.legend(final_handles, final_labels,
                        loc='best', 
                        fontsize=11,
                        frameon=True,
                        fancybox=False,
                        edgecolor='#DDDDDD')
        
        plt.tight_layout()
        self.save_figure(fig, output_file)
        
    
    def plot_scan_length_throughput(self, df, output_file):
        """Generate scan throughput vs scan length plots for different traces."""
        traces = self.get_traces(df)
        trace_titles = [self.trace_map.get(t, t) for t in traces]
        ordered_index_types = self.get_ordered_index_types(df)
        
        # Define scan lengths
        scan_lengths = [100, 500, 1000, 1500, 2000, 2500]
        
        # Create subplots: 1 row, 6 columns
        fig, axes = plt.subplots(nrows=1, ncols=6, figsize=(16, 3.5), sharex=True, sharey=True)
        axes = axes.flatten()
        
        # Add main title at the top of the figure
        fig.suptitle("Scan Throughput vs Scan Length", fontsize=16, y=0.95)
        
        legend_handles = {}
        
        for j, trace in enumerate(traces[:6]):  # Limit to 6 traces
            if j >= len(axes):
                break
                
            ax = axes[j]
            
            for index_name in ordered_index_types:
                throughputs = []
                valid_lengths = []
                
                for scan_length in scan_lengths:
                    # Filter data for specific trace, index type, and scan length                    
                    data = df[
                        (df['trace'] == trace) & 
                        (df['index_type'] == index_name) & 
                        (df['scan_num'] == scan_length) 
                    ]
                    
                    if not data.empty:
                        # Use a specific thread count or average if multiple exist
                        if 40 in data['thread_num'].values:
                            throughput = data[data['thread_num'] == 40]['throughput'].iloc[0] / 1e6
                        else:
                            throughput = data['throughput'].mean() / 1e6
                        throughputs.append(throughput)
                        valid_lengths.append(scan_length)
                
                # Plot if we have valid data points
                if throughputs and valid_lengths:
                    color = self.color_mapping.get(index_name, 'gray')
                    marker = self.marker_mapping.get(index_name, 'o')
                    
                    line, = ax.plot(valid_lengths, throughputs,
                                color=color, marker=marker, linestyle='-',
                                markersize=6, linewidth=1.8, markeredgecolor='black',
                                markeredgewidth=0.1, label=index_name)
                    
                    legend_handles[index_name] = line
            
            # Set labels and formatting
            ax.set_xlabel("Scan Length", fontsize=11)
            if j == 0:  # Only leftmost subplot gets y-label
                ax.set_ylabel("Throughput (Mop/s)", fontsize=11)
            
            # Add trace title below each subplot
            trace_title = trace_titles[j] if j < len(trace_titles) else trace
            ax.annotate(trace_title, xy=(0.5, -0.25), xycoords='axes fraction', 
                    ha='center', va='center', fontsize=11)
            
            # Set x-axis ticks and limits
            ax.set_xticks(scan_lengths)
            ax.set_xlim(0, 2600)
            
            self.style_axis(ax)
            ax.set_ylim(bottom=0)
        
        # Hide unused subplots if fewer than 6 traces
        for k in range(len(traces), 6):
            if k < len(axes):
                axes[k].set_visible(False)
        
        # Add legend at the top
        self.add_legend(fig, legend_handles, ordered_index_types)
        
        plt.tight_layout(rect=[0, 0.02, 1, 0.88])
        plt.subplots_adjust(bottom=0.15, wspace=0.25, hspace=0.35)
        
        self.save_figure(fig, output_file)
    

def main():
    """Main function with argument parsing."""
    parser = argparse.ArgumentParser(description='Generate various index performance plots with PDF/A compliance')
    parser.add_argument('--input', type=str, help='Input CSV file path')
    parser.add_argument('--output', type=str, help='Output file path (without extension)')
    parser.add_argument('--plot_type', type=str, default='basic',
                       choices=['basic', 'read_ratio', 'range_search', 'memory', 'alpha_beta', 'ycsbe', 'scan_length'],
                       help='Type of plot to generate (default: basic)')
    parser.add_argument('--dpi', type=int, default=600, help='DPI for the output figure (default: 600)')
    parser.add_argument('--no-pdfa', action='store_true', help='Disable PDF/A compliance (default: PDF/A enabled)')
    
    args = parser.parse_args()
    
    # Generate default file paths if not provided
    if args.input is None or args.output is None:
        date_tag = datetime.now().strftime("%m%d")
        
    if args.input is None:
        args.input = f"../result/libox/out_libox_{date_tag}.csv"
    
    if args.output is None:
        args.output = f"../result/libox/graph_{args.plot_type}_{date_tag}"
    
    # Initialize plotter with PDF/A compliance setting
    pdfa_compliance = not args.no_pdfa
    plotter = ThroughputPlotter(dpi=args.dpi, pdfa_compliance=pdfa_compliance)
    
    if pdfa_compliance:
        print("PDF/A compliance enabled - output will be PDF/A compliant")
    else:
        print("PDF/A compliance disabled - standard PDF output")
    
    # Load and preprocess data
    df = plotter.load_and_preprocess_data(args.input)
    
    # Generate appropriate plot
    output_file = f"{args.output}.pdf"
    
    if args.plot_type == 'basic':
        plotter.plot_basic_throughput(df, output_file)
    elif args.plot_type == 'read_ratio':
        plotter.plot_read_ratio_throughput(df, output_file)
    elif args.plot_type == 'range_search':
        plotter.plot_range_search_throughput(df, output_file)
    elif args.plot_type == 'memory':
        plotter.plot_memory_vs_traces(df, output_file)
    elif args.plot_type == 'alpha_beta':
        plotter.plot_alpha_beta_analysis(df, output_file)
    elif args.plot_type == 'ycsbe':
        plotter.plot_single_trace_range_search(df, output_file)
    elif args.plot_type == 'scan_length':
        plotter.plot_scan_length_throughput(df, output_file)
    
    print(f"Plot generation completed for type: {args.plot_type}")

if __name__ == "__main__":
    main()