import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import argparse
import os
from datetime import datetime

class ThroughputPlotter:
    """
    A comprehensive plotting class for index structure performance analysis.
    """
    
    def __init__(self, dpi=600):
        self.dpi = dpi
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
            'alex': 'ALEX',
            'lipp': 'LIPP',
            'artunsync': 'ART'
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
            'msr_web': '(a) msr_web', 
            'w048': '(b) w048', 
            'longitudes-200M': '(c) Longitudes',
            'fb': '(d) Facebook', 
            'genome': '(e) Genome', 
            'osm': '(f) OSM'
        }
        
    def setup_matplotlib(self):
        """Configure matplotlib for publication-quality figures."""
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
        })
    
    def load_and_preprocess_data(self, input_file):
        """Load and preprocess the CSV data."""
        print(f"Loading data from: {input_file}")
        df = pd.read_csv(input_file)
        
        # Map index types
        df['index_type'] = df['index_type'].map(self.index_type_map).fillna(df['index_type'])
        
        # Add trace information
        df['trace'] = df['key_path'].apply(lambda x: next((t for t in self.trace_map if t in x), 'unknown'))
        
        # Add workload type
        df['workload'] = df.apply(self.get_workload_type, axis=1)
        
        return df
    
    def get_workload_type(self, row):
        """Determine workload type from read/insert ratios."""
        read_ratio = row.get('read_ratio', 0)
        insert_ratio = row.get('insert_ratio', 0)
        
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
        else:
            return f"R:{read_ratio:.1f}/W:{insert_ratio:.1f}"
    
    def get_ordered_index_types(self, df):
        """Get ordered list of index types present in data."""
        actual_indices = df['index_type'].unique()
        desired_order = ['LiBox', 'ALEX', 'LIPP', 'ART'] # for memory vs traces
        # desired_order = ['LiBox', 'ALEX+', 'LIPP+', 'ART', 'XIndex', 'B+tree', 'LOFT'] # for basic type
        return [idx for idx in desired_order if idx in actual_indices]
    
    def get_traces(self, df):
        """Get ordered list of traces present in data."""
        traces = ['msr_web', 'w048', 'longitudes-200M', 'fb', 'genome', 'osm']
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
        ordered_index_types = self.get_ordered_index_types(df)
        
        # Filter for thread_num = 40 and range search workload
        df_filtered = df[df['thread_num'] == 40]
        
        fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(15, 10), sharex=True, sharey=True)
        axes = axes.flatten()
        
        legend_handles = {}
        
        for j, trace in enumerate(traces):
            if j >= len(axes):
                break
                
            ax = axes[j]
            subset = df_filtered[df_filtered['trace'] == trace]
            
            # Assuming range search data has a specific column or workload type
            # Modify this condition based on your actual data structure
            range_search_data = subset  # You may need to filter further based on range search criteria
            
            for index_name in ordered_index_types:
                data = range_search_data[range_search_data['index_type'] == index_name]
                if not data.empty:
                    color = self.color_mapping.get(index_name, 'gray')
                    marker = self.marker_mapping.get(index_name, 'o')
                    
                    # For range search, you might want to plot against different parameters
                    # This is a placeholder - adjust based on your range search data structure
                    line, = ax.plot(data.index, data['throughput'] / 1e6,
                                  color=color, marker=marker, linestyle='-',
                                  markersize=6, linewidth=1.8, markeredgecolor='black',
                                  markeredgewidth=0.1, label=index_name)
                    
                    legend_handles[index_name] = line
            
            ax.set_title(trace_titles[j], fontsize=12)
            ax.set_xlabel("Range Size", fontsize=11)  # Adjust based on your x-axis
            if j % 3 == 0:
                ax.set_ylabel("Range Search Throughput (Mop/s)", fontsize=11)
            
            self.style_axis(ax)
            ax.set_ylim(bottom=0)
        
        # Hide unused subplots
        for k in range(len(traces), len(axes)):
            axes[k].set_visible(False)
        
        self.add_legend(fig, legend_handles, ordered_index_types)
        plt.tight_layout(rect=[0, 0.02, 1, 0.92])
        plt.subplots_adjust(bottom=0.15, wspace=0.25, hspace=0.35)
        
        self.save_figure(fig, output_file)

    def plot_memory_vs_traces(self, df, output_file):
        """Generate memory consumption vs bulk-loading ratio plots."""
        traces = self.get_traces(df)
        trace_titles = [self.trace_map.get(t, t) for t in traces]
        ordered_index_types = self.get_ordered_index_types(df)
        # print(f"index types: {ordered_index_types}")
        
        bulk_ratios = [0.2, 0.4, 0.6, 0.8, 1.0]
        bulk_labels = ['20%', '40%', '60%', '80%', '100%']
        
        fig, axes = plt.subplots(nrows=1, ncols=6, figsize=(16, 3.5), sharex=True)
        axes = axes.flatten()
        
        legend_handles = {}
        
        for j, trace in enumerate(traces):
            if j >= len(axes):
                break
                
            ax = axes[j]
            
            for index_name in ordered_index_types:
                memories = []
                for bulk_ratio in bulk_ratios:
                    data = df[(df['trace'] == trace) & 
                             (df['init_table_ratio'] == bulk_ratio) & 
                             (df['index_type'] == index_name)]
                    # print(f"index: {index_name} memory: {data['memory_consumption']}")
                    if not data.empty:
                        memory = data['memory_consumption'].mean() / 1e9  # Convert to GB
                        memories.append(memory)
                    else:
                        memories.append(None)
                
                # Filter out None values
                valid_indices = [i for i, m in enumerate(memories) if m is not None]
                valid_ratios = [bulk_ratios[i] * 100 for i in valid_indices]
                valid_memories = [memories[i] for i in valid_indices]
                
                if valid_memories:
                    color = self.color_mapping.get(index_name, 'gray')
                    marker = self.marker_mapping.get(index_name, 'o')
                    
                    line, = ax.plot(valid_ratios, valid_memories,
                                  color=color, marker=marker, linestyle='-',
                                  markersize=6, linewidth=1.8, markeredgecolor='black',
                                  markeredgewidth=0.1, label=index_name)
                    
                    legend_handles[index_name] = line
            
            ax.set_xlabel("Bulk-loading Ratio (%)", fontsize=11)
            if j % 3 == 0:
                ax.set_ylabel("Memory Consumption (GB)", fontsize=11)
                
            ax.annotate(trace_titles[j], xy=(0.5, -0.25), xycoords='axes fraction', 
                      ha='center', va='center', fontsize=11)
            
            ax.set_xticks([20, 40, 60, 80, 100])
            ax.set_xlim(15, 105)
            self.style_axis(ax)
            ax.set_ylim(bottom=0)
        
        # Hide unused subplots
        for k in range(len(traces), len(axes)):
            axes[k].set_visible(False)
        
        self.add_legend(fig, legend_handles, ordered_index_types)
        plt.tight_layout(rect=[0, 0.02, 1, 0.92])
        plt.subplots_adjust(bottom=0.15, wspace=0.25, hspace=0.35)
        
        self.save_figure(fig, output_file)

    def plot_alpha_beta_analysis(self, df, output_file):
        """Generate Alpha and Beta analysis plots for LiBox only."""
        # Filter for w048 trace, thread_num=40, and LiBox only
        df_filtered = df[(df['trace'] == 'w048') & 
                        (df['thread_num'] == 40) & 
                        (df['index_type'] == 'LiBox')]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Alpha analysis (subplot 1)
        alpha_values = [0.1, 0.2, 0.3, 0.4, 0.5]
        alpha_labels = ['10%', '20%', '30%', '40%', '50%']
        
        ro_throughputs = []
        wo_throughputs = []
        memories = []
        
        for alpha in alpha_values:
            # You'll need to add alpha column to your data or filter based on alpha parameter
            # This is a placeholder - adjust based on your actual data structure
            alpha_data = df_filtered  # Filter by alpha value when available
            
            ro_data = alpha_data[alpha_data['workload'] == 'RO']
            wo_data = alpha_data[alpha_data['workload'] == 'WO']
            
            if not ro_data.empty:
                ro_throughputs.append(ro_data['throughput'].iloc[0] / 1e6)
                memories.append(ro_data['memory_consumption'].iloc[0] / 1e9)
            else:
                ro_throughputs.append(None)
                memories.append(None)
                
            if not wo_data.empty:
                wo_throughputs.append(wo_data['throughput'].iloc[0] / 1e6)
            else:
                wo_throughputs.append(None)
        
        # Plot throughput on left y-axis
        ax1_twin = ax1.twinx()
        
        valid_alphas = [i for i, (ro, wo) in enumerate(zip(ro_throughputs, wo_throughputs)) 
                       if ro is not None or wo is not None]
        
        if valid_alphas:
            alpha_pcts = [alpha_values[i] * 100 for i in valid_alphas]
            valid_ro = [ro_throughputs[i] for i in valid_alphas if ro_throughputs[i] is not None]
            valid_wo = [wo_throughputs[i] for i in valid_alphas if wo_throughputs[i] is not None]
            valid_mem = [memories[i] for i in valid_alphas if memories[i] is not None]
            
            if valid_ro:
                ax1.plot(alpha_pcts, valid_ro, 'o-', color='blue', label='RO Throughput', linewidth=2)
            if valid_wo:
                ax1.plot(alpha_pcts, valid_wo, 's-', color='red', label='WO Throughput', linewidth=2)
            if valid_mem:
                ax1_twin.plot(alpha_pcts, valid_mem, '^-', color='green', label='Memory', linewidth=2)
        
        ax1.set_xlabel('Alpha (%)', fontsize=11)
        ax1.set_ylabel('Throughput (Mop/s)', color='black', fontsize=11)
        ax1_twin.set_ylabel('Memory (GB)', color='green', fontsize=11)
        ax1.set_title('Alpha vs Performance & Memory', fontsize=12)
        
        # Beta analysis (subplot 2)
        beta_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
        
        # Similar logic for beta analysis
        # This is a placeholder - implement based on your beta parameter data
        
        ax2.set_xlabel('Beta (%)', fontsize=11)
        ax2.set_ylabel('Throughput (Mop/s)', color='black', fontsize=11)
        ax2.set_title('Beta vs Performance & Memory', fontsize=12)
        
        self.style_axis(ax1)
        self.style_axis(ax2)
        
        # Add legends
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax1_twin.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='best')
        
        plt.tight_layout()
        self.save_figure(fig, output_file)

    def save_figure(self, fig, output_file):
        """Save figure in both PDF and PNG formats."""
        # Create output directory if it doesn't exist
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        print(f"Saving figure to: {output_file}")
        plt.savefig(output_file, bbox_inches='tight', dpi=self.dpi)
        
        # Also save as PNG
        png_output = output_file.replace('.pdf', '.png')
        plt.savefig(png_output, bbox_inches='tight', dpi=300)
        print(f"Also saved as: {png_output}")
        
        plt.close()
        print("Figure generation complete.")

def main():
    """Main function with argument parsing."""
    parser = argparse.ArgumentParser(description='Generate various index performance plots')
    parser.add_argument('--input', type=str, help='Input CSV file path')
    parser.add_argument('--output', type=str, help='Output file path (without extension)')
    parser.add_argument('--plot_type', type=str, default='basic',
                       choices=['basic', 'read_ratio', 'range_search', 'memory', 'alpha_beta'],
                       help='Type of plot to generate (default: basic)')
    parser.add_argument('--dpi', type=int, default=600, help='DPI for the output figure (default: 600)')
    
    args = parser.parse_args()
    
    # Generate default file paths if not provided
    if args.input is None or args.output is None:
        date_tag = datetime.now().strftime("%m%d")
        
    if args.input is None:
        args.input = f"../result/libox/out_libox_{date_tag}.csv"
    
    if args.output is None:
        args.output = f"../result/libox/graph_{args.plot_type}_{date_tag}"
    
    # Initialize plotter
    plotter = ThroughputPlotter(dpi=args.dpi)
    
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
    
    print(f"Plot generation completed for type: {args.plot_type}")

if __name__ == "__main__":
    main()