import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from matplotlib.patches import Patch

# =============================================================================
# 1. CONFIGURATION
# =============================================================================
replicas = ["md1", "md2", "md3", "md4", "md5"]
filename_template = "plot12-{rep}.xvg" # e.g., plot12-md1.xvg

systems = [
    {
        "system_id": "open_bound",
        "title": "Open conformation, bound",
        "data_dir": "../final_data/open/bound/pca/",
        "results_dir": "../final_plots/open/bound/",
    },
    {
        "system_id": "open_apo",
        "title": "Open conformation, apo",
        "data_dir": "../final_data/open/apo/pca/",
        "results_dir": "../final_plots/open/apo/",
    },
    {
        "system_id": "closed_bound",
        "title": "Closed conformation, bound",
        "data_dir": "../final_data/closed/bound/pca/",
        "results_dir": "../final_plots/closed/bound/",
    },
    {
        "system_id": "closed_apo",
        "title": "Closed conformation, apo",
        "data_dir": "../final_data/closed/apo/pca/",
        "results_dir": "../final_plots/closed/apo/",
    },
]

# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================
def read_data(file_path):
    """Reads x and y coordinates from .xvg file, skipping header."""
    x, y = [], []
    try:
        with open(file_path, 'r') as file:
            # Skip first 17 lines (header).
            lines = file.readlines()[17:]
            for line in lines:
                if line.startswith(('@', '#')) or not line.strip():
                    continue
                v = line.split()
                if len(v) >= 2:
                    x.append(float(v[0]))
                    y.append(float(v[1]))
        return x, y
    except FileNotFoundError:
        return [], []

# =============================================================================
# 3. DATA LOADING AND PLOTTING
# =============================================================================
LABEL_SIZE = 14
TITLE_SIZE = 16
LEGEND_SIZE = 12
TICK_SIZE = 12
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']

for system in systems:
    system_id = system["system_id"]
    title = system["title"]
    data_dir = system["data_dir"]
    results_dir = system["results_dir"]

    print(f"\n>>> Generating Scatter Plot for: {title} <<<")
    os.makedirs(results_dir, exist_ok=True)

    plt.figure(figsize=(8, 6))
    legend_elements_scatter = []
    data_found = False

    for i, rep in enumerate(replicas):
        file_name = filename_template.format(rep=rep)
        file_path = os.path.join(data_dir, file_name)
        
        x, y = read_data(file_path)
        
        if not x:
            print(f"  [WARNING] Skipping {rep} (file not found or empty: {file_name})")
            continue
        
        data_found = True
        print(f"  - Loading {rep} data...")
        
        # Plotting the scatter points for the current replica
        plt.scatter(x, y, label=rep, color=colors[i], marker='.', alpha=0.5, s=5, zorder=0)
        
        # Creating patch for the legend
        legend_elements_scatter.append(
            Patch(facecolor=colors[i], edgecolor=colors[i], alpha=0.45, label=rep)
        )

    if not data_found:
        print(f"  [ERROR] No data found for system {system_id}. Skipping plot.")
        plt.close()
        continue

    # Axes lines
    plt.axhline(0, color='black', linewidth=0.6, zorder=1)
    plt.axvline(0, color='black', linewidth=0.6, zorder=1)

    # Grid
    plt.grid(True, linestyle='--', alpha=0.5, zorder=0)

    # Labels and title
    plt.xlabel("trajectory projection on PC1 (nm)", fontsize=LABEL_SIZE)
    plt.ylabel("trajectory projection on PC2 (nm)", fontsize=LABEL_SIZE)
    plt.title(title, fontsize=TITLE_SIZE)
    plt.xticks(fontsize=TICK_SIZE)
    plt.yticks(fontsize=TICK_SIZE)

    # Legend
    plt.legend(handles=legend_elements_scatter, loc='upper right', fontsize=LEGEND_SIZE, framealpha=1.0)

    plt.tight_layout()

    # Save plot
    out_scatter = os.path.join(results_dir, f"pca_scatter_{system_id}_replicas.png")
    plt.savefig(out_scatter, dpi=600, bbox_inches='tight')
    plt.close()

    print(f"  Chart saved in: {out_scatter}")

print("\nAnalysis completed successfully!")

