import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import os

# =============================================================================
# 1. SETTINGS AND FILENAMES
# =============================================================================

# Input files
data_files = {
    "open_apo_pc1": "../final_data/open/apo/pca/eigrmsf-5t-6v3f-no-col-pc1.xvg",
    "open_apo_pc2": "../final_data/open/apo/pca/eigrmsf-5t-6v3f-no-col-pc2.xvg",
    "open_bound_pc1": "../final_data/open/bound/pca/eigrmsf-5t-6v3f-col-pc1.xvg",
    "open_bound_pc2": "../final_data/open/bound/pca/eigrmsf-5t-6v3f-col-pc2.xvg",
    "closed_apo_pc1": "../final_data/closed/apo/pca/eigrmsf-5t-6v3h-no-col-pc1.xvg",
    "closed_apo_pc2": "../final_data/closed/apo/pca/eigrmsf-5t-6v3h-no-col-pc2.xvg",
    "closed_bound_pc1": "../final_data/closed/bound/pca/eigrmsf-5t-6v3h-col-pc1.xvg",
    "closed_bound_pc2": "../final_data/closed/bound/pca/eigrmsf-5t-6v3h-col-pc2.xvg",
}

# Output directories
out_dirs = {
    "open": "../final_plots/open/bound-vs-apo/",
    "closed": "../final_plots/closed/bound-vs-apo/"
}

# Create output directories
for d in out_dirs.values():
    os.makedirs(d, exist_ok=True)

# Colors from 6-rmsf-analysis-ca-protein-cumulative-trajs.py
colors = {
    "open_bound": "blue",
    "open_apo": "red",
    "closed_bound": "green",
    "closed_apo": "purple"
}

# =============================================================================
# 2. FUNCTIONS
# =============================================================================
def read_xvg(file_path):
    """Function to read data from XVG file"""
    x = []
    y = []
    if not os.path.exists(file_path):
        print(f"Warning: File {file_path} not found.")
        return [], []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('@') or line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                x.append(float(parts[0]))
                y.append(float(parts[1]))
    return x, y

def plot_comparison(file1_path, file2_path, label1, label2, color1, color2, title, output_path, y_max=1.6):
    """Generic function to plot two RMSF datasets for comparison"""
    x1, y1 = read_xvg(file1_path)
    x2, y2 = read_xvg(file2_path)

    if not x1 or not x2:
        print(f"Skipping plot {title} due to missing data.")
        return

    # Add +21 to residue numbers (as requested)
    res1 = [xi + 21 for xi in x1]
    res2 = [xi + 21 for xi in x2]

    plt.figure(figsize=(8, 5))
    ax = plt.gca()

    LABEL_SIZE = 14
    TITLE_SIZE = 16
    LEGEND_SIZE = 12
    TICK_SIZE = 12

    # Plot lines
    ax.plot(res1, y1, color=color1, linewidth=1.0, alpha=0.85, label=label1, zorder=3)
    ax.plot(res2, y2, color=color2, linewidth=1.0, alpha=0.85, label=label2, zorder=3)

    # Axes and grid configuration
    ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(50))
    ax.grid(True, which='major', axis='both', linestyle='--', linewidth=0.5, alpha=0.7, zorder=1)
    ax.grid(True, which='minor', axis='x', linestyle=':', linewidth=0.4, alpha=0.5, zorder=1)

    # Labels and title
    plt.xlabel('Residue number', fontsize=LABEL_SIZE)
    plt.ylabel('RMSF (nm)', fontsize=LABEL_SIZE)
    plt.title(title, fontsize=TITLE_SIZE)
    plt.xticks(fontsize=TICK_SIZE)
    plt.yticks(fontsize=TICK_SIZE)
    plt.ylim(0, y_max)

    # Legend
    legend_patches = [
        mpatches.Patch(color=color1, label=label1),
        mpatches.Patch(color=color2, label=label2)
    ]
    plt.legend(handles=legend_patches, fontsize=LEGEND_SIZE, loc='upper right', framealpha=1.0)

    # Coloured background bars to highlight domains (NTD, MLD, CTD)
    ax.axvspan(22, 263,   color='#b2e7fa', alpha=0.7, zorder=0, label='_nolegend_')  # NTD
    ax.axvspan(379, 632,  color='#b2e7ca', alpha=0.7, zorder=0, label='_nolegend_')  # MLD
    ax.axvspan(882, 1105, color='#f3d3f0', alpha=0.7, zorder=0, label='_nolegend_')  # CTD

    plt.tight_layout()
    plt.savefig(output_path, dpi=600, transparent=False, facecolor='white')
    plt.close()
    print(f"Chart saved in: {output_path}")

# =============================================================================
# 3. EXECUTION
# =============================================================================

# 1) confronto open bound - open apo, su PC1
plot_comparison(
    data_files["open_bound_pc1"], data_files["open_apo_pc1"],
    "Open bound", "Open apo",
    colors["open_bound"], colors["open_apo"],
    "Cα RMSF on PC1 ",
    os.path.join(out_dirs["open"], "rmsf-pc1-bound-vs-apo.png"),
    y_max=1.6
)

# 2) confronto open bound - open apo, su PC2
plot_comparison(
    data_files["open_bound_pc2"], data_files["open_apo_pc2"],
    "Open bound", "Open apo",
    colors["open_bound"], colors["open_apo"],
    "Cα RMSF on PC2 ",
    os.path.join(out_dirs["open"], "rmsf-pc2-bound-vs-apo.png"),
    y_max=0.8
)

# 3) confronto closed bound - closed apo, su PC1
plot_comparison(
    data_files["closed_bound_pc1"], data_files["closed_apo_pc1"],
    "Closed bound", "Closed apo",
    colors["closed_bound"], colors["closed_apo"],
    "Cα RMSF on PC1 ",
    os.path.join(out_dirs["closed"], "rmsf-pc1-bound-vs-apo.png"),
    y_max=1.6
)

# 4) confronto closed bound - closed apo, su PC2
plot_comparison(
    data_files["closed_bound_pc2"], data_files["closed_apo_pc2"],
    "Closed bound", "Closed apo",
    colors["closed_bound"], colors["closed_apo"],
    "Cα RMSF on PC2 ",
    os.path.join(out_dirs["closed"], "rmsf-pc2-bound-vs-apo.png"),
    y_max=0.8
)
