import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
from matplotlib.patches import Patch
import os

# =============================================================================
# 1. CONFIGURATION
# =============================================================================
# Define which pairs you want to plot together
# Each tuple contains the IDs of the two systems to compare
plot_pairs = [
    
    #("open_bound", "open_apo"),
    ("closed_bound", "closed_apo"),
    #("open_bound", "closed_bound"),
    #("open_apo", "closed_apo"),
    
]


#title = f"Open conformation, bound vs apo"
title = f"Closed conformation, bound vs apo"
#title = f"Open vs closed conformation, bound"
#title = f"Open vs closed conformation, apo"


###### to compare bound and apo in the same structure: i.e. open or closed
#results_dir = "../final_plots/open/bound-vs-apo/"
results_dir = "../final_plots/closed/bound-vs-apo/"


###### to compare open and closed in the same structure: i.e. bound or apo
#results_dir = "../final_plots/comparison_open-vs-closed/"


os.makedirs(results_dir, exist_ok=True)

# Define the 4 systems with their respective data files and colors
# !!!   here we have 2D projections of the trajectories onto the common PC1-PC2 subspaces [
# obtained from the PCA performed on the concatenated trajectories for the following pairs of systems
# open apo + open bound
# closed apo + closed bound
# open bound + closed bound 
# open apo + closed apo


systems = {

    "open_bound": {
        "label": "Open bound",
        "color": "blue",
        "cmap": "Blues",
        # compare open apo vs open bound
        "file": "../final_data/open/bound-vs-apo/pca-2d-proj/plot-col-5t.xvg" 
        # compare bound open vs closed
        #"file": "../final_data/comparison_open-vs-closed/bound/pca-2d-proj/plot-col-6v3f.xvg"
    },

    "open_apo": {
        "label": "Open apo",
        "color": "red",
        "cmap": "Reds",
        # compare open apo vs open bound
        "file": "../final_data/open/bound-vs-apo/pca-2d-proj/plot-no-col-5t.xvg" 
        # compare apo open vs closed
        #"file": "../final_data/comparison_open-vs-closed/apo/pca-2d-proj/plot-no-col-6v3f.xvg"
    },

    "closed_bound": {
        "label": "Closed bound",
        "color": "green",
        "cmap": "Greens",
        # compare closed apo vs open bound
        "file": "../final_data/closed/bound-vs-apo/pca-2d-proj/plot-col-5t.xvg" 
        # compare bound open vs closed
        #"file": "../final_data/comparison_open-vs-closed/bound/pca-2d-proj/plot-col-6v3h.xvg"
    },

    "closed_apo": {
        "label": "Closed apo",
        "color": "purple",
        "cmap": "Purples",
        # compare closed apo vs open bound
        "file": "../final_data/closed/bound-vs-apo/pca-2d-proj/plot-no-col-5t.xvg" 
        # compare apo open vs closed
        #"file": "../final_data/comparison_open-vs-closed/apo/pca-2d-proj/plot-no-col-6v3h.xvg"
    }
}



# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================
def read_data(file_path):
    """Reads x and y coordinates from file, skipping header and taking every 5th frame."""
    x, y = [], []
    try:
        with open(file_path) as f:
            lines = f.readlines()[17:][::5]
            for line in lines:
                v = line.split()
                x.append(float(v[0]))
                y.append(float(v[1]))
        return np.array(x), np.array(y)
    except FileNotFoundError:
        print(f"Warning: File not found: {file_path}")
        return np.array([]), np.array([])

def calculate_kde(x, y, xx, yy):
    """Calculates the 2D Gaussian KDE and normalizes it."""
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([x, y])
    kernel = gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)
    
    # Normalize to 0-1
    f_norm = (f - f.min()) / (f.max() - f.min())
    return f_norm

# =============================================================================
# 3. PLOTTING LOOP
# =============================================================================

for sys1_id, sys2_id in plot_pairs:
    print(f"\nPlotting {sys1_id} vs {sys2_id}...")
    
    sys1 = systems[sys1_id]
    sys2 = systems[sys2_id]
    
    # Load data
    x1, y1 = read_data(sys1["file"])
    x2, y2 = read_data(sys2["file"])
    
    if len(x1) == 0 or len(x2) == 0:
        print(f"Skipping plot {sys1_id}_vs_{sys2_id} due to missing data.")
        continue

    # Create common grid based on both datasets
    xmin = min(x1.min(), x2.min())
    xmax = max(x1.max(), x2.max())
    ymin = min(y1.min(), y2.min())
    ymax = max(y1.max(), y2.max())

    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]

    # Calculate KDE surfaces
    print("  - Calculating KDEs...")
    z1 = calculate_kde(x1, y1, X, Y)
    z2 = calculate_kde(x2, y2, X, Y)

    # Calculate areas covered by the KDE surfaces (threshold >= 5% of maximum)
    threshold = 0.05
    dx = (xmax - xmin) / 99  # 100 grid points mean 99 intervals
    dy = (ymax - ymin) / 99
    cell_area = dx * dy
    
    # Calculate individual areas
    mask1 = z1 >= threshold
    mask2 = z2 >= threshold
    
    area1 = np.sum(mask1) * cell_area
    area2 = np.sum(mask2) * cell_area
    
    # Calculate intersection area and percentages
    mask_intersect = mask1 & mask2
    area_intersect = np.sum(mask_intersect) * cell_area
    
    perc_sys1_in_sys2 = (area_intersect / area1) * 100 if area1 > 0 else 0
    perc_sys2_in_sys1 = (area_intersect / area2) * 100 if area2 > 0 else 0
    
    print(f"  - Area {sys1_id}: {area1:.2f} nm²")
    print(f"  - Area {sys2_id}: {area2:.2f} nm²")
    print(f"  - Overlap area: {area_intersect:.2f} nm²")
    print(f"  - {perc_sys1_in_sys2:.1f}% of {sys1_id} overlaps with {sys2_id}")
    print(f"  - {perc_sys2_in_sys1:.1f}% of {sys2_id} overlaps with {sys1_id}")


    # --------------------------------------------------
    # Plotting
    # --------------------------------------------------
    LABEL_SIZE = 14
    TITLE_SIZE = 16
    LEGEND_SIZE = 12
    TICK_SIZE = 12
    
    plt.figure(figsize=(8, 6))

    levels = np.linspace(0.1, 1.0, 10)
    outer_level = [0.05]   # 5% of max for the thin contour line

    # Plot System 2 first (background)
    plt.contourf(X, Y, z2, levels=levels, cmap=sys2["cmap"], alpha=0.80)
    plt.contour(X, Y, z2, levels=outer_level, colors=sys2["color"], linewidths=0.6)

    # Plot System 1 second (foreground)
    plt.contourf(X, Y, z1, levels=levels, cmap=sys1["cmap"], alpha=0.80)
    plt.contour(X, Y, z1, levels=outer_level, colors=sys1["color"], linewidths=0.6)

    # Axes lines and grid
    plt.axhline(0, color='black', linewidth=0.6)
    plt.axvline(0, color='black', linewidth=0.6)
    plt.grid(True, linestyle='--', alpha=0.4, zorder=0)
    
    #plt.xlim(-32, 36)
    #plt.ylim(-25, 25)

    # Labels and title
    plt.xlabel("trajectory projection on PC1 (nm)", fontsize=LABEL_SIZE)
    plt.ylabel("trajectory projection on PC2 (nm)", fontsize=LABEL_SIZE)
    plt.xticks(fontsize=TICK_SIZE)
    plt.yticks(fontsize=TICK_SIZE)
    
    #title = f"{sys1['label'].split(' ')[0]} conformation, {sys1['label'].split(' ')[1].strip()} vs {sys2['label'].split(' ')[1].strip()}"
    plt.title(title, fontsize=TITLE_SIZE)

    # Legend
    legend_elements = [
        Patch(facecolor=sys1["color"], edgecolor=sys1["color"], label=sys1["label"], alpha=0.8),
        Patch(facecolor=sys2["color"], edgecolor=sys2["color"], label=sys2["label"], alpha=0.8)
    ]
    plt.legend(handles=legend_elements, loc='upper right', fontsize=LEGEND_SIZE, framealpha=1.0)

    plt.tight_layout()

    # Save plot
    out_filename = f"pca_surface_{sys1_id}_vs_{sys2_id}.png"
    out_path = os.path.join(results_dir, out_filename)
    plt.savefig(out_path, dpi=600, bbox_inches='tight')
    plt.close() # Close the figure to free memory
    print(f"Chart saved in: {out_path}")
