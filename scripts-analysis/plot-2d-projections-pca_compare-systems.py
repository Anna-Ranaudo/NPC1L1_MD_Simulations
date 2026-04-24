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
system_metadata = {
    "open_bound": {"label": "Open bound", "color": "blue", "cmap": "Blues"},
    "open_apo": {"label": "Open apo", "color": "red", "cmap": "Reds"},
    "closed_bound": {"label": "Closed bound", "color": "green", "cmap": "Greens"},
    "closed_apo": {"label": "Closed apo", "color": "purple", "cmap": "Purples"},
}

tasks = [
    {
        "title": "Open conformation, bound vs apo",
        "results_dir": "../final_plots/open/bound-vs-apo/",
        "sys1": {**system_metadata["open_bound"], "file": "../final_data/open/bound-vs-apo/pca-2d-proj/plot-col-5t.xvg"},
        "sys2": {**system_metadata["open_apo"], "file": "../final_data/open/bound-vs-apo/pca-2d-proj/plot-no-col-5t.xvg"},
        "id1": "open_bound", "id2": "open_apo"
    },
    {
        "title": "Closed conformation, bound vs apo",
        "results_dir": "../final_plots/closed/bound-vs-apo/",
        "sys1": {**system_metadata["closed_bound"], "file": "../final_data/closed/bound-vs-apo/pca-2d-proj/plot-col-5t.xvg"},
        "sys2": {**system_metadata["closed_apo"], "file": "../final_data/closed/bound-vs-apo/pca-2d-proj/plot-no-col-5t.xvg"},
        "id1": "closed_bound", "id2": "closed_apo"
    },
    {
        "title": "Open vs closed conformation, bound",
        "results_dir": "../final_plots/comparison_open-vs-closed/",
        "sys1": {**system_metadata["open_bound"], "file": "../final_data/comparison_open-vs-closed/bound/pca-2d-proj/plot-col-6v3f.xvg"},
        "sys2": {**system_metadata["closed_bound"], "file": "../final_data/comparison_open-vs-closed/bound/pca-2d-proj/plot-col-6v3h.xvg"},
        "id1": "open_bound", "id2": "closed_bound"
    },
    {
        "title": "Open vs closed conformation, apo",
        "results_dir": "../final_plots/comparison_open-vs-closed/",
        "sys1": {**system_metadata["open_apo"], "file": "../final_data/comparison_open-vs-closed/apo/pca-2d-proj/plot-no-col-6v3f.xvg"},
        "sys2": {**system_metadata["closed_apo"], "file": "../final_data/comparison_open-vs-closed/apo/pca-2d-proj/plot-no-col-6v3h.xvg"},
        "id1": "open_apo", "id2": "closed_apo"
    },
]

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
                if len(v) >= 2:
                    x.append(float(v[0]))
                    y.append(float(v[1]))
        return np.array(x), np.array(y)
    except FileNotFoundError:
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

for task in tasks:
    title = task["title"]
    results_dir = task["results_dir"]
    sys1 = task["sys1"]
    sys2 = task["sys2"]
    sys1_id = task["id1"]
    sys2_id = task["id2"]

    print(f"\n>>> Processing Comparison: {title} <<<")
    
    # Load data
    x1, y1 = read_data(sys1["file"])
    x2, y2 = read_data(sys2["file"])
    
    if len(x1) == 0 or len(x2) == 0:
        print(f"  [WARNING] Skipping comparison {sys1_id}_vs_{sys2_id} due to missing data.")
        print(f"    Missing: {sys1['file'] if len(x1)==0 else ''} {sys2['file'] if len(x2)==0 else ''}")
        continue

    os.makedirs(results_dir, exist_ok=True)

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
    
    print(f"  - Overlap: {perc_sys1_in_sys2:.1f}% of {sys1_id} overlaps with {sys2_id}")
    print(f"  - Overlap: {perc_sys2_in_sys1:.1f}% of {sys2_id} overlaps with {sys1_id}")


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
    
    # Labels and title
    plt.xlabel("trajectory projection on PC1 (nm)", fontsize=LABEL_SIZE)
    plt.ylabel("trajectory projection on PC2 (nm)", fontsize=LABEL_SIZE)
    plt.xticks(fontsize=TICK_SIZE)
    plt.yticks(fontsize=TICK_SIZE)
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
    print(f"  Chart saved in: {out_path}")

print("\nAll comparison tasks completed successfully!")

