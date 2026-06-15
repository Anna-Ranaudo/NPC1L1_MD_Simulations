import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import numpy as np
import os
from functools import reduce

# =============================================================================
# 1. DEFINE THE 4 SYSTEMS
# =============================================================================
results_dir = "../final_data/comparison_four_systems/"
plots_dir = "../final_plots/comparison_four_systems/"
os.makedirs(results_dir, exist_ok=True)
os.makedirs(plots_dir, exist_ok=True)

replicas = ["run-md1", "run-md2", "run-md3", "run-md4", "run-md5"]

# Base path for all trajectories
base_path = "/mnt/h/Il mio Drive/LAVORO_MD_NPC1L1_nov25/MD_6V3F_6V3H_500ns_sept25/"

systems = [
    {
        "id": "open_bound",
        "label": "Open bound",
        "color": "blue",
        "traj_base": os.path.join(base_path, "6v3f/data/col/"),
        "topology": "prot-lig.prmtop",
        "traj_name": "07-08-prot-lig-pbc.nc"
    },
    {
        "id": "open_apo",
        "label": "Open apo",
        "color": "red",
        "traj_base": os.path.join(base_path, "6v3f/data/no-col/"),
        "topology": "prot.prmtop",
        "traj_name": "07-08-prot-pbc.nc"
    },
    {
        "id": "closed_bound",
        "label": "Closed bound",
        "color": "green",
        "traj_base": os.path.join(base_path, "6v3h/data/col/"),
        "topology": "prot-lig.prmtop",
        "traj_name": "07-08-prot-lig-pbc.nc"
    },
    {
        "id": "closed_apo",
        "label": "Closed apo",
        "color": "purple",
        "traj_base": os.path.join(base_path, "6v3h/data/no-col/"),
        "topology": "prot.prmtop",
        "traj_name": "07-08-prot-pbc.nc"
    }
]


# =============================================================================
# 2. RMSF CALCULATION
# =============================================================================

# define residues for fitting and residues for rmsf calculation
fit_sel = 'name CA and (resnum 612:782 or resnum 1092:1242)'
rmsf_sel = 'resnum 1:1261 and name CA'
dataframes = []

for sys in systems:
    print(f"Analyzing system: {sys['label']}...")
    
    # 1. generate the 5 replicas paths for the current system and check existence
    potential_traj_paths = [os.path.join(sys["traj_base"], rep, sys["traj_name"]) for rep in replicas]
    traj_paths = [p for p in potential_traj_paths if os.path.exists(p)]
    
    top_path = os.path.join(sys["traj_base"], sys["topology"])
    
    if not os.path.exists(top_path):
        print(f"  [ERROR] Topology not found: {top_path}. Skipping system.")
        continue
        
    if not traj_paths:
        print(f"  [WARNING] No replicas found for {sys['label']}. Skipping.")
        continue
        
    if len(traj_paths) < len(replicas):
        print(f"  [INFO] Found {len(traj_paths)}/{len(replicas)} replicas.")

    # 2. load the available replicas in a MDAnalysis universe
    u = mda.Universe(top_path, *traj_paths)
    
    # 3. calculate average structure + fitting
    print(f"  - calculating average structure...")
    average = align.AverageStructure(u, u, select=fit_sel, ref_frame=0).run(step=100)
    ref = average.results.universe
    
    print(f"  - aligning trajectory...")
    align.AlignTraj(u, ref, select=fit_sel, in_memory=True).run(step=100)
    
    # 4. Calculate RMSF
    print(f"  - Calculating RMSF...")
    c_alphas = u.select_atoms(rmsf_sel)
    R = rms.RMSF(c_alphas).run(step=100)
    
    # 5. add +21 to res number (residue numbering in .nc trajectory starts with 1, residue numbering in the original pdb files starts with 22)
    # and save the dataframe
    res_ids = [res.resid + 21 for res in c_alphas.residues]
    col_name = f'RMSF {sys["id"]} (Å)'
    
    df = pd.DataFrame({
        'res id': res_ids,
        col_name: R.results.rmsf.round(2)
    })
    dataframes.append(df)

# =============================================================================
# 3. merge dataframes
# =============================================================================
# merge the 4 dataframes on column 'res id'
merged_df = reduce(lambda left, right: pd.merge(left, right, on='res id'), dataframes)

# save .csv file
csv_path = os.path.join(results_dir, "rmsf-ca-protein-ave-5rep-all-systems.csv")
merged_df.to_csv(csv_path, sep="\t", index=False)
print(f"\ndata saved in: {csv_path}")

# =============================================================================
# 4. PLOTTING
# =============================================================================
LABEL_SIZE = 14
TITLE_SIZE = 16
LEGEND_SIZE = 12
TICK_SIZE = 12

plt.figure(figsize=(8, 5))

# Plot 
for sys in systems:
    col_name = f'RMSF {sys["id"]} (Å)'
    plt.plot(merged_df['res id'], merged_df[col_name], color=sys['color'], linewidth=1.0, alpha=0.85)

# axes and grid
ax = plt.gca()
ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(50))
ax.grid(True, which='major', axis='both', linestyle='--', linewidth=0.5, alpha=0.7)
ax.grid(True, which='minor', axis='x', linestyle=':', linewidth=0.4, alpha=0.5)

plt.xlabel('Residue number', fontsize=LABEL_SIZE)
plt.ylabel('RMSF (Å)', fontsize=LABEL_SIZE)
plt.title('Cα RMSF', fontsize=TITLE_SIZE)
plt.xticks(fontsize=TICK_SIZE)
plt.yticks(fontsize=TICK_SIZE)

# legend
legend_patches = [mpatches.Patch(color=sys['color'], label=sys['label']) for sys in systems]
plt.legend(handles=legend_patches, fontsize=LEGEND_SIZE, loc='upper right')


# coloured background bars to highlight residues belonging to the domains NTD, MLD, CTD
ax.axvspan(22, 263,   color='#b2e7fa', alpha=0.7, zorder=0, label='_nolegend_')  # NTD
ax.axvspan(379, 632,  color='#b2e7ca', alpha=0.7, zorder=0, label='_nolegend_')  # MLD
ax.axvspan(882, 1105, color='#f3d3f0', alpha=0.7, zorder=0, label='_nolegend_')  #CTD

plt.tight_layout()

# save the png
plot_path = os.path.join(plots_dir, "rmsf-ca-protein-ave-5rep-all-systems.png")
plt.savefig(plot_path, dpi=600, transparent=True)
print(f"chart saved in: {plot_path}")
