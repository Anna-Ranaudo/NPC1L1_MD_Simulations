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
fit_sel = 'name CA and resnum 1:242'

# Target residues for RMSF (pdb numbering, need to subtract 21 for 1-indexed trajectory numbering)
target_residues_pdb = [34, 37, 52, 53, 54, 55, 95, 98, 99, 102, 103, 106, 120, 124, 127, 128, 156, 187, 188, 205, 206, 207, 211, 213, 214, 215, 216, 218]
target_residues_traj = [res - 21 for res in target_residues_pdb]
resnum_str = " ".join(map(str, target_residues_traj))

rmsf_sel = f'resnum {resnum_str} and name CA'
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
csv_path = os.path.join(results_dir, "rmsf-ca-protein-NTD-cholesterol-residues.csv")
merged_df.to_csv(csv_path, sep="\t", index=False)
print(f"\ndata saved in: {csv_path}")

# =============================================================================
# 4. PLOTTING
# =============================================================================
import seaborn as sns

LABEL_SIZE = 14
TITLE_SIZE = 16
LEGEND_SIZE = 12
TICK_SIZE = 12

# Filter systems to plot
target_sys_ids = ["open_bound", "open_apo", "closed_bound", "closed_apo"]
cols_to_plot = [f'RMSF {sys_id} (Å)' for sys_id in target_sys_ids]

# Select only relevant columns
plot_df = merged_df[['res id'] + cols_to_plot].copy()

# Rename columns for the legend (using labels from the systems list)
rename_dict = {}
for sys in systems:
    if sys["id"] in target_sys_ids:
        rename_dict[f'RMSF {sys["id"]} (Å)'] = sys["label"]

plot_df = plot_df.rename(columns=rename_dict)

# Melt the dataframe for seaborn grouped bar plot
df_melted = plot_df.melt(id_vars='res id', var_name='System', value_name='RMSF (Å)')

sns.set_theme(style="whitegrid", palette="muted")
plt.figure(figsize=(14, 6))

# Define colors based on original systems definitions
palette = {sys["label"]: sys["color"] for sys in systems if sys["id"] in target_sys_ids}

bar_plot = sns.barplot(
    data=df_melted, 
    x='res id', 
    y='RMSF (Å)', 
    hue='System',
    palette=palette,
    edgecolor='0.2',
    linewidth=0.5
)

# Customizing axes
plt.title('RMSF Comparison: residues near cholesterol', fontsize=TITLE_SIZE, fontweight='bold', pad=20)
plt.xlabel('Residue number', fontsize=LABEL_SIZE, labelpad=10)
plt.ylabel('RMSF (Å)', fontsize=LABEL_SIZE, labelpad=10)
plt.xticks(rotation=45, fontsize=TICK_SIZE)
plt.yticks(fontsize=TICK_SIZE)

# Legend adjustment
plt.legend(title='System', title_fontsize=LABEL_SIZE, fontsize=LEGEND_SIZE, loc='upper right', frameon=True, shadow=True)

plt.tight_layout()

# save the png
plot_path = os.path.join(plots_dir, "rmsf-ca-protein-NTD-cholesterol-residues.png")
plt.savefig(plot_path, dpi=600, transparent=True)
print(f"chart saved in: {plot_path}")
