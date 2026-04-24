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
# 1. DEFINE THE 4 SYSTEMS (uncomment the block you want to analyze)
# =============================================================================
replicas = ["run-md1", "run-md2", "run-md3", "run-md4", "run-md5"]

# --- BLOCK 1A: Open bound ---
system_id = "open_bound"
system_title = "Open bound"
traj_base = "/mnt/h/Il mio Drive/LAVORO_MD_NPC1L1_nov25/MD_6V3F_6V3H_500ns_sept25/6v3f/data/col/"
topology = "prot-lig.prmtop"
traj_name = "07-08-prot-lig-pbc.nc"
results_dir = "../final_data/open/bound/"
plots_dir = "../final_plots/open/bound/"

# --- BLOCK 1B: Open apo ---
# system_id = "open_apo"
# system_title = "Open apo"
# traj_base = "/mnt/h/Il mio Drive/LAVORO_MD_NPC1L1_nov25/MD_6V3F_6V3H_500ns_sept25/6v3f/data/no-col/"
# topology = "prot.prmtop"
# traj_name = "07-08-prot-pbc.nc"
# results_dir = "../final_data/open/apo/"
# plots_dir = "../final_plots/open/apo/"

# --- BLOCK 1C: Closed bound ---
#system_id = "closed_bound"
#system_title = "Closed bound"
#traj_base = "/mnt/h/Il mio Drive/LAVORO_MD_NPC1L1_nov25/MD_6V3F_6V3H_500ns_sept25/6v3h/data/col/"
#topology = "prot-lig.prmtop"
#traj_name = "07-08-prot-lig-pbc.nc"
#results_dir = "../final_data/closed/bound/"
#plots_dir = "../final_plots/closed/bound/"

# --- BLOCK 1D: Closed apo ---
# system_id = "closed_apo"
# system_title = "Closed apo"
# traj_base = "/mnt/h/Il mio Drive/LAVORO_MD_NPC1L1_nov25/MD_6V3F_6V3H_500ns_sept25/6v3h/data/no-col/"
# topology = "prot.prmtop"
# traj_name = "07-08-prot-pbc.nc"
# results_dir = "../final_data/closed/apo/"
# plots_dir = "../final_plots/closed/apo/"


# Ensure directories exist
os.makedirs(results_dir, exist_ok=True)
os.makedirs(plots_dir, exist_ok=True)
# -----------------------------------------------------------------------------


# =============================================================================
# 2. RMSF CALCULATION (LOOP OVER REPLICAS)
# =============================================================================
# define residues for fitting and residues for rmsf calculation
fit_sel = 'name CA and (resnum 612:782 or resnum 1092:1242)'
rmsf_sel = 'resnum 1:1261 and name CA'

dataframes = []

for i, rep in enumerate(replicas, start=1):
    print(f"Analysing {system_title} - Replica {i}: {rep}...")
    
    traj_path = os.path.join(traj_base, rep, traj_name)
    top_path = os.path.join(traj_base, topology)
    u = mda.Universe(top_path, traj_path)
    
    # 1. Calculate average structure
    print(f"  - Calculating average structure...")
    average = align.AverageStructure(u, u, select=fit_sel, ref_frame=0).run(step=100)
    ref = average.results.universe
    
    # 2. Align trajectory to average structure
    print(f"  - Aligning trajectory...")
    align.AlignTraj(u, ref, select=fit_sel, in_memory=True).run(step=100)
    
    # 3. Calculate RMSF
    print(f"  - Calculating RMSF...")
    c_alphas = u.select_atoms(rmsf_sel)
    R = rms.RMSF(c_alphas).run(step=100)
    
    # 4. Adjust residue numbering (+21)  ((residue numbering in .nc trajectory starts with 1, 
    #residue numbering in the original pdb files starts with 22)
    res_ids = [res.resid + 21 for res in c_alphas.residues]
    col_name = f'RMSF_md{i} (Å)'
    
    df = pd.DataFrame({
        'res id': res_ids,
        col_name: R.results.rmsf.round(2)
    })
    dataframes.append(df)

# =============================================================================
# 3. MERGE DATAFRAMES
# =============================================================================
# Merge all replicas on column 'res id'
merged_df = reduce(lambda left, right: pd.merge(left, right, on='res id'), dataframes)

# Save .csv file
csv_filename = f"rmsf-ca-single-replica-{system_id}.csv"
csv_path = os.path.join(results_dir, csv_filename)
merged_df.to_csv(csv_path, sep="\t", index=False)
print(f"\nData saved in: {csv_path}")

# =============================================================================
# 4. PLOTTING
# =============================================================================
LABEL_SIZE = 14
TITLE_SIZE = 16
LEGEND_SIZE = 12
TICK_SIZE = 12

# Colors and labels corresponding to the 5 replicas
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
labels = [f'md{i}' for i in range(1, len(replicas) + 1)]

fig, ax = plt.subplots(figsize=(8, 5))

# Plot lines for each replica
for i, rep in enumerate(replicas, start=1):
    col_name = f'RMSF_md{i} (Å)'
    ax.plot(merged_df['res id'], merged_df[col_name], color=colors[i-1], 
            linewidth=1.0, alpha=0.85, label=labels[i-1], zorder=3)

# Axes and grid configuration
ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(50))
ax.grid(True, which='major', axis='both', linestyle='--', linewidth=0.5, alpha=0.7, zorder=1)
ax.grid(True, which='minor', axis='x', linestyle=':', linewidth=0.4, alpha=0.5, zorder=1)

# Labels and title
plt.xlabel('Residue number', fontsize=LABEL_SIZE)
plt.ylabel('RMSF (Å)', fontsize=LABEL_SIZE)
plt.title(f'Cα RMSF - {system_title}', fontsize=TITLE_SIZE)
plt.xticks(fontsize=TICK_SIZE)
plt.yticks(fontsize=TICK_SIZE)
plt.ylim(0,16.0)

# Legend
legend_patches = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(len(replicas))]
plt.legend(handles=legend_patches, fontsize=LEGEND_SIZE, loc='upper right', framealpha=1.0)


# Coloured background bars to highlight domains (NTD, MLD, CTD)
ax.axvspan(22, 263,   color='#b2e7fa', alpha=0.7, zorder=0, label='_nolegend_')  # NTD
ax.axvspan(379, 632,  color='#b2e7ca', alpha=0.7, zorder=0, label='_nolegend_')  # MLD
ax.axvspan(882, 1105, color='#f3d3f0', alpha=0.7, zorder=0, label='_nolegend_')  # CTD

# Adjust layout
plt.tight_layout()

# Save the plot
plot_filename = f"rmsf-ca-single-replica-{system_id}.png"
plot_path = os.path.join(plots_dir, plot_filename)
plt.savefig(plot_path, dpi=600, transparent=False, facecolor='white')
print(f"Chart saved in: {plot_path}")
