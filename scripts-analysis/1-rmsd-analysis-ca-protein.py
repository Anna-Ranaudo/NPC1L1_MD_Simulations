import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from functools import reduce

# =============================================================================
# 1. DEFINE THE 4 SYSTEMS (uncomment the block you want to analyze)
# =============================================================================
replicas = ["run-md1", "run-md2", "run-md3", "run-md4", "run-md5"]

# --- BLOCK 1A: Open bound ---
#traj_dir = "/mnt/h/Il mio Drive/LAVORO_MD_NPC1L1_nov25/MD_6V3F_6V3H_500ns_sept25/6v3f/"
#system_title = "Open bound"
#data_dir = os.path.join(traj_dir, "data/col/")
#topology = "prot-lig.prmtop"
#traj_name = "07-08-prot-lig-pbc.nc"
#results_dir = "../final_data/open/bound/"
#plots_dir = "../final_plots/open/bound/"

# --- BLOCK 1B: Open apo ---
traj_dir = "/mnt/h/Il mio Drive/LAVORO_MD_NPC1L1_nov25/MD_6V3F_6V3H_500ns_sept25/6v3f/"
system_title = "Open apo"
data_dir = os.path.join(traj_dir, "data/no-col/")
topology = "prot.prmtop"
traj_name = "07-08-prot-pbc.nc"
results_dir = "../final_data/open/apo/"
plots_dir = "../final_plots/open/apo/"

# --- BLOCK 1C: Closed bound ---
#traj_dir = "/mnt/h/Il mio Drive/LAVORO_MD_NPC1L1_nov25/MD_6V3F_6V3H_500ns_sept25/6v3h/"
#system_title = "Closed bound"
#data_dir = os.path.join(traj_dir, "data/col/")
#topology = "prot-lig.prmtop"
#traj_name = "07-08-prot-lig-pbc.nc"
#results_dir = "../final_data/closed/bound/"
#plots_dir = "../final_plots/closed/bound/"

# --- BLOCK 1D: Closed apo ---
#traj_dir = "/mnt/h/Il mio Drive/LAVORO_MD_NPC1L1_nov25/MD_6V3F_6V3H_500ns_sept25/6v3h/"
#system_title = "Closed apo"
#data_dir = os.path.join(traj_dir, "data/no-col/")
#topology = "prot.prmtop"
#traj_name = "07-08-prot-pbc.nc"
#results_dir = "../final_data/closed/apo/"
#plots_dir = "../final_plots/closed/apo/"

os.makedirs(results_dir, exist_ok=True)
os.makedirs(plots_dir, exist_ok=True)
# ------------------------------------------


# =============================================================================
# 2. RMSD CALCULATION
# =============================================================================
# define residues for fitting and residues for rmsf calculation
fit_sel = 'name CA and (resnum 612:774 or resnum 1092:1242)'
rmsd_sel = 'resnum 1:1261 and name CA'

dataframes = []

for i, rep in enumerate(replicas, start=1):
    print(f"Analysing replica {i}: {rep}...")
    top_path = os.path.join(data_dir, topology)
    traj_path = os.path.join(data_dir, rep, traj_name)
    u = mda.Universe(top_path, traj_path)

    # align
    align.AlignTraj(u, u, select=fit_sel, ref_frame=0, in_memory=True).run()
    
    # calculate RMSD
    R = rms.RMSD(u, u, select=rmsd_sel, ref_frame=0, superposition=False)
    R.run()

    # create dataframe
    col_name = f'RMSD_md{i} (Å)'
    df = pd.DataFrame(R.results.rmsd, columns=['Frame', 'Time (ps)', col_name])
    
    if i > 1:
        df = df.drop(columns=['Time (ps)'])
        
    dataframes.append(df)

# =============================================================================
# 3. MERGE DATAFRAMES
# =============================================================================
# Merge all dataframes on column 'Frame'
merged_df = reduce(lambda left, right: pd.merge(left, right, on='Frame'), dataframes)
merged_df = merged_df.drop(columns=['Frame'])

# # 1st column: Time, in ns
cols = merged_df.columns.tolist()
cols.insert(0, cols.pop(cols.index("Time (ps)")))
merged_df = merged_df[cols]
merged_df['Time (ps)'] = merged_df['Time (ps)'] / 1000
merged_df = merged_df.rename(columns={"Time (ps)": "Time (ns)"})

# round RMSD
rmsd_cols = [f'RMSD_md{i} (Å)' for i in range(1, len(replicas) + 1)]
merged_df[rmsd_cols] = merged_df[rmsd_cols].round(2)

# save .csv
csv_path = os.path.join(results_dir, "rmsd-ca-protein.csv")
merged_df.to_csv(csv_path, sep="\t", index=False)
print(f"Data saved in: {csv_path}")

# =============================================================================
# 4. PLOTTING
# =============================================================================
LABEL_SIZE = 14
TITLE_SIZE = 16
LEGEND_SIZE = 12
TICK_SIZE = 12

fig, (ax_kde, ax_time) = plt.subplots(
    1, 2, figsize=(12, 5), sharey=True, gridspec_kw={'width_ratios': [1, 3]}
)

colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
labels = [f'md{i}' for i in range(1, len(replicas) + 1)]

# left: KDE plot
for i, col in enumerate(rmsd_cols):
    sns.kdeplot(data=merged_df, y=col, ax=ax_kde, fill=True, alpha=0.5, color=colors[i], label=labels[i])

ax_kde.invert_xaxis()
ax_kde.set_xlabel('')
ax_kde.set_xticks([])
ax_kde.set_ylabel('RMSD (Å)', fontsize=LABEL_SIZE)
ax_kde.legend(loc='upper left', fontsize=LEGEND_SIZE)

# right: RMSD vs time
for i, col in enumerate(rmsd_cols):
    ax_time.plot(merged_df['Time (ns)'], merged_df[col], label=labels[i], color=colors[i])

ax_time.set_xlabel('Time (ns)', fontsize=LABEL_SIZE)
ax_time.set_title(f'Cα RMSD - {system_title}', fontsize=TITLE_SIZE)
ax_time.set_ylabel('')

# axes limits and ticks size
y_max = 8.0 #common for all systems
ax_kde.set_ylim(0, y_max)
ax_kde.tick_params(axis='y', labelsize=TICK_SIZE)
ax_time.tick_params(axis='x', labelsize=TICK_SIZE)

# save .png
plt.tight_layout()
plot_path = os.path.join(plots_dir, "rmsd-ca-protein.png")
plt.savefig(plot_path, dpi=600)
print(f"Chart saved in: {plot_path}")
