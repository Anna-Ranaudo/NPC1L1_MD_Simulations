import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from functools import reduce

# =============================================================================
# 1. DEFINE THE SYSTEMS
# =============================================================================
systems = [
    {
        "system_name": "Open conformation",
        "traj_dir": "/mnt/h/Il mio Drive/LAVORO_MD_NPC1L1_nov25/MD_6V3F_6V3H_500ns_sept25/6v3f/",
        "results_dir": "../final_data/open/bound/",
        "plots_dir": "../final_plots/open/bound/",
    },
    {
        "system_name": "Closed conformation",
        "traj_dir": "/mnt/h/Il mio Drive/LAVORO_MD_NPC1L1_nov25/MD_6V3F_6V3H_500ns_sept25/6v3h/",
        "results_dir": "../final_data/closed/bound/",
        "plots_dir": "../final_plots/closed/bound/",
    },
]

# common for both systems:
replicas = ["run-md1", "run-md2", "run-md3", "run-md4", "run-md5"]
fit_sel = 'name CA and (resnum 612:782 or resnum 1092:1242)'
rmsd_sel = 'resnum 1262 and not type h*'

for system in systems:
    system_name = system["system_name"]
    traj_dir = system["traj_dir"]
    results_dir = system["results_dir"]
    plots_dir = system["plots_dir"]

    print(f"\n>>> Analyzing System: {system_name} <<<")
    
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)

    data_dir = os.path.join(traj_dir, "data/col/")
    topology = os.path.join(data_dir, "prot-lig.prmtop")
    traj_name = "07-08-prot-lig-pbc.nc"

    # =============================================================================
    # 2. RMSD CALCULATION
    # =============================================================================
    dataframes = []

    for i, rep in enumerate(replicas, start=1):
        print(f"  Analysing replica {i}: {rep}...")
        traj_path = os.path.join(data_dir, rep, traj_name)
        
        if not os.path.exists(topology) or not os.path.exists(traj_path):
            print(f"  [WARNING] Path not found, skipping: {traj_path}")
            continue
            
        u = mda.Universe(topology, traj_path)

        # align
        align.AlignTraj(u, u, select=fit_sel, ref_frame=0, in_memory=True).run()
        
        # calculate RMSD
        R = rms.RMSD(u, u, select=rmsd_sel, ref_frame=0, superposition=False)
        R.run()

        # create dataframe
        col_name = f'RMSD_md{i} (Å)'
        df = pd.DataFrame(R.results.rmsd, columns=['Frame', 'Time (ps)', col_name])
        
        if len(dataframes) > 0:
            df = df.drop(columns=['Time (ps)'])
            
        dataframes.append(df)

    if not dataframes:
        print(f"  [ERROR] No data found for system {system_name}. Skipping...")
        continue

    # =============================================================================
    # 3. MERGE DATAFRAMES
    # =============================================================================
    # Merge all dataframes on column 'Frame'
    merged_df = reduce(lambda left, right: pd.merge(left, right, on='Frame'), dataframes)
    merged_df = merged_df.drop(columns=['Frame'])

    # 1st column: Time, in ns
    cols = merged_df.columns.tolist()
    cols.insert(0, cols.pop(cols.index("Time (ps)")))
    merged_df = merged_df[cols]
    merged_df['Time (ps)'] = merged_df['Time (ps)'] / 1000
    merged_df = merged_df.rename(columns={"Time (ps)": "Time (ns)"})

    # round RMSD
    rmsd_cols = [c for c in merged_df.columns if c.startswith('RMSD_')]
    merged_df[rmsd_cols] = merged_df[rmsd_cols].round(2)

    # save .csv
    csv_path = os.path.join(results_dir, "rmsd-chol.csv")
    merged_df.to_csv(csv_path, sep="\t", index=False)
    print(f"  Data saved in: {csv_path}")

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
    labels = [c.split(' ')[0].replace('RMSD_', '') for c in rmsd_cols]

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
    ax_time.set_title(f'Cholesterol RMSD - {system_name}', fontsize=TITLE_SIZE)
    ax_time.set_ylabel('')

    # axes limits and ticks size
    y_max = 2.0 #common for both systems
    ax_kde.set_ylim(0, y_max)
    ax_kde.tick_params(axis='y', labelsize=TICK_SIZE)
    ax_time.tick_params(axis='x', labelsize=TICK_SIZE)

    # save .png
    plt.tight_layout()
    plot_path = os.path.join(plots_dir, "rmsd-chol.png")
    plt.savefig(plot_path, dpi=600)
    plt.close(fig) # Close plot to free memory
    print(f"  Chart saved in: {plot_path}")

