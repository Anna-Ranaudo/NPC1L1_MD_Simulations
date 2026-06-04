import MDAnalysis as mda
from MDAnalysis.analysis import align
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
chl_sel = 'resnum 1262 and not type h*'
wat_sel = f'resname WAT and around 4 of ({chl_sel})'

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
    # 2. CHECK WATER MOLECULES WITHIN 4 A OF CHL, per replica
    # =============================================================================
    dataframes = []

    for i, rep in enumerate(replicas, start=1):
        print(f"  Analysing replica {i}: {rep}...")
        traj_path = os.path.join(data_dir, rep, traj_name)
        
        if not os.path.exists(topology) or not os.path.exists(traj_path):
            print(f"  [WARNING] Path not found, skipping: {traj_path}")
            continue
            
        u = mda.Universe(topology, traj_path)

        # align the trajectory on the protein segments used for fitting
        align.AlignTraj(u, u, select=fit_sel, ref_frame=0, in_memory=True).run()

        counts = []
        times = []
        for ts in u.trajectory:
            water_sel = u.select_atoms(wat_sel)
            counts.append(len(water_sel.residues))
            times.append(u.trajectory.time / 1000.0)

        col_name = f'Waters_md{i}'
        df = pd.DataFrame({
            'Frame': range(len(counts)),
            'Time (ns)': times,
            col_name: counts,
        })

        if len(dataframes) > 0:
            df = df.drop(columns=['Time (ns)'])
        dataframes.append(df)

    if not dataframes:
        print(f"  [ERROR] No data found for system {system_name}. Skipping...")
        continue

    # =============================================================================
    # 3. MERGE DATAFRAMES
    # =============================================================================
    merged_df = reduce(lambda left, right: pd.merge(left, right, on='Frame'), dataframes)
    merged_df = merged_df.drop(columns=['Frame'])

    cols = merged_df.columns.tolist()
    cols.insert(0, cols.pop(cols.index('Time (ns)')))
    merged_df = merged_df[cols]

    # save .csv
    csv_path = os.path.join(results_dir, "water-counts-chol-4A.csv")
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

    count_cols = [c for c in merged_df.columns if c.startswith('Waters_md')]
    labels = [c.replace('Waters_', '') for c in count_cols]

    for i, col in enumerate(count_cols):
        sns.kdeplot(data=merged_df, y=col, ax=ax_kde, fill=True, alpha=0.5, label=labels[i])

    ax_kde.invert_xaxis()
    ax_kde.set_xlabel('')
    ax_kde.set_xticks([])
    ax_kde.set_ylabel('Number of waters within 4 Å', fontsize=LABEL_SIZE)
    ax_kde.legend(loc='upper left', fontsize=LEGEND_SIZE)

    for i, col in enumerate(count_cols):
        ax_time.plot(merged_df['Time (ns)'], merged_df[col], label=labels[i])

    ax_time.set_xlabel('Time (ns)', fontsize=LABEL_SIZE)
    ax_time.set_title(f'Water count near cholesterol - {system_name}', fontsize=TITLE_SIZE)
    ax_time.set_ylabel('Number of waters', fontsize=LABEL_SIZE)
    ax_time.tick_params(axis='x', labelsize=TICK_SIZE)
    ax_kde.tick_params(axis='y', labelsize=TICK_SIZE)

    plt.tight_layout()
    plot_path = os.path.join(plots_dir, "water-counts-chol-4A.png")
    plt.savefig(plot_path, dpi=600)
    plt.close(fig)
    print(f"  Chart saved in: {plot_path}")

