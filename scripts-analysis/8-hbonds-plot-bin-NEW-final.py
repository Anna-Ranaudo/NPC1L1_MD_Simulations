import os
import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis


#######  6V3F col ##########
data = "data/col/"
results = "results/col/confronti-r-1-5"

replicas = {
    "md1": os.path.join(data, "run-md1/07-08-prot-lig-pbc.nc"),
    "md2": os.path.join(data, "run-md2/07-08-prot-lig-pbc.nc"),
    "md3": os.path.join(data, "run-md3/07-08-prot-lig-pbc.nc"),
    "md4": os.path.join(data, "run-md4/07-08-prot-lig-pbc.nc"),
    "md5": os.path.join(data, "run-md5/07-08-prot-lig-pbc.nc"),
}

topology = os.path.join(data, "prot-lig.prmtop")

# Step 1: Collect all unique H-bond labels across replicas
all_labels = set()
replica_data = {}

for rep, traj in replicas.items():
    u = mda.Universe(topology, traj)
    hbonds = HydrogenBondAnalysis(universe=u, between=["resnum 1262", "protein"])
    hbonds.run()

    hb_array = hbonds.results.hbonds
    labels = []
    for row in hb_array:
        don_idx, acc_idx = int(row[1]), int(row[3])
        don_atom = u.atoms[don_idx]
        acc_atom = u.atoms[acc_idx]
        label = f"{don_atom.resname}{don_atom.resid}-{don_atom.name}...{acc_atom.resname}{acc_atom.resid}-{acc_atom.name}"
        labels.append(label)

    all_labels.update(labels)
    replica_data[rep] = {"hbonds": hbonds, "array": hb_array, "labels": labels}

# Step 2: Determine which H-bonds pass the 10% threshold in at least one replica
label_counts = {label: [] for label in all_labels}

for rep, data in replica_data.items():
    hb_array = data["array"]
    hbonds = data["hbonds"]
    n_frames = int(hb_array[:, 0].max()) + 1
    labels = data["labels"]

    # Count occurrences per label
    for row, label in zip(hb_array, labels):
        label_counts[label].append(1)
    # Fill with zeros for missing occurrences
    for label in all_labels:
        if label not in labels:
            label_counts[label].append(0)

# Fraction per label across frames
filtered_labels = []
for label in all_labels:
    keep = False
    for rep, data in replica_data.items():
        hb_array = data["array"]
        n_frames = int(hb_array[:, 0].max()) + 1 if len(hb_array) > 0 else 1
        count = sum(1 for l in data["labels"] if l == label)
        if count / n_frames >= 0.10:  # at least 10%
            keep = True
            break
    if keep:
        filtered_labels.append(label)

filtered_labels = sorted(filtered_labels)
label_to_col = {label: i for i, label in enumerate(filtered_labels)}

# Step 3: Generate binary matrix for each replica (only filtered labels, fixed order)
for rep, data in replica_data.items():
    hbonds = data["hbonds"]
    hb_array = data["array"]
    labels = data["labels"]

    n_frames = int(hb_array[:, 0].max()) + 1
    matrix = np.zeros((n_frames, len(filtered_labels)), dtype=int)

    for row, label in zip(hb_array, labels):
        if label in label_to_col:
            frame = int(row[0])
            col = label_to_col[label]
            matrix[frame, col] = 1

    time_ns = hbonds.times / 1000

    df_time = pd.DataFrame(time_ns, columns=["Time (ns)"])
    df_matrix = pd.DataFrame(matrix, columns=filtered_labels)
    df = df_time.join(df_matrix)

    out_csv = os.path.join(results, f"{rep}_hbonds_presence_final.csv")
    df.to_csv(out_csv, index=False)

    # Plot
    plt.figure(figsize=(12, 6))
    cmap = plt.cm.get_cmap('Greys', 2)
    im = plt.imshow(matrix.T, aspect='auto', cmap=cmap, interpolation='nearest')

    # X-axis: 0–500 ns, ticks every 50
    plt.xticks(
        ticks=np.linspace(0, n_frames - 1, 11),
        labels=np.arange(0, 501, 50)
    )
    plt.xlabel("Time (ns)")
    plt.yticks(
        ticks=np.arange(len(filtered_labels)),
        labels=filtered_labels,
        rotation=45,
        va='center'
    )
    plt.title(f"Protein-cholesterol hydrogen bonds - 6V3F-{rep}")
    #plt.title(f"{rep}")
    #cbar = plt.colorbar(im, ticks=[0, 1])
    #cbar.ax.set_yticklabels(['Absent', 'Present'])
    #cbar.set_label("H-bond")
    plt.tight_layout()
    plt.savefig(os.path.join(results, f"{rep}_hbonds_heatmap_final.png"), dpi=300)
    plt.close()

