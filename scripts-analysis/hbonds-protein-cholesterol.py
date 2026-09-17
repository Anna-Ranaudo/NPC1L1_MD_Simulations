import argparse
import os
import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis


# replicas and systems following the same conventions as other single-replica scripts
replicas = ["run-md1", "run-md2", "run-md3", "run-md4", "run-md5"]

systems = [
    {
        "system_id": "open_bound",
        "traj_base": "/mnt/h/Il mio Drive/LAVORO_MD_NPC1L1_nov25/MD_6V3F_6V3H_500ns_sept25/6v3f/data/col/",
        "topology": "prot-lig.prmtop",
        "traj_name": "07-08-prot-lig-pbc.nc",
        "results_dir": "../final_data/open/bound/prot_chol_hbond/",
    },
    {
        "system_id": "closed_bound",
        "traj_base": "/mnt/h/Il mio Drive/LAVORO_MD_NPC1L1_nov25/MD_6V3F_6V3H_500ns_sept25/6v3h/data/col/",
        "topology": "prot-lig.prmtop",
        "traj_name": "07-08-prot-lig-pbc.nc",
        "results_dir": "../final_data/closed/bound/prot_chol_hbond/",
    },
]

# selection for hydrogen bond analysis (kept as in original script)
hb_selection = ["resnum 1262", "protein"]
RESIDUE_OFFSET = 21


def atom_label(atom):
    resid = atom.resid if atom.resname == "CHL" else atom.resid + RESIDUE_OFFSET
    return f"{atom.resname}{resid}-{atom.name}"

# threshold to keep labels (fraction of frames)
threshold = 0.10

# Make relative result paths robust: resolve relative paths against this script folder
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser()
parser.add_argument(
    "--system",
    choices=["open", "closed", "all"],
    default="all",
    help="System to process (default: all).",
)
args = parser.parse_args()

selected_systems = systems
if args.system != "all":
    selected_systems = [
        system
        for system in systems
        if system["system_id"].startswith(f"{args.system}_")
    ]

for system in selected_systems:
    system_id = system["system_id"]
    traj_base = system["traj_base"]
    topology = system["topology"]
    traj_name = system["traj_name"]
    results_dir = system["results_dir"]

    # if paths are relative, interpret them relative to the script directory
    if not os.path.isabs(results_dir):
        results_dir = os.path.normpath(os.path.join(SCRIPT_DIR, results_dir))

    os.makedirs(results_dir, exist_ok=True)
    # --- AGGREGATED LOGIC: collect H-bond labels across all replicas ---
    replica_data = {}
    all_labels = set()

    for rep in replicas:
        traj_path = os.path.join(traj_base, rep, traj_name)
        top_path = os.path.join(traj_base, topology)

        print(f"Processing {system_id} replica {rep}...")
        if not os.path.exists(top_path) or not os.path.exists(traj_path):
            print(f"  [WARNING] Path not found, skipping: {traj_path}")
            continue

        u = mda.Universe(top_path, traj_path)
        hba = HydrogenBondAnalysis(universe=u, between=hb_selection)
        hba.run()

        hb_array = hba.results.hbonds
        labels = []
        for row in hb_array:
            don_idx, acc_idx = int(row[1]), int(row[3])
            don_atom = u.atoms[don_idx]
            acc_atom = u.atoms[acc_idx]
            label = f"{atom_label(don_atom)}...{atom_label(acc_atom)}"
            labels.append(label)

        n_frames = int(hb_array[:, 0].max()) + 1 if len(hb_array) > 0 else 0
        times = getattr(hba, 'times', None)
        if times is None and n_frames > 0:
            times = np.arange(n_frames)
        elif times is not None:
            times = times / 1000

        replica_data[rep] = {
            'hbonds': hba,
            'array': hb_array,
            'labels': labels,
            'n_frames': n_frames,
            'times': times,
        }

        all_labels.update(labels)

    if not replica_data:
        print(f"  [ERROR] No replicas processed for {system_id}. Skipping system.")
        continue

    # Determine which labels pass the threshold in at least one replica
    filtered_labels = []
    for label in all_labels:
        keep = False
        for rep, data in replica_data.items():
            n_frames = data['n_frames'] if data['n_frames'] > 0 else 1
            count = sum(1 for l in data['labels'] if l == label)
            if (count / n_frames) >= threshold:
                keep = True
                break
        if keep:
            filtered_labels.append(label)

    filtered_labels = sorted(filtered_labels)
    if not filtered_labels:
        print(f"  [INFO] No H-bond labels pass threshold={threshold} for system {system_id}.")
        continue

    label_to_col = {label: idx for idx, label in enumerate(filtered_labels)}

    # Generate binary matrices and outputs per replica using the common filtered labels
    for rep, data in replica_data.items():
        n_frames = data['n_frames']
        if n_frames == 0:
            print(f"  [INFO] Replica {rep} has 0 frames, skipping.")
            continue

        matrix = np.zeros((n_frames, len(filtered_labels)), dtype=int)
        for row, label in zip(data['array'], data['labels']):
            if label in label_to_col:
                frame = int(row[0])
                col = label_to_col[label]
                matrix[frame, col] = 1

        time_ns = data['times'] if data['times'] is not None else np.arange(n_frames)

        df_time = pd.DataFrame(time_ns, columns=["Time (ns)"])
        df_matrix = pd.DataFrame(matrix, columns=filtered_labels)
        df = df_time.join(df_matrix)

        replica_id = rep.replace("run-", "")
        out_csv = os.path.join(results_dir, f"{replica_id}_hbonds_presence_final.csv")
        df.to_csv(out_csv, index=False)
        print(f"  Saved CSV: {out_csv}")


# Run the plotting script after both systems have been processed.
if __name__ == "__main__":
    import subprocess
    import sys

    plot_script = os.path.join(SCRIPT_DIR, "plot-hbonds-protein-cholesterol.py")
    plot_configs = [(args.system, "bound")] if args.system != "all" else [
        ("open", "bound"),
        ("closed", "bound"),
    ]

    for state, binding in plot_configs:
        input_dir = os.path.join(
            SCRIPT_DIR, "..", "final_data", state, binding, "prot_chol_hbond"
        )
        output_dir = os.path.join(
            SCRIPT_DIR, "..", "final_plots", state, binding
        )
        output_name = f"protein-cholesterol-hbonds_5replicas_{state}_{binding}.png"
        title = f"Protein-cholesterol hydrogen bonds | {state} {binding}"
        subprocess.run(
            [
                sys.executable,
                plot_script,
                "--input-dir",
                input_dir,
                "--output-dir",
                output_dir,
                "--output-name",
                output_name,
                "--title",
                title,
            ],
            check=True,
        )
