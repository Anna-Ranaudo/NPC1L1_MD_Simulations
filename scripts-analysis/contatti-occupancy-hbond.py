# -*- coding: utf-8 -*-
"""
script for inter-domain HYDROGEN BOND occupancy analysis
from Molecular Dynamics simulations using MDAnalysis.

It performs the following steps:

"""

# --- ANALISI DELLO SCRIPT (IN ITALIANO) ---
"""
Questo script Python utilizza la libreria MDAnalysis per eseguire un'analisi
approfondita dei legami a idrogeno (H-Bond) tra diversi domini di una proteina
in simulazioni di dinamica molecolare (MD).

OBIETTIVO PRINCIPALE:
Identificare quali coppie di residui formano legami a idrogeno stabili
tra domini specifici (es. NTD-MLD, NTD-CTD) e confrontare l'occupancy
(la "frequenza") di questi legami in diversi sistemi (es. 'Open Apo',
'Closed Bound').

COME FUNZIONA:
1.  CONFIGURAZIONE: L'utente definisce i percorsi dei file, i criteri
    per i legami a idrogeno, i sistemi da analizzare (4 in questo caso)
    e i domini della proteina.
2.  CALCOLO PARALLELO: Per ogni sistema, lo script analizza le traiettorie
    MD. Utilizza il modulo 'multiprocessing' per parallelizzare il lavoro,
    dividendolo tra i core della CPU per velocizzare l'analisi.
3.  ANALISI H-BOND: Per ogni "chunk" di frame, calcola i legami a idrogeno
    tra i domini specificati usando 'MDAnalysis.analysis.hydrogenbonds'.
4.  CALCOLO OCCUPANCY: Conta quante volte ogni coppia di residui forma un
    legame e lo divide per il numero totale di frame analizzati,
    ottenendo un'occupancy (es. 0.75 = presente nel 75% del tempo).
5.  FILTRAGGIO: Conserva solo i legami "significativi" (sopra una
    soglia, es. 5%).
6.  OUTPUT:
    -   Crea un file .csv per ogni sistema con i legami significativi.
    -   Genera un file .txt di riepilogo che confronta i legami
        comuni a tutti i sistemi.
    -   Crea delle "heatmap" (mappe di calore) in .png che mostrano
        visivamente il confronto dell'occupancy di tutti i legami
        tra i quattro sistemi.
"""
# --- FINE ANALISI ---


import os
import glob
import pandas as pd
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from collections import defaultdict
import multiprocessing

# Ignore common MDAnalysis warnings for a cleaner output
warnings.filterwarnings('ignore', category=UserWarning, message='Failed to guess the mass for')
warnings.filterwarnings('ignore', category=UserWarning, message='Found no information for attr')

# --- CONFIGURATION PARAMETERS (USER TO MODIFY) ---

# 1. File Paths
BASE_DIR = './'
TRAJECTORY_PATTERN = 'run-md*/07-08-prot*.nc'

# 2. Analysis Parameters
# H-Bond geometric criteria
HBOND_DISTANCE_CUTOFF = 3.0
HBOND_ANGLE_CUTOFF = 150.0
# valori standard riportati in https://docs.mdanalysis.org/stable/documentation_pages/analysis/hydrogenbonds.html#MDAnalysis.analysis.hydrogenbonds.hbond_analysis.HydrogenBondAnalysis

RESIDUE_OFFSET = 21  # in .nc i residui partono da 1, nei miei file da 22
OCCUPANCY_THRESHOLD = 0.05 
# Analizza 1 frame ogni N. Imposta a 1 per un'analisi completa.
FRAME_STEP = 1 

# 3. System and Domain Definitions
SYSTEMS = {
    'Open Apo': {
        'trajectories': '6v3f/data/no-col',
        'topology': '6v3f/data/no-col/prot.prmtop'
    },
    'Open Cholesterol': {
        'trajectories': '6v3f/data/col',
        'topology': '6v3f/data/col/prot-lig.prmtop'
    },
    'Closed Cholesterol': {
        'trajectories': '6v3h/data/col',
        'topology': '6v3h/data/col/prot-lig.prmtop'
    },
    'Closed Apo': {
        'trajectories': '6v3h/data/no-col',
        'topology': '6v3h/data/no-col/prot.prmtop'
    },
#    'Closed Sitosterol': {
#        'trajectories': '../Tesi-LT-Garavaglia-Simone/Risultati_MD/data/sit',
#        'topology': '../Tesi-LT-Garavaglia-Simone/Risultati_MD/data/sit/prot-lig.prmtop'
#    }
}


DOMAINS_RESIDS = {
    'NTD': (1, 242),
    'MLD': (358, 611),
    'CTD': (861, 1084)
}
DOMAIN_PAIRS = [
    ('NTD', 'MLD'),
    ('NTD', 'CTD'),
    ('MLD', 'CTD')
]

# 4. Output Parameters
OUTPUT_DIR = './analysis_hbond_occupancy_thre05_step1'
os.makedirs(OUTPUT_DIR, exist_ok=True)


# --- PARALLEL WORKER AND ANALYSIS FUNCTIONS ---

def worker_function(args):
    """
    Function executed by each parallel process. It analyzes a chunk of frames
    for HYDROGEN BONDS.
    """
    topo_path, traj_files, frame_chunk, sel1, sel2, hbond_params = args
    
    u = mda.Universe(topo_path, traj_files)
    
    hba = HydrogenBondAnalysis(
        universe=u,
        between=[sel1, sel2],
        d_a_cutoff=hbond_params['distance'],
        d_h_a_angle_cutoff=hbond_params['angle'],
    )
    
    contact_counts = defaultdict(int)
    
    # Itera solo sui frame assegnati in questo chunk
    for ts in u.trajectory[frame_chunk]:
        hba.results.hbonds = []
        hba.run(start=ts.frame, stop=ts.frame + 1, step=1)
        
        contacting_residue_pairs = set()
        for hb in hba.results.hbonds:
            donor_res = u.atoms[int(hb[1])].residue
            acceptor_res = u.atoms[int(hb[3])].residue
            
            res_tuple1 = (donor_res.resname, donor_res.resid)
            res_tuple2 = (acceptor_res.resname, acceptor_res.resid)
            
            pair = tuple(sorted((res_tuple1, res_tuple2), key=lambda x: x[1]))
            contacting_residue_pairs.add(pair)
            
        for pair in contacting_residue_pairs:
            contact_counts[pair] += 1
            
    return contact_counts

def calculate_hbond_occupancy_parallel(topo_path, traj_files, sel1, sel2, hbond_params, offset, frame_step):
    """
    Orchestrates the parallel calculation of H-Bond occupancies.
    """
    u_for_info = mda.Universe(topo_path, traj_files)
    n_frames_total = len(u_for_info.trajectory)
    del u_for_info
    
    #  Seleziona i frame da analizzare ***
    frames_to_analyze = list(range(0, n_frames_total, frame_step))
    n_frames_sampled = len(frames_to_analyze)

    if n_frames_sampled == 0:
        print(f"  WARNING: No frames to analyze with frame_step={frame_step}. n_frames_total={n_frames_total}.")
        return {}

    n_procs = min(multiprocessing.cpu_count(), n_frames_sampled)
    if n_procs == 0: n_procs = 1

    # Divide la *lista* di frame (non un range)
    frame_chunks = np.array_split(frames_to_analyze, n_procs)
    
    tasks = [(topo_path, traj_files, chunk, sel1, sel2, hbond_params)
             for chunk in frame_chunks if len(chunk) > 0]

    print(f"  Total frames: {n_frames_total}. Analyzing {n_frames_sampled} frames (1 every {frame_step}).")
    print(f"  Distributing H-Bond analysis over {len(tasks)} processes...")

    with multiprocessing.Pool(processes=n_procs) as pool:
        results = pool.map(worker_function, tasks)
        
    merged_counts = defaultdict(int)
    for count_dict in results:
        for pair, count in count_dict.items():
            merged_counts[pair] += count
            
    # L'occupancy è basata sul numero di frame *campionati*
    occupancies = {
        f"{res1[0]}{res1[1] + offset}-{res2[0]}{res2[1] + offset}": count / n_frames_sampled
        for (res1, res2), count in merged_counts.items()
    }
    
    return occupancies

def get_domain_selection_string(domain_name, resids):
    """
    Creates the selection string for MDAnalysis.
    *** Uses ORIGINAL resids (no offset) ***
    """
    start_res, end_res = resids
    return f"protein and resid {start_res}-{end_res}"


# --- MAIN SCRIPT ---

def main():
    print("Starting PARALLEL H-Bond occupancy analysis...")

    hbond_params = {
        'distance': HBOND_DISTANCE_CUTOFF,
        'angle': HBOND_ANGLE_CUTOFF
    }
    
    all_occupancy_data = defaultdict(dict)

    for system_label, system_info in SYSTEMS.items():
        print(f"\n--- Analyzing system: {system_label} ---")
        
        trajectories_folder_path = os.path.join(BASE_DIR, system_info['trajectories'])
        if not os.path.isdir(trajectories_folder_path):
            print(f"WARNING: Trajectory directory not found: {trajectories_folder_path}. Skipping.")
            continue
            
        topo_path = os.path.join(BASE_DIR, system_info['topology'])
        traj_files = sorted(glob.glob(os.path.join(trajectories_folder_path, TRAJECTORY_PATTERN)))
        
        if not traj_files or not os.path.exists(topo_path):
            print(f"WARNING: Files missing for {system_label}. Trajectories found: {len(traj_files)}, Topology exists: {os.path.exists(topo_path)} at '{topo_path}'. Skipping.")
            continue

        # Definiamo 'safe_system_label' qui
        safe_system_label = system_label.replace(' ', '_')
        
        # Rimuoviamo questa lista, non accumuliamo più i risultati qui
        # system_results_df = [] # <-- RIMOSSO

        for domain1_name, domain2_name in DOMAIN_PAIRS:
            pair_label = f'{domain1_name}-{domain2_name}'
            print(f"  Calculating H-Bond occupancy for: {pair_label}...")

            sel1 = get_domain_selection_string(domain1_name, DOMAINS_RESIDS[domain1_name])
            sel2 = get_domain_selection_string(domain2_name, DOMAINS_RESIDS[domain2_name])
            
            occupancies = calculate_hbond_occupancy_parallel(
                topo_path, traj_files, sel1, sel2, hbond_params, RESIDUE_OFFSET, FRAME_STEP
            )
            
            significant_contacts = {
                pair: occ for pair, occ in occupancies.items() if occ >= OCCUPANCY_THRESHOLD
            }
            all_occupancy_data[system_label][pair_label] = significant_contacts
            
            # Creiamo una lista *specifica* per questa coppia di domini
            pair_results_list = []
            for pair, occ in significant_contacts.items():
                pair_results_list.append({'Domain Pair': pair_label, 'H-Bond': pair, 'Occupancy': occ})

            # Salviamo il file CSV per questa specifica coppia, se ci sono dati
            if pair_results_list:
                df_pair = pd.DataFrame(pair_results_list).sort_values(by='Occupancy', ascending=False)
                
                # Creiamo il nuovo nome file come richiesto
                csv_filename = f'hbond_occupancy_{safe_system_label}_{pair_label}.csv'
                csv_filepath = os.path.join(OUTPUT_DIR, csv_filename)
                
                # Salviamo il file
                df_pair.to_csv(csv_filepath, index=False)
                print(f"  Significant H-Bonds for {pair_label} saved to: {csv_filepath}")

        # Rimuoviamo il vecchio blocco di salvataggio che raggruppava tutto
        # if system_results_df:
        #     df = pd.DataFrame(system_results_df).sort_values(by=['Domain Pair', 'Occupancy'], ascending=[True, False])
        #     safe_system_label = system_label.replace(' ', '_')
        #     csv_filepath = os.path.join(OUTPUT_DIR, f'hbond_occupancy_{safe_system_label}.csv')
        #     df.to_csv(csv_filepath, index=False)
        #     print(f"  Significant H-Bonds saved to: {csv_filepath}")

    # --- COMPARISON AND PLOTTING ---
    if not all_occupancy_data:
        print("\nNo data was analyzed. Cannot perform comparison.")
        return

    print("\n--- Comparing H-Bonds across systems ---")
    
    summary_file_path = os.path.join(OUTPUT_DIR, 'common_hbonds_summary.txt')
    with open(summary_file_path, 'w') as f:
        f.write(f"Summary of Common H-Bonds (Occupancy >= {OCCUPANCY_THRESHOLD:.0%})\n")
        f.write("="*60 + "\n")

        for domain1, domain2 in DOMAIN_PAIRS:
            pair_label = f'{domain1}-{domain2}'
            f.write(f"\nAnalysis for Domain Pair: {pair_label}\n")
            
            contact_sets = [set(all_occupancy_data[sys_label].get(pair_label, {}).keys()) for sys_label in SYSTEMS.keys()]
            
            if not any(contact_sets):
                f.write("  No significant H-bonds found for this pair in any system.\n")
                continue

            common_contacts = set.intersection(*contact_sets)
            
            if common_contacts:
                f.write(f"\n  Found {len(common_contacts)} H-bonds common to ALL 4 systems:\n")
                for contact in sorted(list(common_contacts)):
                    f.write(f"    - {contact}\n")
            else:
                f.write("\n  No H-bonds were found to be common across ALL 4 systems.\n")

            all_contacts_for_pair = sorted(list(set.union(*contact_sets)))
            
            if not all_contacts_for_pair:
                continue

            heatmap_data = []
            for contact in all_contacts_for_pair:
                row = {'H-Bond': contact}
                for sys_label in SYSTEMS.keys():
                    occupancy = all_occupancy_data[sys_label].get(pair_label, {}).get(contact, 0)
                    row[sys_label] = occupancy
                heatmap_data.append(row)
            
            heatmap_df = pd.DataFrame(heatmap_data).set_index('H-Bond')
            
            # *** NUOVO: Limite di altezza per stabilità ***
            height = min(100, max(8, len(heatmap_df) * 0.3)) # Max 100 pollici
            
            plt.figure(figsize=(10, height))
            sns.heatmap(heatmap_df, cmap="viridis", annot=True, fmt=".2f", linewidths=.5,
                        vmin=0, vmax=1) # Fissiamo la scala colori 0-1)
            plt.title(f'H-Bond Occupancy Heatmap for {pair_label}', fontsize=16)
            plt.xlabel('System', fontsize=12)
            plt.ylabel('Residue Pair', fontsize=12)
            plt.xticks(rotation=0, ha='center')
            plt.tight_layout()
            
            heatmap_path = os.path.join(OUTPUT_DIR, f'heatmap_hbond_occupancy_{domain1}_{domain2}.png')
            plt.savefig(heatmap_path, dpi=300)
            print(f"  Heatmap saved to: {heatmap_path}")
            plt.close()
    
    print(f"\nComparison summary saved to: {summary_file_path}")
    print("\nAdvanced H-Bond analysis successfully completed!")

if __name__ == '__main__':
    multiprocessing.freeze_support()
    main()
