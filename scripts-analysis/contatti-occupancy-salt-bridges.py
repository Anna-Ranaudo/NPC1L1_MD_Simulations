# -*- coding: utf-8 -*-
"""
script for inter-domain SALT BRIDGES occupancy analysis
from Molecular Dynamics simulations using MDAnalysis.

It performs the following steps:
1.  Configuration: User defines paths, salt bridge distance cutoff, 
    systems, and domains.
2.  Parallel Calculation: For each system, it parallelizes the analysis
    across CPU cores.
3.  Salt Bridge Analysis: For each frame chunk, it calculates distances 
    between charged groups (Basic-Acidic and Acidic-Basic) of the 
    specified domains. It uses MDAnalysis.analysis.distances.
4.  Occupancy Calculation: Counts how many frames each residue pair is 
    within the distance cutoff.
5.  Filtering: Keeps only significant interactions (above threshold).
6.  Output: Saves .csv files for each domain pair, a summary .txt, 
    and comparative heatmaps.
"""

# --- ANALISI DELLO SCRIPT (IN ITALIANO) ---
"""
Questo script Python utilizza la libreria MDAnalysis per eseguire un'analisi
approfondita dei salt bridges (ponti salini) tra diversi domini di una proteina
in simulazioni di dinamica molecolare (MD).

OBIETTIVO PRINCIPALE:
Identificare quali coppie di residui formano ponti salini stabili
tra domini specifici (es. NTD-MLD, NTD-CTD) e confrontare l'occupancy
(la "frequenza") di questi legami in diversi sistemi.

COME FUNZIONA:
1.  CONFIGURAZIONE: L'utente definisce i percorsi, il cutoff di distanza
    per i ponti salini (es. 4.5 Å), i sistemi da analizzare (5 in questo caso)
    e i domini della proteina.
2.  CALCOLO PARALLELO: Per ogni sistema, lo script analizza le traiettorie
    MD. Utilizza il modulo 'multiprocessing' per parallelizzare il lavoro.
3.  ANALISI PONTI SALINI: Per ogni "chunk" di frame, calcola le distanze
    tra i gruppi carichi dei domini. Esempio per NTD-MLD:
    a) (ARG/LYS di NTD) <-> (ASP/GLU di MLD)
    b) (ASP/GLU di NTD) <-> (ARG/LYS di MLD)
    Usa la libreria 'MDAnalysis.analysis.distances'.
4.  CALCOLO OCCUPANCY: Conta quante volte ogni coppia di residui forma un
    ponte (distanza < cutoff) e lo divide per il numero di frame analizzati.
5.  FILTRAGGIO: Conserva solo i legami "significativi" (sopra la
    soglia, es. 5%).
6.  OUTPUT:
    -   Crea file .csv separati per ogni coppia di domini (es. ..._NTD-MLD.csv).
    -   Genera un file .txt di riepilogo che confronta i legami.
    -   Crea delle "heatmap" (mappe di calore) in .png per il confronto.
"""
# --- FINE ANALISI ---


import os
import glob
import pandas as pd
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
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
# Salt bridges (SB) geometric criteria
SB_DISTANCE_CUTOFF = 4.5  # Distanza in Å tra i gruppi carichi

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
    'Closed Apo': {
        'trajectories': '6v3h/data/no-col',
        'topology': '6v3h/data/no-col/prot.prmtop'
    },
    'Closed Cholesterol': {
        'trajectories': '6v3h/data/col',
        'topology': '6v3h/data/col/prot-lig.prmtop'
    },
    'Closed Sitosterol': {
        'trajectories': '../Tesi-LT-Garavaglia-Simone/Risultati_MD/data/sit',
        'topology': '../Tesi-LT-Garavaglia-Simone/Risultati_MD/data/sit/prot-lig.prmtop'
    }
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


# Definiamo i residui che possono formare ponti salini
# Atomi carichi per i residui basici (cationici)
sel_basic = "(resname ARG LYS) and (name NH* NZ)"
# Atomi carichi per i residui acidi (anionici)
sel_acidic = "(resname GLU ASP) and (name OE* OD*)"

# 4. Output Parameters
OUTPUT_DIR = './analysis_salt_bridges_occupancy_thre05_step1'
os.makedirs(OUTPUT_DIR, exist_ok=True)


# --- PARALLEL WORKER AND ANALYSIS FUNCTIONS ---

def worker_function(args):
    """
    Function executed by each parallel process. It analyzes a chunk of frames
    for SALT BRIDGES based on distance.
    """
    # Ora riceviamo 4 selezioni (Basic1, Acidic2, Acidic1, Basic2)
    topo_path, traj_files, frame_chunk, sel_b1, sel_a2, sel_a1, sel_b2, sb_params = args
    
    u = mda.Universe(topo_path, traj_files)
    
    # Selezioniamo i gruppi di atomi *una sola volta*
    g_b1 = u.select_atoms(sel_b1) # Basic Domain 1
    g_a2 = u.select_atoms(sel_a2) # Acidic Domain 2
    g_a1 = u.select_atoms(sel_a1) # Acidic Domain 1
    g_b2 = u.select_atoms(sel_b2) # Basic Domain 2
    
    cutoff = sb_params['distance']
    contact_counts = defaultdict(int)

    # Controlliamo se i gruppi sono vuoti per evitare calcoli inutili
    calc_b1_a2 = g_b1.n_atoms > 0 and g_a2.n_atoms > 0
    calc_a1_b2 = g_a1.n_atoms > 0 and g_b2.n_atoms > 0

    if not calc_b1_a2 and not calc_a1_b2:
        return contact_counts # Niente da calcolare

    # Itera solo sui frame assegnati in questo chunk
    for ts in u.trajectory[frame_chunk]:
        contacting_residue_pairs = set()
        
        # Combinazione 1: Basic(D1) - Acidic(D2)
        if calc_b1_a2:
            # Calcola la matrice delle distanze
            dist_array = distances.distance_array(g_b1.positions, g_a2.positions)
            # Trova gli indici dove la distanza è <= cutoff
            contacts = np.where(dist_array <= cutoff)
            
            # 'contacts' è una tupla (array_indici_g_b1, array_indici_g_a2)
            for i_b1, i_a2 in zip(*contacts):
                # Risaliamo dall'atomo al residuo
                res_b1 = g_b1[i_b1].residue
                res_a2 = g_a2[i_a2].residue
                
                res_tuple1 = (res_b1.resname, res_b1.resid)
                res_tuple2 = (res_a2.resname, res_a2.resid)
                
                # Creiamo una tupla ordinata per rappresentare la coppia
                pair = tuple(sorted((res_tuple1, res_tuple2), key=lambda x: x[1]))
                contacting_residue_pairs.add(pair)

        # Combinazione 2: Acidic(D1) - Basic(D2)
        if calc_a1_b2:
            dist_array = distances.distance_array(g_a1.positions, g_b2.positions)
            contacts = np.where(dist_array <= cutoff)
            
            for i_a1, i_b2 in zip(*contacts):
                res_a1 = g_a1[i_a1].residue
                res_b2 = g_b2[i_b2].residue
                
                res_tuple1 = (res_a1.resname, res_a1.resid)
                res_tuple2 = (res_b2.resname, res_b2.resid)
                
                pair = tuple(sorted((res_tuple1, res_tuple2), key=lambda x: x[1]))
                contacting_residue_pairs.add(pair)
        
        # Aggiungiamo +1 al conteggio per ogni coppia *unica* trovata in questo frame
        for pair in contacting_residue_pairs:
            contact_counts[pair] += 1
            
    return contact_counts

# Rinominiamo la funzione
def calculate_salt_bridge_occupancy_parallel(topo_path, traj_files, sel_b1, sel_a2, sel_a1, sel_b2, sb_params, offset, frame_step):
    """
    Orchestrates the parallel calculation of Salt Bridge occupancies.
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
    
    # Passiamo le 4 selezioni al task
    tasks = [(topo_path, traj_files, chunk, sel_b1, sel_a2, sel_a1, sel_b2, sb_params)
             for chunk in frame_chunks if len(chunk) > 0]

    print(f"  Total frames: {n_frames_total}. Analyzing {n_frames_sampled} frames (1 every {frame_step}).")
    print(f"  Distributing Salt Bridge analysis over {len(tasks)} processes...")

    with multiprocessing.Pool(processes=n_procs) as pool:
        results = pool.map(worker_function, tasks)
        
    merged_counts = defaultdict(int)
    for count_dict in results:
        for pair, count in count_dict.items():
            merged_counts[pair] += count
            
    # L'occupancy è basata sul numero di frame *campionati*
    # Applichiamo l'offset qui per la stampa
    occupancies = {
        f"{res1[0]}{res1[1] + offset}-{res2[0]}{res2[1] + offset}": count / n_frames_sampled
        for (res1, res2), count in merged_counts.items()
    }
    
    return occupancies

# Questa funzione non è più necessaria, la sua logica è stata spostata in main()
# def get_domain_selection_string(...):


# --- MAIN SCRIPT ---

def main():
    print("Starting PARALLEL Salt Bridge occupancy analysis...")

    # Parametri per i ponti salini
    sb_params = {
        'distance': SB_DISTANCE_CUTOFF,
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

        safe_system_label = system_label.replace(' ', '_')

        for domain1_name, domain2_name in DOMAIN_PAIRS:
            pair_label = f'{domain1_name}-{domain2_name}'
            print(f"  Calculating Salt Bridge occupancy for: {pair_label}...")

            # Prendiamo i residui *originali* (senza offset)
            resids_1 = DOMAINS_RESIDS[domain1_name]
            resids_2 = DOMAINS_RESIDS[domain2_name]
            
            # Creiamo la stringa di selezione per il dominio
            domain_sel_1 = f"(protein and resid {resids_1[0]}-{resids_1[1]})"
            domain_sel_2 = f"(protein and resid {resids_2[0]}-{resids_2[1]})"

            # Costruiamo le 4 stringhe di selezione come da tuo commento
            # sel_b1: Basic(Domain1)
            sel_b1 = f"({domain_sel_1}) and ({sel_basic})"
            # sel_a2: Acidic(Domain2)
            sel_a2 = f"({domain_sel_2}) and ({sel_acidic})"
            
            # sel_a1: Acidic(Domain1)
            sel_a1 = f"({domain_sel_1}) and ({sel_acidic})"
            # sel_b2: Basic(Domain2)
            sel_b2 = f"({domain_sel_2}) and ({sel_basic})"
            
            # Chiamiamo la nuova funzione
            occupancies = calculate_salt_bridge_occupancy_parallel(
                topo_path, traj_files, sel_b1, sel_a2, sel_a1, sel_b2, 
                sb_params, RESIDUE_OFFSET, FRAME_STEP
            )
            
            significant_contacts = {
                pair: occ for pair, occ in occupancies.items() if occ >= OCCUPANCY_THRESHOLD
            }
            all_occupancy_data[system_label][pair_label] = significant_contacts
            
            pair_results_list = []
            for pair, occ in significant_contacts.items():
                # Aggiorniamo l'intestazione
                pair_results_list.append({'Domain Pair': pair_label, 'Salt Bridge': pair, 'Occupancy': occ})

            if pair_results_list:
                df_pair = pd.DataFrame(pair_results_list).sort_values(by='Occupancy', ascending=False)
                
                # Aggiorniamo il nome del file
                csv_filename = f'salt_bridge_occupancy_{safe_system_label}_{pair_label}.csv'
                csv_filepath = os.path.join(OUTPUT_DIR, csv_filename)
                
                df_pair.to_csv(csv_filepath, index=False)
                print(f"  Significant Salt Bridges for {pair_label} saved to: {csv_filepath}")

    # --- COMPARISON AND PLOTTING ---
    if not all_occupancy_data:
        print("\nNo data was analyzed. Cannot perform comparison.")
        return

    print("\n--- Comparing Salt Bridges across systems ---")
    
    # Aggiorniamo i nomi dei file e i titoli
    summary_file_path = os.path.join(OUTPUT_DIR, 'common_salt_bridges_summary.txt')
    with open(summary_file_path, 'w') as f:
        f.write(f"Summary of Common Salt Bridges (Occupancy >= {OCCUPANCY_THRESHOLD:.0%})\n")
        f.write("="*60 + "\n")

        for domain1, domain2 in DOMAIN_PAIRS:
            pair_label = f'{domain1}-{domain2}'
            f.write(f"\nAnalysis for Domain Pair: {pair_label}\n")
            
            # Usiamo i dati dei ponti salini
            contact_sets = [set(all_occupancy_data[sys_label].get(pair_label, {}).keys()) 
                            for sys_label in SYSTEMS.keys()]
            
            if not any(contact_sets):
                f.write("  No significant salt bridges found for this pair in any system.\n")
                continue

            common_contacts = set.intersection(*contact_sets)
            
            if common_contacts:
                f.write(f"\n  Found {len(common_contacts)} Salt Bridges common to ALL {len(SYSTEMS)} systems:\n")
                for contact in sorted(list(common_contacts)):
                    f.write(f"    - {contact}\n")
            else:
                f.write(f"\n  No Salt Bridges were found to be common across ALL {len(SYSTEMS)} systems.\n")

            all_contacts_for_pair = sorted(list(set.union(*contact_sets)))
            
            if not all_contacts_for_pair:
                continue

            heatmap_data = []
            for contact in all_contacts_for_pair:
                row = {'Salt Bridge': contact} # Aggiorniamo etichetta
                for sys_label in SYSTEMS.keys():
                    occupancy = all_occupancy_data[sys_label].get(pair_label, {}).get(contact, 0)
                    row[sys_label] = occupancy
                heatmap_data.append(row)
            
            heatmap_df = pd.DataFrame(heatmap_data).set_index('Salt Bridge') # Aggiorniamo etichetta
            
            height = min(100, max(8, len(heatmap_df) * 0.3)) 
            
            plt.figure(figsize=(10, height))
            sns.heatmap(heatmap_df, cmap="viridis", annot=True, fmt=".2f", linewidths=.5)
            # Aggiorniamo titoli e etichette
            plt.title(f'Salt Bridge Occupancy Heatmap for {pair_label}', fontsize=16)
            plt.xlabel('System', fontsize=12)
            plt.ylabel('Residue Pair', fontsize=12)
            plt.xticks(rotation=0, ha='center') # Rimesso a 0 per i 5 sistemi
            plt.tight_layout()
            
            heatmap_path = os.path.join(OUTPUT_DIR, f'heatmap_salt_bridge_occupancy_{domain1}_{domain2}.png')
            plt.savefig(heatmap_path, dpi=300)
            print(f"  Heatmap saved to: {heatmap_path}")
            plt.close()
    
    print(f"\nComparison summary saved to: {summary_file_path}")
    print("\nAdvanced Salt Bridge analysis successfully completed!")

if __name__ == '__main__':
    multiprocessing.freeze_support()
    main()
