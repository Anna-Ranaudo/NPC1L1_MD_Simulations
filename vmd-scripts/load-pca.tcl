# ===============================
# LISTA DEI PDB DA CARICARE
# ===============================
set pdb_list {
    5t-no-col-extr-eigenvect1.pdb
    5t-no-col-extr-eigenvect2.pdb
    ../../../col/all-pca-gromacs-07-08/pca-fit-transmemb-calc-all/5t-col-extr-eigenvect1.pdb
    ../../../col/all-pca-gromacs-07-08/pca-fit-transmemb-calc-all/5t-col-extr-eigenvect2.pdb
}

# ===============================
# CICLO SU OGNI PDB
# ===============================
foreach pdb $pdb_list {

    mol new $pdb type pdb waitfor all
    set molid [molinfo top]

    # -------------------------------
    # RAPPRESENTAZIONE BASE: TUTTO GRIGIO
    # -------------------------------
    mol delrep 0 $molid

    mol selection "all"
    mol representation Tube 0.5
    mol color ColorID 2      ;# grigio
    mol addrep $molid

    # -------------------------------
    # GRUPPI DI RESIDUI COLORATI
    # -------------------------------

    # resid 0–240 → ColorID 0 (blu)
    mol selection "resid 0 to 240"
    mol representation Tube 0.5
    mol color ColorID 0
    mol addrep $molid

    # resid 357–610 → ColorID 7 (verde)
    mol selection "resid 357 to 610"
    mol representation Tube 0.5
    mol color ColorID 7
    mol addrep $molid

    # resid 860–1083 → ColorID 27 (magenta)
    mol selection "resid 860 to 1083"
    mol representation Tube 0.5
    mol color ColorID 27
    mol addrep $molid

}


