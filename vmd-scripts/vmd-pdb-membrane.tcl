# ===============================
# LIST PDB
# ===============================
set pdb_list {
    input_structures/6v3h/bound/6v3h-bound-input-no-sol.pdb
    input_structures/6v3f/bound/6v3f-bound-input-no-sol.pdb
}

# ===============================
# CICLO SU OGNI PDB
# ===============================
foreach pdb $pdb_list {

# ===============================
# LOAD
# ===============================
mol new $pdb type pdb waitfor all
set molid [molinfo top]
mol delrep 0 $molid

color Display Background silver

# ===============================
# PROTEIN
# ===============================
mol selection "protein"
mol representation NewCartoon
mol color ColorID 2
mol addrep $molid

# ===============================
# LIPID TAILS
# ===============================
mol selection "resname OL PA SA ST and not hydrogen"
mol representation Lines 1.6
mol color ColorID 2
mol addrep $molid

# ===============================
# HEADGROUPS
# ===============================
proc add_lipid {molid resname colorid} {
    mol selection "resname $resname and not hydrogen"
    mol representation Lines 1.6
    mol color ColorID $colorid
    mol addrep $molid
}

add_lipid $molid PC   0
add_lipid $molid PE   10
add_lipid $molid PS   1
add_lipid $molid PI   6
add_lipid $molid SPM  7
add_lipid $molid H2C  3

# ===============================
# CHOLESTEROL
# ===============================
mol selection "resname CHL and not hydrogen"
mol representation Lines 1.6
mol color ColorID 4
mol addrep $molid

}

