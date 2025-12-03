#!/bin/bash
# this is run in /home/rame/NPC1L1_MD_Simulations/scripts/set-up-simulations
# requires AmberTools24

# 0. Start with already prepared .pdb files (here they were prepared with the protein preparation workflow in Maestro)
# they are the bound forms of the two conformations of NPC1L1 (6v3f = open conformation, 6v3h = closed conformation)

# 1. Prepare .pdb files for Amber
perl PDB4amber.pl ../../input_structures/6v3f/bound/6V3F-x-MD.pdb
perl PDB4amber.pl ../../input_structures/6v3h/bound/6V3H-x-MD.pdb
# the files generated are 6V3F-x-MD.amber.pdb and 6V3H-x-MD.amber.pdb

# 2. Run packmol-memgen to build the membrane and .top and .crd files (input for MD simulations with Amber)
cd ../../input_structures/6v3f/bound

packmol-memgen  --pdb 6V3F-x-MD.amber.pdb  --pdb2pqr  --lipids  \
POPC:DOPC:POPE:DOPE:POPS:PSM:SSM:OSM:POPI:POPI45H:CHL1//POPC:DOPC:POPE:DOPE:PSM:SSM:OSM:CHL1 \ 
ratio 0.16:0.31:0.21:0.68:0.38:0.08:0.08:0.13:0.16:0.05:1.00//0.29:0.57:0.04:0.13:0.28:0.28:0.45:0.98 \
--ppm --keep --parametrize  --ffwat opc --ffprot ff19SB --fflip lipid21  --keepligs --salt
# run this command twice, because the first time it stops at pdb2pqr
# it first generates this file (.pdb format, with the membrane)
#      bilayer_6V3F-x-MD.amber.pdb
# then the --parametrize flag makes tleap run, reading leap.in, with the ff specified
# and it generates the following files:
#      bilayer_6V3F-x-MD.amber_lipid.pdb
#      bilayer_6V3F-x-MD.amber_lipid.crd  # coordinates for MD
#      bilayer_6V3F-x-MD.amber_lipid.top  # topology, can be renamed to .prmtop 

# 2.1 do the same for 6v3h
cd ../../6v3h/bound
packmol-memgen  --pdb 6V3H-x-MD.amber.pdb  --pdb2pqr  --lipids  \
POPC:DOPC:POPE:DOPE:POPS:PSM:SSM:OSM:POPI:POPI45H:CHL1//POPC:DOPC:POPE:DOPE:PSM:SSM:OSM:CHL1 \
ratio 0.16:0.31:0.21:0.68:0.38:0.08:0.08:0.13:0.16:0.05:1.00//0.29:0.57:0.04:0.13:0.28:0.28:0.45:0.98 \
--ppm --keep --parametrize  --ffwat opc --ffprot ff19SB --fflip lipid21  --keepligs --salt

cd ../../

# 3. For the apo structures, we need to remove the ligands (in order to preserve the membrane we just built)
# we copy e.g. the bilayer_6V3F-x-MD.amber.pdb to the apo structure directory
cp 6v3f/bound/bilayer_6V3F-x-MD.amber.pdb 6v3f/apo/
cp 6v3h/bound/bilayer_6V3H-x-MD.amber.pdb 6v3h/apo/
# we rename them and we remove the ligands (the ligands are the residues n 1262)
mv 6v3f/apo/bilayer_6V3F-x-MD.amber.pdb 6v3f/apo/bilayer_6V3F-x-MD-apo.amber.pdb
mv 6v3h/apo/bilayer_6V3H-x-MD.amber.pdb 6v3h/apo/bilayer_6V3H-x-MD-apo.amber.pdb
# we copy the leap.in files to the apo structure directory and modify them accordingly
cp 6v3f/bound/leap.in 6v3f/apo/
cp 6v3h/bound/leap.in 6v3h/apo/
# we go the apo directories and run tleap to generate the .prmtop and .crd files
cd 6v3f/apo
tleap -f leap.in
cd ../../6v3h/apo
tleap -f leap.in

cd ../../../
