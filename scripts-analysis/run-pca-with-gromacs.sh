#!/bin/bash
# 1. convert nc to xtc with MDAnalysis

# 2. if need, concatenate the trajectory files
# gmx trjcat -f ... # list the xtc files of the single replica

# 3. run PCA performing a fitting on the transmembrane region
echo 1 0 | gmx covar -f ../5t07-08-prot-lig-pbc.xtc -s ../../run-md1/f0-07-08-prot-lig-pbc.gro -n index.ndx
# (in the above command, 1 is the index of the transmembrane region, 0 is the index of the whole protein)

# 3.1 project the individual trajectories on the cumulative PC1-PC2 space
echo 1 0 | gmx anaeig -f ../../run-md1/07-08-prot-lig-pbc.xtc -s ../../run-md1/pca-gromacs-07-08/average.pdb -v eigenvec.trr -eig eigenval.xvg -first 1 -last 2 -2d plot12-md1.xvg -n index.ndx
echo 1 0 | gmx anaeig -f ../../run-md2/07-08-prot-lig-pbc.xtc -s ../../run-md2/pca-gromacs-07-08/average.pdb -v eigenvec.trr -eig eigenval.xvg -first 1 -last 2 -2d plot12-md2.xvg -n index.ndx
echo 1 0 | gmx anaeig -f ../../run-md3/07-08-prot-lig-pbc.xtc -s ../../run-md3/pca-gromacs-07-08/average.pdb -v eigenvec.trr -eig eigenval.xvg -first 1 -last 2 -2d plot12-md3.xvg -n index.ndx
echo 1 0 | gmx anaeig -f ../../run-md4/07-08-prot-lig-pbc.xtc -s ../../run-md4/pca-gromacs-07-08/average.pdb -v eigenvec.trr -eig eigenval.xvg -first 1 -last 2 -2d plot12-md4.xvg -n index.ndx
echo 1 0 | gmx anaeig -f ../../run-md5/07-08-prot-lig-pbc.xtc -s ../../run-md5/pca-gromacs-07-08/average.pdb -v eigenvec.trr -eig eigenval.xvg -first 1 -last 2 -2d plot12-md5.xvg -n index.ndx

# 3.2 calculate the RMSF of the first two eigenvectors
echo 0 | gmx anaeig -f ../5t07-08-prot-lig-pbc.xtc -s ../../run-md1/f0-07-08-prot-lig-pbc.gro -first 1 -last 1 -rmsf eigrmsf-5t-6v3f-col-pc1.xvg -n index.ndx
echo 0 | gmx anaeig -f ../5t07-08-prot-lig-pbc.xtc -s ../../run-md1/f0-07-08-prot-lig-pbc.gro -first 2 -last 2 -rmsf eigrmsf-5t-6v3f-col-pc2.xvg -n index.ndx

# 3.3 calculation on the average structure of the two extreme projections along the trajectory
# and interpolation of 100 frames in between them
# as explained in https://manual.gromacs.org/2024.2/onlinehelp/gmx-anaeig.html
echo 1 0 |  gmx anaeig -f ../5t07-08-prot-lig-pbc.xtc -s ../../run-md1/f0-07-08-prot-lig-pbc.gro -first 1 -last 1 -extr 5t-col-extr-eigenvect1.pdb -nframes 100 -n index.ndx
echo 1 0 |  gmx anaeig -f ../5t07-08-prot-lig-pbc.xtc -s ../../run-md1/f0-07-08-prot-lig-pbc.gro -first 2 -last 2 -extr 5t-col-extr-eigenvect2.pdb -nframes 100 -n index.ndx
