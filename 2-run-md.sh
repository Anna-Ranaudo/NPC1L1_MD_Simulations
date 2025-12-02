#!/bin/bash
# here add cluster specific options and slurm commands, e.g.:
#SBATCH --job-name job-name
#SBATCH -N1 --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --gres=gpu:4
#SBATCH --account=
#SBATCH --partition=
##SBATCH --qos=boost_qos_dbg
#SBATCH --mail-type=ALL
#SBATCH --mail-user=

# load amber (the following lines are specific to the cluster that was used 
# to run the simulations of this study)
module load profile/lifesc
module load  amber/22--openmpi--4.1.4--gcc--11.3.0-tk-cuda-11.8

. $AMBERHOME/amber.sh

export OMP_NUM_THREADS=8

# define topology, starting coordinates, path to the md config files and working directory
export prmtop="./input_structures/6v3f/bound/bilayer_6v3f-x-MD.amber_lipid.top"
export coord="./input_structures/6v3f/bound/bilayer_6v3f-x-MD.amber_lipid.crd"
export md_config_files="./amber_md_config_files/6v3f/bound/"
export working_dir="./run/6v3f/bound/"

# create working directory for replica n
export rep="rep1"
mkdir -p $working_dir/$rep

### Equilibration ##
mpirun -np 4 $AMBERHOME/bin/pmemd.cuda_SPFP.MPI \
-O \
-i $md_config_files/03_Heat.in \
-o $working_dir/$rep/03_Heat.out \
-p $prmtop \
-c $working_dir/02_Min2.rst \
-r $working_dir/$rep/03_Heat.rst \
-x $working_dir/$rep/03_Heat.nc \
-ref $working_dir/02_Min2.rst \
-inf $working_dir/$rep/03_Heat.mdinfo

mpirun -np 4 $AMBERHOME/bin/pmemd.cuda_SPFP.MPI \
-O \
-i $md_config_files/04_Heat2.in \
-o $working_dir/$rep/04_Heat2.out \
-p $prmtop \
-c $working_dir/$rep/03_Heat.rst \
-r $working_dir/$rep/04_Heat2.rst \
-x $working_dir/$rep/04_Heat2.nc \
-ref $working_dir/$rep/03_Heat.rst \
-inf $working_dir/$rep/04_Heat2.mdinfo

## Protein backbone and ligand restrained
mpirun -np 4 $AMBERHOME/bin/pmemd.cuda_SPFP.MPI \
-O \
-i $md_config_files/05_Back.in \
-o $working_dir/$rep/05_Back.out \
-p $prmtop \
-c $working_dir/$rep/04_Heat2.rst \
-r $working_dir/$rep/05_Back.rst \
-x $working_dir/$rep/05_Back.nc \
-ref $working_dir/$rep/04_Heat2.rst \
-inf $working_dir/$rep/05_Back.mdinfo

## Protein C-alpha atoms and ligand restrained
mpirun -np 4 $AMBERHOME/bin/pmemd.cuda_SPFP.MPI \
-O \
-i $md_config_files/06_Calpha.in \
-o $working_dir/$rep/06_Calpha.out \
-p $prmtop \
-c $working_dir/$rep/05_Back.rst \
-r $working_dir/$rep/06_Calpha.rst \
-x $working_dir/$rep/06_Calpha.nc \
-ref $working_dir/$rep/05_Back.rst \
-inf $working_dir/$rep/06_Calpha.mdinfo

## 100ns NPT run ##  equil all
mpirun -np 4 $AMBERHOME/bin/pmemd.cuda_SPFP.MPI \
-O \
-i $md_config_files/07_Prod.in \
-o $working_dir/$rep/07_Prod.out \
-p $prmtop \
-c $working_dir/$rep/06_Calpha.rst \
-r $working_dir/$rep/07_Prod.rst \
-x $working_dir/$rep/07_Prod.nc \
-ref $working_dir/$rep/06_Calpha.rst \
-inf $working_dir/$rep/07_Prod.mdinfo
    
## 400ns prod ##
mpirun -np 4 $AMBERHOME/bin/pmemd.cuda_SPFP.MPI \
-O \
-i $md_config_files/08_Long.in \
-o $working_dir/$rep/08_Long.out \
-p $prmtop \
-c $working_dir/$rep/07_Prod.rst \
-r $working_dir/$rep/08_Long.rst \
-x $working_dir/$rep/08_Long.nc \
-ref $working_dir/$rep/07_Prod.rst \
-inf $working_dir/$rep/08_Long.mdinfo
# needs probably to be restarted if a walltime is present on the cluster
# and the simulation is not finished
