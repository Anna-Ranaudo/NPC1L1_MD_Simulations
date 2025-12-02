#!/bin/bash
# here add cluster specific options and slurm commands, e.g.:
#SBATCH --job-name job-name
#SBATCH -N1 --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --gres=gpu:1
#SBATCH --account=
#SBATCH --partition=
#SBATCH --qos=boost_qos_dbg
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
mkdir -p $working_dir/


## 1st minimisation ##
$AMBERHOME/bin/pmemd \
-O \
-i $md_config_files/01_Min.in \
-o $working_dir/01_Min.out \
-p $prmtop \
-c $coord \
-r $working_dir/01_Min.rst \
-inf $working_dir/01_Min.mdinfo

## 2nd minimisation ##
$AMBERHOME/bin/pmemd.cuda_SPFP \
-O \
-i $md_config_files/02_Min2.in \
-o $working_dir/02_Min2.out \
-p $prmtop \
-c $working_dir/01_Min.rst \
-r $working_dir/02_Min2.rst \
-inf $working_dir/02_Min2.mdinfo
