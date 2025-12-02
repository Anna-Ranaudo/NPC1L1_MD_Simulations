# NPC1L1 MD Simulations

Repository containing files and scripts for molecular dynamics simulations of the transmembrane protein NPC1L1.

## The repository is structured as follows:                
- `input_structures`:
    - `6v3f`: Contains the input structures for the NPC1L1 open conformation (PDB ID: [6V3F](https://www.rcsb.org/structure/6V3F)).
    - `6v3h`: Contains the input structures for the NPC1L1 closed conformation (PDB ID: [6V3H](https://www.rcsb.org/structure/6V3H)).
- `scripts`:
    - `set-up-simulations`: Contains the scripts and protocols to prepare the files for running the simulations.
    - `analyze-simulations`: Contains the scripts for analyzing the simulations.
- `amber_md_config_files`: Contains the amber md configuration files.



## To build the membrane and topology files:
- `scripts/set-up-simulations/0-build-membrane-and-top-files.sh`.

## To run the energy minimization + MD simulation:
- `1-run-energy-min.sh`
- `2-run-md.sh`

## To run the analysis:
- `scripts/analyze-simulations/...`
