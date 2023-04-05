#!/bin/bash

# Request memory
#$ -l h_vmem=16G

# Request runtime
#$ -l h_rt=00:30:00

# Specify output file destination
#$ -o /stanley/WangLab/kamal/outputs/spin_tests/

# Join output and error files
#$ -j y


# Set up Python environment
source /broad/software/scripts/useuse
reuse Anaconda3
source activate /stanley/WangLab/envs/spin

# Set up arguments
n_nbrs=30
n_samples=12
resolution="0.6"
adata_path=/stanley/WangLab/kamal/data/spin_tests/xspecies.h5ad
write_path=/stanley/WangLab/kamal/data/spin_tests/xspecies_spin.h5ad

# Run script
python /stanley/WangLab/kamal/code/integration/spatial/spin/src/spin/spin.py --adata_paths $adata_path --write_path $write_path --batch_key dataset --batch_labels marmoset mouse --n_nbrs $n_nbrs --n_samples $n_samples --resolution $resolution