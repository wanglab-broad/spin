#!/bin/bash

# Request memory
#$ -l h_vmem=4G

# Request runtime
#$ -l h_rt=00:10:00

# Specify output file destination
#$ -o /stanley/WangLab/kamal/outputs/spin_tests/

# Join output and error files
#$ -j y


# Set up Python environment
source /broad/software/scripts/useuse
reuse Anaconda3
source activate /stanley/WangLab/envs/spin

# Set up arguments
adata_path=/stanley/WangLab/kamal/data/spin_tests/single.h5ad
write_path=/stanley/WangLab/kamal/data/spin_tests/single_spin.h5ad
resolution="0.5"
n_pcs=50

# Run script
python /stanley/WangLab/kamal/code/integration/spatial/spin/src/spin/spin.py --adata_paths $adata_path --write_path $write_path --resolution $resolution --n_pcs $n_pcs