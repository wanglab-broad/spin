#!/bin/bash

# Request memory
#$ -l h_vmem=8G

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
adata_path_1=/stanley/WangLab/kamal/code/remote_notebooks/mpfc_xtech/adata_starmap.h5ad
adata_path_2=/stanley/WangLab/kamal/code/remote_notebooks/mpfc_xtech/adata_slideseq.h5ad
adata_path_3=/stanley/WangLab/kamal/code/remote_notebooks/mpfc_xtech/adata_visium.h5ad
write_path=/stanley/WangLab/kamal/data/spin_tests/xtech_spin.h5ad
resolution="0.4"
n_pcs=20

# Run script
python /stanley/WangLab/kamal/code/integration/spatial/spin/src/spin/spin.py --adata_paths $adata_path_1  $adata_path_2  $adata_path_3 --write_path $write_path --batch_key tech --batch_labels starmap slideseq visium --n_nbrs 30 400 20 --n_samples 12 133 7 --resolution $resolution --n_pcs $n_pcs