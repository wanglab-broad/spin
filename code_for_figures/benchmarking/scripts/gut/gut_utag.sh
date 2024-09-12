#!/bin/bash

#$ -l h_vmem=256G
#$ -l h_rt=01:30:00
#$ -o /stanley/WangLab/kamal/outputs/spin/reviews/gut/utag/
#$ -j y
#$ -t 1-110

source /broad/software/scripts/useuse
reuse Anaconda3
source activate /stanley/WangLab/bridget/envs/utag/

sample_rates=($(seq 0 0.1 1))
random_seeds=($(seq 0 1 10))
# random_seeds=($(seq 0 0)) # test using a single random seed

n_sample_rates=${#sample_rates[@]}
sample_rate=${sample_rates[($((SGE_TASK_ID-1))%$n_sample_rates)]}
random_seed=${random_seeds[($((SGE_TASK_ID-1))/$n_sample_rates)]}

radius=260 # ~ 50 nbrs
k_latent=15

python /stanley/WangLab/kamal/code/projects/spin/reviews/scripts/utag_.py \
--read_path /stanley/WangLab/kamal/data/projects/spin/reviews/gut/adata.h5ad \
--sample_rate $sample_rate \
--radius $radius \
--k_latent $k_latent \
--n_pcs 50 \
--leiden_res "0.2" \
--write_path /stanley/WangLab/kamal/data/projects/spin/reviews/results/gut/utag/reconstructions_${sample_rate}_${radius}radius_${k_latent}latent/ \
--random_seed $random_seed