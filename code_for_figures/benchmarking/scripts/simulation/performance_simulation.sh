#!/bin/bash

tissue=simulation

simulation_methods=(knn radius utag stagate graphst)
simulation_jobids=(42429477 42429759 42429768 42430818 42430894)
n_methods="${#simulation_methods[@]}"

# Write job resource metrics to disk
for i in $(seq 0 $((n_methods-1))); do
    method=${simulation_methods[i]}
    jobid=${simulation_jobids[i]}
    write_path=/stanley/WangLab/kamal/data/projects/spin/reviews/results/${tissue}/${method}/
    performance=${write_path}performance_${tissue}_${method}_${jobid}.txt
    echo "Writing $performance"
    touch $performance
    nonsub_tasks=($(seq 11 11 110)) # only compare using non-subsampling runs
    for task in ${nonsub_tasks[@]}; do
        qacct -j $jobid -t $task >> $performance;
    done;
done
echo "Done"