#!/bin/bash

tissue=gut

methods=(knn radius utag stagate graphst)
jobids=(42437430 42437787 42437788 42437789 42437790)
n_methods="${#methods[@]}"

# Write job resource metrics to disk
for i in $(seq 0 $((n_methods-1))); do
    method=${methods[i]}
    jobid=${jobids[i]}
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