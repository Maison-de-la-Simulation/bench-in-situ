#!/bin/bash
#SBATCH --account=cad14985
#sbatch --output=one-64NSizeBench.out
#SBATCH --job-name=I-bis_nd
#SBATCH --constraint=GENOA
#SBATCH --nodes=1
##SBATCH --exclusive
#SBATCH --time=20:00:00
##SBATCH --nodelist=c1155


for iteration in {1..1}; do
    sbatch --wait -o resTest.out launcherStrongScaling_64_compEOY.sh
    echo iteration
    grep -rw /lus/home/CT6/cad14985/jauriac/bench-in-situ/results_EOY24/NoDeisa/64 -e "RESULT" >> one-64-output_EOY24.txt
done
