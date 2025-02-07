#!/bin/bash

for iteration in {1..1}; do
    /bin/bash --wait -o resTest.out launcherStrongScaling_64_withGenoa.sh
    echo iteration
    grep -rw /lus/home/CT6/cad14985/jauriac/bench-in-situ/results_withGenoa/NoDeisa/64 -e "RESULT" >> one-64-output_withGenoa.txt
done
