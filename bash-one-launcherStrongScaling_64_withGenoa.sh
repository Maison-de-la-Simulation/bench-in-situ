#!/bin/bash

for iteration in {1..1}; do
    /bin/bash --wait -o resTest.out launcherStrongScaling_64_withGenoa.sh
    echo iteration
    grep -rw ${PWD}/results_withGenoa/NoDeisa/64 -e "RESULT" >> one-64-output_withGenoa.txt
done
