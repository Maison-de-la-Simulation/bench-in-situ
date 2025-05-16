## compilation
The environnement file for compilation "env.sh" are in the directory scripts_cluster
Configure this file before compilation for example:
- GPU_ARCH
- CASE_BATCH

cd scripts_cluster
bash build_all.sh

## launch script
The file "param_launcher.ini" give the launcher parameter.

bash launcherJeanZay.sh

The results of the bench is writting on the directory 

strong/bench_result_directory

or

weak/bench_result_directory

where an example of "bench_result_directory" is "jeanzay_20250409_131805_V100_16G"

## post processing

### strong case:
In strong directory,

bash post_proccessing.sh bench_result_directory

output:
- A temporary .txt file with the different execution time are generated:
    strong/bench_result_directory/bench_output_bench_result_directory.txt
- The plot are generated in: 
    strong/bench_result_directory/Strong_scaling_plot_bench_result_directory.pdf

### weak case:
In weak directory,

bash post_proccessing.sh bench_result_directory

output:
- A temporary .txt file with the different execution time are generated:
    weak/bench_result_directory/bench_output_bench_result_directory.txt
- The plot are generated in:
    weak/bench_result_directory/Weak_scaling_plot_bench_result_directory.pdf
