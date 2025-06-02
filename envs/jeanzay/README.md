## compilation
The environnement file for compilation "env.sh" are in the directory "envs/jeanzay/scripts_cluster"
Configure this file before compilation for example:
- GPU_ARCH: V100, A100, H100
- CASE_BATCH: strong or weak

``` bash
cd envs/jeanzay/
cd scripts_cluster
bash build_all.sh
```

## launch script
The file "env/jeanzay/param_launcher.ini" set the launcher parameter. 
Configure this file before launching the script:
``` bash
cd envs/jeanzay/
bash launcherJeanZay.sh
```

The results of the bench is writting in the directory 

strong/bench_result_directory

or

weak/bench_result_directory

An example of "bench_result_directory" is "jeanzay_20250409_131805_V100_16G"

## post processing

### strong case:
In strong directory (envs/jeanzay/strong),

``` bash
cd envs/jeanzay/
cd strong
bash post_proccessing.sh bench_result_directory
```

The output files generated are:
- A temporary .txt file with the different execution time are generated:
    strong/bench_result_directory/bench_output_bench_result_directory.txt
- The plot are generated in: 
    strong/bench_result_directory/Strong_scaling_plot_bench_result_directory.pdf

### weak case:
In weak directory (envs/jeanzay/weak),
``` bash
cd envs/jeanzay
cd weak
bash post_proccessing.sh bench_result_directory
```

The output files generated are:
- A temporary .txt file with the different execution time are generated:
    weak/bench_result_directory/bench_output_bench_result_directory.txt
- The plot are generated in:
    weak/bench_result_directory/Weak_scaling_plot_bench_result_directory.pdf
