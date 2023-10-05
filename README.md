# eclip quality control

### Please make sure to read the information below before running the pipeline.

Snakemake has a ton of additional parameters, we are only using a few here.

`--snakefile`: snakemake pipeline file to use, usually you don't make changes to that file unless needed. It's a file in python language.

`--configfile`: specifies the config.yaml file that contains some parameters. Adjust this config.yaml file accordingly to your input files path. Importantly, specify the number of reads you want to consider for your debugging purpose. This number represents the number of sequences that you want to take from your unmapped reads. Higher this number is, more information, maybe more accurate, your results will be. However, it does make the pipeline to run longer time. Below 1,000,000 is recommended, but set it to any values as needed.

`--use-singularity`: Use Docker containers for running pipelines.

`--cores 8`: makes the run to use more sources and run faster. 

`-j 8`: makes the run to process 8 jobs at the same time, also make the run faster.

`-latency-wait 30`: allows some extra 30 seconds waiting time for the pipeline to wait for some slow processed files, instead of throwing an error right away.

`-p` helps to print the current running commands which is helpful for debugging.

`--singularity-args`: lets snakemake pipeline know where to find your files. For example, your input bam file's absolute path is `/project/test/test.bam`. Snakemake doesn't know to start with the `project` folder in your computer's root directory to look for your bam file. Therefore, you can set this parameter to `--singularity-args "--bind /project"`, now `/project` becomes one of the options for snakemake to look for resource files.

### Please make sure to update the path of your config file accordingly in the beginning of the script file `scripts/blastresults_piechart.py`. The path format should be similar to `/absolute/path/to/eclip-qc/config.yaml`.

### Extra notes:
- For testing and debugging purposes, you may also look at the examples in the unmapped_counts or pieChart folder to have an idea what the output files will look like. 
- For most important result information, pieChart and unmapped_counts folder should be where you want to check out for sure.

# Start

To run the pipeline:

```bash
# load the latest snakemake package (This has been pre-installed for you on TSCC and will provide the path to the snakemake command)
module load snakemake/7.17.1
# use the command format below to start a fresh run or continue an unfinished run
snakemake --snakefile Snakefile --configfile config.yaml -j 8 --use-singularity --cores 8 --latency-wait 30 --singularity-args "--bind /oasis --bind /projects --bind /home"
```
