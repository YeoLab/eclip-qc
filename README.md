# Eclip Quality Control

### Description
From the RBP eClip pipeline, a sample may contain sequences that are unmapped to the genome of interest. We want to blast these unmapped sequences to interpret what a majority of these sequences are. 

#### Workflow Steps

1. Take the input sample bam file and count the number of sequences overall.
2. Randomly downsample the original bam file to the number of sequences specified in the config.yaml file.
3. Convert downsampled bam file to an unfiltered fasta file.
4. Run BLAST against sequences using a database of synthetic/artificial sequences.
5. Filter out synthetic sequences in unfiltered fasta file from `unmapped_counts/*_blast_synthetic.tsv`.
- We select the highest `pident` (percent of identical matches) for each `qseqid` from the blast_synthetic.tsv file, and use `qstart` and `qend` to remove those synthetic sequences in the unfiltered fasta file. 
- For example, one `qseqid` is “AAA” and its sequence in the unfiltered fasta file is “ABCDEF”. The blast_synthetic.tsv has `qstart` = 2 and `qend` = 5. The filtering will be done by removing “BCDE” from the sequence. Now your filtered fasta file will still have a sequence identifier to be “AAA”, but its sequence will only have “AF”.
6. Run filtered fasta file against BLASTn and BLASTx
7. Generate pie chart summarizing top BLASTn, BLASTx outputs
- In each of the BLAST outputs, each `qseqid` can have at most five hits. We have included additional priorities based on categories of interest. For the results of each `qseqid`, `pident` is prioritized over the Category priority. In the case that all the results of a `qseqid` have the same `pident` value and same Category priority, we take the alphabetically `qseqid` first hit. 

#### Category Identity:
- priority['Homo sapiens'] = Largest_count_of_matched_species + 100
- priority['Mus musculus'] = priority['Homo sapiens'] - 1
- priority['bacteria'] = priority['Mus musculus'] - 1
- priority['bacterium'] = priority['bacteria']
- priority['virus'] = priority['bacteria']
- priority['strain'] = priority['bacteria']

### How to interpret output files and troubleshoot if encounter errors.

1. `unmapped_counts/*unmappedblast_downsampled_blastn.tsv` and `unmapped_counts/*unmappedblast_downsampled_blastx.tsv` are the key files to look at. These files contain the blast results that are used to build the figures which are information-limited. In general, people would want to look at these blast files in different aspects.

2. Piechart file: `pieChart/ [sample].png` contains all the figures that summarize unmapped sequences for ths sample.
- Blastx/blastn mapping percentage: indicates what percentage of the total sequences blasted had a result
- Filtered blastn/blastx piechart: indicates what percentage of the total sequences are of the top six hits 
- Top Three Hits for Blastn/Blastx: the top three hits for blastn/blastx and the frequency of these hits. There might appear "No_Species" in the figures. It means our databases cannot find a name for the `sseqid` for that specific hit. When larger database is employed in the future, the issue should be solved.
- Blastn/Blastx Breakdown of ‘Other’ of Interest Section: From the results within the ‘Other’ category, a breakdown of the number of hits of the results that are in categories of Homo sapiens, Bacteria/Virus, or Mus Musculus 

3. Sometimes, the pipeline can fail due to unstable internet, simply rerun the pipeline with `--rerun-incomplete` parameter in the command line might solve the problem. If the pipeline failure is obviously cause by other errors, if using cluster, `stderr` and `stdout` folders will record pipeline progress for each rule and each sample. Therefore, `stderr` will record the error messages if the pipeline fails. Hope this is helpful. Otherwise, users can refer to snakemake message logs at `cd absolute/path/eclip-qc/.snakemake/log/` to find the error messages.

### Please make sure to read the information below before running the pipeline.

Snakemake has a ton of additional parameters, we are only using a few here.

`--snakefile` or `s`: snakemake pipeline file to use, usually you don't make changes to that file unless needed. It's a file in python language.

`--configfile`: specifies the config.yaml file that contains some parameters. Adjust this config.yaml file accordingly to your input files path. Importantly, specify the number of reads you want to consider for your debugging purpose. This number represents the number of sequences that you want to take from your unmapped reads. Higher this number is, more information, maybe more accurate, your results will be. However, it does make the pipeline to run longer time. Below 1,000,000 is recommended, but set it to any values as needed.

`--use-singularity`: Use Docker containers for running pipelines.

`--cores 8`: makes the run to use more sources and run faster. 

`-j 8`: makes the run to process 8 jobs at the same time, also make the run faster.

`--latency-wait 30` or `-w 30`: allows some extra 30 seconds waiting time for the pipeline to wait for some slow processed files, instead of throwing an error right away.

`-p`: helps to print the current running commands which is helpful for debugging.

`--keeping-going` or `-k`: Normally, when a job in a Snakemake workflow fails, the entire workflow is halted immediately. However, with this option, Snakemake will continue to execute other independent jobs that are not affected by the failure of the current job. 

`--singularity-args`: lets snakemake pipeline know where to find your files. For example, your input bam file's absolute path is `/project/test/test.bam`. Snakemake doesn't know to start with the `project` folder in your computer's root directory to look for your bam file. Therefore, you can set this parameter to `--singularity-args "--bind /project"`, now `/project` becomes one of the options for snakemake to look for resource files.

`--cluster`: If applicant, users may want to use online clusters to run the pipeline to have more computational power.


### Extra notes:
- For testing and debugging purposes, you may also look at the examples in the unmapped_counts or pieChart folder to have an idea what the output files will look like. 
- For most important result information, pieChart and unmapped_counts folder should be where you want to check out for sure.
- The default config.yaml file contains paths to sample files (which can be found in the samplefiles folder). These are smaller example bam files and can be used for a test run to see how the outputs will look like.

# Start

To run the pipeline:

```bash
# load the latest snakemake package (This has been pre-installed for you on TSCC and will provide the path to the snakemake command)
module load snakemake/7.17.1
# use the command format below to start a fresh run or continue an unfinished run
snakemake -kps Snakefile --configfile config.yaml -j 8 --use-singularity --cores 8 -w 30 --singularity-args "--bind /oasis --bind /projects --bind /scratch"
# Use the clusters qsub
snakemake -kps Snakefile --configfile config.yaml -j 8 --use-singularity --cores 8 -w 30 --singularity-args "--bind /oasis --bind /projects --bind /scratch" --cluster "qsub -e {params.error_file} -o {params.out_file} -l walltime={params.run_time} -l nodes=1:ppn={threads} -q home-yeo"
```

