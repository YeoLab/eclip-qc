# eclip-qc

To run a snakemake pipeline:

```bash
module load snakemake/7.17.1;  # load the latest snakemake package (This has been pre-installed for you on TSCC and will provide the path to the snakemake command)
snakemake --snakefile /path/to/Snakefile -j 8 --use-conda -p  # use the --snakefile parameter to specify which Snakefile to run, use -j to specify how many cores to use. -p prints the command out, useful for debugging.
```

Snakemake has a ton of additional options (such as allowing the use of full cluster resources), but you can use the above command to start.

# Extra notes:
- For testing and debugging purposes, you may also 