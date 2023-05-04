configfile: "config.yaml"
import os
OUTDIR = config["output_folder"]
SAMPLESPATH = config["SAMPLESPATH"]
SOURCE = os.path.dirname(SAMPLESPATH[0])
SCRIPTS = os.path.dirname(config["SCRIPTS"])

import sys

d = {}
SAMPLES = []
for eachSample in SAMPLESPATH:
    basename = os.path.basename(eachSample)
    newBase = basename.replace('.genome-mapped.bam', '')
    SAMPLES.append(newBase)
    d[newBase] = eachSample
    
rule all:
    input:
        expand(os.path.join(OUTDIR, "pieChart", "{sample}" + ".png"), sample=list(d.keys()))
        
rule unmapped_count:
    input:
        bam=lambda wildcards: d[wildcards.sample]
    output:
        readnum=os.path.join(OUTDIR, "unmapped_counts", "{sample}" + ".txt")
    conda:
        "envs/samtools.yaml"
    params:
        num_threads = 1,
        run_time = "4:00:00"
    shell:
        """
        module load samtools;
        samtools view -cf 4 {input.bam} > {output.readnum}
        echo {SOURCE} >> {output.readnum}
        """
        
rule unmapped_bam:
    params:
        N_downsample_reads=config["N_downsample_reads"],
        num_threads = 1,
        run_time = "4:00:00"
    input:
        txt=os.path.join(OUTDIR, "unmapped_counts", "{SAMPLES}.txt")
    output:
        bam=os.path.join(OUTDIR, "unmapped_counts/", "{SAMPLES}_unmapped_downsampled.bam")
    conda:
        "envs/python3.yaml"
    shell:
        """
        module load python3essential;
        python3 {SCRIPTS}/downsampled_seqs.py -input {input.txt} -N_downsample {params.N_downsample_reads} -output {output.bam}
        """

rule unmapped_fasta:
    input:
        bam=os.path.join(OUTDIR, "unmapped_counts", "{SAMPLES}_unmapped_downsampled.bam")
    output:
        fasta=os.path.join(OUTDIR, "unmapped_counts", "{SAMPLES}_unmapped_downsampled.fasta")
    conda:
        "envs/samtools.yaml"
    params:
        num_threads = 1,
        run_time = "4:00:00"
    shell:
        """
        module load samtools;
        samtools fasta {input.bam} > {output.fasta}
        """

rule unmapped_blastn:
    threads: 8
    params:
        error_file = "unmapped_count/blast.err",
        out_file = "unmapped_count/blast.out",
        run_time = "24:00:00",
        memory = "200",
        job_name = "blast",
        DB = config["DB_N"],
        outfmt = 6,
        max_target_seqs = 5,
        max_hsps = 1,
        num_threads = 8
    input:
        fasta=os.path.join(OUTDIR, "unmapped_counts", "{SAMPLES}_unmapped_downsampled.fasta")
    output:
        tsv=os.path.join(OUTDIR, "unmapped_counts", "{SAMPLES}_unmappedblast_downsampled_blastn.tsv")
    conda:
        "envs/blast.yaml"
    shell:
        """
        module load blast;
        blastn -db {params.DB} -query {input.fasta} -out {output.tsv} -outfmt {params.outfmt} -max_target_seqs {params.max_target_seqs} -max_hsps {params.max_hsps} -num_threads {params.num_threads}
        """

rule unmapped_blastx:
    threads: 8
    params:
        error_file = "unmapped_count/blast.err",
        out_file = "unmapped_count/blast.out",
        run_time = "24:00:00",
        memory = "200",
        job_name = "blast",
        DB = config["DB_X"],
        outfmt = 6,
        max_target_seqs = 5,
        max_hsps = 1,
        num_threads = 8
    input:
        fasta=os.path.join(OUTDIR, "unmapped_counts", "{SAMPLES}_unmapped_downsampled.fasta")
    output:
        tsv=os.path.join(OUTDIR, "unmapped_counts", "{SAMPLES}_unmappedblast_downsampled_blastx.tsv")
    conda:
        "envs/diamond.yaml"
    shell:
        """
        module load diamond;
        diamond blastx -d {params.DB} -q {input.fasta} -o {output.tsv} -k {params.max_target_seqs} --threads {threads} --max-hsps {params.max_hsps}
        """

rule unmapped_pie:
    input:
       input1=os.path.join(OUTDIR, "unmapped_counts", "{SAMPLES}_unmappedblast_downsampled_blastn.tsv"),
       input2=os.path.join(OUTDIR, "unmapped_counts", "{SAMPLES}_unmappedblast_downsampled_blastx.tsv")
    output:
        pie=os.path.join(OUTDIR, "pieChart", "{SAMPLES}.png")
    params:
        num_threads = 1,
        run_time = "1:00:00"
    conda:
        "envs/python3.yaml"
    shell:
        """
        module load python3essential;
        python3 {SCRIPTS}/blastresults_piechart.py {input.input1} {input.input2} {output.pie}
        """
