configfile: "config.yaml"

SAMPLES = config["SAMPLES"]
SOURCE = config["SOURCE"]

rule all:
    input:
        expand("pieChart/{sample}.png", sample=SAMPLES)
        
rule unmapped_count:
    input:
        "/oasis/tscc/scratch/eczhang/larp6/larp6_GRCh38/results/{SAMPLES}.genome-mapped.bam"
    output:
        "unmapped_counts/{SAMPLES}.txt"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools view -c -f 4 {input} > {output}
        echo {SOURCE} >> {output}
        """
        
rule unmapped_bam:
    params:
        N_downsample_reads=config["N_downsample_reads"]
    input:
        "unmapped_counts/{SAMPLES}.txt"
    output:
        "unmapped_counts/{SAMPLES}_unmapped_downsampled.bam"
    conda:
        "envs/python3.yaml"
    shell:
        "python3 script/downsampled_seqs.py -input {input} -N_downsample {params.N_downsample_reads} -output {output}"

rule unmapped_fasta:
    input:
        "unmapped_counts/{SAMPLES}_unmapped_downsampled.bam"
    output:
        "unmapped_counts/{SAMPLES}_unmapped_downsampled.fasta"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools fasta {input} > {output}"
        
rule unmapped_blast:
    threads: 8
    params:
        error_file = "unmapped_count/blast.err",
        out_file = "unmapped_count/blast.out",
        run_time = "24:00:00",
        memory = "200",
        job_name = "blast",
        DB = config["DB"],
        outfmt = 6,
        max_target_seqs = 5,
        max_hsps = 1,
        num_threads = 8
    input:
        "unmapped_counts/{SAMPLES}_unmapped_downsampled.fasta"
    output:
        "unmapped_counts/{SAMPLES}_unmappedblast_downsampled.tsv"
    conda:
        "envs/blast.yaml"
    shell:
        "blastn -db {params.DB} -query {input} -out {output} -outfmt {params.num_threads} -max_target_seqs {params.max_target_seqs} -max_hsps {params.max_hsps} -num_threads {params.num_threads}"
        
rule unmapped_pie:
    input:
        "unmapped_counts/{SAMPLES}_unmappedblast_downsampled.tsv"
    output:
        "pieChart/{SAMPLES}.png"
    conda:
        "envs/python3.yaml"
    shell:
        "python3 script/blastresults_piechart.py {input} {output}"
