SAMPLES = ["LARP6.CTRL_IN1.umi.r1.fq", "LARP6.CTRL_IN2.umi.r1.fq", "LARP6.CTRL_IP1.umi.r1.fq", "LARP6.CTRL_IP2.umi.r1.fq", "LARP6.TGFb_IN1.umi.r1.fq", "LARP6.TGFb_IN2.umi.r1.fq", "LARP6.TGFb_IP1.umi.r1.fq", "LARP6.TGFb_IP2.umi.r1.fq"]

configfile: "config.yaml"

rule all:
    input:
        expand("pieChart/{sample}.png", sample=SAMPLES)
        
rule unmapped_count:
    input:
        "config["SOURCE"]/{SAMPLES}.genome-mapped.bam"
    output:
        "unmapped_counts/{SAMPLES}.txt"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools view -c -f 4 {input} > {output}
        echo config["SOURCE"] >> {output}
        """
        
rule unmapped_bam:
    input:
        "unmapped_counts/{SAMPLES}.txt"
    output:
        "unmapped_counts/{SAMPLES}_unmapped_downsampled.bam"
    conda:
        "envs/python3.yaml"
    shell:
        "python3 script/one_hundredk_seqs.py {input} {output}"
 
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
        num_threads=8
    input:
        "unmapped_counts/{SAMPLES}_unmapped_downsampled.fasta"
    output:
        "unmapped_counts/{SAMPLES}_unmappedblast_downsampled.tsv"
    conda:
        "envs/blast.yaml"
    shell:
        "blastn -db config["DB"] -query {input} -out {output} -outfmt 6 -max_target_seqs 5 -max_hsps 1 -num_threads {params.num_threads}"
        
rule unmapped_pie:
    input:
        "unmapped_counts/{SAMPLES}_unmappedblast_downsampled.tsv"
    output:
        "pieChart/{SAMPLES}.png"
    conda:
        "envs/python3.yaml"
    shell:
        "python3 script/blastresults_piechart.py {input} {output}"
