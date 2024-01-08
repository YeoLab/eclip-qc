configfile: "config.yaml"
import os

if not os.path.exists("stderr"): os.makedirs("stderr")
if not os.path.exists("stdout"): os.makedirs("stdout")

SAMPLESPATH = config["SAMPLESPATH"]
SOURCE = os.path.dirname(SAMPLESPATH[0])
d = {}
SAMPLES = []
for eachSample in SAMPLESPATH:
    basename = os.path.basename(eachSample)
    SAMPLES.append(basename)
    d[basename] = eachSample
    
rule all:
    params:
        error_file = "stderr/all.err",
        out_file = "stdout/all.out",
        run_time = "00:04:00",
        memory = "40",
        job_name = "all"
    threads: 1
    input:
        expand(os.path.join("pieChart", "{sample}" + ".png"), sample=list(d.keys()))

rule unmapped_count:
    input:
        bam=lambda wildcards: d[wildcards.sample]
    output:
        readnum=os.path.join("unmapped_counts", "{sample}" + ".txt")
    threads: 1
    params:
        error_file = lambda wildcards: "stderr/{}.count.err".format(wildcards.sample),
        out_file = lambda wildcards: "stdout/{}.count.out".format(wildcards.sample),
        run_time = "00:10:00",
        memory = "4000",
        job_name = "count"
    container:
        "docker://howardxu520/skipper:samtools_1.17"
    shell:
        """
        samtools view -cf 4 {input.bam} > {output.readnum}
        echo {SOURCE} >> {output}
        """
        
rule unmapped_bam:
    input:
        info="unmapped_counts/{SAMPLES}.txt"
    output:
        "unmapped_counts/{SAMPLES}_unmapped_downsampled.bam"
    threads: 1
    params:
        N_downsample_reads=config["N_downsample_reads"],
        error_file = "stderr/{SAMPLES}.bam.err",
        out_file = "stdout/{SAMPLES}.bam.out",
        run_time = "00:20:00",
        memory = "4000",
        job_name = "bam"
    container:
        "docker://howardxu520/eclip-qc:python"
    shell:
        "python3 scripts/downsampled_seqs.py -input {input.info} -N_downsample {params.N_downsample_reads} -output {output}"

rule unmapped_fasta:
    input:
        bam="unmapped_counts/{SAMPLES}_unmapped_downsampled.bam"
    output:
        "unmapped_counts/{SAMPLES}_unmapped_downsampled.fasta"
    threads: 1
    params:
        error_file = "stderr/{SAMPLES}.fasta.err",
        out_file = "stdout/{SAMPLES}.fasta.out",
        run_time = "00:20:00",
        memory = "40000",
        job_name = "fasta"
    container:
        "docker://howardxu520/skipper:samtools_1.17"
    shell:
        "samtools fasta {input.bam} > {output}"

rule unmapped_blast_synthetic:
    threads: 4
    params:
        error_file = "stderr/{SAMPLES}.blast_synthetic.err",
        out_file = "stdout/{SAMPLES}.blast_synthetic.out",
        run_time = "10:00:00",
        memory = "40000",
        job_name = "blast_synthetic",
        DB = config["DB_synthetic"],
        outfmt = 6,
        max_target_seqs = 5,
        max_hsps = 1
    input:
        fasta="unmapped_counts/{SAMPLES}_unmapped_downsampled.fasta"
    output:
        "unmapped_counts/{SAMPLES}_unmappedblast_downsampled_blast_synthetic.tsv"
    container:
        "docker://howardxu520/eclip-qc:blast_2.13.0"
    shell:
        "blastn -db {params.DB} -query {input.fasta} -out {output} -outfmt {params.outfmt} -max_target_seqs {params.max_target_seqs} -max_hsps {params.max_hsps} -num_threads {threads}"

rule filter_synthetic:
    input:
       blast="unmapped_counts/{SAMPLES}_unmappedblast_downsampled_blast_synthetic.tsv",
       fasta="unmapped_counts/{SAMPLES}_unmapped_downsampled.fasta"
    output:
        "unmapped_counts/{SAMPLES}_unmapped_downsampled_filtered.fasta"
    threads: 1
    params:
        error_file = "stderr/{SAMPLES}.synthetic.err",
        out_file = "stdout/{SAMPLES}.synthetic.out",
        run_time = "10:00:00",
        memory = "400",
        job_name = "synthetic"  
    container:
        "docker://howardxu520/eclip-qc:python"
    shell:
        "python3 scripts/filter_synthetic.py {input.blast} {input.fasta} {output}"

rule unmapped_blastn:
    threads: 4
    input:
        fasta="unmapped_counts/{SAMPLES}_unmapped_downsampled_filtered.fasta"
    params:
        error_file = "stderr/{SAMPLES}.blastn.err",
        out_file = "stdout/{SAMPLES}.blastn.out",
        run_time = "10:00:00",
        memory = "40000",
        job_name = "blastn",
        DB = config["DB_N"],
        outfmt = 6,
        max_target_seqs = 5,
        max_hsps = 1
    output:
        "unmapped_counts/{SAMPLES}_unmappedblast_downsampled_blastn.tsv"
    container:
        "docker://howardxu520/eclip-qc:blast_2.13.0"
    shell:
        "blastn -db {params.DB} -query {input.fasta} -out {output} -outfmt {params.outfmt} -max_target_seqs {params.max_target_seqs} -max_hsps {params.max_hsps} -num_threads {threads}"

rule unmapped_blastx:
    threads: 4
    input:
        fasta="unmapped_counts/{SAMPLES}_unmapped_downsampled_filtered.fasta"
    params:
        error_file = "stderr/{SAMPLES}.blastx.err",
        out_file = "stdout/{SAMPLES}.blastx.out",
        run_time = "10:00:00",
        memory = "40000",
        job_name = "blastx",
        DB = config["DB_X"],
        outfmt = 6,
        max_target_seqs = 5,
        max_hsps = 1
    output:
        "unmapped_counts/{SAMPLES}_unmappedblast_downsampled_blastx.tsv"
    container:
        "docker://howardxu520/eclip-qc:diamond"
    shell:
        "diamond blastx -d {params.DB} -q {input.fasta} -o {output} -k {params.max_target_seqs} --threads {threads} --max-hsps {params.max_hsps}"

rule unmapped_pie:
    threads: 4
    input:
       input1="unmapped_counts/{SAMPLES}_unmappedblast_downsampled_blastn.tsv",
       input2="unmapped_counts/{SAMPLES}_unmappedblast_downsampled_blastx.tsv",
       input3="unmapped_counts/{SAMPLES}_unmappedblast_downsampled_blast_synthetic.tsv",
    output:
        "pieChart/{SAMPLES}.png"
    params:
        error_file = "stderr/{SAMPLES}.pie.err",
        out_file = "stdout/{SAMPLES}.pie.out",
        memory = "4000",
        run_time = "10:00:00",
        job_name = "pie" 
    container:
        "docker://howardxu520/eclip-qc:python"
    shell:
        "python3 scripts/blastresults_piechart.py {input.input1} {input.input2} {output}"
