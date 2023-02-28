configfile: "config.yaml"
import os

SAMPLESPATH = config["SAMPLESPATH"]
SOURCE = os.path.dirname(SAMPLESPATH[0])
d = {}
SAMPLES = []
for eachSample in SAMPLESPATH:
    basename = os.path.basename(eachSample)
    newBase = basename.replace('.genome-mapped.bam', '')
    SAMPLES.append(newBase)
    d[newBase] = eachSample
    
rule all:
    input:
        expand(os.path.join("pieChart", "{sample}" + ".png"), sample=list(d.keys()))
        
rule unmapped_count:
    input:
        bam=lambda wildcards: d[wildcards.sample]
    output:
        readnum=os.path.join("unmapped_counts", "{sample}" + ".txt")
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools view -cf 4 {input.bam} > {output.readnum}
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
        "python3 scripts/downsampled_seqs.py -input {input} -N_downsample {params.N_downsample_reads} -output {output}"

rule unmapped_fasta:
    input:
        "unmapped_counts/{SAMPLES}_unmapped_downsampled.bam"
    output:
        "unmapped_counts/{SAMPLES}_unmapped_downsampled.fasta"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools fasta {input} > {output}"
       

rule unmapped_blastx:
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
	"unmapped_counts/{SAMPLES}_unmappedblast_downsampled_blastx.tsv"
    conda:
	"envs/diamond.yaml"
    shell:
        "diamond blastx -d /projects/ps-yeolab3/bay001/annotations/nr/nr -q {input} -o {output} -k 5"

rule unmapped_pie_blastx:
    input:
	"unmapped_counts/{SAMPLES}_unmappedblast_downsampled_blastx.tsv"
    output:
	"pieChart/{SAMPLES}.png"
    conda:
	"envs/python3.yaml"
    shell:
	"python3 script/blastresults_piechart.py {input} {output}"




