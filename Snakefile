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
        expand(os.path.join("unmapped_counts", "{sample}" + "_blastx.tsv"), sample=list(d.keys()))

rule diamond_blastx:
    input:
        fname_fastq = "{sample}_unmapped_downsampled.fasta",
        fname_db = "/projects/ps-yeolab3/bay001/annotations/nr/nr.dmnd"
    output:
        fname = "{sample}_blastx.tsv"
    threads: 8
    wrapper:
        "v1.23.4-17-gebb541ee/bio/diamond/blastx"
    
