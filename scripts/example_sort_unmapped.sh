#!/bin/bash

#need to fix all the bash locations because these are all currently saved to my personal scratch directories
cd /oasis/tscc/scratch/eczhang/larp6/larp6_GRCh38/results

#extracting ~1million genome-mapped unmapped reads from bam files to convert to fasta file
module load samtools

file_ext1="_unmapped_1mil.bam"
declare -i unmapped_count

loc="/oasis/tscc/scratch/eczhang/larp6/larp6_GRCh38/results/"

for bam in *genome-mapped.bam; do
        #to check sequences before with samtools view -c -f 4 $bam
        echo $bam
        unmapped_count=$(samtools view -c -f 4 $loc$bam)
        echo $unmapped_count

        outfile=/oasis/tscc/scratch/eczhang/larp6_unmapped/${bam:0:14}$file_ext1

        #realistically when converted to python, float division can be easily done
        cd /oasis/tscc/scratch/eczhang
        python3 one_mil_seqs.py $unmapped_count $bam $outfile

        samtools view -c /oasis/tscc/scratch/eczhang/larp6_unmapped/${bam:0:14}$file_ext1;

done

#converting new bam files with ~1million reads into fasta file
file_ext2="_unmapped_1mil.fasta"

cd /oasis/tscc/scratch/eczhang/larp6_unmapped

for bam in *_unmapped_1mil.bam; do
        samtools fasta $bam > /oasis/tscc/scratch/eczhang/larp6_unmapped/${bam:0:14}$file_ext2
done

#use fasta files as query to run against blastn
#set searchpath for blast, use echo $BLASTDB to check
export BLASTDB=/projects/ps-yeolab3/bay001/annotations/nr

export INPUT_TARGET=nt

module load blast

file_ext3="_1mil_unmappedblast.tsv"
for unmapped in *_unmapped_1mil.fasta; do
        export OUTPUT_FASTA=$unmapped
        export OUTPUT_BLAST_RESULT=/oasis/tscc/scratch/eczhang/larp6_unmapped_blast/${unmapped:0:14}$file_ext3

        #run blast
        blastn -db $INPUT_TARGET -query $OUTPUT_FASTA -out $OUTPUT_BLAST_RESULT -outfmt 6 -max_target_seqs 5;
done

#convert each of the blastn result tsv files into their specific pie charts, stored under larp6_unmapped_blast dir
cd /oasis/tscc/scratch/eczhang/larp6_unmapped_blast

for blast_file in *_1mil_unmappedblast.tsv; do
        python3 blastresults_piechart.py $blast_file;
done


echo "Runtime: $SECONDS"
