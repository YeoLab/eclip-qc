import sys
import os
import argparse

parser = argparse.ArgumentParser(description='DownSample samtools')
parser.add_argument("-N_downsample", "--N_downsample_reads", help="Number of reads to downsample", type=int)
parser.add_argument("-input", "--input_file", help="input Number of reads file", type=str)
parser.add_argument("-output", "--output_file", help="output file", type=str)

args = parser.parse_args()

#number of unmapped sequences in each file
unmapped_count = args.input_file
N_downsample_reads = args.N_downsample_reads
out_file = args.output_file

file1 = open(unmapped_count, "r+")
sample = os.path.splitext(unmapped_count)[0]
sample2 = os.path.split(sample)
bam_file = sample2[1] + ".genome-mapped.bam"
Lines = file1.readlines()
count = 0
unmappedCount = ""
sourcePath = ""
for line in Lines:
    count += 1
    if count == 1:
        unmappedCount = int(line.replace("\n", ""))
    else:
        sourcePath = line.replace("\n", "")

#get float value to get N downsample reads, if file has less than N downsample reads, just use all of them
if unmappedCount < N_downsample_reads:
        N_downsample_unmapped=1
else:
        N_downsample_unmapped = N_downsample_reads/unmappedCount

#output is samtools command
cmd = "samtools view -f 4 -s " + str(N_downsample_unmapped) + " " + sourcePath + "/" +  bam_file + ">" + out_file

print(cmd)
os.system(cmd)
