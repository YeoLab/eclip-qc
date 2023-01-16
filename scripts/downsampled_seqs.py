import sys
import os

#number of unmapped sequences in each file
unmapped_count = sys.argv[1]
N_downsample_reads = int(sys.argv[2])
out_file = sys.argv[3]

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
