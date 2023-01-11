import sys
import os

#number of unmapped sequences in each file
unmapped_count = sys.argv[1]
#bam_file = sys.argv[2]
out_file = sys.argv[2]

file1 = open(unmapped_count, "r+")
sample = os.path.splitext(unmapped_count)[0]
sample2 = os.path.split(sample)
bam_file = sample2[1] + ".genome-mapped.bam"
unmappedCount = int(file1.read())


hundred_thousand = 100000

#get float value to get 1 mil reads, if file has less than a million reads, just use all of them
if unmappedCount < hundred_thousand:
        hundred_thousand_unmapped=1
else:
        hundred_thousand_unmapped = hundred_thousand/unmappedCount

#output is samtools command
cmd = "samtools view -f 4 -s " + str(hundred_thousand_unmapped) + " /oasis/tscc/scratch/eczhang/larp6/larp6_GRCh38/results/" + bam_file + ">" + out_file

print(cmd)
os.system(cmd)