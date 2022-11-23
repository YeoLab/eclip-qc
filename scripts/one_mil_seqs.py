import sys
import os

#number of unmapped sequences in each file
unmapped_count = int(sys.argv[1])
bam_file = sys.argv[2]
out_file = sys.argv[3]

one_mil = 1000000

#get float value to get 1 mil reads, if file has less than a million reads, just use all of them
if unmapped_count < one_mil:
        one_mil_unmapped=1
else:
     	one_mil_unmapped = one_mil/unmapped_count

#output is samtools command
cmd = "samtools view -f 4 -s " + str(one_mil_unmapped) + " /oasis/tscc/scratch/eczhang/larp6/larp6_GRCh38/results/" + bam_file + ">" + out_file

print(cmd)
os.system(cmd)