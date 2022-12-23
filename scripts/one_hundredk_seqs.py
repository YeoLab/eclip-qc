import sys
import os

#number of unmapped sequences in each file
unmapped_count = sys.argv[1]
out_file = sys.argv[2]

file1 = open(unmapped_count, "r+")
sample = os.path.splitext(unmapped_count)[0]
sample2 = os.path.split(sample)
bam_file = sample2[1] + ".genome-mapped.bam"
unmappedCount = int(file1.read())
sourcePath = file1.read()


one_hundredk = 100000

#get float value to get 100,000 reads, if file has less than 100,000 reads, just use all of them
if unmappedCount < one_hundredk:
        one_hundredk_unmapped=1
else:
        one_hundredk_unmapped = one_hundredk/unmappedCount

#output is samtools command
cmd = "samtools view -f 4 -s " + str(one_hundredk_unmapped) + " " + sourcePath + "/" +  bam_file + ">" + out_file

print(cmd)
os.system(cmd)
