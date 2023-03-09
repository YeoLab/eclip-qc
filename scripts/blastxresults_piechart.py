#select intepreter and change environment to anaconda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import urllib
import xmltodict
from xml.dom import minidom


#handles command line input parameter
import sys
import os

blast_tsv_file = sys.argv[1]
piePath = sys.argv[2]
piePathSample = os.path.splitext(piePath)[0]
piePathSample2 = os.path.split(piePathSample)
pieName = piePathSample2[1]

#matplotlib to build the piechart
fig, ax = plt.subplots()

df = pd.read_csv(blast_tsv_file, header=None, sep='\t')
num_seqs = df.size
df.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']

df2 = df[['qseqid','sseqid','evalue']].copy()

#blast output default is by best hit/lowest e value score, therefore add the first sseqid result for each qseqid
sseqid_list = []
qseqid_list = []
for idx in df.index:
        if df['qseqid'][idx] not in qseqid_list:
                qseqid_list.append(df['qseqid'][idx])
                sseqid_list.append(df['sseqid'][idx])
sseqid_np = np.array(sseqid_list)

df3 = pd.DataFrame(data = sseqid_np, columns=['sseqid'])
print(df3)

#get the frequency of each blast result
sseq_count_series = df3['sseqid'].value_counts()
print(sseq_count_series)

#loop through series to determine which elements to remove and add into "other" column
to_remove = []
other_count = 0

print(num_seqs)

#taking 1% of number of sequences
for index,values in sseq_count_series.iteritems():
        if(values < (0.001)*(num_seqs)):
                to_remove.append(index)
                other_count += values


#remove these elements from series
sseq_count_series2 = sseq_count_series.drop(to_remove)

#generate new pandas series with new element to concatenate with old series
d = {'Other':other_count}
ser = pd.Series(data=d, index=['Other'])

#append 'Other' element
sseq_count = sseq_count_series2.append(ser)
#print(sseq_count)

#value_counts returns a pandas series so convert to a data frame
sseq_count_df = pd.DataFrame({'sseqid':sseq_count.index, 'count':sseq_count.values})
print(sseq_count_df)

#add respective columns
count = sseq_count_df['count']
sseq = sseq_count_df['sseqid']

ax.pie(count, labels = sseq, colors=None,autopct='%1.1f%%',startangle=45,
        wedgeprops={"linewidth": 1, "edgecolor": "white"})

plt.title('Unmapped Sequences Summary ' + pieName)
plt.show(block=True)
plt.savefig(piePath,format='png',bbox_inches='tight')
