#select intepreter and change environment to anaconda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
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

#ncbi query functions
def esearch(term, db='gds'):
    """
    Queries NCBI using the esearch utility. GEO ('gds') database is used as default for search term.
    """
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db={db}&term={term}&retmax=5000&usehistory=y'
    response = urllib.request.urlopen(url)
    return response.read()


def get_esummary(esearch_string, db='gds'):
    """
    Parses a http response in XML format to obtain the webenv and querykey tokens.
    Uses NCBI eutils to transform these tokens into web summaries of GEO (db='gds') datasets.
    """
    xmldoc = minidom.parseString(esearch_string)
    try:
        webenv = xmldoc.getElementsByTagName('WebEnv')[0].firstChild.data
        querykey = xmldoc.getElementsByTagName('QueryKey')[0].firstChild.data
        host = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
        params = f'?db={db}&version=2.0&query_key={querykey}&WebEnv={webenv}'
        url = host + params
        response = urllib.request.urlopen(url)
        return response.read()
    except IndexError as e:
        print(f"Unparsable publication string ({e}, search={esearch_string}")
        return ""

#matplotlib to build the piechart
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 5))
fig.subplots_adjust(wspace=0)

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

sseq_count_series = df3['sseqid'].value_counts()
print(sseq_count_series)

#loop through series to determine which elements to remove and add into "other" column
to_remove = []
otherDict = {}
other_count = 0

print(num_seqs)

#taking 1% of number of sequences
for index,values in sseq_count_series.iteritems():
        if(values < (0.001)*(num_seqs)):
                to_remove.append(index)
                otherDict[values] = index
                other_count += values

otherAnswer = []
increment = 0
for key,value in sorted(otherDict.items(),reverse=True):
        if increment <= 2:
                otherAnswer.append((key,value))
                increment+=1
        else:
             	break

#bar chart parameters
gene_ratios = []
gene_labels = []
gene_labels1 = []
bottom = 1
width = .2

for idx in range(len(otherAnswer)):
        #otherAnswer[idx][0] = otherAnswer[idx][0]/other_count
        gene_ratios.append(otherAnswer[idx][0]/other_count)
        gene_labels1.append(otherAnswer[idx][1])

#replace index sseqids with ncbi name
for idx in range(len(gene_labels1)):
        term = str(gene_labels1[idx])
        print(term)
        esearch_string = esearch(term=term, db='protein')
        result = get_esummary(esearch_string=esearch_string, db='protein')
        result = xmltodict.parse(result)
        sseq_name = result['eSummaryResult']['DocumentSummarySet']['DocumentSummary']['Title']
        gene_labels.append(sseq_name)


#remove these elements from series
sseq_count_series2 = sseq_count_series.drop(to_remove)

sseq_name_list = []
sseqid = []
#replace index sseqids with ncbi name
for index,values in sseq_count_series2.iteritems():
        sseqid.append(index)
        term = str(index)
        esearch_string = esearch(term=term, db='protein')
        result = get_esummary(esearch_string=esearch_string, db='protein')
        result = xmltodict.parse(result)
        sseq_name = result['eSummaryResult']['DocumentSummarySet']['DocumentSummary']['Title']
        sseq_name_list.append(sseq_name)
        #sseq_count_series2.rename(index={index:sseq_name})
#replace the sseqids with ncbi query names
replacements = {sseqid:sseq_name_list for sseqid, sseq_name_list in zip(sseqid, sseq_name_list)}
sseq_count_series3 = sseq_count_series2.rename(replacements)
#sseq_count_series2 = sseq_count_series2.rename(index=dict(zip(sseq_name_list,sseqid)))

#generate new pandas series with new element to concatenate with old series
d = {'Other':other_count}
ser = pd.Series(data=d, index=['Other'])

#append new element
sseq_count = sseq_count_series3.append(ser)
#print(sseq_count)

#value_counts returns a pandas series so convert to a data frame
sseq_count_df = pd.DataFrame({'sseqid':sseq_count.index, 'count':sseq_count.values})
print(sseq_count_df)

count = sseq_count_df['count']
sseq = sseq_count_df['sseqid']

wedges, *_ = ax1.pie(count, labels = sseq, colors=None,autopct='%1.1f%%',startangle=45,
        wedgeprops={"linewidth": 1, "edgecolor": "white"})

# Adding from the top matches the legend.
for j, (height, label) in enumerate(reversed([*zip(gene_ratios, gene_labels)])):
    bottom -= height
    bc = ax2.bar(0, height, width, bottom=bottom, color='C0', label=label,
                 alpha=0.1 + 0.25 * j)
    ax2.bar_label(bc, labels=[f"{height:.0%}"], label_type='center')

ax2.set_title('Other unmapped reads')
ax2.legend()
ax2.axis('off')
ax2.set_xlim(- 2.5 * width, 2.5 * width)

# use ConnectionPatch to draw lines between the two plots
theta1, theta2 = wedges[0].theta1, wedges[0].theta2
center, r = wedges[0].center, wedges[0].r
bar_height = sum(gene_ratios)

# draw top connecting line
x = r * np.cos(np.pi / 180 * theta2) + center[0]
y = r * np.sin(np.pi / 180 * theta2) + center[1]
con = ConnectionPatch(xyA=(-width / 2, bar_height), coordsA=ax2.transData,
                      xyB=(x, y), coordsB=ax1.transData)
con.set_color([0, 0, 0])
con.set_linewidth(4)
ax2.add_artist(con)

# draw bottom connecting line
x = r * np.cos(np.pi / 180 * theta1) + center[0]
y = r * np.sin(np.pi / 180 * theta1) + center[1]
con = ConnectionPatch(xyA=(-width / 2, 0), coordsA=ax2.transData,
                      xyB=(x, y), coordsB=ax1.transData)
con.set_color([0, 0, 0])
ax2.add_artist(con)
con.set_linewidth(4)

plt.title('Blastx Unmapped Sequences Summary ' + pieName)
plt.show(block=True)
plt.savefig(piePath,format='png',bbox_inches='tight')
