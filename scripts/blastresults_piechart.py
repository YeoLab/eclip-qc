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

larp6_file = sys.argv[1]
pieName = sys.argv[2]

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

#matplotlib to build the piechart, use pandas to create dataframe from blast output tsv
fig, ax = plt.subplots()

df = pd.read_csv(larp6_file, header=None, sep='\t')
num_seqs = df.size
df.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']

#matplotlib to build the piechart
fig, ax = plt.subplots()

df = pd.read_csv(larp6_file, header=None, sep='\t')
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
other_count = 0

print(num_seqs)

#taking 1% of number of sequences
for index,values in sseq_count_series.iteritems():
        if(values < (0.01)*(num_seqs)):
                to_remove.append(index)
                other_count += values


#remove these elements from series
sseq_count_series2 = sseq_count_series.drop(to_remove)

sseq_name_list = []
sseqid = []
#replace index sseqids with ncbi name
for index,values in sseq_count_series2.iteritems():
        sseqid.append(index)
        term = str(index)
        esearch_string = esearch(term=term, db='nucleotide')
        result = get_esummary(esearch_string=esearch_string, db='nucleotide')
        sseq_count_series = df['sseqid'].value_counts()
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

ax.pie(count, labels = sseq, colors=None,autopct='%1.1f%%',startangle=45,
wedgeprops={"linewidth": 1, "edgecolor": "white"})

plt.title('unmapped sequences summary LARP6')
plt.show(block=True)
plt.savefig(pieName,format='png',bbox_inches='tight')
