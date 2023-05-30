#select intepreter and change environment to anaconda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import urllib
import xmltodict
from xml.dom import minidom
import time


#handles command line input parameter
import sys
import os

# testing file is from 10,000 reads, so mapping percentage is expected to be out of 10,000
blast_tsv_file = sys.argv[1]
#blast_tsv_file2 = sys.argv[2]


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
fig, ((ax1, ax2, ax5), (ax3, ax4, ax6)) = plt.subplots(2, 3, figsize=(9, 5))
fig.tight_layout(h_pad=2)
plt.subplots_adjust(bottom=0, right=1.5, top=1.5)
ax1.set_title('blastn piechart')
ax2.set_title('blastn other')
ax3.set_title('blastx piechart')
ax4.set_title('blastx other')
ax5.set_title('blastn percentage')
ax6.set_title('blastx percentage')


filterlist = ['bacter', 'bacterium', 'strain', 'Homo sapiens', 'Mus musculus', 'virus']


# Read the blast n file
df = pd.read_csv(blast_tsv_file, header=None, sep='\t')
df.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']

df2 = df[['qseqid','sseqid','evalue']].copy()
print(df2)
# Read the blast x file
#df4 = pd.read_csv(blast_tsv_file2, header=None, sep='\t')
#df4.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']

#df3 = df4[['qseqid','sseqid','evalue']].copy()
#print(df3)

to_delete1 = []
qseqid_checker1 = []
for index, row in df2.iterrows():
    checker = 0
    term = str(row['sseqid'])
    print(term)
    esearch_string = esearch(term=term, db='nucleotide')
    time.sleep(0.1)
    result = get_esummary(esearch_string=esearch_string, db='nucleotide')
    result = xmltodict.parse(result)
    sseq_name = result['eSummaryResult']['DocumentSummarySet']['DocumentSummary']['Title']
    for i in range(len(filterlist)):
        if filterlist[i] in sseq_name:
            checker += 1
            qseqid_checker1.append(row['qseqid']) 
    if checker == 0:
        to_delete1.append(index)
    #if none of them have species of interest, keep the top one
    if ((index+1)%5 == 0) and (row['qseqid'] not in qseqid_checker1):
        to_delete1.remove(index-4)

for idx in range(len(to_delete1)):
    df2 = df2.drop(to_delete1[idx])
print(df2)

