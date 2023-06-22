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
#piePath = sys.argv[3]
#piePathSample = os.path.splitext(piePath)[0]
#piePathSample2 = os.path.split(piePathSample)
#pieName = piePathSample2[1]
filterlist = ['bacter', 'bacterium', 'strain', 'Homo sapiens', 'Mus musculus', 'virus']

def esearch(term, db='gds', api_key='9e09d5d38c680a8358426f7fac6d154b4f08'):
    """
    Queries NCBI using the esearch utility. GEO ('gds') database is used as default for search term.
    """
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db={db}&term={term}&retmax=5000&usehistory=y&api_key={api_key}'
    try:
        response = urllib.request.urlopen(url)
    except urllib.error.HTTPError:
        return 'No_Species'        
    return response.read()


def get_esummary(esearch_string, db='gds', api_key='9e09d5d38c680a8358426f7fac6d154b4f08'):
    """
    Parses a http response in XML format to obtain the webenv and querykey tokens.
    Uses NCBI eutils to transform these tokens into web summaries of GEO (db='gds') datasets.
    """
    xmldoc = minidom.parseString(esearch_string)
    try:
        webenv = xmldoc.getElementsByTagName('WebEnv')[0].firstChild.data
        querykey = xmldoc.getElementsByTagName('QueryKey')[0].firstChild.data
        host = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
        params = f'?db={db}&version=2.0&query_key={querykey}&WebEnv={webenv}&api_key={api_key}'
        url = host + params
        response = urllib.request.urlopen(url)
        return response.read()
    except IndexError as e:
        print(f"Unparsable publication string ({e}, search={esearch_string}")
        return ""


def ncbi_parse(result):
    try:
        return result['eSummaryResult']['DocumentSummarySet']['DocumentSummary']['Title']
    except TypeError:
        return 'No_Species'
    except KeyError:
        return 'No_Species'

# Read the blast n file
df = pd.read_csv(blast_tsv_file, header=None, sep='\t')
df.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']
original_numseqs1 = df.size
df2 = df[['qseqid','sseqid','pident']].copy()

# Read the blast x file
#df3 = pd.read_csv(blast_tsv_file2, header=None, sep='\t')
#df3.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']
#df4 = df3[['qseqid','sseqid','pident']].copy()

dictionary1 = {}
for eachid,eachp in zip(df2['sseqid'], df2['pident']):
    if eachid not in dictionary1:
        dictionary1[eachid] = [eachp, 1]
    elif float(eachp) > dictionary1[eachid][0]:
        dictionary1[eachid] = [float(eachp), 1]
    elif float(eachp) == dictionary1[eachid][0]:
        dictionary1[eachid][1] += 1

dictionary1 = sorted(dictionary1.items(), key=lambda x:x[1], reverse=True)

checkqseqid = []
dict1 = {}
for eachqid,eachsid,eachp in zip(df2['qseqid'],df2['sseqid'],df2['pident']):
    if eachqid not in checkqseqid:
        #first instance of qseqid encountered, add top hit
        checkqseqid.append(eachqid)
        if eachsid not in dict1:
            dict1[eachsid] = 1
        else:
           dict1[eachsid] += 1
        continue
    #not the first hit for the qseqid
    if eachp == 100:
        if eachsid not in dict1:
            dict1[eachsid] = 1
        else:
            dict1[eachsid] += 1
    else:
        #ignore if the score is less than 100 until next qseqid is reached
        continue

qseqidDf = pd.DataFrame(list(dict1.items()), columns=['qseqid','frequency'])

#filter by frequency if frequency is too small TODO

toDrop = []
#filter by keywords!
for index,row in qseqidDf.iterrows():
    term = str(row['qseqid'])
    #print(term)
    esearch_string = esearch(term=term, db='nucleotide')
    if (esearch_string == 'No_Species'):
        continue
    time.sleep(0.5)
    result = get_esummary(esearch_string=esearch_string, db='nucleotide')
    result = xmltodict.parse(result)
    #print(result)
    sseq_name = ncbi_parse(result)
    #print(sseq_name)
    for index,row in qseqidDf.iterrows():
    found = False
    for i in range(len(filterlist)):
        if filterlist[i] in row['qseqid']:
            found = True
    if found == False:
        toDrop.append(index)

for i in range(len(toDrop)):
    qseqidDf = qseqidDf.drop(toDrop[i])
#qseqidDf

num_seqs = 0 
for index,row in qseqidDf.iterrows():
    num_seqs += int(row['frequency'])

to_remove = []
other_count = 0
for index,row in qseqidDf.iterrows():
    if(row['frequency'] < (0.001)*(num_seqs)):
        to_remove.append(index)
        other_count += int(row['frequency'])

ncbiDf = qseqidDf.drop(index=to_remove)
ncbiDf = ncbiDf.rename(columns{'qseqid':'name'})

ncbiDf.loc[len(ncbiDf)] = ['Other', other_count]
ncbiDf.loc[len(ncbiDf)] = ['Nonrelevant species', original_numseqs - num_seqs]

#blastnSort1 = pd.DataFrame(ncbilist1, columns=['Name'])
#blastnSort2 = blastnSort1['Name'].value_counts()
#print(blastnSort2)

