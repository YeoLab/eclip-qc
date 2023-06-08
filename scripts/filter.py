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
blast_tsv_file2 = sys.argv[2]
#piePath = sys.argv[3]
#piePathSample = os.path.splitext(piePath)[0]
#piePathSample2 = os.path.split(piePathSample)
#pieName = piePathSample2[1]
filterlist = ['bacter', 'bacterium', 'strain', 'Homo sapiens', 'Mus musculus', 'virus']

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


# Read the blast n file
df = pd.read_csv(blast_tsv_file, header=None, sep='\t')
df.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']

df2 = df[['qseqid','sseqid','pident']].copy()

# Read the blast x file
df3 = pd.read_csv(blast_tsv_file2, header=None, sep='\t')
df3.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']
df4 = df3[['qseqid','sseqid','pident']].copy()

dictionary1 = {}
for eachid,eachp in zip(df2['sseqid'], df2['pident']):
    if eachid not in dictionary1:
        dictionary1[eachid] = [eachp, 1]
    elif float(eachp) > dictionary1[eachid][0]:
        dictionary1[eachid] = [float(eachp), 1]
    elif float(eachp) == dictionary1[eachid][0]:
        dictionary1[eachid][1] += 1
        
dictionary1 = sorted(dictionary1.items(), key=lambda x:x[1], reverse=True)
ncbilist1 = []

for each in dictionary1:
    if float(each[1][0]) == float(100):
        checker = 0
        term = str(each[0])
        #print(term)
        esearch_string = esearch(term=term, db='nucleotide')
        time.sleep(0.1)
        result = get_esummary(esearch_string=esearch_string, db='nucleotide')
        #print(result)
        result = xmltodict.parse(result)
        for i in range(len(filterlist)):
            if filterlist[i] in result:
                ncbilist1.append(result)

blastnSort1 = pd.DataFrame(ncbilist1, columns=['Name'])
blastnSort2 = blastnSort1['Name'].value_counts()
print(blastnSort2)

dictionary2 = {}
for eachid,eachp in zip(df4['sseqid'], df4['pident']):
    if eachid not in dictionary2:
        dictionary2[eachid] = [eachp, 1]
    elif float(eachp) > dictionary2[eachid][0]:
        dictionary2[eachid] = [float(eachp), 1]
    elif float(eachp) == dictionary2[eachid][0]:
        dictionary2[eachid][1] += 1
        
dictionary2 = sorted(dictionary2.items(), key=lambda x:x[1], reverse=True)
ncbilist2 = []

for each in dictionary2:
    if float(each[1][0]) == float(100):
        checker = 0
        term = str(each[0])
        #print(term)
        esearch_string = esearch(term=term, db='protein')
        time.sleep(0.1)
        result = get_esummary(esearch_string=esearch_string, db='protein')
        #print(result)
        result = xmltodict.parse(result)
        for i in range(len(filterlist)):
            if filterlist[i] in result:
                ncbilist2.append(result)

blastxSort1 = pd.DataFrame(ncbilist2, columns=['Name'])
blastxSort2 = blastxSort1['Name'].value_counts()
print(blastxSort2)
