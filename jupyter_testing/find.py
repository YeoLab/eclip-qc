#select intepreter and change environment to anaconda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import urllib
import xmltodict
from xml.dom import minidom
import time
import yaml
import traceback
import urllib.error

#handles command line input parameter
import sys
import os

with open("/home/s5xu/projects/eclip-qc/config.yaml", 'r') as stream:
    config = yaml.safe_load(stream)
    

# testing file is from 10,000 reads, so mapping percentage is expected to be out of 10,000
blastn_file = 'LARP6.TGFb_IN1.umi.r1.fq_unmappedblast_downsampled_blastn.tsv'
blastx_file = 'LARP6.TGFb_IN1.umi.r1.fq_unmappedblast_downsampled_blastx.tsv'
# piePath = sys.argv[3]
# piePathSample = os.path.splitext(piePath)[0]
# piePathSample2 = os.path.split(piePathSample)
# pieName = piePathSample2[1]


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
        try:
            response = urllib.request.urlopen(url)
        except urllib.error.HTTPError as e:
            print("http error")
            return ""
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


# Creating plots for blastn

# Read the blast n file
origin_blastn = pd.read_csv(blastn_file, header=None, sep='\t')
origin_blastn.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']

origin_blastn2 = origin_blastn[['qseqid','sseqid','pident']].copy()
original_blastn_numseqs = len(origin_blastn)

all_blastn_freq_dict = {}
for eachid,eachp in zip(origin_blastn2['sseqid'], origin_blastn2['pident']):
    if eachid not in all_blastn_freq_dict:
        all_blastn_freq_dict[eachid] = [eachp, 1]
    elif float(eachp) > all_blastn_freq_dict[eachid][0]:
        all_blastn_freq_dict[eachid] = [float(eachp), 1]
    elif float(eachp) == all_blastn_freq_dict[eachid][0]:
        all_blastn_freq_dict[eachid][1] += 1

all_blastn_freq_dict = sorted(all_blastn_freq_dict.items(), key=lambda x:x[1], reverse=True)

checkqseqid_n = []
top_blastn_freq_dict = {}
for eachqid,eachsid,eachp in zip(origin_blastn2['qseqid'],origin_blastn2['sseqid'],origin_blastn2['pident']):
    if eachqid not in checkqseqid_n:
        #first instance of qseqid encountered, add top hit
        checkqseqid_n.append(eachqid)
        if eachsid not in top_blastn_freq_dict:
            top_blastn_freq_dict[eachsid] = 1
        else:
            top_blastn_freq_dict[eachsid] += 1
        continue
    #not the first hit for the qseqid
    if eachp == 100:
        if eachsid not in top_blastn_freq_dict:
            top_blastn_freq_dict[eachsid] = 1
        else:
            top_blastn_freq_dict[eachsid] += 1
    else:
        #ignore if the score is less than 100 until next qseqid is reached
        continue

top_blastn_freq_table = pd.DataFrame(list(top_blastn_freq_dict.items()), columns=['qseqid','frequency'])

try:  
    #filter by frequency if frequency is too small TODO
    small_hits_interests_of_blastn = []
    #filter by keywords!
    for index,row in top_blastn_freq_table.iterrows():
        term = str(row['qseqid'])
        print(str(index))
        #print(term)
        esearch_string = esearch(term=term, db='nucleotide')
        if (esearch_string == 'No_Species'):
            continue
        time.sleep(0.5)
        result = get_esummary(esearch_string=esearch_string, db='nucleotide')
        if result == "":
            continue
        result = xmltodict.parse(result)
        #print(result)
        sseq_name = ncbi_parse(result)
        #print(sseq_name)
        found = False
        for i in range(len(filterlist)):
            if filterlist[i] in row['qseqid']:
                top_blastn_freq_table.loc[index, 'qseqid'] = sseq_name
                found = True
        if found == False:
            small_hits_interests_of_blastn.append(index)
    for i in range(len(small_hits_interests_of_blastn)):
        top_blastn_freq_table = top_blastn_freq_table.drop(small_hits_interests_of_blastn[i])
    #top_blastn_freq_table

except Exception as e:
    error_message = traceback.format_exc()  # Get the error message as a string
    with open('error_log.txt', 'w') as f:
        f.write(error_message)