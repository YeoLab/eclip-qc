#select intepreter and change environment to anaconda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import urllib
import xmltodict
from xml.dom import minidom
from xml.parsers.expat import ExpatError
import time
import yaml


#handles command line input parameter
import sys
import os

with open("/home/s5xu/projects/eclip-qc/config.yaml", 'r') as stream:
    config = yaml.safe_load(stream)

# testing file is from 10,000 reads, so mapping percentage is expected to be out of 10,000
blastn_file = sys.argv[1]
blastx_file = sys.argv[2]
piePath = sys.argv[3]
piePathSample = os.path.splitext(piePath)[0]
piePathSample2 = os.path.split(piePathSample)
pieName = piePathSample2[1]

filterlist = ['bacteria', 'bacterium', 'strain', 'Homo sapiens', 'Mus musculus', 'virus']

def esearch(term, db='gds', api_key='9e09d5d38c680a8358426f7fac6d154b4f08'):
    """
    Queries NCBI using the esearch utility. GEO ('gds') database is used as default for search term.
    """
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db={db}&term={term}&retmax=5000&usehistory=y&api_key={api_key}'
    try:
        response = urllib.request.urlopen(url)
    except (urllib.error.HTTPError, urllib.error.URLError):
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
        except (urllib.error.HTTPError, urllib.error.URLError) as e:
            print("http or URL error")
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


#matplotlib to build the piechart
fig, ((ax1, ax2, ax7, ax5), (ax3, ax4, ax8, ax6)) = plt.subplots(2, 4, figsize=(9, 5))
fig.tight_layout(h_pad=2)
plt.subplots_adjust(bottom=0, right=1.5, top=1.5)
ax1.set_title('filtered blastn piechart')
ax2.set_title('blastn breakdown of other of interest section')
ax3.set_title('filtered blastx piechart')
ax4.set_title('blastx breakdown of other of interest section')
ax5.set_title('blastn mapping percentage')
ax6.set_title('blastx mapping percentage')





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

#filter by frequency if frequency is too small TODO
unrelated_hits_of_blastn = []
rerun_list_blastn = []

#filter by keywords!
for index,row in top_blastn_freq_table.iterrows():
    term = str(row['qseqid'])
    print(str(index))
    esearch_string = esearch(term=term, db='nucleotide')
    if (esearch_string == 'No_Species'):
        rerun_list_blastn.append((index,term))
        continue
    time.sleep(0.5)
    result = get_esummary(esearch_string=esearch_string, db='nucleotide')
    if result == "":
        rerun_list_blastn.append((index,term))
        continue
    try:
        result = xmltodict.parse(result)
    except xmltodict.expat.ExpatError:
        result = ""
    if result == "":
        rerun_list_blastn.append((index,term))
        continue
    sseq_name = ncbi_parse(result)
    found = False
    for i in range(len(filterlist)):
        if filterlist[i] in sseq_name:
            top_blastn_freq_table.loc[index, 'qseqid'] = sseq_name
            found = True
    if found == False:
        unrelated_hits_of_blastn.append(index)

#refilter the rerun list that encountered errors before!
for index,term in rerun_list_blastn:
    print(str(index))
    esearch_string = esearch(term=term, db='nucleotide')
    if (esearch_string == 'No_Species'):
        continue
    time.sleep(0.5)
    result = get_esummary(esearch_string=esearch_string, db='nucleotide')
    if result == "":
        continue
    try:
        result = xmltodict.parse(result)
    except xmltodict.expat.ExpatError:
        result = ""
    if result == "":
        continue
    sseq_name = ncbi_parse(result)
    found = False
    for i in range(len(filterlist)):
        if filterlist[i] in sseq_name:
            top_blastn_freq_table.loc[index, 'qseqid'] = sseq_name
            found = True
    if found == False:
        unrelated_hits_of_blastn.append(index)
        
for i in range(len(unrelated_hits_of_blastn)):
    top_blastn_freq_table = top_blastn_freq_table.drop(unrelated_hits_of_blastn[i])
# top_blastn_freq_table

output_n_csv = piePathSample + '_n.csv'
top_blastn_freq_table.to_csv(output_n_csv, index=False)
# top_blastn_freq_table = pd.read_csv(output_n_csv)

num_seqs_of_interests_blastn = 0 
for index,row in top_blastn_freq_table.iterrows():
    num_seqs_of_interests_blastn += int(row['frequency'])

to_remove_small_hits_blastn = []
other_count_interest_blastn = 0
top_three_blastn = []
top_blastn_freq_table = top_blastn_freq_table.sort_values(by="frequency",ascending=False)
for index,row in top_blastn_freq_table.iterrows():
    if len(top_three_blastn) < 3:
        top_three_blastn.append((row['frequency'],row['qseqid']))
    if(row['frequency'] < (0.01)*(original_blastn_numseqs)):
        to_remove_small_hits_blastn.append(index)
        other_count_interest_blastn += int(row['frequency'])

filtered_top_blastn_of_interest = top_blastn_freq_table.drop(index=to_remove_small_hits_blastn)
filtered_top_blastn_of_interest.rename(columns = {'qseqid':'name'}, inplace = True)
filtered_top_blastn_of_interest = filtered_top_blastn_of_interest.reset_index()
filtered_top_blastn_of_interest = filtered_top_blastn_of_interest.drop(columns = ('index'))
filtered_top_blastn_of_interest.loc[len(filtered_top_blastn_of_interest)] = ['Other of interest', other_count_interest_blastn]
filtered_top_blastn_of_interest.loc[len(filtered_top_blastn_of_interest)] = ['Nonrelevant species', original_blastn_numseqs - num_seqs_of_interests_blastn]

output_n_interest_csv = piePathSample + '_n_interest.csv'
filtered_top_blastn_of_interest.to_csv(output_n_interest_csv, index=False)
# filtered_top_blastn_of_interest = pd.read_csv(output_n_interest_csv)

blastn_pie1 = filtered_top_blastn_of_interest['frequency']
blastn_pie1_label = filtered_top_blastn_of_interest['name']
blastn_pie_explode = [0 for i in range(len(blastn_pie1)-1)]
blastn_pie_explode.append(0.1)
wedges_blastn, *_ = ax1.pie(blastn_pie1, labels = blastn_pie1_label,explode=blastn_pie_explode,colors=None,autopct='%1.1f%%',startangle=45,
        wedgeprops={"linewidth": 1, "edgecolor": "white"})

filtered_other_blastn_of_interest= top_blastn_freq_table.filter(items = to_remove_small_hits_blastn, axis=0)
blastn_bar_category = ['Homo sapiens','bacteria/virus', 'Mus musculus']
blastn_bar_sum = [0,0,0]
for index,row in filtered_other_blastn_of_interest.iterrows():
    if ('Homo sapiens') in row['qseqid']:
        blastn_bar_sum[0] += row['frequency']
    elif ('bacteria') in row['qseqid']:
        blastn_bar_sum[1] += row['frequency']
    elif ('bacterium') in row['qseqid']:
        blastn_bar_sum[1] += row['frequency']
    elif ('strain') in row['qseqid']:
        blastn_bar_sum[1] += row['frequency']
    elif ('virus') in row['qseqid']:
        blastn_bar_sum[1] += row['frequency']
    elif ('Mus musculus') in row['qseqid']:
        blastn_bar_sum[2] += row['frequency']
        
ax2.bar(blastn_bar_category,blastn_bar_sum)
ax2.set_ylabel('num of hits')
ax2.set_xlabel('Categories')

### Show the top three hits for blastn by horizontal bar
ax7_three_seqs_freq = [x[0] for x in top_three_blastn]
ax7_three_seqs_freq.reverse()
ax7_three_seqs_label = [x[1] for x in top_three_blastn]
ax7_three_seqs_label.reverse()
ax7.barh(ax7_three_seqs_label,ax7_three_seqs_freq, align='center')
ax7.set_xlabel('Frequency')
ax7.set_title('Top Three Hits For Blast_n')

### Get the percentage of the the blast n/x mapped results
# for ax5
ax5_percentages = [config['N_downsample_reads']*5 - original_blastn_numseqs, original_blastn_numseqs]
ax5_labels = ['Unmapped', 'Mapped']
ax5.pie(ax5_percentages, labels=ax5_labels, autopct='%1.1f%%',
       colors=['skyblue', 'gray'], labeldistance=1.1)
# Adding Circle in Pie chart
circle1 = plt.Circle((0, 0), radius=0.6, color='white')
ax5.add_patch(circle1)
text1 = str(original_blastn_numseqs) + " out of " + str(config['N_downsample_reads']*5) + " sequences mapped"
ax5.text(0.95, 0.95, text1, transform=ax5.transAxes, fontsize=14,
        verticalalignment='top', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))

print('finish blastn')











# Creating plots for blastx

# Read the blast x file
origin_blastx = pd.read_csv(blastx_file, header=None, sep='\t')
origin_blastx.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']

origin_blastx2 = origin_blastx[['qseqid','sseqid','pident']].copy()
original_blastx_numseqs = len(origin_blastx)

all_blastx_freq_dict = {}
for eachid,eachp in zip(origin_blastx2['sseqid'], origin_blastx2['pident']):
    if eachid not in all_blastx_freq_dict:
        all_blastx_freq_dict[eachid] = [eachp, 1]
    elif float(eachp) > all_blastx_freq_dict[eachid][0]:
        all_blastx_freq_dict[eachid] = [float(eachp), 1]
    elif float(eachp) == all_blastx_freq_dict[eachid][0]:
        all_blastx_freq_dict[eachid][1] += 1

all_blastx_freq_dict = sorted(all_blastx_freq_dict.items(), key=lambda x:x[1], reverse=True)

checkqseqid_x = []
top_blastx_freq_dict = {}
for eachqid,eachsid,eachp in zip(origin_blastx2['qseqid'],origin_blastx2['sseqid'],origin_blastx2['pident']):
    if eachqid not in checkqseqid_x:
        #first instance of qseqid encountered, add top hit
        checkqseqid_x.append(eachqid)
        if eachsid not in top_blastx_freq_dict:
            top_blastx_freq_dict[eachsid] = 1
        else:
            top_blastx_freq_dict[eachsid] += 1
        continue
    #not the first hit for the qseqid
    if eachp == 100:
        if eachsid not in top_blastx_freq_dict:
            top_blastx_freq_dict[eachsid] = 1
        else:
            top_blastx_freq_dict[eachsid] += 1
    else:
        #ignore if the score is less than 100 until next qseqid is reached
        continue

top_blastx_freq_table = pd.DataFrame(list(top_blastx_freq_dict.items()), columns=['qseqid','frequency'])

# filter by frequency if frequency is too small TODO
unrelated_hits_of_blastx = []
rerun_list_blastx = []

#filter by keywords!
for index,row in top_blastx_freq_table.iterrows():
    term = str(row['qseqid'])
    print(str(index))
    esearch_string = esearch(term=term, db='protein')
    if (esearch_string == 'No_Species'):
        rerun_list_blastx.append((index,term))
        continue
    time.sleep(0.5)
    result = get_esummary(esearch_string=esearch_string, db='protein')
    if result == "":
        rerun_list_blastx.append((index,term))
        continue
    try:
        result = xmltodict.parse(result)
    except xmltodict.expat.ExpatError:
        result = ""
    if result == "":
        rerun_list_blastx.append((index,term))
        continue
    sseq_name = ncbi_parse(result)
    found = False
    for i in range(len(filterlist)):
        if filterlist[i] in sseq_name:
            top_blastx_freq_table.loc[index, 'qseqid'] = sseq_name
            found = True
    if found == False:
        unrelated_hits_of_blastx.append(index)

#refilter the rerun list that encountered errors before!
for index,term in rerun_list_blastx:
    print(str(index))
    esearch_string = esearch(term=term, db='protein')
    if (esearch_string == 'No_Species'):
        continue
    time.sleep(0.5)
    result = get_esummary(esearch_string=esearch_string, db='protein')
    if result == "":
        continue
    try:
        result = xmltodict.parse(result)
    except xmltodict.expat.ExpatError:
        result = ""
    if result == "":
        continue
    sseq_name = ncbi_parse(result)
    found = False
    for i in range(len(filterlist)):
        if filterlist[i] in sseq_name:
            top_blastn_freq_table.loc[index, 'qseqid'] = sseq_name
            found = True
    if found == False:
        unrelated_hits_of_blastx.append(index)

for i in range(len(unrelated_hits_of_blastx)):
    top_blastx_freq_table = top_blastx_freq_table.drop(unrelated_hits_of_blastx[i])
#top_blastx_freq_table

output_x_csv = piePathSample + '_x.csv'
top_blastx_freq_table.to_csv(output_x_csv, index=False)
# top_blastx_freq_table = pd.read_csv('qseqidDf_x.csv')

num_seqs_of_interests_blastx = 0 
for index,row in top_blastx_freq_table.iterrows():
    num_seqs_of_interests_blastx += int(row['frequency'])

to_remove_small_hits_blastx = []
other_count_interest_blastx = 0
top_three_blastx = []
top_blastx_freq_table = top_blastx_freq_table.sort_values(by="frequency",ascending=False)
for index,row in top_blastx_freq_table.iterrows():
    if len(top_three_blastx) < 3:
        top_three_blastx.append((row['frequency'],row['qseqid']))
    if(row['frequency'] < (0.001)*(original_blastx_numseqs)):
        to_remove_small_hits_blastx.append(index)
        other_count_interest_blastx += int(row['frequency'])

filtered_top_blastx_of_interest = top_blastx_freq_table.drop(index=to_remove_small_hits_blastx)
filtered_top_blastx_of_interest.rename(columns = {'qseqid':'name'}, inplace = True)
filtered_top_blastx_of_interest = filtered_top_blastx_of_interest.reset_index()
filtered_top_blastx_of_interest = filtered_top_blastx_of_interest.drop(columns = ('index'))
filtered_top_blastx_of_interest.loc[len(filtered_top_blastx_of_interest)] = ['Other of interest', other_count_interest_blastx]
filtered_top_blastx_of_interest.loc[len(filtered_top_blastx_of_interest)] = ['Nonrelevant species', original_blastx_numseqs - num_seqs_of_interests_blastx]

output_x_interest_csv = piePathSample + '_x_interest.csv'
filtered_top_blastx_of_interest.to_csv(output_x_interest_csv, index=False)
# filtered_top_blastx_of_interest = pd.read_csv('qseqidDf_x_interest.csv')

blastx_pie1 = filtered_top_blastx_of_interest['frequency']
blastx_pie1_label = filtered_top_blastx_of_interest['name']
blastx_pie_explode = [0 for i in range(len(blastx_pie1)-1)]
blastx_pie_explode.append(0.1)
wedges_blastx, *_ = ax3.pie(blastx_pie1, labels = blastx_pie1_label,explode=blastx_pie_explode,colors=None,autopct='%1.1f%%',startangle=45,
        wedgeprops={"linewidth": 1, "edgecolor": "white"})

filtered_other_blastx_of_interest= top_blastx_freq_table.filter(items = to_remove_small_hits_blastx, axis=0)
blastx_bar_category = ['Homo sapiens','bacteria/virus', 'Mus musculus']
blastx_bar_sum = [0,0,0]
for index,row in filtered_other_blastx_of_interest.iterrows():
    if ('Homo sapiens') in row['qseqid']:
        blastx_bar_sum[0] += row['frequency']
    elif ('bacteria') in row['qseqid']:
        blastx_bar_sum[1] += row['frequency']
    elif ('bacterium') in row['qseqid']:
        blastx_bar_sum[1] += row['frequency']
    elif ('strain') in row['qseqid']:
        blastx_bar_sum[1] += row['frequency']
    elif ('virus') in row['qseqid']:
        blastx_bar_sum[1] += row['frequency']
    elif ('Mus musculus') in row['qseqid']:
        blastx_bar_sum[2] += row['frequency']
        
ax4.bar(blastx_bar_category,blastx_bar_sum)
ax4.set_ylabel('num of hits')
ax4.set_xlabel('Categories')

### Show the top three hits for blastx by horizontal bar
ax8_three_seqs_freq = [x[0] for x in top_three_blastx]
ax8_three_seqs_freq.reverse()
ax8_three_seqs_label = [x[1] for x in top_three_blastx]
ax8_three_seqs_label.reverse()
ax8.barh(ax8_three_seqs_label,ax8_three_seqs_freq, align='center')
ax8.set_xlabel('Frequency')
ax8.set_title('Top Three Hits For Blast_x')

### Get the percentage of the the blast n/x mapped results
# for ax5
ax6_percentages = [config['N_downsample_reads']*5 - original_blastx_numseqs, original_blastx_numseqs]
ax6_labels = ['Unmapped', 'Mapped']
ax6.pie(ax6_percentages, labels=ax6_labels, autopct='%1.1f%%',
       colors=['skyblue', 'gray'], labeldistance=1.1)
# Adding Circle in Pie chart
circle2 = plt.Circle((0, 0), radius=0.6, color='white')
ax6.add_patch(circle2)
text2 = str(original_blastx_numseqs) + " out of " + str(config['N_downsample_reads']*5) + " sequences mapped"
ax6.text(0.95, 0.95, text2, transform=ax6.transAxes, fontsize=14,
        verticalalignment='top', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))


print('finish blastx')


# plt.subplots_adjust(hspace=0.4, wspace=0.9)
fig.text(0.3, 0.85, '⬇Unrelevant', ha='center', va='center', fontsize=18)
fig.text(0.3, 0.05, '⬇Unrelevant', ha='center', va='center', fontsize=18)

plt.title('Blastn and Blastx Unmapped Sequences Summary ' + pieName)
fig.tight_layout(pad=5.0)
plt.show(block=True)
plt.savefig(piePath,format='png',bbox_inches='tight')
