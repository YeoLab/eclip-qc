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
import math


#handles command line input parameter
import sys
import os

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######

# Make sure to modify the path to your config.yaml below!!!

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
with open("/home/s5xu/scratch/temp_eclip/eclip-qc/config.yaml", 'r') as stream:
    config = yaml.safe_load(stream)

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


# Data sorting and filtering
# Read the blast n file
origin_blastn = pd.read_csv(blastn_file, header=None, sep='\t')
origin_blastn.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']

origin_blastn2 = origin_blastn[['qseqid','sseqid','pident']].copy()
original_blastn_numseqs = len(origin_blastn)

priority = origin_blastn2['sseqid'].value_counts().to_dict()
priority['Homo sapiens'] = 1000000000
priority['Mus musculus'] = priority['Homo sapiens'] - 1
priority['bacteria'] = priority['Mus musculus'] - 1
priority['bacterium'] = priority['bacteria']
priority['virus'] = priority['bacteria']
priority['strain'] = priority['bacteria']

# Generate sseqids blastn dict if applicable
rerun_list_blastn = []
unconverted_blast_n_set = set(origin_blastn2['sseqid'])
default_none = None
blast_n_dict = {sseqid: default_none for sseqid in unconverted_blast_n_set}
blastn_dict_tsv_readin = pd.read_csv('other_source/blastn_dict.tsv', sep='\t')
blast_n_readin_dict = dict(zip(blastn_dict_tsv_readin['sseqid'], blastn_dict_tsv_readin['sseq_name']))

for each_sseqid in blast_n_dict:
    if each_sseqid not in blast_n_readin_dict.keys():
        term = str(each_sseqid)
        print(term)
        esearch_string = esearch(term=term, db='nucleotide')
        if (esearch_string == 'No_Species'):
            rerun_list_blastn.append(term)
            continue
        time.sleep(0.5)
        result = get_esummary(esearch_string=esearch_string, db='nucleotide')
        if result == "":
            rerun_list_blastn.append(term)
            continue
        try:
            result = xmltodict.parse(result)
        except xmltodict.expat.ExpatError:
            result = ""
        if result == "":
            rerun_list_blastn.append(term)
            continue
        sseq_name = ncbi_parse(result)
        blast_n_dict[each_sseqid] = sseq_name
        

#refilter the rerun list that encountered errors before!
for term in rerun_list_blastn:
    print(term)
    esearch_string = esearch(term=term, db='nucleotide')
    if (esearch_string == 'No_Species'):
        blast_n_dict[term] = term
        continue
    time.sleep(0.5)
    result = get_esummary(esearch_string=esearch_string, db='nucleotide')
    if result == "":
        blast_n_dict[term] = term
        continue
    try:
        result = xmltodict.parse(result)
    except xmltodict.expat.ExpatError:
        result = ""
    if result == "":
        blast_n_dict[term] = term
        continue
    sseq_name = ncbi_parse(result)
    blast_n_dict[term] = sseq_name

blast_n_dict_combine = {**blast_n_dict, **blast_n_readin_dict}
save_n_dict = pd.DataFrame(blast_n_dict_combine.items(), columns=['sseqid', 'sseq_name'])
save_n_dict.to_csv('other_source/blastn_dict.tsv', sep='\t', index=False)

print('output the blast n dict file')

# fill in the sseq_names into blastn tsv table
for index,row in origin_blastn2.iterrows():
    origin_blastn2.loc[index, 'sseq_name'] = blast_n_dict_combine[row['sseqid']]
    
output_n_csv = piePathSample + '_n.csv'
origin_blastn2.to_csv(output_n_csv, index=False)
# origin_blastn2 = pd.read_csv(output_n_csv)
print('output the', output_n_csv, 'n.csv file')

top_blastn_freq_dict = {}
curr_count_n = 0
potential_n = []
curr_n = None
next_n = None

for index,row in origin_blastn2.iterrows():
    curr_n = row['qseqid']
    if index == original_blastn_numseqs - 1:
        next_n = 'end'
    else:
        next_n = origin_blastn2.loc[index+1, 'qseqid']
    found = False
    species_n_type = None
    if len(potential_n) == 0:
        for i in range(len(filterlist)):
            if filterlist[i] in row['sseq_name']:
                found = True
                species_n_type = filterlist[i]
                potential_n.append((row['sseq_name'], row['pident'], priority[species_n_type]))
                break
        if found == False:
            potential_n.append((row['sseq_name'], row['pident'], priority[row['sseqid']]))
    elif curr_n == next_n:
        for i in range(len(filterlist)):
            if filterlist[i] in row['sseq_name']:
                found = True
                species_n_type = filterlist[i]
                potential_n.append((row['sseq_name'], row['pident'], priority[species_n_type]))
                break
        if found == False:
            potential_n.append((row['sseq_name'], row['pident'], priority[row['sseqid']]))
    if curr_n != next_n:
        if len(potential_n) == 0:
            continue
        sorted_n_potential = sorted(potential_n, key=lambda x: (x[1], x[2]), reverse=True)
        if sorted_n_potential[0][0] not in top_blastn_freq_dict:
            top_blastn_freq_dict[sorted_n_potential[0][0]] = 1
        else:
            top_blastn_freq_dict[sorted_n_potential[0][0]] += 1
        potential_n = []

top_blastn_freq_table = pd.DataFrame(top_blastn_freq_dict.items(), columns=['qseqid', 'frequency'])
output_n_csv = piePathSample + '_filtered_n.csv'
top_blastn_freq_table.to_csv(output_n_csv, index=False)

print('output the filtered', output_n_csv, 'n.csv file')


# Read the blast x file
origin_blastx = pd.read_csv(blastx_file, header=None, sep='\t')
origin_blastx.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']

origin_blastx2 = origin_blastx[['qseqid','sseqid','pident']].copy()
original_blastx_numseqs = len(origin_blastx)

priority = origin_blastx2['sseqid'].value_counts().to_dict()
priority['Homo sapiens'] = 1000000000
priority['Mus musculus'] = priority['Homo sapiens'] - 1
priority['bacteria'] = priority['Mus musculus'] - 1
priority['bacterium'] = priority['bacteria']
priority['virus'] = priority['bacteria']
priority['strain'] = priority['bacteria']

# Generate sseqids blastx dict if applicable
rerun_list_blastx = []
unconverted_blast_x_set = set(origin_blastx2['sseqid'])
default_none = None
blast_x_dict = {sseqid: default_none for sseqid in unconverted_blast_x_set}
blastx_dict_tsv_readin = pd.read_csv('other_source/blastx_dict.tsv', sep='\t')
blast_x_readin_dict = dict(zip(blastx_dict_tsv_readin['sseqid'], blastx_dict_tsv_readin['sseq_name']))

for each_sseqid in blast_x_dict:
    if each_sseqid not in blast_x_readin_dict.keys():
        term = str(each_sseqid)
        print(term)
        esearch_string = esearch(term=term, db='protein')
        if (esearch_string == 'No_Species'):
            rerun_list_blastx.append(term)
            continue
        time.sleep(0.5)
        result = get_esummary(esearch_string=esearch_string, db='protein')
        if result == "":
            rerun_list_blastx.append(term)
            continue
        try:
            result = xmltodict.parse(result)
        except xmltodict.expat.ExpatError:
            result = ""
        if result == "":
            rerun_list_blastx.append(term)
            continue
        sseq_name = ncbi_parse(result)
        blast_x_dict[each_sseqid] = sseq_name

#refilter the rerun list that encountered errors before!
for term in rerun_list_blastx:
    print(term)
    esearch_string = esearch(term=term, db='protein')
    if (esearch_string == 'No_Species'):
        blast_x_dict[term] = term
        continue
    time.sleep(0.5)
    result = get_esummary(esearch_string=esearch_string, db='protein')
    if result == "":
        blast_x_dict[term] = term
        continue
    try:
        result = xmltodict.parse(result)
    except xmltodict.expat.ExpatError:
        result = ""
    if result == "":
        blast_x_dict[term] = term
        continue
    sseq_name = ncbi_parse(result)
    blast_x_dict[term] = sseq_name

blast_x_dict_combine = {**blast_x_dict, **blast_x_readin_dict}
save_x_dict = pd.DataFrame(blast_x_dict_combine.items(), columns=['sseqid', 'sseq_name'])
save_x_dict.to_csv('other_source/blastx_dict.tsv', sep='\t', index=False)
  
print('output the blast x dict file')

# fill in the sseq_names into blastn tsv table
for index,row in origin_blastx2.iterrows():
    origin_blastx2.loc[index, 'sseq_name'] = blast_x_dict_combine[row['sseqid']]

output_x_csv = piePathSample + '_x.csv'
origin_blastx2.to_csv(output_x_csv, index=False)
# origin_blastx2 = pd.read_csv('test_n.csv')
print('output the', piePathSample, 'x.csv file')

top_blastx_freq_dict = {}
curr_count_x = 0
unrelated_hits_of_blastx = []
potential_x = []
curr_x = None
next_x = None

for index,row in origin_blastx2.iterrows():
    curr_x = row['qseqid']
    if index == original_blastx_numseqs - 1:
        next_x = 'end'
    else:
        next_x = origin_blastx2.loc[index+1, 'qseqid']
    found = False
    species_x_type = None
    if len(potential_x) == 0:
        for i in range(len(filterlist)):
            if filterlist[i] in row['sseq_name']:
                found = True
                species_x_type = filterlist[i]
                potential_x.append((row['sseq_name'], row['pident'], priority[species_x_type]))
                break
        if found == False:
            potential_x.append((row['sseq_name'], row['pident'], priority[row['sseqid']]))
    elif curr_x == next_x:
        for i in range(len(filterlist)):
            if filterlist[i] in row['sseq_name']:
                found = True
                species_x_type = filterlist[i]
                potential_x.append((row['sseq_name'], row['pident'], priority[species_x_type]))
                break
        if found == False:
            potential_x.append((row['sseq_name'], row['pident'], priority[row['sseqid']]))
    if curr_x != next_x:
        if len(potential_x) == 0:
            continue
        sorted_x_potential = sorted(potential_x, key=lambda x: (x[1], x[2]), reverse=True)
        if sorted_x_potential[0][0] not in top_blastx_freq_dict:
            top_blastx_freq_dict[sorted_x_potential[0][0]] = 1
        else:
            top_blastx_freq_dict[sorted_x_potential[0][0]] += 1
        potential_x = []

top_blastx_freq_table = pd.DataFrame(top_blastx_freq_dict.items(), columns=['qseqid', 'frequency'])
output_x_csv = piePathSample + '_filtered_x.csv'
top_blastx_freq_table.to_csv(output_x_csv, index=False)

print('output the filtered', piePathSample, 'x.csv file')

# Data visualization/ploting

#matplotlib to build the piechart
fig, ((ax5, ax1, ax7, ax2), (ax6, ax3, ax8, ax4)) = plt.subplots(2, 4, figsize=(9, 5))
fig.tight_layout(h_pad=2)
plt.subplots_adjust(bottom=0, right=1.5, top=1.5)
ax1.set_title('filtered blastn piechart')
ax2.set_title('blastn breakdown of other of interest section')
ax3.set_title('filtered blastx piechart')
ax4.set_title('blastx breakdown of other of interest section')
ax5.set_title('blastn mapping percentage')
ax6.set_title('blastx mapping percentage')


# ploting for blast n

### Maybe rename this to filtered selected blastn representatives
num_seqs_of_interests_blastn = 0 
for index,row in top_blastn_freq_table.iterrows():
    num_seqs_of_interests_blastn += int(row['frequency'])

top_three_blastn = []
top_blastn_freq_table = top_blastn_freq_table.sort_values(by="frequency",ascending=False)
for index,row in top_blastn_freq_table.iterrows():
    if len(top_three_blastn) < 3:
        top_three_blastn.append((row['frequency'],row['qseqid'])) 
    
blastn_pie1 = [num_seqs_of_interests_blastn, original_blastn_numseqs - num_seqs_of_interests_blastn]
blastn_pie1_total_sum = sum(blastn_pie1)

# Custom pie chart labelling function
def pie_label_n(pct):
    absolute = math.ceil(pct/100.*blastn_pie1_total_sum)
    return f'{pct:.1f}% ({absolute}/{blastn_pie1_total_sum})'

# blastn_pie1_label = filtered_top_blastn_of_interest['name']
blastn_pie1_label = ['Selected Representatives','Unselected']
blastn_pie_explode = [0 for i in range(len(blastn_pie1)-1)]
blastn_pie_explode.append(0.1)
wedges_blastn, *_ = ax1.pie(blastn_pie1, labels = blastn_pie1_label,explode=blastn_pie_explode,colors=None,autopct=pie_label_n,startangle=45,
        wedgeprops={"linewidth": 1, "edgecolor": "white"})

filtered_other_blastn_of_interest= top_blastn_freq_table
blastn_bar_category = ['Homo sapiens','bacteria/virus', 'Mus musculus', 'Other']
blastn_bar_sum = [0,0,0,0]
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
    else:
        blastn_bar_sum[3] += row['frequency']

colors = ['red', 'blue', 'green', 'orange'] 
colors2 = ['red', 'blue', 'green']   

#blasn barchat for other interest categories
bars_n = ax2.bar(blastn_bar_category,blastn_bar_sum, color = colors)
ax2.set_ylabel('num of hits')
ax2legend_handles = [plt.Rectangle((0, 0), 1, 1, color=color) for color in colors]
ax2.legend(ax2legend_handles, blastn_bar_category, title="Categories", loc="upper left")
ax2.xticks([])
#ax2.set_xlabel('Categories')

for bar in bars_n:
    yval = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2.0, yval, int(yval), va='bottom')

### Show the top three hits for blastn by horizontal bar
ax7_three_seqs_freq = [x[0] for x in top_three_blastn]
ax7_three_seqs_freq.reverse()
ax7_three_seqs_label = [x[1] for x in top_three_blastn]
ax7_three_seqs_label.reverse()
ax7.barh(ax7_three_seqs_label,ax7_three_seqs_freq, align='center', color = colors2)
ax7.set_xlabel('Frequency')
ax7.set_title('Top Three Hits For Blast_n')

ax7legend_handles = [plt.Rectangle((0, 0), 1, 1, color=color) for color in colors2]
ax7.legend(ax7legend_handles, ax7_three_seqs_label, title="Categories", loc="lower right")
ax7.xticks([])

### Get the percentage of the the blast n/x mapped results
# for ax5
# Fix to out of N_downsample_reads rather than out of the total hits. Rename variables will be needed later.
ax5_percentages = [config['N_downsample_reads'] - len(set(origin_blastn2['qseqid'])), len(set(origin_blastn2['qseqid']))]
ax5_percentages_sum = config['N_downsample_reads']

def donut_n_label(pct):
    absolute = math.ceil(pct/100.*ax5_percentages_sum)
    return f'{pct:.1f}% ({absolute}/{ax5_percentages_sum})'

ax5_labels = ['Unmapped', 'Mapped']
ax5.pie(ax5_percentages, labels=ax5_labels, autopct=donut_n_label,
       colors=['skyblue', 'gray'], labeldistance=1.1)
# Adding Circle in Pie chart
circle1 = plt.Circle((0, 0), radius=0.6, color='white')
ax5.add_patch(circle1)
text1 = str(len(set(origin_blastn2['qseqid']))) + " out of " + str(config['N_downsample_reads']) + " sequences mapped"
ax5.text(0.00, 0.9, text1, transform=ax5.transAxes, fontsize=14,
        verticalalignment='top', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))

print('Finished drawing blast n')

# ploting for blast x

### Maybe rename this to filtered selected blastx representatives
num_seqs_of_interests_blastx = 0 
for index,row in top_blastx_freq_table.iterrows():
    num_seqs_of_interests_blastx += int(row['frequency'])

top_three_blastx = []
top_blastx_freq_table = top_blastx_freq_table.sort_values(by="frequency",ascending=False)
for index,row in top_blastx_freq_table.iterrows():
    if len(top_three_blastx) < 3:
        top_three_blastx.append((row['frequency'],row['qseqid']))
    
blastx_pie1 = [num_seqs_of_interests_blastx, original_blastx_numseqs - num_seqs_of_interests_blastx]
blastx_pie1_total_sum = sum(blastx_pie1)

# Custom pie chart labelling function
def pie_label_x(pct):
    absolute = math.ceil(pct/100.*blastx_pie1_total_sum)
    return f'{pct:.1f}% ({absolute}/{blastx_pie1_total_sum})'

blastx_pie1_label = ['Selected Representatives','Unselected']
blastx_pie_explode = [0 for i in range(len(blastx_pie1)-1)]
blastx_pie_explode.append(0.1)
wedges_blastx, *_ = ax3.pie(blastx_pie1, labels = blastx_pie1_label,explode=blastx_pie_explode,colors=None,autopct=pie_label_x,startangle=45,
        wedgeprops={"linewidth": 1, "edgecolor": "white"})
    

filtered_other_blastx_of_interest= top_blastx_freq_table
blastx_bar_category = ['Homo sapiens','bacteria/virus', 'Mus musculus', 'Other']
blastx_bar_sum = [0,0,0,0]
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
    else:
        blastx_bar_sum[3] += row['frequency']

#blastx barchart for other interest categories       
bars_x = ax4.bar(blastx_bar_category,blastx_bar_sum, color = colors)
ax4.set_ylabel('num of hits')
ax4legend_handles = [plt.Rectangle((0, 0), 1, 1, color=color) for color in colors]
ax4.legend(ax4legend_handles, blastx_bar_category, title="Categories", loc="upper left")
ax4.xticks([])
#ax4.set_xlabel('Categories')


for bar in bars_x:
    yval = bar.get_height()
    ax4.text(bar.get_x() + bar.get_width()/2.0, yval, int(yval), va='bottom')

### Show the top three hits for blastx by horizontal bar
ax8_three_seqs_freq = [x[0] for x in top_three_blastx]
ax8_three_seqs_freq.reverse()
ax8_three_seqs_label = [x[1] for x in top_three_blastx]
ax8_three_seqs_label.reverse()
ax8.barh(ax8_three_seqs_label,ax8_three_seqs_freq, align='center', color = colors2)
ax8.set_xlabel('Frequency')
ax8.set_title('Top Three Hits For Blast_x')

ax8legend_handles = [plt.Rectangle((0, 0), 1, 1, color=color) for color in colors2]
ax8.legend(ax8legend_handles, ax8_three_seqs_label, title="Categories", loc="lower right")
ax8.xticks([])

### Get the percentage of the the blast n/x mapped results
# for ax6
# Fix to out of N_downsample_reads rather than out of the total hits. Rename variables will be needed later.
ax6_percentages = [config['N_downsample_reads'] - len(set(origin_blastx2['qseqid'])), len(set(origin_blastx2['qseqid']))]
ax6_percentages_sum = config['N_downsample_reads']

def donut_x_label(pct):
    absolute = math.ceil(pct/100.*ax6_percentages_sum)
    return f'{pct:.1f}% ({absolute}/{ax6_percentages_sum})'

ax6_labels = ['Unmapped', 'Mapped']
ax6.pie(ax6_percentages, labels=ax6_labels, autopct=donut_x_label,
       colors=['skyblue', 'gray'], labeldistance=1.1)
# Adding Circle in Pie chart
circle2 = plt.Circle((0, 0), radius=0.6, color='white')
ax6.add_patch(circle2)
text2 = str(len(set(origin_blastx2['qseqid']))) + " out of " + str(config['N_downsample_reads']) + " sequences mapped"
ax6.text(0.00, 0.95, text2, transform=ax6.transAxes, fontsize=14,
        verticalalignment='top', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))


plt.title('Blastn and Blastx Unmapped Sequences Summary ' + pieName)
fig.tight_layout(pad=5.0)
plt.show(block=True)
plt.savefig(piePath,format='png',bbox_inches='tight')

print('All Done!!!')
