import pandas as pd
import numpy as np
import networkx as nx
from sklearn.preprocessing import normalize
from rwr import rwr_vec
import matplotlib.pyplot as plt
from scipy.special import softmax
import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument('-n', '--network', default='data/KnowEnG/9606.hn_IntNet.edge', help='Name of the network file')
parser.add_argument('-m', '--mapping', default='data/KnowEnG/9606.hn_IntNet.node_map', help='Name of the gene mapping file')
parser.add_argument('-a', '--assoc', default='data/DisGeNet/all_gene_disease_associations.tsv', help='Name of the disease-gene assoc file')
parser.add_argument('-d', '--disease', default='data/disease_of_interest.txt', help='Name of the file containing the list of diseases of interest')
parser.add_argument('-g', '--gene', default='data/genes_of_interest.txt', help='Name of the file containing the list of genes of interest')
parser.add_argument('-w', '--weight', default='gda', help='[gda] | ones')
args = parser.parse_args() 

assoc_file = args.assoc
network_file = args.network
mapping_file = args.mapping
genes_file = args.gene 
diseases_file = args.disease
weight = args.weight
threshold = 0.4


# create results folder
results_dir = 'results/'
if not os.path.exists(results_dir):
	os.mkdir(results_dir)

# load diseases of interest
doi = []
with open(diseases_file) as f:
	for line in f:
		doi.append(line.strip())

# load the disease-gene association file
assoc = pd.read_csv(assoc_file, sep='\t')
assoc = assoc.loc[assoc['diseaseName'].isin(doi)]
assoc = assoc.loc[assoc['score'] >= threshold]

# load the network
network = pd.read_csv(network_file, sep='\t', header=None)
network = network[[0,1,2]]
network.columns = ['gene1', 'gene2', 'weight']

# load the gene name mapping 
mapping = pd.read_csv(mapping_file, sep='\t', header=None)

# load genes of interest
goi = []
with open(genes_file) as f:
	for line in f:
		goi.append(line.strip())

# check for mappings
goi_mapping = mapping.loc[mapping[3].isin(goi)]
goi_mapping = goi_mapping.set_index(3)
assoc_mapping = mapping.loc[mapping[3].isin(assoc['geneSymbol'].unique())]
assoc_mapping = assoc_mapping.set_index(3)

not_mapped  = []
for gene in goi:
	if gene not in goi_mapping.index:
		not_mapped.append(gene)

if len(not_mapped) > 0:
	print("Genes not in the network:", not_mapped)
goi = list(goi_mapping.index)

# scale weights of the gene-gene network
network['weight'] = network['weight']/network['weight'].max()
G = nx.from_pandas_edgelist(network, 'gene1', 'gene2', 'weight')

# add disease assoc
for disease in doi:
	disease_assoc = assoc.loc[assoc['diseaseName'] == disease]
	disease_assoc = disease_assoc.loc[disease_assoc['geneSymbol'].isin(assoc_mapping.index)]
	genes = disease_assoc['geneSymbol'].unique()

	if disease_assoc['geneSymbol'].nunique() != len(disease_assoc):
		print('ensure one entry per gene')
		exit()

	if weight == 'gda':
		weights = disease_assoc['score'].values
	elif weight == 'ones':
		weights = np.ones(len(genes))
	else:
		print('weight parameter invalid')

	genes = assoc_mapping.loc[genes][0]

	# add disease node and edges
	G.add_node(disease)
	for i, gene in enumerate(genes):
		G.add_edge(disease, gene, weight=weights[i])


node_names = list(G.nodes)
network_matrix = nx.to_numpy_matrix(G, nodelist=node_names, weight="weight")
network_matrix = normalize(network_matrix, norm='l1', axis=1)

# create restart vector
restart_vec = np.zeros(len(node_names))
for g in goi:
	x = goi_mapping.loc[g][0]
	restart_vec[node_names.index(x)]  = 1

restart_vec = restart_vec/np.sum(restart_vec)

# run RWR
(n_iter, residual, steady_prob) = rwr_vec(node_names, network_matrix, \
            restart_vec, restart_prob=0.5, max_iter=100, tolerance=1e-8)


results = pd.DataFrame(index=node_names)
results['steady_prob_hlh'] = steady_prob

# all genes as restart vector
restart_vec = np.ones(len(node_names))
for d in doi: # zero prob of restarting to a disease
	restart_vec[node_names.index(d)]  = 0
restart_vec = restart_vec/np.sum(restart_vec)

# run RWR
(n_iter, residual, steady_prob) = rwr_vec(node_names, network_matrix, \
            restart_vec, restart_prob=0.5, max_iter=100, tolerance=1e-8)
results['steady_prob_all'] = steady_prob
results = results.loc[doi]

results['difference'] = results['steady_prob_hlh'] - results['steady_prob_all']
results['ratio'] = results['steady_prob_hlh']/results['steady_prob_all']
results['norm_diff'] = (results['difference']/(np.abs(results['difference']).max()))/2 + 0.5
results.sort_values('norm_diff', ascending=False).to_csv(results_dir + 'rwr_steady_prob.csv')

ranking = pd.DataFrame(index=range(1, len(results)+1))
ranking['hlh'] = results.sort_values('steady_prob_hlh', ascending=False).index
ranking['all'] = results.sort_values('steady_prob_all', ascending=False).index

ranking.to_csv(results_dir + 'rwr_ranking.csv')
results['gene_assoc_count'] = assoc.groupby('diseaseName')['geneId'].count()


vmax = results['steady_prob_hlh'].max() + results['steady_prob_hlh'].std()
plt.scatter(results['gene_assoc_count'].values, results['steady_prob_hlh'].values)
plt.ylim(top=vmax, bottom=-results['steady_prob_hlh'].std())
plt.xlabel('number of associated genes')
plt.ylabel('steady state probability')
plt.title('RWR using HLH')
for d in results.index:
	plt.annotate(d[:10], (results.loc[d]['gene_assoc_count'], results.loc[d]['steady_prob_hlh']))
plt.savefig(results_dir + 'RWR_HLH.png')
# plt.show()

plt.clf()
vmax = results['steady_prob_all'].max() + results['steady_prob_all'].std()
plt.scatter(results['gene_assoc_count'].values, results['steady_prob_all'].values)
plt.ylim(top=vmax, bottom=-results['steady_prob_all'].std())
plt.xlabel('number of associated genes')
plt.ylabel('steady state probability')
plt.title('RWR using all genes')
for d in results.index:
	plt.annotate(d[:10], (results.loc[d]['gene_assoc_count'], results.loc[d]['steady_prob_all']))
plt.savefig(results_dir + 'RWR_all.png')
# plt.show()




