# foRWaRD (Informed Random Walk for Ranking of Diseases)

![Pipeline]
(https://github.com/ddhostallero/foRWaRD/forRWaRD_pipe.png)

# Setup
The following dependencies are required:
 - Python 3.7
 - Numpy
 - Pandas
 - NetworkX
 - SciKit Learn
 - Matplotlib
 
# Input Files
## Network file
Gene-gene network files and gene-mapping files can be downloaded directly from the [KnowEnG KN-Fetcher](https://github.com/KnowEnG/KN_Fetcher/blob/master/Contents.md). 

### Custom Network Files
For custom network files, the file should be tab-separted and contains no header row. The first and second column should be the Ensembl ID of the genes and the third column is the normalized weight of the edges such that the total weight sums up to 1. A gene-mapping file should also be provided as a headerless, tab-separated file. The first column should be the Ensembl ID and the third column should be the gene name.

##  Gene-Disease Association
The gene-disease association file can be downloaded from [DisGeNet](https://www.disgenet.org/)

## Disease of interests
List the diseases of interest in the `disease_of_interest.txt`. Make sure that all these diseases exist in the gene-disease association file.

# Running the Program
To run the program specify the following arguments:
 - `--network` | `-n`: Name of the network file (e.g. "data/KnowEnG/9606.hn_IntNet.edge")
 - `--mapping` |`-m` : Name of the gene-mapping file (e.g. "data/KnowEnG/9606.hn_IntNet.node_map")
 - `--assoc` | `-a`: Name of the disease-gene association file (e.g. "data/DisGeNet/all_gene_disease_associations.tsv")
 - `--disease` | `-d`: Name of the file containing the list of diseases of interest (e.g. "data/disease_of_interest.txt")
 - `--gene` | `-g`: Name of the file containing the list of genes of interest (e.g. "data/genes_of_interest.txt")
 - `--weight` | `-w`: Edge weight scheme of the disease-gene association edges (either "gda" or "ones")

```
python main.py -w gda -g data/genes_of_interest.txt -d data/disease_of_interest.txt \
	-a data/DisGeNet/all_gene_disease_associations.tsv \
	-m data/KnowEnG/9606.hn_IntNet.node_map \
	-n data/KnowEnG/9606.hn_IntNet.edge 
```