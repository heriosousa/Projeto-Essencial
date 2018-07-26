# Gene Essentiality Prediction Using Topological Features From Metabolic Networks

This is the code repository from the paper "Gene Essentiality Prediction Using Topological Features From Metabolic Networks" to be presented at 7th Brazilian Conference on Intelligent Systems (BRACIS) in São Paulo, Brazil.

In this work, we evaluate the influence and the contribution of network topological features for the essencial gene prediction task in metabolic networks.

## Data sources

Selected organisms:

- *Escherichia coli*
- *Mycoplasma genitalium*
- *Pseudomonas aeruginosa*
- *Saccharomyces cerevisiae*

All metabolic networks were construted from the Kyoto Encyclopedia of Genes and Genomes (KEGG) metabolic pathway database, while essentiality information for each organism was collected from diferent sources:
- Profiling of *Escherichia coli* Chromosome Database (PEC) for *E. coli* 
- Database of Essential Genes (DEG) for *M. genitalium* and *P. aerguinosa*
- Saccharomyces Genome Deletion Project for *S. cerevisiae*. 

## Scripts and steps

DataVisualization.ipynb contains how data was collected for both networks' construction and gene essentiality information, and how these data are integrated in a graph (for each organism) to be used later in the feature extraction step. Some graph information is also avaliable for all organisms, including its visualization according to node essentiality and frequency.

PreProcessing.ipynb shows feature extraction is performed from the constructed graphs, using several topology features from different domains, to be later discussed.

*InitialTest.ipynb comento melhor assim que decidirmos se vamos cortar ou não essa parte*

MachineLearningApproach.ipynb displays both experiments scenarios: pairwise and leave-one-out, including its results (assessed by ROC curve) and runtime, as well as the selected classifiers.

## Selecting the feature set

We have opted to build graphs using an undirected representation, increasing the number of potential topology features to be selected into the model, including:

### Graph centrality measures

In graph theory, centrality is a measure that evaluates how important a node is in a network, according to different criteria depending on how importance is characterized in the centrality. Graph centrality measures are commonly used in gene essentiality prediction in different biological networks, the most common ones being:

- Degree centrality, the number of neighbors a node has
- Betweeness centrality, the number of shortest paths that pass through a node
- Closeness centrality, the average length of the shortest paths between a node and all other nodes 
- Eigenvector centrality, the sum of the centrality values of the node neighbors

In additional to these measures, we decided to include other graph centralities:



### Other topology features

## Citation

Available after publishing.
