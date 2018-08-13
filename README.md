# Gene Essentiality Prediction Using Topological Features From Metabolic Networks

This is the code repository for the paper "Gene Essentiality Prediction Using Topological Features From Metabolic Networks" to be presented at 7th Brazilian Conference on Intelligent Systems (BRACIS) in São Paulo, Brazil.

In this work, we evaluate the influence and the contribution of network topological features for the essencial gene prediction task in metabolic networks.

## Data sources

Selected organisms:

- *Escherichia coli*
- *Mycoplasma genitalium*
- *Pseudomonas aeruginosa*
- *Saccharomyces cerevisiae*

All metabolic networks were construted from the [Kyoto Encyclopedia of Genes and Genomes (KEGG)](https://www.genome.jp/kegg/) metabolic pathway database, while essentiality information for each organism was collected from diferent sources:
- [Profiling of *Escherichia coli* Chromosome Database (PEC)](https://shigen.nig.ac.jp/ecoli/pec/) for *E. coli* 
- [Database of Essential Genes (DEG)](http://www.essentialgene.org/) for *M. genitalium* and *P. aerguinosa*
- [Saccharomyces Genome Deletion Project](http://www-sequence.stanford.edu/group/yeast_deletion_project/deletions3.html) for *S. cerevisiae*. 

## Scripts and steps

DataVisualization.ipynb contains how data was collected for both networks' construction and gene essentiality information, and how these data are integrated in a graph (for each organism) to be used later in the feature extraction step. Some graph information is also avaliable for all organisms, including its visualization according to node essentiality and frequency.

PreProcessing.ipynb shows feature extraction is performed from the constructed graphs, using several topology features from different domains, to be later discussed.

MachineLearningApproach.ipynb displays both experiments scenarios: pairwise and leave-one-out, including its results (assessed by ROC curve) and runtime, as well as the selected classifiers.

## Selecting the feature set

We have opted to build graphs using an undirected representation, increasing the number of potential topology features to be selected into the model, including:

### Graph centrality measures

In graph theory, centrality is a measure that evaluates how important a node is in a network, according to different criteria depending on how importance is characterized in the centrality. Graph centrality measures are commonly used in gene essentiality prediction in different biological networks, the most common ones being:

- Degree centrality, the number of neighbors a node has [More info](https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.degree_centrality.html#networkx.algorithms.centrality.degree_centrality)
- Betweeness centrality, the number of shortest paths that pass through a node [More info](https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.betweenness_centrality.html#networkx.algorithms.centrality.betweenness_centrality)
- Closeness centrality, the average length of the shortest paths between a node and all other nodes [More info](https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.closeness_centrality.html#networkx.algorithms.centrality.closeness_centrality)
- Eigenvector centrality, the sum of the centrality values of the node neighbors [More info](https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.eigenvector_centrality.html#networkx.algorithms.centrality.eigenvector_centrality)

In addition to these measures, we decided to include other four graph centralities:

- Load centrality, a variant of betweeness centrality [More info](https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.load_centrality.html#networkx.algorithms.centrality.load_centrality)
- Local reaching centrality, the number of nodes that can be reached from a node [More info](https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.local_reaching_centrality.html#networkx.algorithms.centrality.local_reaching_centrality)
- Harmonic centrality, the harmonic mean of the shortest paths between a node and all other nodes [More info](https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.harmonic_centrality.html#networkx.algorithms.centrality.harmonic_centrality)
- Subgraph centrality, the number of subgraphs in the graph a node is part of [More info](https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.subgraph_centrality.html#networkx.algorithms.centrality.subgraph_centrality)

### Other topology features

Besides centrality, another topology measure is widely used in essentiality prediction in biological networks:

- Clustering coefficient, the amount of triplets/triangles (nodes connected by two or three undirected links) in the graph a node is part of [More info](https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.cluster.clustering.html#networkx.algorithms.cluster.clustering)

Since this measure is closely related to social network analysis and information retrieval, we decided to use two link analysis measures from this domain:

- HITS algorithm, which estimates a node value based on incoming and outgoing link scores [More info](https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.link_analysis.hits_alg.hits.html#networkx.algorithms.link_analysis.hits_alg.hits)
- PageRank, which ranks nodes based on the quality of incoming links [More info](https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.link_analysis.pagerank_alg.pagerank.html#networkx.algorithms.link_analysis.pagerank_alg.pagerank)

Two more metrics were implemented, the first one from graph theory:

- Length of a random maximal independent (stable) set, an independent set that is not a subset of any other independent set

And the last one designed specifically for gene essenciality prediction in metabolic networks:

- Damage, proposed by [Lemke et al. 2004](https://pdfs.semanticscholar.org/6530/2a07106438acb8e6d59891e3907f3915d8db.pdf)

## Citation

NAGAI, J. S.; SOUSA, H.; AONO, A. H.; LORENA, A. C.; KUROSHU, R. M. Gene Essentiality Prediction Using Topological Features From Metabolic Networks. In: 7th Brazilian Conference on Intelligent Systems (BRACIS), São Paulo, 2018 (No prelo).
