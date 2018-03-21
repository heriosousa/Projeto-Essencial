#! usr/bin/python3

from Bio.KEGG import REST as kegg
from Bio.KEGG.KGML import KGML_parser as parse
import networkx as nx
import itertools as it

class KeggParser():
    # parser constructor needs a pathway name
    def __init__(self, path=None, ref='ec'):
        # generator of KEGG Pathway network
        if path is not None:
            self.genes = parse.read(kegg.kegg_get(path, 'kgml'))
            self.genes_default = parse.read(kegg.kegg_get(
                                            "path:" + ref+path[-5:], 'kgml'))
        else:
            self.genes = None
            self.genes_default = None
        # ec and org_name
        self.ec_org = {}
        # ec and org_name
        self.ec_org_target = {}
        # path name
        self.path = path
        # Digraph structure for centrality
        self.grafo = nx.DiGraph()
        # dict for reactions substrates
        self.reactsub = {}
        # dict for reactions product
        self.reactpro = {}
        # dict for ec reactions
        self.ecreact = {}
        # dict for ec products
        self.ecproducts = {}
        # dict for ec substrates
        self.ecsubstrates = {}
        # dict for ec_pathway
        self.ecpathway = {}

    # using path to get a kgml pathway form and retrieve a relations
    def get_relations(self):
        '''
            O método é responsável por separar as relações por EC number
            fazendo o mapeamento com a via de referencia.
        '''
        for i in self.genes_default.relations:
            if i.entry1.name[0:4] != 'path' and i.entry1.name[0:4] != 'PATH':
                for ec in i.entry1.name.split():
                    self.ec_org[i.entry1.id] = ec
            if i.entry2.name[0:4] != 'path' and i.entry2.name[0:4] != 'PATH':
                for ec in i.entry2.name.split():
                    self.ec_org[i.entry2.id] = ec
        for i in self.genes.relations:
            if (i.entry1.name[0:4] != 'path' and
                    i.entry1.id in self.ec_org.keys()):
                for ec in self.ec_org[i.entry1.id].split():
                    self.ec_org_target[ec] = i.entry1.name
                    if ec not in self.ecpathway.keys():
                        self.ecpathway[ec] = [self.path[0]]
                    elif self.path[0] not in self.ecpathway[ec]:
                        self.ecpathway[ec].append(self.path[0])
            if (i.entry2.name[0:4] != 'path' and
                    i.entry2.id in self.ec_org.keys()):
                for ec in self.ec_org[i.entry2.id].split():
                    self.ec_org_target[ec] = i.entry2.name
                    if ec not in self.ecpathway.keys():
                        self.ecpathway[ec] = [self.path[0]]
                    elif self.path[0] not in self.ecpathway[ec]:
                        self.ecpathway[ec].append(self.path[0])

        # In relations if rel isn't a path add to ecreact dict
        for rel in self.genes.relations:
            if (rel.entry1.name[0:4] != 'path' and
                    rel.entry1.id in self.ec_org.keys()):
                if self.ec_org[rel.entry1.id] not in self.ecreact:
                    if (self.ec_org[rel.entry1.id] is not 'undefined' and
                            rel.entry1.reaction.split(" ") is not " "):
                        self.ecreact[self.ec_org[rel.entry1.id]] = str(rel.entry1.reaction).split(" ")
                else:
                    myelements1 = list(self.ecreact[self.ec_org[rel.entry1.id]])
                    myreacts1 = str(rel.entry1.reaction).split(" ")
                    for i in myreacts1:
                        if i not in myelements1:
                            myelements1.append(i)
                    self.ecreact[self.ec_org[rel.entry1.id]] = myelements1
            if (rel.entry2.name[0:4] != 'path' and
                    rel.entry2.id in self.ec_org.keys()):
                if self.ec_org[rel.entry2.id] not in self.ecreact:
                    if (self.ec_org[rel.entry2.id] is not 'undefined' and
                            rel.entry1.reaction.split(" ") is not " "):
                        self.ecreact[self.ec_org[rel.entry2.id]] = str(rel.entry2.reaction).split(" ")
                else:
                    myelements2 = list(self.ecreact[self.ec_org[rel.entry2.id]])
                    myreacts2 = str(rel.entry2.reaction).split(" ")
                    for i in myreacts2:
                        if i not in myelements2:
                            myelements2.append(i)
                    self.ecreact[self.ec_org[rel.entry2.id]] = myelements2

    # using retrieved relations to get a reactions
    def get_reactions(self):
        '''
        Utiliza de cada relação para obter informações das reações e
        assim obter substratos e produtos de acordo com o tipo de reação.
        '''
        # getting reactions products and substrates using reaction
        for rea in self.genes.reactions:
            reaction = rea.name
            reaction = reaction.split(" ")
            if rea.type == 'irreversible':
                substrates = []
                products = []
                for i in rea.substrates:
                    substrates.append(i.name)
                for j in rea.products:
                    products.append(j.name)
                for k in reaction:
                    self.reactsub[k] = substrates
                    self.reactpro[k] = products
            elif rea.type == 'reversible':
                subpro = []
                for i in rea.substrates:
                    subpro.append(i.name)
                for j in rea.products:
                    subpro.append(j.name)
                for k in reaction:
                    self.reactsub[k] = subpro
                    self.reactpro[k] = subpro

    # matching ec numbers and making a graph
    def matching_ec_reac(self):
        '''
        Faz a correspondencia entre as reações e as respectivas EC numbers
        '''
        for h in self.ecreact.items():
            enzime = h[0]
            reaction = h[1]
            p, s = [], []
            for i in reaction:
                if i in self.reactpro:
                    if isinstance(self.reactpro[i], str):
                        p.append(self.reactpro[i])
                    else:
                        for j in self.reactpro[i]:
                            p.append(j)
                if i in self.reactsub:
                    if isinstance(self.reactsub[i], str):
                        s.append(self.reactsub[i])
                    else:
                        for j in self.reactsub[i]:
                            s.append(j)
            self.ecsubstrates[enzime] = s
            self.ecproducts[enzime] = p

    # buiding a graph using kegg informations
    def building_graph(self, opt=1):
        '''
        Monta o grafo a partir dos produtos e reagentes
        '''
        if opt == 1:
            self.get_relations()
            self.get_reactions()
    #    self.selection_len()
    #    self.reac_len()
        self.matching_ec_reac()
        for i in self.ecreact.keys():
            self.grafo.add_node(i)
        reactions = set()
        for i in self.ecreact.values():
            for j in i:
                reactions.add(str(j))
        for i in self.reactpro.keys():
            if i not in reactions:
                self.ecreact[i] = i
        pairs = []
        for pair in it.combinations(list(self.ecreact.keys()), 2):
            pairs.append(pair)
        cnt = 0
        print("Building graphs.")
        for i in pairs:
            gen1 = i[0]
            gen2 = i[1]
            dir1 = 0
            dir2 = 0
            if('ec' in gen1 and 'ec' in gen2):
                for j in self.ecproducts[gen1]:
                    if j in self.ecsubstrates[gen2]:
                        dir1 += 1
                for j in self.ecproducts[gen2]:
                    if j in self.ecsubstrates[gen1]:
                        dir2 += 1
            elif('ec' in gen1 and 'rn' in gen2):
                for j in self.ecproducts[gen1]:
                    if j in self.reactsub[gen2]:
                        dir1 += 1
                for j in self.reactpro[gen2]:
                    if j in self.ecsubstrates[gen1]:
                        dir2 += 1
            elif('rn' in gen1 and 'ec' in gen2):
                for j in self.reactpro[gen1]:
                    if j in self.ecsubstrates[gen2]:
                        dir1 += 1
                for j in self.ecproducts[gen2]:
                    if j in self.reactsub[gen1]:
                        dir2 += 1
            elif('rn' in gen1 and 'rn' in gen2):
                for j in self.reactpro[gen1]:
                    if j in self.reactsub[gen2]:
                        dir1 += 1
                for j in self.reactpro[gen2]:
                    if j in self.reactsub[gen1]:
                        dir2 += 1
            if dir1 >= 1 and dir2 == 0:
                cnt += 1
                self.grafo.add_edge(gen1, gen2, color='blue')
            if dir1 == 0 and dir2 >= 1:
                cnt += 1
                self.grafo.add_edge(gen1, gen2, color='blue')
            if dir1 >= 1 and dir2 >= 1:
                cnt += 1
                self.grafo.add_edge(gen1, gen2, color='green4')
                self.grafo.add_edge(gen2, gen1, color='green4')

