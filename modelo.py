import kegg_parser as k
import networkx as nx
import sys
import os
from Bio.KEGG import REST as r
from Bio.KEGG.KGML import KGML_parser as parse
import time


def make_directories(name):
    os.mkdir('./'+str(name)[5:])
    os.mkdir('./'+str(name)[5:]+'/img')
    os.mkdir('./'+str(name)[5:]+'/centrality')


def img_ploter(grafo, opt, name=None, i=None):
    if opt == 1:
        img1 = nx.nx_agraph.to_agraph(grafo.grafo)
        nx.nx_agraph.write_dot(grafo.grafo, "completo"+name+".dot")
        img1.layout(prog='dot')
        if grafo.path is not None:
            img1.draw('./'+str(grafo.path)[5:] +
                      '/img/'+str(grafo.path[5:])+"_complete.png")
        else:
            img1.draw('final'+name+'.pdf')
    elif opt == 2:
        img1 = nx.nx_agraph.to_agraph(grafo)
        img1.layout(prog='dot')
        img1.draw('./'+str(name)[5:]+'/img/'+str(name)+'_no_'+str(i)+'.png')


def centrality(grafo, nodes=None):
    '''
    Cálculo da centralidade: Opção de plot
    '''
    d_centrality = nx.degree_centrality(grafo)
    c_centrality = nx.closeness_centrality(grafo)
    b_centrality = nx.betweenness_centrality(grafo)
    if nodes is None:
        return ([d_centrality, c_centrality, b_centrality])
    else:
        d_aux = {}
        c_aux = {}
        b_aux = {}
        for i in nodes:
            d_aux[i] = d_centrality[i]
            c_aux[i] = c_centrality[i]
            b_aux[i] = b_centrality[i]
        return(d_aux, c_aux, b_aux)


def articulation(grafo):
    return (list(nx.articulation_points(grafo.grafo.to_undirected())))


def clust_coefficient(grafo):
    return (nx.average_clustering(grafo.grafo.to_undirected()))


def MBB(grafo):
    MBB = nx.strongly_connected_components(grafo.grafo)
    MBB1 = []
    count = 0
    corresp = {}
    for blocks in MBB:
        corresp[count] = blocks
        MBB1.append(blocks)
        count = count + 1
    grafo.grafo = nx.condensation(grafo.grafo, scc=MBB1)
    return (grafo, corresp)


def get_network(org, opt='ec'):
    # Creating a Parser Object
    graph = k.KeggParser()
    # Store pathways that doesn't have EC numbers
    error = []
    # Getting organism
    list1 = r.kegg_list('pathway', org).read()
    list1 = list1.split('\n')
    list1.remove('')
    print('Retrieving data from KEGG PATHWAY database. '+str(time.ctime()))

    # For each path getting enzymes and reactions
    for path in list1:
        try:
            path = path.split('\t')
            # print (path[0])
            graph.genes = parse.read(r.kegg_get(path[0], 'kgml'))
            graph.genes_default = parse.read(r.kegg_get("path:" +
                                             opt+path[0][-5:], 'kgml'))
            graph.path = path
        except Exception:
            error.append(path[0])
            continue
        # print ("getting relations")
        graph.get_relations()
        # print ("getting reaction")
        graph.get_reactions()
    # print ('Unretrieved data',error)
    genes = 0
    for i in graph.ec_org_target.items():
        genes += len(i[1].split())

    # print (graph.ec_org_target.keys())
    # Building Graph
    graph.building_graph(2)
    return (graph)


def get_pathway(path):
    try:
        n_path = k.KeggParser(path)
    except Exception as e:
        print('Erro: Favor inserir o nome da via path:via')
        sys.exit(-1)
    n_path.building_graph()
    return (n_path)


if __name__ == '__main__':
    main()
