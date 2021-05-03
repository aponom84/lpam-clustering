
#!/usr/bin/env python

import networkx as nx
import cdlib
from cdlib import ensemble
from cdlib import algorithms
from cdlib import evaluation
from cdlib import viz
from cdlib import NodeClustering
from sklearn.metrics import adjusted_mutual_info_score
import typing
from typing import List
from multiprocessing import Pool
from pyclustering.cluster.kmedoids import kmedoids
from pyclustering.cluster.center_initializer import kmeans_plusplus_initializer
from pyclustering.cluster import cluster_visualizer
from pyclustering.utils import read_sample
from pyclustering.samples.definitions import FCPS_SAMPLES
import numpy as np
from tqdm import tqdm
from collections import defaultdict
import sys
import traceback
import warnings
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from cdlib import evaluation
import itertools
import copy
import time

warnings.filterwarnings("ignore")
NUM_THREADS = 24

def isOverlap(community)->bool:
    return any([ set(a).intersection(b) for a,b in itertools.combinations(community, 2)])

def getOverlappingNumber(community)->int:
    overpaping = set()
    for a,b in itertools.combinations(community, 2):
        overpaping = set.union(overpaping,  set(a).intersection(b))
    return len(overpaping)


def copra(graph, param)->List[List[int]]:
    inputFile = 'copra_input.txt'
    z = get_ipython().getoutput(f'rm -f {inputFile} || true')
    nx.write_edgelist(graph, inputFile,  data=False)
    outputDir = "COPRA"
    outputFile = outputDir + "/" + "clusters-" + inputFile.split('/')[-1]
    get_ipython().system(f'rm -rf {outputDir}')
    get_ipython().system(f'mkdir {outputDir}')
    cmd = f'(cd {outputDir} && java -cp ../../related_methods/OSLOM2/copra.jar COPRA ../{inputFile} -v {param} -repeat 100 -mo -nosplit)'
    # print(f'running: [{cmd}]')
    tmp = get_ipython().getoutput(cmd)
    com = []
    with open(outputFile) as f:
        lines = f.readlines()
        com = [[int(a) for a in line.split()] for line in lines]

    return NodeClustering(communities=com, graph=graph, method_name='copra', method_parameters=param, overlap=isOverlap(com))


def oslom(graph, t=0.5, cp=0.5, seed=13)->List[List[int]]:
    '''
    graph - NetworkX
    t - threashould ?
    cp  - xxx?
    seed - random seed
    '''
    inputFile = 'oslom_tmp_input.txt'
    z = get_ipython().getoutput(f'rm -f {inputFile} || true')
    z = get_ipython().getoutput('rm -f tp || true')
    z = get_ipython().getoutput('rf -rf oslom_tmp_input.txt_oslo_files || true')
    nx.write_edgelist(graph, inputFile, data=False)
    param = f'-t {t} -cp {cp} -seed {seed}'
    cmd = f'../related_methods/OSLOM2/oslom_undir -f {inputFile} -uw {param}'
    tmp = get_ipython().getoutput(cmd)
    lines = get_ipython().getoutput('cat oslom_tmp_input.txt_oslo_files/partitions_level_0')
    com = [ [int(v) for v in line.split()] for line in lines if not line.startswith('#')]
    return NodeClustering(communities=com, graph=graph, method_name='oslom', method_parameters=param, overlap=isOverlap(com))
def GCE(graph, param):
    inputFile = 'GCE_input.txt'
    z = get_ipython().getoutput(f'rm -f {inputFile} || true')
    z = get_ipython().getoutput('rm -f tp || true')
    nx.write_edgelist(graph, inputFile, data=False)
    out = get_ipython().getoutput(f'../related_methods/GCECommunityFinder/build/GCECommunityFinder {inputFile} 4 0.9 {param} .75')
    com = [[int(a) for a in line.split()] for line in out[out.index('Finished')+1:]]
    return NodeClustering(communities=com, graph=graph, method_name='GCE', method_parameters=param, overlap=isOverlap(com))


disntaces_cache = {}
def lpam_python_wrong_amp(graph, k=7, threshold=0.5, seed=0):
    def getAmp(G):
        """
        Returns the matrix of amplified commute distance
        """
        verts = list(G.nodes)
        n = len(verts)

        vol = nx.volume(G,verts)

        # use NetworkX to get Laplacian
        L = nx.laplacian_matrix(G)
        L = L.todense()
        C_AMP = np.zeros([n,n])

        #get Moore-Penrose pseudo inverse
        L_MP_pseudo_inv = np.linalg.pinv(L, rcond=1e-5, hermitian=True)
        for i in tqdm(range(n-1)) :
            e_i = np.zeros([n,1])
            e_i[i, 0] = 1
            for j in range(i+1, n):
                e_j = np.zeros([n,1])
                e_j[j, 0] = 1
                a = (e_i - e_j)
                b = L_MP_pseudo_inv @ a
                c_ij = ( vol * (np.transpose(a) @ b) )[0,0]
                c_ij_amp = (c_ij / vol) - (1/G.degree(verts[i])) - (1/G.degree(verts[j])) + ( 2/( G.degree(verts[i]) *G.degree(verts[j]) ) )
                C_AMP[i,j] = c_ij_amp
                C_AMP[j,i] = c_ij_amp

        return C_AMP


    line_graph = nx.line_graph(graph)

    global disntaces_cache
    D = disntaces_cache[graph] if graph in disntaces_cache else getAmp(line_graph)
    disntaces_cache[graph] = D

    _n = len(line_graph.nodes())
    np.random.seed(0)
    initial_medoids = np.random.choice(_n, k, replace=False)
    kmedoids_instance = kmedoids(D, initial_medoids, data_type='distance_matrix')
    # run cluster analysis and obtain results
    kmedoids_instance.process()

    clusters = kmedoids_instance.get_clusters()
    medoids = kmedoids_instance.get_medoids()

    final_clusters = {}
    for c_i, c in enumerate(clusters):
        for line_vertex in c:
            source, target = list(line_graph.nodes()) [line_vertex]
            if source not in final_clusters:
                final_clusters[source] = []
            final_clusters[source].append(c_i)
            if target not in final_clusters:
                final_clusters[target] = []

            final_clusters[target].append(c_i)

    res_clusters = {}
    for v, l in final_clusters.items():
        degree = len(l)
        res = defaultdict(list)
        for x in l: res[x].append(x)
        covering = np.zeros(k)
        for c_i, l in res.items():
            covering[c_i] = len(l)/degree

        res_clusters[v] = covering

    _res_clusters = [ [] for i in range(k) ]

    for v, l in res_clusters.items():
        for i in range(k):
            if l[i] >= threshold:
                _res_clusters[i].append(v)

    return NodeClustering(communities=_res_clusters, graph=graph, method_name='lpam_amp_wrong', overlap = True)


disntaces_cache_amp = {}
res_clusters_cache_amp = {}
def lpam_python_amp(graph, k=7, threshold=0.5, seed=0):
    def getAmp(G):
        """
        Returns the matrix of amplified commute distance
        """

        verts = list(G.nodes)
        n = len(verts)

        #get adj matrix
        A = nx.adjacency_matrix(G)
        A = A.todense()

        # use NetworkX to get Laplacian
        L = nx.laplacian_matrix(G)
        L = L.todense()
        Gamma = L + (1/n) * np.ones([n,n])
        C_AMP = np.zeros([n,n])

        #get Moore-Penrose pseudo inverse
        Gamma_pinv = np.linalg.pinv(Gamma, rcond=1e-4)
        # for i in tqdm(range(n-1)) :
        for i in range(n):
            for j in range(i+1,n):
                r_ij = Gamma_pinv[i,i] + Gamma_pinv[j,j] - 2 * Gamma_pinv[i,j] #res dist
                d_i = G.degree(list(G.nodes())[i])
                d_j = G.degree(list(G.nodes())[j])
                if d_i !=0 and d_j !=0:
                    s_ij = r_ij - (1 / d_i) - (1 / d_j)
                    w_ij = A[i,j]
                    w_ii = A[i,i]
                    w_jj = A[j,j]
                    u_ij = ( ( 2 * w_ij) / (d_i * d_j) ) - ( w_ii / (d_i**2) ) - ( w_jj / (d_j**2) )
                    C_AMP[i,j] = s_ij + u_ij
                    C_AMP[j,i] = s_ij + u_ij

                else:
                    C_AMP[i,j] =  np.NaN
                    C_AMP[j,i] =  np.NaN
            #end for j
        #end for i
        return C_AMP


    line_graph = nx.line_graph(graph)


    global disntaces_cache_amp
    global res_clusters_cache_amp
    if (graph,k) in res_clusters_cache_amp:
        res_clusters = res_clusters_cache_amp[(graph,k)]
    else:
        D = disntaces_cache_amp[graph] if graph in disntaces_cache_amp else getAmp(line_graph)
        disntaces_cache_amp[graph] = D

        _n = len(line_graph.nodes())
        np.random.seed(0)
        initial_medoids = np.random.choice(_n, k, replace=False)
        kmedoids_instance = kmedoids(D, initial_medoids, data_type='distance_matrix')
        # run cluster analysis and obtain results
        kmedoids_instance.process()

        clusters = kmedoids_instance.get_clusters()
        medoids = kmedoids_instance.get_medoids()

        final_clusters = {}
        for c_i, c in enumerate(clusters):
            for line_vertex in c:
                source, target = list(line_graph.nodes())[line_vertex]
                if source not in final_clusters:
                    final_clusters[source] = []
                final_clusters[source].append(c_i)
                if target not in final_clusters:
                    final_clusters[target] = []

                final_clusters[target].append(c_i)

        res_clusters = {}
        for v, l in final_clusters.items():
            degree = len(l)
            res = defaultdict(list)
            for x in l: res[x].append(x)
            covering = np.zeros(k)
            for c_i, l in res.items():
                covering[c_i] = len(l)/degree

            res_clusters[v] = covering
        res_clusters_cache_amp[(graph,k)] = res_clusters

    _res_clusters = [ [] for i in range(k) ]

    for v, l in res_clusters.items():
        for i in range(k):
            if l[i] >= threshold:
                _res_clusters[i].append(v)

    return NodeClustering(communities=[c for c in _res_clusters if len(c)>0], graph=graph, method_name='lpam_amp', overlap = True)

disntaces_cache_cm = {}
res_clusters_cache_cm = {}
def lpam_python_cm(graph, k=7, threshold=0.5, seed=0):
    def getCommuteDistace(G):
        """
        Returns the matrix of commute distance
        """
        verts = list(G.nodes)
        n = len(verts)
        vol = nx.volume(G,verts)

        #get adj matrix
        A = nx.adjacency_matrix(G)
        A = A.todense()

        # use NetworkX to get Laplacian
        L = nx.laplacian_matrix(G)
        L = L.todense()
        Gamma = L + (1/n) * np.ones([n,n])
        CM = np.zeros([n,n])

        #get Moore-Penrose pseudo inverse
        Gamma_pinv = np.linalg.pinv(Gamma, rcond=1e-4)
        # for i in tqdm(range(n-1)) :
        for i in range(n):
            for j in range(i+1,n):
                CM[i,j] = vol*(Gamma_pinv[i,i] + Gamma_pinv[j,j] - 2 * Gamma_pinv[i,j])
                CM[j,i] = CM[i,j]
        return CM


    line_graph = nx.line_graph(graph)

    global disntaces_cache_cm
    global res_clusters_cache_cm
    res_clusters = {}
    if (graph,k) in res_clusters_cache_cm:
        res_clusters = res_clusters_cache_cm[(graph,k)]
    else:
        # print('calculating distance')
        D = disntaces_cache_cm[graph] if graph in disntaces_cache_cm else getCommuteDistace(line_graph)
        disntaces_cache_cm[graph] = D
        # print('solve p-median problems')

        _n = len(line_graph.nodes())
        np.random.seed(0)
        initial_medoids = np.random.choice(_n, k, replace=False)
        kmedoids_instance = kmedoids(D, initial_medoids, data_type='distance_matrix')
        # run cluster analysis and obtain results
        kmedoids_instance.process()

        clusters = kmedoids_instance.get_clusters()
        medoids = kmedoids_instance.get_medoids()

        final_clusters = {}
        for c_i, c in enumerate(clusters):
            for line_vertex in c:
                source, target = list(line_graph.nodes()) [line_vertex]
                if source not in final_clusters:
                    final_clusters[source] = []
                final_clusters[source].append(c_i)
                if target not in final_clusters:
                    final_clusters[target] = []

                final_clusters[target].append(c_i)


        for v, l in final_clusters.items():
            degree = len(l)
            res = defaultdict(list)
            for x in l: res[x].append(x)
            covering = np.zeros(k)
            for c_i, l in res.items():
                covering[c_i] = len(l)/degree

            res_clusters[v] = covering
        res_clusters_cache_cm[(graph,k)] = res_clusters

    _res_clusters = [ [] for i in range(k) ]
    for v, l in res_clusters.items():
        for i in range(k):
            if l[i] >= threshold:
                _res_clusters[i].append(v)

    return NodeClustering(communities=[c for c in _res_clusters if len(c)>0], graph=graph, method_name='lpam_cm', overlap = True)

def complete_partition(partition, g, mode='new_cluster'):
    pc = partition
    if type(partition) != list:
        pc = partition.communities
    nodes_cluster = {node:nid for nid, cluster in enumerate(pc) for node in cluster}
    def getFixedClusterId(i, nodes_cluster, addition_cluster_id):
        if i in nodes_cluster:
            return nodes_cluster[i]

        if mode == 'new_cluster':
            addition_cluster_id = addition_cluster_id + 1

        return addition_cluster_id-1
    clusters = [ []  for i in range(len(pc)+1) ]

    for i in range( len(g.nodes()) ):
        clusters[getFixedClusterId(i, nodes_cluster, len(pc))].append(i)

    return [c for c in clusters if c]

def save_community(reference, path):
    with open(path, 'w') as f:
        f.write('\n'.join([' '.join( [str(i) for i in c]) for c in reference]))


def read_com(groundTruth):
    file_content = get_ipython().getoutput('cat {groundTruth}')
    com = [[int(i) for i in line.split()] for line in file_content]
    return com
def read_graph(inputFile, groundTruth, shift=0):
    file_content = get_ipython().getoutput('cat {inputFile}')
    edges_list = [[int(i)+shift for i in line.split()] for line in file_content]
    g = nx.Graph()
    g.add_edges_from(edges_list)

    datasetName = inputFile.split('/')[-2]
    g.name = datasetName
    com = read_com(groundTruth)
    return g, com


def generateFARZ(n=200, m=5, k=5, beta = 0.9, alpha=0.2, gamma=0.5):
    while True:
        get_ipython().system('rm -f network.gml')
        get_ipython().system(f'python ../FARZ/src/FARZ.py -n {n} -m {m} -k {k} --beta {beta} --alpha {alpha} --gamma {gamma} -r {k} -q 0.05 --format comb')
        g = nx.read_gml("network.gml")
        com = read_com('network.lgt')
        get_ipython().system(f'cp network.gml ../datasets/FARZ/farz_{beta:.2f}.gml')
        get_ipython().system(f'cp network.lgt ../datasets/FARZ/farz_{beta:.2f}.lgt')

        if nx.is_connected(g):
            return g, com


def getBench(N=200, k = 15, maxk=50, mu=0.1, minc=5, maxc = 50, on= 20, om=2):
    '''
    Generate Lancichinetti's benchmark graphs
    Flags:
    -N		number of nodes
    -k		average degree
    -maxk		maximum degree
    -mu		mixing parameter
    -t1		minus exponent for the degree sequence
    -t2		minus exponent for the community size distribution
    -minc		minimum for the community sizes
    -maxc		maximum for the community sizes
    -on		number of overlapping nodes
    -om		number of memberships of the overlapping nodes
    -C              [average clustering coefficient]
    '''
    while True:
        print(f"Generate Lancichinetti'benchmark graphs N={N}, k = {k}, maxk={maxk}, mu={mu}, minc={minc}, maxc = {maxc}")
        out = get_ipython().getoutput('rm -f network.dat community.dat')
        out = get_ipython().getoutput(f'../binary_networks_benchmark/benchmark -N {N} -k {k} -maxk {maxk} -mu {mu} -minc {minc} -maxc {maxc} -on {on} -om {om} -t1 2 -t2 1')
        print(out)
        # out = get_ipython().getoutput('../binary_networks_benchmark/benchmark -N 30 -k 4 -maxk 12 -mu 0.1 -minc 5 -maxc 20 -on 2 -om 2 -t1 2 -t2 1')
        get_ipython().system(f'cp network.dat ../datasets/bench/bench_{mu:.2f}.data')
        get_ipython().system(f'cp community.dat ../datasets/bench/bench_truth_{mu:.2f}.dat')

        g = nx.Graph()
        with open('network.dat') as f:
            for line in f.readlines():
                g.add_edge(int(line.split()[0])-1, int(line.split()[1]) -1)

        com = defaultdict(list)
        with open('community.dat') as f:
            for line in f.readlines():
                splits = line.split()
                p = int(splits[0])
                for i in splits[1:]:
                    com[p].append(int(i) - 1)
                # com[ int(line.split()[1])-1].append(int(line.split()[0])-1)

        if nx.is_connected(g):
            return g, list(com.values())

def getAllScoresDict(g, _reference, _communities, executionTime):
    scores = {}
    scores['time'] = executionTime
    reference = copy.deepcopy(_reference)
    reference.communities = complete_partition(reference.communities, g, mode='new_cluster')
    communities = copy.deepcopy(_communities)
    communities.communities = complete_partition(communities.communities, g, mode='new_cluster')



    # scores['adjusted_mutual_information'] = evaluation.adjusted_mutual_information(reference,communities).score

    # returns MatchingResult object
    # scores['adjusted_rand_index'] = evaluation.adjusted_rand_index(reference,communities).score
    # Compute the average F1 score of the optimal algorithms matches among the partitions in input.
    try:
        scores['f1'] = evaluation.f1(reference, communities).score
    except:
        scores['f1'] = np.nan
    # Compute the Normalized F1 score of the optimal algorithms matches among the partitions in input.
    try:
        scores['nf1'] = evaluation.nf1(reference, communities).score
    except:
        scores['nf1'] = np.nan
    # Normalized Mutual Information between two clusterings.
    # scores['normalized_mutual_information'] = evaluation.normalized_mutual_information(reference, communities)[0]
    # Index of resemblance for overlapping, complete coverage, network clusterings.
    try:
        scores['omega'] = evaluation.omega(reference, communities).score
    except:
        scores['omega'] = np.nan
    # Overlapping Normalized Mutual Information between two clusterings.
    try:
        scores['overlapping_normalized_mutual_information_LFK'] = evaluation.overlapping_normalized_mutual_information_LFK(reference, communities)[0]
    except:
        scores['overlapping_normalized_mutual_information_LFK']  = np.nan
    # Overlapping Normalized Mutual Information between two clusterings.
    # scores['overlapping_normalized_mutual_information_MGH'] =  evaluation.overlapping_normalized_mutual_information_MGH(reference, communities)[0]
    # Variation of Information among two nodes partitions.
    # scores['variation_of_information'] =  evaluation.variation_of_information(reference, communities)[0]

    # scores['avg_distance'] = evaluation.avg_distance(g,communities, summary=True)
    try:
        scores['avg_embeddedness'] = evaluation.avg_embeddedness(g,communities, summary=True).score
    except:
        scores['avg_embeddedness'] = np.nan
    try:
        scores['average_internal_degree'] = evaluation.average_internal_degree(g,communities, summary=True).score
    except:
        scores['average_internal_degree'] = np.nan
    # scores['avg_transitivity']  = evaluation.avg_transitivity(g,communities, summary=True)
    # Fraction of total edge volume that points outside the community.
    try:
        scores['conductance']  = evaluation.conductance(g,communities, summary=True).score
    except:
        scores['conductance'] = np.nan
    # Fraction of existing edges (out of all possible edges) leaving the community.
    try:
        scores['cut_ratio']  = evaluation.cut_ratio(g,communities, summary=True).score
    except:
        scores['cut_ratio'] = np.nan

    # Number of edges internal to the community
    try:
        scores['edges_inside']  = evaluation.edges_inside(g,communities, summary=True).score
    except:
        scores['edges_inside'] = np.nan
    # Number of edges per community node that point outside the cluster
    try:
        scores['expansion']  = evaluation.expansion(g,communities, summary=True).score
    except:
        scores['expansion'] = np.nan
    # Fraction of community nodes of having internal degree higher than the median degree value.
    try:
        scores['fraction_over_median_degree']  = evaluation.fraction_over_median_degree(g,communities, summary=True).score
    except:
        scores['fraction_over_median_degree'] = np.nan
    # The hub dominance of a community is defined as the ratio of the degree of its most connected node w.r.t. the theoretically maximal degree within the community.
    # scores['hub_dominance']  = evaluation.hub_dominance(g,communities, summary=True)
    # The internal density of the community set.
    try:
        scores['internal_edge_density'] = evaluation.internal_edge_density(g,communities, summary=True).score
    except:
        scores['internal_edge_density'] = np.nan
    # Normalized variant of the Cut-Ratio
    try:
        scores['normalized_cut']  = evaluation.normalized_cut(g,communities, summary=True).score
    except:
        scores['normalized_cut'] = np.nan
    # Maximum fraction of edges of a node of a community that point outside the community itself.
    # scores['max_odf']  = evaluation.max_odf(g,communities, summary=True)
    # Average fraction of edges of a node of a community that point outside the community itself.
    # scores['avg_odf']  = evaluation.avg_odf(g,communities, summary=True)
    # Fraction of nodes in S that have fewer edges pointing inside than to the outside of the community.
    # scores['flake_odf']  = evaluation.flake_odf(g,communities, summary=True)
    # The scaled density of a community is defined as the ratio of the community density w.r.t. the complete graph density.
    try:
        scores['scaled_density']  = evaluation.scaled_density(g,communities, summary=True).score
    except:
        scores['scaled_density'] = np.nan
    # Significance estimates how likely a partition of dense communities appear in a random graph.
    try:
        scores['significance'] = evaluation.significance(g,communities).score
    except:
        scores['significance'] = np.nan
    # Size is the number of nodes in the community
    try:
        scores['size']  = evaluation.size(g,communities, summary=True).score
    except:
        scores['size'] = np.nan
    # Surprise is statistical approach proposes a quality metric assuming that edges between vertices emerge randomly according to a hyper-geometric distribution.
    # According to the Surprise metric, the higher the score of a partition, the less likely it is resulted from a random realization, the better the quality of the community structure.
    try:
        scores['surprise'] = evaluation.surprise(g,communities).score
    except:
        scores['surprise'] = np.nan

    try:
        scores['modularity_density'] = evaluation.modularity_density(g,communities).score
    except:
        scores['modularity_density'] = np.nan

    # Fraction of community nodes that belong to a triad.
    # scores['triangle_participation_ratio']  = evaluation.triangle_participation_ratio(g,communities, summary=True)
    # Purity is the product of the frequencies of the most frequent labels carried by the nodes within the communities
    # scores['purity'] = evaluation.purity(communities)
    return scores

def getPoint(args):
    g, name, reference, clustering_method, parameters = args
    results = []
    try:
        startTime = time.time()
        for communities in ensemble.grid_execution(graph=g, method=clustering_method, parameters=parameters):
            executionTime = time.time() - startTime
            results.append((name, getAllScoresDict(g, reference, communities, executionTime)))
            startTime = time.time()
        return results
    except:
        print('graph:', name)
        print('parameters:', parameters)
        print(traceback.format_exc())
        raise

def getResultsParallel(clustering_method, graphs, names, references, parameters_list, num_threads=NUM_THREADS):
    print(clustering_method)
    results = []
    # tasks = zip(graphs, names, references, [clustering_method]*len(graphs), [scoring_method]*len(graphs),  parameters_list)
    tasks = zip(graphs, names, references, [clustering_method for i in range(len(graphs))],   parameters_list)
    # print('--------------------------- tasks: ----------------')
    # print(list(tasks))
    with Pool(processes=num_threads) as pool:
        with tqdm(total=len(graphs), position=0, leave=True) as pbar:
            for res in pool.imap_unordered(getPoint, tasks):
                results.extend(res)
                pbar.update()

    return results

import sys
import warnings
warnings.filterwarnings("ignore")
def getResults(clustering_method, graphs, names, references, parameters_list):
    print(f'Testing method: {clustering_method}')
    results = []
    for g, name, reference, parameters in tqdm(list(zip(graphs, names, references, parameters_list))):
        print(f'graph: {name}')
        try:
            startTime = time.time()
            for communities in ensemble.grid_execution(graph=g, method=clustering_method, parameters=parameters):
                executionTime =  time.time() - startTime
                results.append((name, getAllScoresDict(g, reference, communities, executionTime)))
                startTime = time.time()
        except:
            print('graph:', name)
            print('parameters:', parameters)
            print("Unexpected error:", sys.exc_info()[0])
            raise

    return results


if __name__ == '__main__':
    school_g, school_gt = read_graph('../datasets/school_friendship/school-2.dat',  '../datasets/school_friendship/truth-school_overlap.dat', shift=-1)
    karate_g, karate_gt = read_graph("../datasets/karate/karate.dat",  "../datasets/karate/truth_karate.dat", shift=0)
    adj_g, adj_gt = read_graph("../datasets/adjnoun/adjnoun.dat",  "../datasets/adjnoun/truth_adjnoun.dat", shift=0)
    # football_g, football_gt = read_graph("../datasets/football/footballTSEinput_original.dat",  "../datasets/football/truth_footballTSEinput.dat", shift=-1)
    football_g, football_gt = read_graph("../datasets/football/footballTSEinput_original.dat",  "../datasets/football/truth_footballTSEinput_overlap.dat", shift=-1)
    poolbooks_g, poolbooks_gt = read_graph("../datasets/polbooks/polbooks.dat",  "../datasets/polbooks/truth_polbooks.dat", shift=0)
    lattice8x8_g, lattice8x8_gt = read_graph( "../datasets/lattice_8x8/lattice_8x8.dat",  "../datasets/lattice_8x8/truth_4_b.dat", shift=0)

    graphs = [ school_g, karate_g, adj_g, football_g, poolbooks_g, lattice8x8_g ]
    names = 'school karate adj foolball poolbooks lattice8x8'.split()
    references = [school_gt, karate_gt, adj_gt, football_gt, poolbooks_gt, lattice8x8_gt]

    bench_30_g, bench_30_gt = read_graph('../datasets/bench_30/bench_30_network.dat',  '../datasets/bench_30/bench_30_truth.dat', shift=0)
    bench_40_g, bench_40_gt = read_graph('../datasets/bench_40/bench_40_network.dat',  '../datasets/bench_40/bench_40_truth.dat', shift=0)
    bench_50_g, bench_50_gt = read_graph('../datasets/bench_50/bench_50_network.dat',  '../datasets/bench_50/bench_50_truth.dat', shift=0)
    bench_60_g, bench_60_gt = read_graph('../datasets/bench_60/bench_60_network.dat',  '../datasets/bench_60/bench_60_truth.dat', shift=0)
    bench_60_dense_g, bench_60_dense_gt = read_graph('../datasets/bench_60_dense/bench_60_dense_network.dat',  '../datasets/bench_60_dense/bench_60_dense_truth.dat', shift=0)

    graphs.extend([ bench_30_g, bench_40_g, bench_50_g, bench_60_g, bench_60_dense_g])
    names.extedn( 'bench_30 bench_40 bench_50 bench_60, bench_60_dense'.split() )
    references.extend = ([ bench_30_gt, bench_40_gt, bench_50_gt, bench_60_gt, bench_60_dense_gt])




    results_partitioning = {}
    results_scoring = {} #here we will store all scoring/covering results?


    # generate Lancichinetti's benchmark graphs
    for mu in np.arange(0.1, 1.0, 0.05):
        g, com = getBench(N=200, mu = mu )
        graphs.append(g)
        names.append(f'bench_{mu:.2f}')
        references.append(com)


    # generate fraz graph
    for beta in np.arange(0.1, 1.0, 0.05):
        g, com = generateFARZ(beta = beta)
        graphs.append(g)
        names.append(f'farz_{beta:.2f}')
        references.append(com)



    # generate planted partition random graphs
    numClust = 6
    nodesPerClust = 5
    numNodesTot = numClust * nodesPerClust

    p_in_array = np.array([0.2, 0.4, 0.6, 0.8, 1])
    p_out_array = np.array([0.10, 0.15, 0.2, 0.25])

    for p_out in p_out_array:
        for p_in  in p_in_array:
            seed = 0
            while True:
                G = nx. _graph(numClust, nodesPerClust, p_in, p_out, seed=seed, directed=False)
                seed = seed + 1
                if nx.is_connected(G):
                    break
            C = [ [  (i + c*nodesPerClust) for i in range(nodesPerClust) ]  for c in range(numClust)  ]
            graphs.append(G)
            references.append(C)
            names.append(f'PP_{p_in:.2f}_{p_out:.2f}')
            # nx.write_gml(G, f'../datasets/pp/PP_{p_in}_{p_out}.gml')
            nx.write_edgelist(G, f'../datasets/pp/PP_{p_in:.2f}_{p_out:.2f}.dat')
            save_community(C, f'../datasets/pp/PP_{p_in:.2f}_{p_out:.2f}_truth.dat')



    # Generate Stohastic Block Model
    numClust = 7 # num clusters
    avgClusterSize = 10 #avarage nodes per cluster
    nodesPerClust = [int(i) for i in ((np.random.rand(numClust)-0.5)*avgClusterSize*0.1 + avgClusterSize)]
    numNodes = np.sum(nodesPerClust)
    p_in_array = np.array([0.2, 0.4, 0.6, 0.8, 1])
    p_out_array = np.array([0.10, 0.15, 0.2, 0.25])

    for p_in in p_in_array:
        for p_out in p_out_array:
            probs = (np.random.rand(len(nodesPerClust), len(nodesPerClust)) - 0.5)*p_out*0.1 + p_out
            for i in range(len(nodesPerClust)):
                probs[i][i] = p_in
                for j in range(i+1, len(nodesPerClust)):
                    probs[j][i] = probs[i][j]
            G, C = None, None
            while True:
                G = nx.stochastic_block_model( sizes = nodesPerClust, p=probs,  seed=0, directed=False)
                C = []
                i = 0
                for c in nodesPerClust:
                    C.append( [  (i + j) for j in range(c) ] )
                    i = i + c
                if nx.is_connected(G):
                    break
            graphs.append(G)
            references.append(C)
            names.append(f'SBM_{p_in}_{p_out}')
            # nx.write_gml(G, f'../datasets/sbm/SBM_{p_in}_{p_out}.gml')
            nx.write_edgelist(G, f'../datasets/sbm/SBM_{p_in:.2f}_{p_out:.2f}.dat')
            save_community(C, f'../datasets/sbm/SBM_{p_in:.2f}_{p_out:.2f}_truth.dat')


    # check all graph for connectivity
    for g, name in zip(graphs, names):
        if not nx.is_connected(g):
            print(name, 'is not connected')

    references = [NodeClustering(communities=r, graph=g, method_name='reference', method_parameters=None, overlap=isOverlap(r)) for r,g in zip(references, graphs) ]

    ###########################################################################
    ### Testing Methods
    ###########################################################################

    # graphs = [school_g]
    # names = ['school']
    # references = [school_gt]


    results_scoring['lpam_python_cm'] = getResults(
        clustering_method=lpam_python_cm,
        graphs=graphs,
        names=names,
        references=references,
        parameters_list = [ [ensemble.Parameter(name="k", start=max(1,len(r.communities)-3), end=len(r.communities)+1, step=1), \
                            ensemble.Parameter(name="threshold", start=0.25, end=0.75, step=0.05),\
                            ensemble.Parameter(name="seed", start=0, end=10, step=1)
                            ]  for r in references ]
    )
    # right implementaion of amplified commute distance
    results_scoring['lpam_python_amp'] = getResults(
        clustering_method=lpam_python_amp,
        graphs=graphs,
        names=names,
        references=references,
        parameters_list = [ [ensemble.Parameter(name="k", start=max(1,len(r.communities)-3), end=len(r.communities)+1, step=1), \
                            ensemble.Parameter(name="threshold", start=0.25, end=0.75, step=0.05),\
                            ensemble.Parameter(name="seed", start=0, end=10, step=1)
                            ]  for r in references ]
    )



    # The PercoMVC approach composes of two steps. In the first step, the algorithm attempts to determine
    # all communities that the clique percolation algorithm may find. In the second step,
    # the algorithm computes the Eigenvector Centrality method on the output of the first step to measure
    # the influence of network nodes and reduce the rate of the unclassified nodes
    # ok 100%
    results_scoring['percomvc'] = getResultsParallel(
        clustering_method=algorithms.percomvc,
        graphs=graphs,
        names=names,
        references=references,
        parameters_list = [[]]*len(graphs))

    # Parameters:
    # g_original – a networkx/igraph object
    # dimensions – Number of embedding dimensions. Default 8.
    # iterations – Number of training iterations. Default 50.
    # learning_rate – Gradient ascent learning rate. Default is 0.005.
    # ok 100%
    results_scoring['danmf'] = getResultsParallel(
        clustering_method=algorithms.danmf,
        graphs=graphs,
        names=names,
        references=references,
        parameters_list =  [[ensemble.Parameter(name="seed", start=1, end=13, step=1)]]*len(graphs))

    # forDf = {'graph':[], 'method':[], 'score':[]}
    # for method in results_scoring:
    #     for graph, score in results_scoring[method]:
    #         forDf['graph'].append(graph)
    #         forDf['method'].append(method)
    #         forDf['score'].append(score)

    # df = pd.DataFrame(forDf)
    # df.to_csv('cdlib_results.csv')
    # sys.exit(0)

    # SLPA is an overlapping community discovery that extends tha LPA. SLPA consists of the following three stages: 1) the initialization 2) the evolution 3) the post-processing

    # Parameters:
    # g_original – a networkx/igraph object
    # t – maximum number of iterations, default 20
    # r – threshold ∈ [0, 1]. It is used in the post-processing stage: if the probability of seeing a particular label during the whole process is less than r, this label is deleted from a node’s memory. Default 0.1
    #ok 100%
    results_scoring['slpa'] = getResultsParallel(
        clustering_method=algorithms.slpa,
        graphs=graphs,
        names=names,
        references=references,
        parameters_list =  [[ensemble.Parameter(name="r", start=0.1, end=1.0, step=0.1)]]*len(graphs))



    results_scoring['egonet_splitter'] = getResultsParallel(
        clustering_method=algorithms.egonet_splitter,
        graphs=graphs,
        names=names,
        references=references,
        parameters_list = [[ensemble.Parameter(name="resolution", start=0.5, end=1.0, step=0.1)]]*len(graphs))



    # Demon is a node-centric bottom-up overlapping community discovery algorithm.
    # It leverages ego-network structures and overlapping label propagation to identify micro-scale
    # communities that are subsequently merged in mesoscale ones.
    #
    # Parameters:
    # - g_original – a networkx/igraph object
    # - epsilon – merging threshold in [0,1], default 0.25.
    # - min_com_size – minimum community size, default 3.
    # ok 100%
    results_scoring['demon'] = getResultsParallel(
        clustering_method=algorithms.demon,
        graphs=graphs,
        names=names,
        references=references,
        parameters_list =  [[ensemble.Parameter(name="epsilon", start=0.0, end=1.0, step=0.1)]]*len(graphs))


    # Find k-clique communities in graph using the percolation method. A k-clique community is the union of all cliques of size k that can be reached through adjacent (sharing k-1 nodes) k-cliques.

    # Parameters:
    # g_original – a networkx/igraph object
    # k – Size of smallest clique
    # ok 100%
    results_scoring['kclique'] = getResultsParallel(
        clustering_method=algorithms.kclique,
        graphs=graphs,
        names=names,
        references=references,
        parameters_list = [[ensemble.Parameter(name="k", start=3, end=5, step=1)]]*len(graphs))


    # LAIS2 is an overlapping community discovery algorithm based on the density function.
    # In the algorithm considers the density of a group is defined as the average density of the communication
    # exchanges between the actors of the group. LAIS2 IS composed of two procedures LA (Link Aggregate Algorithm)
    # and IS2 (Iterative Scan Algorithm).
    # ok 100%
    results_scoring['lais2'] = getResultsParallel(
        clustering_method=algorithms.lais2,
        graphs=graphs,
        names=names,
        references=references,
        parameters_list = [[]]*len(graphs))

    # Parameters:
    # g_original – a networkx/igraph object
    # seeds – Node list
    # min_com_size – the minimum size of a single community in the network, default 20
    # max_com_size – the maximum size of a single community in the network, default 50
    # expand_step – the step of seed set increasement during expansion process, default 6
    # subspace_dim – dimension of the subspace; choosing a large dimension is undesirable because it would increase the computation cost of generating local spectra default 3
    # walk_steps – the number of step for the random walk, default 3
    # biased – boolean; set if the random walk starting from seed nodes, default False
    # all below are ok
    results_scoring['angel'] = getResultsParallel(
        clustering_method=algorithms.angel,
        graphs=graphs,
        names=names,
        references=references,
        parameters_list =  [[ensemble.Parameter(name="threshold", start=0.3, end=0.7, step=0.1)]]*len(graphs))

    results_scoring['leiden'] = getResultsParallel(
        clustering_method=algorithms.leiden,
        graphs=graphs,
        names=names,
        references=references,
        parameters_list = [[]]*len(graphs))

    results_scoring['label_propagation'] = getResultsParallel(
        clustering_method=algorithms.label_propagation,
        graphs=graphs,
        names=names,
        references=references,
        parameters_list = [[]]*len(graphs))

    results_scoring['oslom'] = getResults(
        clustering_method=oslom,
        graphs=graphs,
        names=names,
        references=references,
        parameters_list = [ [ensemble.Parameter(name="cp", start=0.25, end=0.6, step=0.05), \
                            ensemble.Parameter(name="t", start=0.25, end=0.6, step=0.05),\
                            ensemble.Parameter(name="seed", start=1, end=11, step=1)
                            ]  for r in references ]
    )


    results_scoring['copra'] = getResults(
        clustering_method=copra,
        graphs=graphs,
        names=names,
        references=references,
        parameters_list = [ [ensemble.Parameter(name="param", start=1, end=10, step=1), \
                            ]  for r in references ]
    )


    results_scoring['GCE'] = getResults(
        clustering_method=GCE,
        graphs=graphs,
        names=names,
        references=references,
        parameters_list = [ [ensemble.Parameter(name="param", start=0.05, end=1.0, step=0.05), \
                            ]  for r in references ]
    )

    forDf = {'graph':[], 'method':[], 'score':[]}
    for method in results_scoring:
        for graph, score in results_scoring[method]:
            forDf['graph'].append(graph)
            forDf['method'].append(method)
            forDf['score'].append(score)

    df = pd.DataFrame(forDf)
    df.to_csv('cdlib_results_legacy_part.csv')

    graphStatistics = {'graph': [], 'statistic': [], 'value':[]}
    statistics = {
        '$|V|$': lambda graph, reference: len(graph.nodes),
        '$|E|$': lambda graph, reference: len(graph.edges),
        '$\\hat{k}$': lambda graph, reference: len(reference.communities),
        '$\\hat{o}$': lambda graph, reference: getOverlappingNumber(reference.communities),
        '$\\hat{c}$': lambda graph, reference: nx.average_clustering(graph),
        '\\shortstack[l]{normalized \\\\ cut}': lambda graph, reference: evaluation.normalized_cut(graph,reference, summary=True).score,
        '\\shortstack[l]{internal \\\\ edge \\\\density}': lambda graph, reference: evaluation.internal_edge_density(graph,reference, summary=True).score,
        'significance': lambda graph, reference: evaluation.significance(graph,reference, summary=True).score,
        '\\shortstack[l]{avarage \\\\ internal \\\\degree}': lambda graph, reference: evaluation.average_internal_degree(graph,reference, summary=True).score,
        '\\shortstack[l]{modularity \\\\ density}': lambda graph, reference: evaluation.modularity_density(graph,reference).score,
    }

    for graph, name, reference in zip(graphs, names, references):
        for statistic_name, statistic in statistics.items():
            graphStatistics['graph'].append(name)
            graphStatistics['statistic'].append(statistic_name)
            graphStatistics['value'].append(statistic(graph, reference) )

    graphStatisticsDF = pd.DataFrame(graphStatistics)
    graphStatisticsDF.to_csv('../Results/graphs_stats_legacy.csv')

    # method that were not used, or don't work by some reason

    # wrong implementaion of amplified commute distance
    # results_scoring['lpam_amp'] = getResultsParallel(
    #     clustering_method=lpam_python,
    #     graphs=graphs,
    #     names=names,
    #     references=references,
    #     parameters_list = [ [ensemble.Parameter(name="k", start=max(1,len(r.communities)-3), end=len(r.communities)+1, step=1), \
    #                         ensemble.Parameter(name="threshold", start=0.25, end=0.75, step=0.05),\
    #                         ensemble.Parameter(name="seed", start=0, end=10, step=1)
    #                         ]  for r in references ],
    #                         num_threads=1
    # )

      # ok?
    # CONGO (CONGA Optimized) is an optimization of the CONGA algortithm.
    # The CONGO algorithm is the same as CONGA but using local betweenness.
    # The complete CONGO algorithm is as follows:

    # Calculate edge betweenness of edges and split betweenness of vertices.
    # Find edge with maximum edge betweenness or vertex with maximum split betweenness, if greater.
    # Recalculate edge betweenness and split betweenness:
    # Subtract betweenness of h-region centred on the removed edge or split vertex.
    # Remove the edge or split the vertex.
    # Add betweenness for the same region.
    # Repeat from step 2 until no edges remain.
    # Parameters:
    # - g_original – a networkx/igraph object
    # - number_communities – the number of communities desired
    # - height – The lengh of the longest shortest paths that CONGO considers, default 2
    #
    # results_scoring['congo'] = getResultsParallel(
    #     clustering_method=algorithms.congo,
    #     scoring_method=adjusted_mutual_info_score,
    #     graphs=graphs,
    #     names=names,
    #     references=references,
    #     parameters_list = [[ensemble.Parameter(name="number_communities", start=max(1,len(r)-3), end=len(r)+1, step=1),
    #                     ensemble.Parameter(name="height", start=2, end=5, step=1)] for r in references] )



    # нужен node seeds
    # Lemon is a large scale overlapping community detection method based on local expansion via minimum one norm.

    # The algorithm adopts a local expansion method in order to identify the community members from a few exemplary seed members.
    # The algorithm finds the community by seeking a sparse vector in the span of the local spectra such that the seeds are in its support.
    # LEMON can achieve the highest detection accuracy among state-of-the-art proposals.
    # The running time depends on the size of the community rather than that of the entire graph.



    # slow
    # LFM is based on the local optimization of a fitness function. It finds both overlapping communities and the hierarchical structure.

    # Parameters:
    # g_original – a networkx/igraph object
    # alpha – parameter to controll the size of the communities: Large values of alpha yield very small communities,
    # small values instead deliver large modules. If alpha is small enough, all nodes end up in the same cluster,
    # the network itself. In most cases, for alpha < 0.5 there is only one community,
    # for alpha > 2 one recovers the smallest communities.
    # A natural choise is alpha =1.
    # results_scoring['lfm'] = getResults(
    #     clustering_method=algorithms.lfm,
    #     scoring_method=adjusted_mutual_info_score,
    #     graphs=graphs,
    #     names=names,
    #     references=references,
    #     parameters_list = [[ensemble.Parameter(name="alpha", start=0.7, end=1.2, step=0.1)]]*len(graphs))


    #OSSE is an overlapping community detection algorithm optimizing the conductance community score
    # The algorithm uses a seed set expansion approach; the key idea is to find good seeds,
    # and then expand these seed sets using the personalized PageRank clustering procedure.
    #
    # Parameters:
    # g_original – a networkx/igraph object
    # seeds – Node list
    # ninf – Neighbourhood Inflation parameter (boolean)
    # expansion – Seed expansion: ppr or vppr
    # stopping – Stopping criteria: cond
    # nworkers – Number of Workers: default 1
    # nruns – Number of runs: default 13
    # alpha – alpha value for Personalized PageRank expansion: default 0.99
    # maxexpand – Maximum expansion allowed for approximate ppr: default INF
    # delta – Minimum distance parameter for near duplicate communities: default 0.2
    # results_scoring['OSSE'] = getResults(
    #     clustering_method=algorithms.overlapping_seed_set_expansion,
    #     scoring_method=adjusted_mutual_info_score,
    #     graphs=graphs,
    #     names=names,
    #     references=references,
    #     parameters_list =  [[ensemble.Parameter(name="delta", start=0.0, end=0.2, step=0.05)]]*len(graphs))



    ##threshold, overlap_threshold, min_comm_size=3
    # FileNotFoundError: [Errno 2] No such file or directory: 'tmp/part1_'
    # results_scoring['node_perception'] = getResultsParallel(
    #     clustering_method=algorithms.node_perception,
    #     scoring_method=adjusted_mutual_info_score,
    #     graphs=graphs,
    #     names=names,
    #     references=references,
    #     parameters_list =  [[ensemble.Parameter(name="threshold", start=0.3, end=0.7, step=0.1),
    #                         ensemble.Parameter(name="overlap_threshold", start=0.3, end=0.7, step=0.1)]]*len(graphs))


  # Get an Error "The node indexing is wrong.
    # BigClam is an overlapping community detection method that scales to large networks. The procedure uses gradient ascent to create an embedding which is used for deciding the node-cluster affiliations.

    # Parameters:
    # g_original – a networkx/igraph object
    # dimensions – Number of embedding dimensions. Default 8.
    # iterations – Number of training iterations. Default 50.
    # learning_rate – Gradient ascent learning rate. Default is 0.005.
    # results_scoring['big_clam'] = getResults(
    #     clustering_method=algorithms.big_clam,
    #     scoring_method=adjusted_mutual_info_score,
    #     graphs=graphs,
    #     names=names,
    #     references=references,
    #     parameters_list =  [[ensemble.Parameter(name="dimensions", start=8, end=9, step=1)]]*len(graphs))



    # ASLPAw can be used for disjoint and overlapping community detection and works on weighted/unweighted and directed/undirected networks. ASLPAw is adaptive with virtually no configuration parameters.
    # results_scoring['aslpaw'] = getResultsParallel(
    #     clustering_method=algorithms.aslpaw,
    #     graphs=graphs,
    #     names=names,
    #     references=references,
    #     parameters_list = [[]]*len(graphs))



    # Get an Error "The node indexing is wrong.
    # The procedure uses non-negative matrix factorization in order to learn an unnormalized cluster membership distribution over nodes. The method can be used in an overlapping and non-overlapping way.

    # Parameters:
    # g_original – a networkx/igraph object
    # dimensions – Embedding layer size. Default is 32.
    # iterations – Number of training epochs. Default 10.
    # seed – Random seed for weight initializations. Default 42.
    # learning_rate – Gradient ascent learning rate. Default is 0.005.
    # results_scoring['nnsed'] = getResults(
    #     clustering_method=algorithms.nnsed,
    #     scoring_method=adjusted_mutual_info_score,
    #     graphs=graphs,
    #     names=names,
    #     references=references,
    #     parameters_list = [[ensemble.Parameter(name="seed", start=1, end=13, step=1)]]*len(graphs))



# uncomment to plot the graph
# df = pd.DataFrame(columns=["graph", "score", "method"], data=forDf)
# sns.set(rc={'figure.figsize':(15,10)}) #set the figure size
# ax = sns.lineplot(x="graph", y="score", hue="method", data=df, legend="brief")
# ax.legend(loc='best')
# for tick in ax.get_xticklabels():
#     tick.set_rotation(90)
#     plt.tight_layout()

# %%
