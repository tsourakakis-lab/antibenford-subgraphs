import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from fibheap import FibonacciHeap as FibHeap
import random
import pandas as pd
from scipy.stats import chisquare

class Node:
    def __init__(self, node_id, motif_types):
        self.nodeId = node_id
        self.motif_list = [list() for i in range(motif_types)]
        

def init_heap_and_dict(motif_lists):
    '''
    motif_lists: list of tuples. Each tuple contains (motif_list, motif_value) for certain motif.
    '''
    node_d = dict()
    fibheap = FibHeap()
    motif_types = len(motif_lists)
    for motif_type_index in range(len(motif_lists)):
        motif_tuple = motif_lists[motif_type_index]
        motif_list = motif_tuple[0]
        motif_value = motif_tuple[1]
        for m_index in range(len(motif_list)):
            each_m = motif_list[m_index]
            for n in each_m:
                if n not in node_d:
                    node_d[n] = (Node(n, motif_types), fibheap.insert(motif_value, n))
                else:
                    fibheap.decrease_key(node_d[n][1], node_d[n][1].key+motif_value)
                node_d[n][0].motif_list[motif_type_index].append(m_index)

    return node_d, fibheap


def peel_by_motif(node_dict, fib_heap, motif_lists, least_nodes, log=True, log_interval=1):
    total_degree = sum([len(motif_tuple[0])*motif_tuple[1] for motif_tuple in motif_lists])
    n = node_dict.__len__()
    avg_degree = (float)(total_degree) / n
    node_list = list(node_dict.keys())
    max_avg = avg_degree
    S_size = n
    
    loss_list = [0 for i in range(len(motif_lists))]
    result_motif_nums = [len(motif_lists[i][0])-loss_list[i] for i in range(len(motif_lists))]

    for i in range(n - least_nodes):

        # find min node from graph (remove from heap)
        node_to_remove = fib_heap.extract_min().value
        # if log:
        # 	print(node_to_remove)
        for motif_type_index in range(len(node_dict[node_to_remove][0].motif_list)):
            for motif_index in node_dict[node_to_remove][0].motif_list[motif_type_index]:
                loss_list[motif_type_index] += 1
                total_degree -= motif_lists[motif_type_index][1]
                for each_node in motif_lists[motif_type_index][0][motif_index]:
                    if each_node!= node_to_remove:
                        fib_heap.decrease_key(node_dict[each_node][1], node_dict[each_node][1].key-motif_lists[motif_type_index][1])
                        node_dict[each_node][0].motif_list[motif_type_index].remove(motif_index)

        del node_dict[node_to_remove]
    
        avg_degree = (float)(total_degree) / (n - i - 1)
        
        if max_avg < avg_degree:
            max_avg = avg_degree
            node_list = list(node_dict.keys())
            S_size = n - i - 1
            result_motif_nums = [len(motif_lists[i][0])-loss_list[i] for i in range(len(motif_lists))]
        if log and i%log_interval==0:
            print('iteration', i, 'to remove node', node_to_remove)
            print(fib_heap.total_nodes, avg_degree, max_avg, [len(motif_lists[i][0])-loss_list[i] for i in range(len(motif_lists))])
    return max_avg, node_list, S_size, result_motif_nums

def count_occ(l, dist_len, adjust_idx=0):
    result = [0 for i in range(dist_len)]
    for e in l:
        result[e-adjust_idx] += 1
    return result


# calculate the kl divergence
def kl_divergence(p, q):
    return sum(p[i] * log2(p[i]/q[i]) for i in range(len(p)))

def node_chisquare(edgelist, node_num, dist, adjust_idx=1, directed='both', diff_func='chi'):
    '''
    edgelist: list of edges. Each item is a tuple contains (node_from, node_to, edge_weight) representing a weighted edge.
    node_num: number of nodes.
    dist: list of float sum to 1, describe the distribution used to calculate the chi square statistic. 
    '''
    node_induced_dist = [[] for i in range(node_num)]
    node_chis = []
    for edge in edgelist:
        if directed=='both' or directed=='out':
            node_induced_dist[edge[0]].append(edge[2])
        if directed=='both' or directed=='in':
            node_induced_dist[edge[1]].append(edge[2])
#             G[node][neighbor]['weight']
    for node_dist in node_induced_dist:
        count_dist = count_occ(node_dist, len(dist), adjust_idx)
#        get the chi square statistic, the higher it is, more abnormal the node is
        if diff_func=='chi':
            node_chis.append(chisquare(count_dist, sum(count_dist)*np.array(dist))[0])
        elif diff_func=='kl':
            node_chis.append(kl_divergence(np.array(count_dist)/sum(count_dist), np.array(dist)))
    return node_chis

# peeling algorithm similar to peel_by_motif, but take weighted edges as input.
def peel_by_edge(node_dict, fib_heap, edge_list, least_nodes, log=True, log_interval=1):
    total_degree = sum([edge[-1] for edge in edge_list])
    n = node_dict.__len__()
    avg_degree = (float)(total_degree) / n
    node_list = list(node_dict.keys())
    max_avg = avg_degree
    S_size = n

    for i in range(n - least_nodes):

        # find min node from graph (remove from heap)
        node_to_remove = fib_heap.extract_min().value

        for motif_type_index in range(len(node_dict[node_to_remove][0].motif_list)):
            for edge_index in node_dict[node_to_remove][0].motif_list[motif_type_index]:
                total_degree -= edge_list[edge_index][-1]
                for each_node in edge_list[edge_index][:-1]:
                    if each_node!= node_to_remove:
                        fib_heap.decrease_key(node_dict[each_node][1], node_dict[each_node][1].key-edge_list[edge_index][-1])
                        try:
                            node_dict[each_node][0].motif_list[motif_type_index].remove(edge_index)
                        except ValueError:
                            pass

        del node_dict[node_to_remove]
    
        avg_degree = (float)(total_degree) / (n - i - 1)
        
        if max_avg < avg_degree:
            max_avg = avg_degree
            node_list = list(node_dict.keys())
            S_size = n - i - 1
        if log and i%log_interval==0:
            print('iteration', i, 'to remove node', node_to_remove)
            print(fib_heap.total_nodes, avg_degree, max_avg)
    return max_avg, node_list, S_size

# init function for peel_by_edge.
def init_heap_dict_weighted(edgelist):
    node_d = dict()
    fibheap = FibHeap()
    motif_types = 1
    for m_index in range(len(edgelist)):
        edge = edgelist[m_index]
    #     G.add_edge(edge[0], edge[1], weight=edge[2])
        for n in edge[:-1]:
            if n not in node_d:
                node_d[n] = (Node(n, motif_types), fibheap.insert(edge[-1], n))
            else:
                fibheap.decrease_key(node_d[n][1], node_d[n][1].key+edge[-1])
            node_d[n][0].motif_list[0].append(m_index)

    return node_d, fibheap 

def several_peel_weighted(edges, node_chisquares, peel_times=3, log=True, weight_by='average'):
    pairs = []
    removed_nodes = []
    for edge in edges:
        if weight_by=='average':
            weight = np.average([node_chisquares[i] for i in edge])
        elif weight_by=='product':
            weight = np.prod([node_chisquares[i] for i in edge])
        elif weight_by=='sum':
            weight = sum([node_chisquares[i] for i in edge])
        temp_e = list(edge)
        temp_e.append(weight)
        pairs.append(temp_e)
    for i in range(peel_times):
    #     print('!!', len(pairs[0][0]))
        node_d, fibheap = init_heap_dict_weighted(pairs)
        results = peel_by_edge(node_d, fibheap, pairs, 1, log=log, log_interval=10000)
#         results = peel_by_motif(node_d, fibheap, pairs, 1, False)
#         H = G.subgraph(results[1])
        if log:
            print(results[1])
        removed_nodes.append(results[1])

        new_pairs = []
        for pair in pairs:
            if len(set(pair) & set(results[1]))>0:
                continue
            else:
                new_pairs.append(pair)
        #         print(len(temp_motif_list))
        pairs = new_pairs
    return removed_nodes

def multipass_greedy_peel(G, pos_list_generator, pos_value, neg_list_generator, neg_value):
    density_list = []
    H = G.copy()
#     rates = nx.get_edge_attributes(H, "rate")
#     rate_list.append(rates)
    while True:
        to_remove_list = []

        tris = nx.triangles(H)
        tri_count = sum([tris[i] for i in tris])
        c4s = squares(H)
        squ_count = sum([c4s[i] for i in c4s])
        print('triangle count:', squ_count)
        density = (squ_count/4*squ_value + tri_value*tri_count/3) / len(H.nodes)
        density_list.append(density)
        if density==max(density_list):
            max_graph = list(H.nodes)
        for n in H.nodes:
            if c4s[n]*squ_value + tris[n]*tri_value <= 2 * density:
                to_remove_list.append(n)

        print('to remove', len(to_remove_list), 'nodes, current density:', density, 'current # of nodes', len(H.nodes))

        H.remove_nodes_from(to_remove_list)
#         rates = nx.get_edge_attributes(H, "rate")
#         rate_list.append(rates)

        if (len(set(H.nodes))==0) or (len(to_remove_list)==0):
            print('over')
            break

    print(density_list)
    return max_graph
# edge_list = list(G.edges)
# edge_value = 1
# tri_value = -1
# node_d, fibheap = init_heap_and_dict([(edge_list,1), (cycls_3,-1)])
# total_degree = len(cycls_3)*tri_value + len(cycls_4)*squ_value
# results = peel_by_motif(node_d, fibheap, [(edge_list,1), (cycls_3,-1)], least_nodes=1)
