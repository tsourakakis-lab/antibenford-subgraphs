import networkx as nx

def to_motivo_edgelist(filename, G):
    map_d = {}
    index = 0
    for i in G.nodes:
        map_d[i] = index
        index+= 1
    num_labeled_G = nx.relabel_nodes(G, map_d)
    di_G = nx.to_directed(num_labeled_G)
    cs = str(len(G.nodes))+' '+str(len(G.edges))
    fh = open(filename + ".txt", "w")
    fh.write(cs+'\n')
    fh.close()
    fh = open(filename + ".txt", "ab")
    nx.write_edgelist(di_G,fh, data=False)
    fh.close()
    return map_d
    
def get_start_digit(v):
    if v==0:
        return 0
    if v<0:
    	v = -v
    while v<1:
        v = v*10
    return int(str(v)[:1])

def attribute_to_weight(G, target_digit=1, attribute_name='weight'):
    pos_list = []
    neg_list = []
    for e in G.edges:
        start_digit = get_start_digit(G[e[0]][e[1]][attribute_name])
        if start_digit==0:
            continue
        if start_digit==target_digit:
            pos_list.append(e)
        else:
            neg_list.append(e)
    return pos_list, neg_list
    
def get_key(my_dict, vals):
    result = []
    for key, value in my_dict.items():
        if value in vals:
            result.append(key)
    if len(result) != len(vals):
        print('some vals not in dict')
        return []
    return result
