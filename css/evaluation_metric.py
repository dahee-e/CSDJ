import networkx as nx
from itertools import chain
import time
import argparse
import common_utility
import os
def get_base(file_path) :
    path = os.path.dirname(file_path)
    return path+"/"

def get_average_diameter(G, Cgraph):
    sum = 0
    for cg in Cgraph :
        sum = sum + nx.diameter(cg)
    return sum/len(Cgraph)
def degreeSum(G, nodes):
    sum = 0
    for u in nodes:
        sum = sum + G.degree(u)
    return sum

def volume(G, S, weight=None):
    degree = G.out_degree if G.is_directed() else G.degree
    return sum(d for v, d in degree(S, weight=weight))

def get_average_local_modularity(G, C, Cgraph):

    sum = 0
    for idx, S in enumerate(C):
        T = set(G) - set(S)
        cgraph = Cgraph[idx]
        out_degree = cut_size(G, S, T, None)
        in_degree = cgraph.number_of_edges()
        if out_degree == 0 :
            sum += 0
        else :
            sum += in_degree/out_degree

    return sum/len(C)
def cut_size(G, S, T=None, weight=None):

    edges = nx.edge_boundary(G, S, T, data=weight, default=1)
    if G.is_directed():
        edges = chain(edges, nx.edge_boundary(G, T, S, data=weight, default=1))
    return sum(weight for u, v, weight in edges)

def get_avg_graph_density(Cgraph):
    sum = 0
    for cg in Cgraph:
        n = cg.number_of_nodes()
        e = cg.number_of_edges()
        sum += e/n
    return sum/len(Cgraph)
def get_avg_edge_density(Cgraph):
    sum = 0
    for cg in Cgraph:
        n = cg.number_of_nodes()
        e = cg.number_of_edges()
        #print('n', n, 'e', e)
        sum += 2*e/(n*(n-1))
    return sum/len(Cgraph)
def get_average_modularity(G, Cgraph):
    E = G.number_of_edges()
    sum = 0
    for cg in Cgraph :
        lc = cg.number_of_edges()
        var1 = lc/E
        var2 = (degreeSum(G, cg.nodes)/(2.0*E))**2
        mod = var1 - var2
        sum = sum + mod

    return sum/len(Cgraph)
def get_average_size(C):
    cnum = len(C)
    csum = 0
    for c in C:
        csum += len(c)
    return csum/cnum
def get_num_of_sub(C):
    return len(C)


#new evaluation metric
def get_avg_clustering_coefficient(C):
    return nx.average_clustering(C)
def get_conductance_local(G, C):
    sum = 0
    if len(C) == 1 :
        if len(C[0]) == len(G) :
            return 1
    for c in C :
        cd = nx.conductance(G, c)
        sum = sum + cd
    return sum / len(C)
def get_conductance_global(G, C,V):
    if G.number_of_nodes() == V:
        return 0
    return nx.conductance(G,C)
def get_average_degree(V,E): #average degree of vertices in C
    return (2*E)/V
def get_internal_density(V,E):
    return (2*E)/(V*(V-1))
def get_cut_ratio(G,C,V):
    boundary = nx.cut_size(G,C)
    print(boundary)
    n = G.number_of_nodes()
    if n == V:
        return 0
    return boundary/(V*(n-V))



#metric
def metric_global(G, C):
    if len(C) == 0 :
        return (None, None, None, None,None,0, 0)

    C = set(sum(C,[]))
    Cgraph = G.subgraph(C)


    # average_local_modularity = get_average_local_modularity(G, C, Cgraph)
    # average_diameter = get_average_diameter(G, Cgraph)
    # average_graph_density = get_avg_graph_density(Cgraph)
    # average_edge_density = get_avg_edge_density(Cgraph)

    V = Cgraph.number_of_nodes()
    E = Cgraph.number_of_edges()
    average_degree = get_average_degree(V,E)
    internal_density = get_internal_density(V,E)
    cut_ratio = 1 - get_cut_ratio(G,Cgraph,V)
    inverse_conductance = 1.0 - get_conductance_global(G, C, V)
    average_coefficient = get_avg_clustering_coefficient(Cgraph)




    return (average_degree,
            internal_density,
            cut_ratio,
            inverse_conductance,
            average_coefficient,
            V,
            E)


def metric_local(G, C):
    if len(C) == 0 :
        return (None, None, None, None, 0, 0)
    Cgraph = []

    for subset in C :
        G0 = G.subgraph(subset)
        Cgraph.append(G0)

    sum_modularity = get_average_modularity(G, Cgraph)
    average_graph_density = get_avg_graph_density(Cgraph)
    average_edge_density = get_avg_edge_density(Cgraph)
    inverse_conductance = 1.0 - get_conductance_local(G, C)
    average_size = get_average_size(C)
    num_of_sub = get_num_of_sub(C)

    return (sum_modularity,
            average_graph_density,
            average_edge_density,
            inverse_conductance,
            average_size,
            num_of_sub)

# mod, local_mod, v_density, e_density, inv_cond, diam, size


parser = argparse.ArgumentParser(description='measure')

parser.add_argument('--dataset', default="../dataset/karate/network.dat",
                    help='algorithm result file name (community C)')
parser.add_argument('--network', default="../dataset/karate/network.dat",
                    help='input graph (graph G)')

parser.add_argument('--local', type=bool, default=False,
                    help='local measure')

parser.add_argument('--global', type=bool, default=True,
                    help='global measure')


args = parser.parse_args()
print("network ", args.network)
print("dataset ", args.dataset)
user_params = args

output = get_base(args.network)
alg_name = args.dataset.split("/")[8].split(".dat")[0]
print(alg_name, args.dataset)
output = output  + alg_name
if args.local == True:
    output = output+".local"
else :
    output = output+".global"
print("output", output)

# read network
G = common_utility.readEdgeList(args.network)
G.remove_edges_from(nx.selfloop_edges(G))

f = open(args.dataset, 'r')
readC = f.readlines()[1:]
f.close()
C = list()
if readC != None:
    for lines in readC:
        line = lines.strip()
        C.append(line.split())


if args.local == True: #local measure

    mod, v_density, e_density, inv_cond,size, num_of_sub = metric_local(G,C)
    print("resultant_statistic ", mod, v_density, e_density, inv_cond,size, num_of_sub)
    with open(output, 'w') as f:
        f.write(alg_name+ '\n')
        f.write("modularity" + "\t" + str(mod) + '\n')
        f.write("vertex_density" + "\t" + str(v_density) + '\n')
        f.write("edge_density" + "\t" + str(e_density) + '\n')
        f.write("inv_cond" + "\t" + str(inv_cond) + '\n')
        f.write("size" + "\t" + str(size) + '\n')
        f.write("number_of_subgraph" + "\t" + str(num_of_sub) + '\n')


    f.close()
else : #global
    avg_degree, inter_density, cut_ratio, inv_cond, avg_coeff, V, E = metric_global(G, C)
    print("resultant_statistic ", avg_degree, inter_density, cut_ratio, inv_cond, avg_coeff, V,E)
    with open(output, 'w') as f:
        f.write(alg_name + '\n')
        # f.write("modularity" + "\t" + str(mod) + '\n')
        f.write("average_degree" + "\t" + str(avg_degree) + '\n')
        f.write("internal_density" + "\t" + str(inter_density) + '\n')
        f.write("cut_ratio" + "\t" + str(cut_ratio) + '\n')
        f.write("inv_cond" + "\t" + str(inv_cond) + '\n')
        f.write("clustering coefficient" + "\t" + str(avg_coeff) + '\n')
        f.write("num of node\t" + str(V) +"\n")

    f.close()
