import networkx as nx
import numpy as np
import pandas as pd
import math

# alphaCore ranking of nodes in a complex, directed network
#
# Iteratively computes a node ranking based on a feature set derived from
# edge attributes and optionally static node features using the
# mahalanobis data depth function at the origin.
#
# @param graph A networkx directed graph
# @param stepSize Defines the stepsize of each iteration as percentage of node count
# @param startEpsi The epsilon to start with. Removes all nodes with depth>epsilon at start
# @param expoDecay Dynamically reduces the step size, to have high cores with few nodes if true
# @return A dataframe of columns nodeID, alpha value, and batchID
def alphaCore(G, stepSize, startEpsi, expoDecay = True):


    weighted = dict()
    for (u, v) in G.edges:
        weighted[(u, v)] = 1

    #1
    data = computeNodeFeatures(G,weighted=weighted)
    #2 compute cov matrix to be used for all remainder of depth calculations
    matrix = data.drop("node", axis=1)  # convert dataframe to numeric matrix by removing first column containing nodeID
    cov = np.cov(matrix.values.T.astype(float))
    #3 calculate the Mahalanobis depth and add it to the respective row of the dataframe
    data['mahal'] = calculateMahalFromCenter(data, 0, cov)
    #4
    epsi = startEpsi
    #5
    node = []
    alphaVals = []
    #6
    batch = []
    #7
    alpha = 1 - epsi
    #8
    alphaPrev = alpha
    #9
    batchID = 0
    #10
    while G.number_of_nodes() > 0:
        #11
        while True:
            depthFound = False  # to simulate do-while loop; used to check if there exists a node with depth >= epsi on current iteration
            #12
            for row in data.itertuples():
                if row.mahal >= epsi:
                    depthFound = True
                    #13
                    node.append(row.node)  # set node core
                    alphaVals.append(alphaPrev)
                    #14
                    batch.append(batchID)
                    #15
                    G.remove_node(row.node)
            #16
            batchID += 1
            #19 while condition of do-while loop of #11
            if G.number_of_nodes() == 0 or not depthFound:
                break
            #17
            data = computeNodeFeatures(G, weighted=weighted)  # recompute node properties
            #18
            data['mahal'] = calculateMahalFromCenter(data, 0, cov)  # recompute depth
        #20
        alphaPrev = alpha
        #21
        if expoDecay and G.number_of_nodes() > 0:  # exponential decay
            localStepSize = math.ceil(G.number_of_nodes() * stepSize)
            data = data.sort_values(ascending=False, by=['mahal'])
            epsi = data.iloc[localStepSize - 1]['mahal']
        else:  # step decay
            epsi -= stepSize
        #22
        alpha = 1 - epsi
    #23
    return pd.DataFrame({'node': node, 'alpha': alphaVals, 'batch': batch})



# Computes the node features of a given directed graph and returns a dataframe containing the features of each node
#
# @param graph A networkx directed graph
# @return A dataframe containing the computed node features with each row as a new entry and columns as different features
def computeNodeFeatures(G,weighted):
    '''
        :param G: graph
        :param weighted: dictionary type,,
        :return:strength : sum of weight
        '''
    df = pd.DataFrame(columns=['node', 'degree', 'strength'])
    cnt = 0

    # strength = defaultdict()
    for u in G.nodes():
        weight = 0
        for v in G[u]:
            try:
                weight += weighted[(u, v)]
            except KeyError:
                weight += weighted[(v, u)]
        df.loc[cnt] = [u, len(G[u]), weight]
        cnt += 1

        # strength[u] = weight
    return df



# Computes the mahalanobis depth from a given center of each row of a given set of data and returns it as an array
#
# @param data Dataframe where each row is a new entry and each column after the first (nodeID) is a type of data
# @param center A center value calculated with respect to when computing mahalanobis depth
# @param cov The covariance matrix of the data matrix
# @return An array containing the mahalanobis depth of each row entry of a given set of data
def calculateMahalFromCenter(data, center, cov):
    matrix = data.drop("node", axis=1)  # convert dataframe to numeric matrix by removing first column containing nodeID
    x_minus_center = matrix.values - center
    x_minus_center_transposed = (matrix.values - center).T
    inv_cov = np.linalg.pinv(cov) # inv-> pinv
    left = np.dot(x_minus_center, inv_cov)
    mahal = np.dot(left, x_minus_center_transposed)
    return np.diagonal(np.reciprocal(1+mahal))  # diagonal contains the depth vaulues corresponding to each row from matrix


def run(G, alpha, stepSize=0.1, startEpsilon=1) :
    G1 = G.copy()
    T = alphaCore(G1, stepSize, startEpsilon)
    #print(T)
    ret = []
    for index, each in T.iterrows() :
        if each['alpha'] >= alpha :
            ret.append(each['node'])

    G0 = G.subgraph(ret)
    return nx.connected_components(G0)