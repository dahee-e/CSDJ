from cdlib import algorithms


def run(G, k, epsilon):
    ret = algorithms.scan(G, epsilon, k).communities
    ret = [C for C in ret if len(C)>1]
    return ret
