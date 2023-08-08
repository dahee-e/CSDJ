
import networkx as nx


#############################################################################
#############################################################################

def check_integrity(args):
    assert type(args.k) is int, "k type error"
    assert type(args.h) is int, "h type error"
    assert type(args.s) is int, "s type error"


    if args.algorithm == 'kcore' or args.algorithm == 'kecc' or args.algorithm == 'kvcc' or args.algorithm=='kpeak' \
            or args.algorithm == 'kscore' or args.algorithm == 'kpcore' or args.algorithm == 'MkMFD':
        assert args.k >= 1

    if args.algorithm == 'ktruss':
        assert args.k >= 2

    if args.algorithm == 'alphacore' :
        args.alpha <= 1.0 and args.alpha >= 0.0

    return True


def readEdgeList(fileName):
    #print(fileName)
    g1 = nx.read_edgelist(fileName)
    print("V=", g1.number_of_nodes(), "\tE=", g1.number_of_edges())
    return g1

