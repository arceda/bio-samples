from ete3 import PhyloTree

from skbio import DistanceMatrix
from skbio.tree import nj



if __name__ == "__main__":
    print("Tree using only a distance matrix")
    data = [[0, 8, 4, 6],
            [8, 0, 8, 8],
            [4, 8, 0, 6],         
            [6, 8, 6, 0]]

    '''data = [[0,     0.4,    0.35,   0.6],
            [0.4,   0,      0.45,   0.7],
            [0.35,  0.45,   0,      0.55],         
            [0.6,   0.7,    0.55,   0]]        
        '''
    ids = list('ABCD')
    dm = DistanceMatrix(data, ids)
    tree = nj(dm)
    print(tree.ascii_art())
    newick_str = nj(dm, result_constructor=str)
    print(newick_str)
    #print(newick_str[:55], "...")
    t = PhyloTree(newick_str)
    t.show()

   