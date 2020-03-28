import numpy as np

from ete3 import PhyloTree, TreeStyle

from skbio import DistanceMatrix
from skbio.tree import nj




#########################################################################
class Node:
        def __init__(self, id, parent_dist, children = None):
                self.id = id
                self.parent_dist = parent_dist
                self.children = children
                self.parent = None

nodes = []

newick_str_gbl = ''
def print_tree(tree_nn):
        global newick_str_gbl
        if tree_nn.children == None: # leaf
                newick_str_gbl += tree_nn.id 

        else:
                newick_str_gbl += "("                
                if tree_nn.children != None:
                        for i in range(len(tree_nn.children)):
                                print_tree(tree_nn.children[i])
                                if i < len(tree_nn.children) - 1:
                                        newick_str_gbl += ","

                newick_str_gbl += ")"

        if tree_nn.parent_dist != None:
                newick_str_gbl += ":" + str(tree_nn.parent_dist)


def compute_q(dist):
        N = dist.shape[0]
        u = dist.sum(1)
        Q = np.zeros(dist.shape)

        for i in range(N):
                for j in range(N):
                        Q[i,j] = (N-2)*dist[i,j] - u[i] - u[j]
                        if i == j:
                               Q[i,j] = 0 

        return u, Q


def find_node(nodes, id):
        for node in nodes:
                if node.id == id:
                        return node
        return None


D =     [[0, 8, 4, 6],
        [8, 0, 8, 8],
        [4, 8, 0, 6],         
        [6, 8, 6, 0]]

D = [[0,  5,  9,  9,  8],
         [5,  0, 10, 10,  9],
         [9, 10,  0,  8,  7],
         [9, 10,  8,  0,  3],
         [8,  9,  7,  3,  0]]

D =     [[0, 8, 4, 6],
        [8, 0, 8, 8],
        [4, 8, 0, 6],         
        [6, 8, 6, 0]]

ides = ['A', 'B', 'C', 'D']

D = np.matrix((D))

for ide in ides:
        nodes.append(Node(ide, None))



N = D.shape[0]
iteration = 1
while N > 2:
        D_old = np.copy(D)
        rows, cols = D.shape
        print(D)
        print(ides)
        u, Q = compute_q(D)
        print(u)
        print(Q)
        
        # get minimun value        
        min_indexes = np.where(Q == np.amin(Q))
        cordinates = list(zip(min_indexes[0], min_indexes[1]))
        print(cordinates)
        #print(np.matrix(min_indexes))

        # choose the mininmun value and compute distances to new node
        i_min = cordinates[0][0] # choose the first one
        j_min = cordinates[0][1]
        d_1 = int(D[i_min, j_min]/2 + (u[i_min] - u[j_min])/(2*(N-2)))
        d_2 = int(D[i_min, j_min]/2 + (u[j_min] - u[i_min])/(2*(N-2)))
        print(i_min, d_1)
        print(j_min, d_2)

        ######################################################################
        # create nodes
        node_1 = find_node(nodes, ides[i_min])
        node_1.parent_dist = d_1
        node_2 = find_node(nodes, ides[j_min])
        node_2.parent_dist = d_2

        #node_1 = Node(ides[i_min], d_1)
        #node_2 = Node(ides[j_min], d_2)
        node_3 = Node('U' + str(iteration), None, [node_1, node_2]) # new node U1
        nodes.append(node_3)
        node_1.parent = node_3
        node_2.parent = node_3
        ides.pop(max(i_min, j_min))
        ides.pop(min(i_min, j_min))
        ides.insert(0, 'U' + str(iteration))
        ######################################################################3

        # delete rows and colums from D in order to insert the new node  U1    
        tmp = np.delete(D, min(i_min, j_min), 0)        # first delete minor i
        tmp = np.delete(tmp, max(i_min, j_min) -1, 0)   # then eliminate max j, -1 because just elimnate one row
        tmp = np.delete(tmp, min(i_min, j_min), 1)
        tmp = np.delete(tmp, max(i_min, j_min) -1, 1)
        
        # build new D
        D = np.zeros((rows -1, cols -1))
        D[1:rows-1, 1:cols-1] = tmp
        #print(D)

        tmp_index = 1
        for i in range(1, D_old.shape[0]): # the fist value is zero, so avoid it                
                if i != i_min and i != j_min:
                        print("here", i)
                        D[tmp_index, 0] =  D_old[i,i_min] -  d_1
                        D[0, tmp_index] =  D_old[i,i_min] -  d_1

                        tmp_index += 1

        #print(D)
        N = D.shape[0]
        
        iteration += 1


# insert the last node
print(ides)
print(D)
node_1 = find_node(nodes, ides[1])
node_2 = find_node(nodes, ides[0])
node_1.parent_dist = D[1,0]
node_2.children.append(node_1)
node_1.parent = node_2


for nod in nodes:
        if nod.children != None:
                str_tmp = ""
                for child in nod.children:
                        str_tmp = str_tmp + " " + child.id
                        
                print(nod.id, nod.parent_dist, "parent_dist:" ,nod.parent_dist, "parent:" ,nod.parent, "children:" ,str_tmp)
        else:
                print(nod.id, nod.parent_dist, "parent_dist:" ,nod.parent_dist, "parent:" ,nod.parent, "children: None")


# find the root:: the node without parent
tree_nn = None
for nod in nodes:
        if nod.parent == None:
                tree_nn = nod

print("\nTree\n")
newick_str_gbl = ''
print_tree(tree_nn)
print(newick_str_gbl)
t = PhyloTree(newick_str_gbl + ";")
t.show()

########################################################################################################3
########################################################################################################
# with skbio




data = [[0,  5,  9,  9,  8],
         [5,  0, 10, 10,  9],
         [9, 10,  0,  8,  7],
         [9, 10,  8,  0,  3],
         [8,  9,  7,  3,  0]]

data = [[0, 8, 4, 6],
        [8, 0, 8, 8],
        [4, 8, 0, 6],         
        [6, 8, 6, 0]]

ids = list('abcd')
dm = DistanceMatrix(data, ids)
tree = nj(dm)
print(tree.ascii_art())
newick_str = nj(dm, result_constructor=str)
print(newick_str)
#print(newick_str[:55], "...")
t = PhyloTree(newick_str)
t.show()