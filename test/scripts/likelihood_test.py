import numpy as np
import random
import numpy.linalg as linalg

class Node():
    
    def __init__(self, id, active, branch_length):
        ''' Initializes a node with given parameters.

            :param id: id of node (int)
            :param branch_length: length of branch that leads to this node (float)
            :param active: state of node and branch above, active or dormant (boolean)
                
            Other initialized fields:
            probs: probability of observed bases beneath this node
                [list of 4 probs for 'ACGT'] (initialized to None)
            children: list of children nodes
            parent: parent node
        '''
        self.id = id
        self.branch_length = branch_length
        self.active = active
        
        self.probs = [None for _ in range(4)]
        self.children = []
        self.parent = None


def likelihood(ordering, seq_list, model="JC69", dormant_scaling=1):
    '''
    Calculates the log-likelihood of a provided tree

    :param ordering: list of nodes in post-traverse order
    :param seq_list: list of sequences on leaves
    :return: log-likelihood
    '''
    if (model == "JC69"):
        pi = np.array([1/4,1/4,1/4,1/4])
        rates = np.array([1, 1, 1, 1, 1, 1])
    elif (model == "GTR"):
        pi = np.array([1/2,1/6,1/4,1/12])
        rates = np.array([0.5,0.75,1.25,0.75,1.25,1.5])
    else:
        pass

        #arbitrary GTR
        # GTR_params_active = np.array([0.5,0.75,1.25,0.75,1.25,1.5]) #human random
        # GTR_params_dorm = np.array([0.25,0.2,0.9,0.6,0.7,0.35]) #human random

        #arbitrary TN93
        # GTR_params_active = np.array([0.5,1.2,1.2,1.2,1.2,0.7]) #human random
        # GTR_params_dorm = np.array([0.2,0.5,0.5,0.5,0.5,0.8]) #human random

    #*** Code now normalizes and rescales mutation rate in GTR_model with 'scale' parameter, in line with Beast2 methods***
    #*** Need to re-rewrite this method implementation if one wishes to again account for arbitrary rate matrices***
    P_params_active = GTR_model(pi, rates)
    P_params_dorm = GTR_model(pi, rates, dormant_scaling)
    Lpi = np.log(pi)

    alph = 'ACGT'
    N = len(seq_list) #number of leaves
    n = len(seq_list[0])   #length of sequences

    LL = 0
    for i in range(n):  #for each element i in sequence
        for j in range(len(ordering)):  #for each node in tree (first N nodes are leaves

            node = ordering[j]
            # print(node.id)
            node.probs = [-np.inf for foo in range(4)]
            children = node.children

            if j < N:
                node.probs[alph.index(seq_list[j][i])] = 0 # = log(1)

            elif len(children) == 1:
                P_params = P_params_active if children[0].active else P_params_dorm  #check branch type
                logP_mid = np.log(GTR_P(children[0].branch_length,P_params))
                for k in range(4):      #ancestor base
                    sum_mid = -np.inf
                    for p in range(4):   #descendant base
                        sum_mid = sumLogProbs(sum_mid,children[0].probs[p]+logP_mid[k,p])
                    node.probs[k] = sum_mid

            elif len(children) == 2:
                logP_left = np.log(GTR_P(children[0].branch_length,P_params_active))
                logP_right = np.log(GTR_P(children[1].branch_length,P_params_active))
                # print(GTR_P(children[0].branch_length,P_params_active))
                # print(GTR_P(children[1].branch_length,P_params_active))

                for k in range(4):      #ancestor base
                    sum_left = -np.inf
                    sum_right = -np.inf
                    for p in range(4):  #descendant base
                        sum_left = sumLogProbs(sum_left, children[0].probs[p]+logP_left[k,p])
                        sum_right = sumLogProbs(sum_right, children[1].probs[p]+logP_right[k,p])
                    node.probs[k] = sum_left+sum_right

            else: 
                raise Exception(f"Node {node.id} has more than 2 children")
            
            # print(node.probs)

        LLi = -np.inf
        root = ordering[-1]
        for k in range(4):
            LLi = sumLogProbs(LLi,Lpi[k]+root.probs[k])
        LL += LLi

    return LL


def sumLogProbs(a,b):
    if a == -np.inf and b == -np.inf:
        return -np.inf
    elif a > b:
        return a + np.log1p(np.exp(b-a))
    else:
        return b + np.log1p(np.exp(a-b))


def GTR_model(pi, rates, scale=1, normalize=True):
    '''
    calculates the eigenvalues and eigenvectors of the rate matrix of a substitution model

    :param pi: base frequencies of the substitution model
    :param rates: mutation rates of the substitution model
    :param scale: scaling factor after the rate matrix is constructed
    :param normalize: whether to normalize the rate matrix so diagonal entries = -1
    :return: eigenvals, eigenvectors, inverse eigenvectors
    '''
    # print("DORM SCALING = " + str(scale))
    a,b,c,d,e,f = rates
    Q = np.array([[0,a*pi[1],b*pi[2],c*pi[3]],[ a*pi[0],0,d*pi[2],e*pi[3]],[b*pi[0],d*pi[1],0,f*pi[3]],[c*pi[0],e*pi[1],f*pi[2],0]])

    if normalize:
        denominators = np.array([sum(q) for q in Q])
        # Q = Q / denominators[:, None]
    
        beta = 1/(2 * (a*pi[0]*pi[1] + b*pi[0]*pi[2] + c*pi[0]*pi[3] + d*pi[1]*pi[2] + e*pi[1]*pi[3] + f*pi[2]*pi[3]))
        Q = Q * beta

    for i in range(4):
        Q[i,i] = -np.sum(Q[i])
    
    Q = Q * scale;
    # print(Q)

    Evals, T = np.linalg.eig(Q)
    Ti = linalg.inv(T)

    P_params = Evals, T, Ti
    return P_params


def GTR_P(t,P_params):
    '''
    :param t: branch length
    :param P_params: triplet of eigenvalues, eigenvectors, inverse eigenvectors of rate matrix
    :return: probability matrix
    '''
    Evals, T, Ti = P_params
    return T@np.diag(np.exp(Evals*t))@Ti


def create_test_tree_topology(type=1):
    '''
    Takes test sequences list
    Makes test tree, returns nodes in postorder
    '''

    def add_coalesce(a, b, length=1, _id=None):
        '''
        :param a: node1 index
        :param b: node2 index
        :param length: length of branch above new coalescent parent node
        '''
        if _id == None: _id = len(nodes)

        parent = Node(_id, True, length)
        node_map[_id] = parent

        parent.children.append(node_map[a])
        parent.children.append(node_map[b])
        node_map[a].parent = parent
        node_map[b].parent = parent
        nodes.append(parent)

    def add_migration(a, dormant_to_active, length=1, _id=None):
        '''
        :param a: node index
        :param dormant_to_active: type of migration
        :param length: length of branch above new migration node
        '''
        if _id == None: _id = len(nodes)

        parent = Node(_id, dormant_to_active, length)
        node_map[_id] = parent

        parent.children.append(node_map[a])
        node_map[a].parent = parent
        nodes.append(parent)

    nodes = []
    node_map = {}

    if type == 1: 
        for i in range(3):
            node_map[i] = Node(i, True, 1)
            nodes.append(node_map[i])

        add_migration(1, False, 1)
        add_migration(3, True, 1)
        add_coalesce(0, 4, 1)
        add_migration(2, False, 1)
        add_migration(6, True, 1)
        add_coalesce(5, 7, 0)

    if type == 2: 
        # same as type 1, with subsumed dormant branches and 0.5 scaling
        node_map[0] = Node(0, True, 1)    
        nodes.append(node_map[0])
        node_map[1] = Node(0, True, 2.5)  
        nodes.append(node_map[1])
        node_map[2] = Node(0, True, 2.5)  
        nodes.append(node_map[2])

        add_coalesce(0, 1, 1)
        add_coalesce(3, 2, 0)

    if type == 3:
        # dormant tips, different branch lengths
        node_map[1] = Node(1, False, 0.3)
        nodes.append(node_map[1])
        node_map[2] = Node(2, True, 0.4)
        nodes.append(node_map[2])
        node_map[3] = Node(3, False, 2.4)
        nodes.append(node_map[3])
        node_map[4] = Node(4, True, 5.0)
        nodes.append(node_map[4])
        
        add_migration(1, True, 0.7, 101)
        add_migration(2, False, 0.2, 201)
        add_migration(201, True, 0.4, 202)
        add_coalesce(101, 202, 0.3, 5)

        add_migration(5, False, 0.3, 501)
        add_migration(501, True, 0.9, 502)
        add_migration(3, True, 0.6, 301)
        add_coalesce(502, 301, 1.1, 6)

        add_coalesce(6, 4, 0.0, 7)
    
    if type == 4:
        for i in range(10):
            node_map[i] = Node(i, True, i+1)
            nodes.append(node_map[i])
        for i in range(10, 20):
            node_map[i] = Node(i, False, i)
            nodes.append(node_map[i])
        
        for i in range(10): #creates migrations 20-29 on 0-9
            add_migration(i, False, i+1)
        for i in range(20, 30): #creates migrations 30-39 on 20-29
            add_migration(i, True, i-15)
        last = 30
        for i in range(31, 40): #creates coalescents 40-48 on 30-39
            add_coalesce(last, i, i-30)
            last = len(nodes)-1
        
        for i in range(10, 20): #creates migrations 49-58 on 10-19
            add_migration(i, True, i)
        last=49 
        for i in range(50, 59): #creates coalescents 59-67 on 49-58
            add_coalesce(last, i, i-49)
            last = len(nodes)-1
        
        add_coalesce(48, 67, 0) #join last remaining lineages

    return nodes

def tree_to_str(root):
    out = f"[&state={str(1 if root.active else 0)}]:{root.branch_length}"

    if len(root.children) == 0:
        out = "sample" + str(root.id) + out
    elif len(root.children) == 1:
        out = f"({tree_to_str(root.children[0])})" + out
    elif len(root.children) == 2:
        out = f"({tree_to_str(root.children[0])},{tree_to_str(root.children[1])})" + out
    else:
        raise Exception("More than two children found")
    return out


"""takes a number of desired sequences and desired length 
and returns a list of sequences of desired length"""
def make_sample_sequences(num_sequences,seq_length):
    alph = "ACGT"
    seq_list = []
    random.seed(83)
    for i in range(num_sequences):
        seq_list.append('')
        for j in range(seq_length):
            seq_list[i] = seq_list[i] + alph[np.random.randint(0,4)]
    return seq_list

def test_seq_3x1():
    return ['A', 'C', 'G']

def test_seq_3x10():
    return ['AATCGGAGTT', 'ATTACCTCAT', 'CACTGTCATC']

def test_seq_4x20():
    return ['TAAGAGGGGAGAGCTCTCAC', 'CTCGGAAGAGCGAGTCCCCC', 'CCCCGGTAGGACGCACCCGG', 'CCCTTATTGCCCGGGATCCA']

def test_seq_20x20():
    return ['ATAACCATGAGTTCTCGACA', 'AGGTGTGTAGCGCACCGCTT', 'GGGTACATGCCAACTAGAGA', 'ACATGCTACGATGATGTCTC', 'TGGCCCCTCTTGTTGACAGT', 'GTACTGAAACTTGTGGGTCA', 'ACACGAATGCAAGCGTAACT', 'CCTTAGCAAGTAAATGGGTC', 'CAGAGGGCAAATAACTCCGG', 'GGTTGGGGCCTAGGCAAAAG', 'AGTTAGTGTATGCGAGATTG', 'GTCCCGAGCAGACATCCTAT', 'CTCTGCCGATTCTTCACCGG', 'AGTTGTCAACTCAAGCGTGT', 'TGTGCTTTCATAACCATTTT', 'ACACAATATTTTACATGTAC', 'TCAAACCTCCGGCATCGATT', 'TGGACTCTGCAATAGCTGCG', 'GAACCACGTTCTCTGGCAAG', 'ATTATCGGTCCTATGGCAAC']

def main():
    # test_list = test_seq_3x10()
    # # test_list = test_seq_3x1()
    # ordering = create_test_tree_topology(1)
    # LL = likelihood(ordering, test_list, "JC69", 0.5)
    # print("LL for full calculation: ", LL)

    # ordering = create_test_tree_topology(2)
    # LL = likelihood(ordering, test_list, "JC69", 1)
    # print("LL for scaled duration:  ",LL)

    # LL = likelihood(create_test_tree_topology(3), test_seq_4x20(), "JC69", 0.7)
    # print("LL for full calculation: ", LL)

    # LL = likelihood(create_test_tree_topology(1), test_seq_3x10(), "GTR", 0.5)
    # print("GTR explicit: ", LL)

    # LL = likelihood(create_test_tree_topology(2), test_seq_3x10(), "GTR")
    # print("GTR implicit: ", LL)

    # LL = likelihood(create_test_tree_topology(3), test_seq_4x20(), "GTR", 0.7)
    # print("LL for full calculation: ", LL)

    # x = create_test_tree_topology(4)
    # print(tree_to_str(x[-1]))

    LL = likelihood(create_test_tree_topology(4), test_seq_20x20(), "GTR", 0.09)
    print("LL for full calculation: ", LL)


if __name__ == "__main__":
    main()