import numpy as np
import random
import numpy.linalg as linalg

class Node():
    ''' Initializes a node with given parameters.
    Arguments:
        active: state of node, active or dormant (refers to branch above)
        branch_length: length of branch that leads to this node (float)
        branch_id: id of branch that leads to this node (int)
    Other param
        probs: probability of observed bases beneath this node
                [list of 4 probs for 'ACGT'] (initialized to None)
        childs: list of children nodes
        parent: parent node
    '''
    def __init__(self, branch_id, active, branch_length):
        self.childs = []
        self.parent = None
        self.branch_length = branch_length
        self.active = active
        self.branch_id = branch_id
        self.probs = [None for _ in range(4)]


def likelihood(ordering, seq_list):
    '''
    Calculates the log-likelihood of a provided tree using the
    TideTree model.
    NOTE:  This is not yet used in this project
    :param ordering: list of nodes in post-traverse order
    :param seq_list: list of sequences on leaves
    :return: log-likelihood
    '''
    alph = 'ACGT'
    N = len(seq_list) #number of leaves
    n = len(seq_list[0])   #length of sequences
    #pi = np.array([1/2,1/6,1/4,1/12])
    pi = np.array([1/4,1/4,1/4,1/4])
    Lpi = np.log(pi)

    #arbitrary GTR
    # GTR_params_active = np.array([0.5,0.75,1.25,0.75,1.25,1.5]) #human random
    # GTR_params_dorm = np.array([0.25,0.2,0.9,0.6,0.7,0.35]) #human random

    #downscaled dormant mutation rate
    GTR_params_active = np.array([0.5,0.75,1.25,0.75,1.25,1.5]) #human random
    GTR_params_dorm = 0.5*np.array([0.5,0.75,1.25,0.75,1.25,1.5]) #downscaled mutation rate

    #arbitrary TN93
    # GTR_params_active = np.array([0.5,1.2,1.2,1.2,1.2,0.7]) #human random
    # GTR_params_dorm = np.array([0.2,0.5,0.5,0.5,0.5,0.8]) #human random


    #*** Code now normalizes and rescales mutation rate in GTR_model with 'scale' parameter, in line with Beast2 methods***
    #*** Need to re-rewrite this method implementation if one wishes to again account for arbitrary rate matrices***
    P_params_active = GTR_model(pi, GTR_params_active, 1)
    P_params_dorm = GTR_model(pi, GTR_params_active, 0.5)



    LL = 0
    for i in range(n):  #for each element i in sequence
        for j in range(len(ordering)):  #for each node in tree (first N nodes are leaves
            node = ordering[j]
            node.probs = [-np.inf for foo in range(4)]
            if j < N:
                node.probs[alph.index(seq_list[j][i])] = 1.0
            else:
                childs = node.childs
                if len(childs) == 1:
                    P_params = P_params_active if childs[0].active else P_params_dorm  #check branch type
                    for k in range(4):      #ancestor base
                        sum_mid = -np.inf
                        logP_mid = np.log(GTR_P(childs[0].branch_length,P_params))
                        for p in range(4):   #descendant base
                            sum_mid = sumLogProbs(sum_mid,childs[0].probs[p]+logP_mid[k,p])
                        node.probs[k] = sum_mid
                elif len(childs) == 2:
                    for k in range(4):      #ancestor base
                        sum_left = -np.inf
                        logP_left = np.log(GTR_P(childs[0].branch_length,P_params_active))
                        sum_right = -np.inf
                        logP_right = np.log(GTR_P(childs[1].branch_length,P_params_active))
                        for p in range(4):  #descendant base
                            sum_left = sumLogProbs(sum_left, childs[0].probs[p]+logP_left[k,p])
                            sum_right = sumLogProbs(sum_right, childs[1].probs[p]+logP_right[k,p])
                        node.probs[k] = sum_left+sum_right
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

def GTR_model(pi,params,scale=1):
    a,b,c,d,e,f = params
    Q = np.array([[0,a*pi[1],b*pi[2],c*pi[3]],[a*pi[0],0,d*pi[2],e*pi[3]],[b*pi[0],d*pi[1],0,f*pi[3]],[c*pi[0],e*pi[1],f*pi[2],0]])
    for i in range(4):
        Q[i,i] = -np.sum(Q[i])
        Q[i] = -scale*Q[i]/Q[i,i]      #rescaled to expected mutation rate of 1 (per unit time)
    Evals, T = np.linalg.eig(Q)
    Ti = linalg.inv(T)  #points off

    P_params = Evals, T, Ti
    return P_params


def GTR_P(t,P_params):
    Evals, T, Ti = P_params
    return T@np.diag(np.exp(Evals*t))@Ti


def create_test_tree_topology(type=1):
    '''
    Takes test sequences list
    Makes test tree, returns nodes in postorder
    '''
    def coalesce(a, b, len_a=1, len_b=1):
        '''
        :param nodes: list of nodes
        :param a: node1 index
        :param b: node2 index
        :return: list of nodes with a and b coalesced into parent node
        '''
        next_num = len(nodes)
        parent = Node(next_num, True, 1)
        parent.childs.append(nodes[a])
        parent.childs.append(nodes[b])
        nodes[a].parent = parent
        nodes[b].parent = parent
        nodes[a].branch_length = len_a
        nodes[b].branch_length = len_b
        nodes.append(parent)


    def state_transition(a, dormant_to_active, len1 = 1):
        next_num = len(nodes)
        parent = Node(next_num, dormant_to_active, 1)
        parent.childs.append(nodes[a])
        parent.branch_length = len1
        nodes[a].parent = parent
        nodes.append(parent)

    nodes = []
    for i in range(3):
        nodes.append(Node(i, True, 1))

    if type == 1:   #code in different sample trees later
        #dormant_active(1)
        state_transition(1,False,1)
        state_transition(3,True,1)
        coalesce(0, 4)
        state_transition(2,False,1)
        state_transition(6,True,1)
        coalesce(5, 7)

    if type == 2:   #code in different sample trees later
        coalesce(0, 1, 1, 2.5)
        coalesce(3, 2, 1, 2.5)

    return nodes


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

def test_seq_3x10():
    return ['AATCGGAGTT', 'ATTACCTCAT', 'CACTGTCATC']

def main():
    test_list = test_seq_3x10()
    ordering = create_test_tree_topology(1)
    LL = likelihood(ordering,test_list)
    print("LL for full calculation: ", LL)
    ordering = create_test_tree_topology(2)
    LL = likelihood(ordering,test_list)
    print("LL for scaled duration:  ",LL)

if __name__ == "__main__":
    main()