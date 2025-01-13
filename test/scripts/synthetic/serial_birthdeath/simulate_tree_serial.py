"""
Modified from https://github.com/evolbioinfo/phylodeep/tree/main/simulators/bd_models

"""

import sys, json
import numpy as np

STOP_UNKNOWN = 0
STOP_DIVERGENCE = 1
STOP_MIGRATION = 2
STOP_REMOVAL_WOS = 3
STOP_SAMPLING = 4
STOP_TIME = 5

class SBTree(object):

    def __init__(self, branch_length):
        self.branch_length = branch_length;
        self.up = None
        self.name = ""
        self.type = ""
        self.children = []

    def __str__(self):
        if self.is_leaf():
            return f"{self.name}[&type={self.type}]:{self.branch_length}"

        newickStr = "("
        if len(self.get_children()) == 1:
            newickStr += str(self.get_children()[0])
        else:
            newickStr += ",".join(str(child) for child in self.get_children())

        newickStr += f"){self.name}[&type={self.type}]:{self.branch_length}"
        return newickStr
    
    def is_root(self):
        return self.up == None
    
    def is_leaf(self):
        return self.children == []
    
    def add_child(self, branch_length):
        nd = SBTree(branch_length)
        nd.up = self
        self.children.append(nd)
        return nd

    def attach_child(self, child):
        self.children.append(child)
    
    def get_children(self):
        return self.children
    
    def remove_child(self, nd):
        nd.up == None
        self.children.remove(nd) 
    
    def traverse(self, order):
        if order == 'postorder':
            ret = []
            for child in self.children:
                ret.append(child.traverse('postorder'))
            return [nd for lst in ret for nd in lst] + [self]
        elif order == 'levelorder':
            ret = []
            layer = 0
            curr = [self]
            while (curr and layer < 1000):
                ret.append(curr)
                nxt = []
                for nd in curr:
                    nxt += nd.get_children()
                curr = nxt
                layer += 1

            if (layer == 1000):
                raise Exception("Error: cycle found or tree depth >= 1000")
            return [nd for lst in ret for nd in lst]
    
    def get_leaves(self):
        ret = []

        if self.is_leaf():
            return [self]
        for child in self.children:
            ret += child.get_leaves()

        return ret
    
    def get_size(self):
        return len(self.traverse('postorder'))


def simulate_bd_tree(divergence_r, migration_r, K, removal_r, sampling_p, size, 
                     max_t=float('inf'), trial_n=100):
    """
    Simulates the tree evolution with infectious hosts based on the given 
    divergence rate, removal rate, sampling probabilities and number of tips
    :param divergence_r: float of divergence rate
    :param removal_r: float of removal rate
    :param migration_r: float of active-dormant migration rate
    :param K: relative Seedbank size, ratio between active and dormant population sizes
    :param sampling_p: float, between 0 and 1, probability for removed leave to be immediately sampled
    :param size: number of sampled leaves in a tree
    :param max_t: maximum time from root simulation
    :param trial_n: number of attempts to get desired tree
    :return: the simulated tree (ete3.Tree).
    
    :param 
    """

    def update_rates(rates, metrics_dc):
        """
        updates rates dictionary
        :param rates: dict, all rate values at previous step
        :param metrics_dc: dict, counts of different individuals, list(s) of different branch types
        :return:
        """
        rates['divergence_rate_i'] = divergence_r * metrics_dc['number_active_lineages']
        rates['ad_migration_rate_i'] = migration_r * K * metrics_dc['number_active_lineages']
        rates['da_migration_rate_i'] = migration_r * metrics_dc['number_dormant_lineages']
        rates['removal_rate_i'] = removal_r * metrics_dc['number_total_lineages']
        
        rates['sum_rates_i'] = rates['divergence_rate_i'] + rates['ad_migration_rate_i'] 
        rates['sum_rates_i'] += rates['da_migration_rate_i'] + rates['removal_rate_i']
        return None

    def divergence_event(which_lf):
        """
        updates the tree, the metrics and leaves_dict following a tranmission event that affects which_lf
        :param which_lf: ete3.Tree, a leaf affected by the divergence event
        :return: void, modifies which_lf, leaves_dict and metrics
        """
        metrics['total_branches'] += 1
        metrics['number_total_lineages'] += 1
        metrics['number_active_lineages'] += 1
        
        # the leaf becomes a divergence node
        which_lf.stop_reason = STOP_DIVERGENCE

        # the leaf branches into two new branches
        left, right = which_lf.add_child(branch_length=0), which_lf.add_child(branch_length=0)
        left.dist_to_start = which_lf.dist_to_start
        right.dist_to_start = which_lf.dist_to_start
        left.type = 'active'
        right.type = 'active'

        # two new infectious leaves
        lineages_dict['active'].append(left)
        lineages_dict['active'].append(right)
        return None
        
    def migration_event(which_lf):
        """
        updates the tree, the metrics and leaves_dict following a migration event that affects which_lf
        :param which_lf: ete3.Tree, a leaf affected by the divergence event
        :return: void, modifies which_lf, leaves_dict and metrics
        """
        if which_lf.type == 'active':
            metrics['number_active_lineages'] -= 1
            metrics['number_dormant_lineages'] += 1
        else:
            metrics['number_dormant_lineages'] -= 1
            metrics['number_active_lineages'] += 1
        
        which_lf.stop_reason = STOP_MIGRATION

        left = which_lf.add_child(branch_length=0)
        left.dist_to_start = which_lf.dist_to_start

        if which_lf.type == 'active':
            left.type = 'dormant'
            lineages_dict['dormant'].append(left)
        else:
            left.type = 'active'
            lineages_dict['active'].append(left)

        return None
        
    def removal_event(which_lf):
        """
        updates the tree, the metrics and leaves_dict following a removal event that affects which_lf
        :param which_lf: ete3.Tree, a leaf affected by the divergence event
        :return: void, modifies which_lf, leaves_dict and metrics
        """
        metrics['number_total_lineages'] -= 1
        if which_lf.type == "active":
            metrics['number_active_lineages'] -= 1
        else:
            metrics['number_dormant_lineages'] -= 1
            
        # sampling upon removal?
        if np.random.rand() < sampling_p:
            metrics['number_sampled'] += 1
            which_lf.stop_reason = STOP_SAMPLING
            metrics['total_removed'] += 1
        # removal without sampling
        else:
            which_lf.stop_reason = STOP_REMOVAL_WOS
            metrics['total_removed'] += 1
        return None

    def which_lf_helper(nb_which_lf, active, time):
        """
        helper function
        """
        if active:
            which_lf = lineages_dict['active'][nb_which_lf]
            del lineages_dict['active'][nb_which_lf]
        else:
            which_lf = lineages_dict['dormant'][nb_which_lf]
            del lineages_dict['dormant'][nb_which_lf]

        which_lf.branch_length = abs(time - which_lf.dist_to_start)
        which_lf.dist_to_start = time
        return which_lf


    # up to 100 times retrial of simulation until reaching correct size
    right_size = False
    trial = 0
    while not right_size and trial < 7:
        # start a tree
        root = SBTree(branch_length = 0)
        root.dist_to_start = 0
        root.type = "active"

        # INITIATE: metrics counting leaves and branches of different types, leaves_dict storing all leaves alive, and
        # rates_i with all rates at given time, for Gillespie algorithm
        metrics = {'total_branches': 1, 'total_removed': 0, 
                   'number_total_lineages': 1, 'number_active_lineages': 1, 
                   'number_dormant_lineages': 0,'number_sampled': 0}
        lineages_dict = {'active': [root], 'dormant': []}

        rates_i = {'removal_rate_i': 0, 'ad_migration_rate_i': 0, 'da_migration_rate_i': 0, 
                   'divergence_rate_i': 0, 'sum_rates_i': 0}
        time = 0

        # simulate while [1] the epidemics do not go extinct, [2] given number of patients were not sampled,
        # [3] maximum time of simulation was not reached
        while metrics['number_total_lineages'] > 0 and metrics['number_sampled'] < size and time < max_t:
            # first we need to re-calculate the rates and take its sum
            update_rates(rates_i, metrics)
            # when does next event take place?
            time_to_next = np.random.exponential(1 / rates_i['sum_rates_i'], 1)[0]
            time = time + time_to_next
            # what event?
            random_event = (np.random.uniform(0, 1, 1) * rates_i['sum_rates_i'])[0]

            # print("++++")
            # print(rates_i)
            # print(random_event)
            
            # divergence event
            if random_event < rates_i['divergence_rate_i']:
                nb_which_lf = int(np.floor(np.random.uniform(0, metrics['number_active_lineages'], 1)[0]))
                which_lf = which_lf_helper(nb_which_lf, True, time)
                divergence_event(which_lf)
            elif random_event < rates_i['ad_migration_rate_i'] + rates_i['divergence_rate_i']:
                nb_which_lf = int(np.floor(np.random.uniform(0, metrics['number_active_lineages'], 1)[0]))
                which_lf = which_lf_helper(nb_which_lf, True, time)
                migration_event(which_lf)
            elif random_event < rates_i['da_migration_rate_i'] + rates_i['ad_migration_rate_i'] + rates_i['divergence_rate_i']:
                nb_which_lf = int(np.floor(np.random.uniform(0, metrics['number_dormant_lineages'], 1)[0]))
                which_lf = which_lf_helper(nb_which_lf, False, time)
                migration_event(which_lf)
            # removal event
            else:
                nb_which_lf = int(np.floor(np.random.uniform(0, metrics['number_total_lineages'], 1)[0]))
                if nb_which_lf < metrics['number_active_lineages']:
                    which_lf = which_lf_helper(nb_which_lf, True, time)
                else:
                    nb_which_lf -= metrics['number_active_lineages']
                    which_lf = which_lf_helper(nb_which_lf, False, time)
                removal_event(which_lf)

            # print(root)
            # print("----")
        print(root.get_size(), len(root.get_leaves()))
        # at the end of simulation, tag the non-removed tips
        for leaflet in root.get_leaves():
            if hasattr(leaflet, "stop_reason") and leaflet.stop_reason in  [STOP_REMOVAL_WOS, STOP_SAMPLING]:
                continue
            leaflet.branch_length = abs(time - leaflet.dist_to_start)
            leaflet.dist_to_start = time
            leaflet.stop_reason = STOP_TIME

        # we sampled correct number of patients
        if metrics['number_sampled'] == size:
            right_size = True
        # we retry again a simulation
        else:
            trial += 1

    i=0
    for node in root.traverse("levelorder"):
        node.name = "n" + str(i)
        i += 1

    # statistics on the number of branches, removed tips, sampled tips, time of simulation and number of sim trials
    metrics['time'] = time
    metrics['trial'] = trial
    return root, metrics

def _merge_node_with_its_child(nd, child=None):
    if not child:
        child = nd.get_children()[0]
    # nd_hist = getattr(nd, HISTORY, [(getattr(nd, state_feature, ''), 0)])
    # nd_hist += [('!', nd.dist - sum(it[1] for it in nd_hist))] \
    #            + getattr(child, HISTORY, [(getattr(child, state_feature, ''), 0)])
    # child.add_features(**{HISTORY: nd_hist})
    child.branch_length += nd.branch_length
    if nd.is_root():
        child.up = None
    else:
        parent = nd.up
        parent.remove_child(nd)
        parent.attach_child(child)
    return child


def remove_certain_leaves(tree):
    """
    Removes all the branches leading to naive leaves from the given tree.
    :param tree: the tree of interest
    :return: the tree with naive branches removed
            or None if all the leaves were naive in the initial tree.
    """

    for nd in tree.traverse("postorder"):
        # If this node is a divergence has only one child branch
        # it means that the other child branch used to lead to a naive leaf and was removed.
        # We can merge this node with its child
        # (the child was already processed and either is a leaf or has 2 children).
        if nd.stop_reason == STOP_DIVERGENCE and len(nd.get_children()) == 1:
            merged_node = _merge_node_with_its_child(nd)
            if merged_node.is_root():
                tree = merged_node
        # migration nodes with 1 child are safe, with 0 child are leafs
        elif nd.is_leaf() and nd.stop_reason != STOP_SAMPLING:
            if nd.is_root():
                return None
            nd.up.remove_child(nd)
    return tree

def main():
    if len(sys.argv) != 2:
        print("Call", sys.argv[0], "<path to config file>")
        sys.exit(1)

    with open(sys.argv[1]) as f:
        config_data = json.load(f)

    sample_size = config_data["sample_size"]
    sample_size_1 = config_data["sample_size_1"]
    sample_size_2 = config_data["sample_size_2"]
    migration_rate = config_data["transition_rate"]
    K = config_data["relative_seedbank_size"]

    migration_rate = 0.5
    divergence_rate = 0.2
    removal_rate = 0.5
    sampling_prob = 1.0

    np.random.seed(2)

    root, vector_count = simulate_bd_tree(divergence_rate, migration_rate, K, removal_rate, sampling_prob, sample_size)
    print(vector_count)
    print("Simulated BD tree")
    print(root)

    # remove unsampled tips
    root = remove_certain_leaves(root)

    print("\nSimulated BD tree after removing unsampled tips")
    print(root)

    with open("simulate_tree_serial_output.txt", "w") as file:
        newick_str = str(root)
        file.write(newick_str + ";\n")

if __name__ == "__main__":
    main()
