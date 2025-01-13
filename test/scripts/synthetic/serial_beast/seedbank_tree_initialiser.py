###################################
# Seedbank coalescent simulation converted from 
# SeedbankTreeInitialiser simulate() function.
#
# python3 seedbank_tree_initialiser.py <config_file>
#
###################################

from typing import List
import numpy as np
import sys, json

class SeedbankNode:
    def __init__(self):
        self.nr = None
        self.ID = None
        self.height = None
        self.node_type = 0 # 0: dormant, 1: active
        self.n_type_changes = 0
        self.change_types = []
        self.change_times = []
        self.parent = None
        self.left = None
        self.right = None
    
    def add_change(self, new_type: int, time: float):
        self.change_types.append(new_type)
        self.change_times.append(time)
        self.n_type_changes += 1

class SBEvent:
    def __init__(self):
        self.time = None
        self.from_type = None
        self.to_type = None

class NullEvent(SBEvent):
    def __init__(self):
        super().__init__()
        self.time = float('inf')

class CoalescenceEvent(SBEvent):
    def __init__(self, type, time):
        self.from_type = type
        self.to_type = type
        self.time = time

class MigrationEvent(SBEvent):
    def __init__(self, from_type, to_type, time):
        self.from_type = from_type
        self.to_type = to_type
        self.time = time

def simulate_tree(n_leaves: int, leaf_names: List[str], leaf_types: List[str], leaf_times: List[float], params) -> SeedbankNode:
    """
    Generates tree using the specified list of active leaf nodes 
    using the seedbank coalescent.

    Returns:
        Root node of generated tree.
    """

    # Initialize node creation counter
    next_node_nr = 0

    # Initialize node
    live_nodes = [[] for _ in range(2)]
    dead_nodes = [[] for _ in range(2)]

    # Add nodes to dead nodes list
    for l in range(n_leaves):
        node = SeedbankNode()
        node.nr = next_node_nr
        node.ID = leaf_names[l]
        dead_nodes[leaf_types[l]].append(node)
        node.height = leaf_times[l]
        node.node_type = leaf_types[l]
        next_node_nr += 1

    # Sort nodes in dead nodes lists in order of increasing age
    for i in range(2):
        dead_nodes[i].sort(key = lambda node: node.height)

    # Allocate prepensity lists
    migration_prop = [0.0, 0.0]
    coalesce_prop = [0.0]
    t = 0

    while total_nodes_remaining(live_nodes) > 1 or total_nodes_remaining(dead_nodes):
        # Step 1: Calculate propensities
        total_prop = update_propensities(migration_prop, coalesce_prop, live_nodes, params)

        # Step 2: Determine next event
        event = get_next_event(migration_prop, coalesce_prop, total_prop, t)

        # Step 3: Handle activation of nodes
        next_node = None
        next_node_type = -1
        next_time = float('inf')
        for i in range(2):
            if not dead_nodes[i]:
                continue
            if dead_nodes[i][0].height < next_time:
                next_node = dead_nodes[i][0]
                next_time = next_node.height
                next_node_type = i
        
        if next_time < event.time:
            t = next_time
            live_nodes[next_node_type].append(next_node)
            dead_nodes[next_node_type].pop(0)
            continue

        # Step 4: Place event on tree
        next_node_nr = update_tree(live_nodes, event, next_node_nr)

        # Step 5: Keep track of time increment
        t = event.time 

    # Return sole remaining live node as root
    for i, node_list in enumerate(live_nodes):
        if node_list:
            return node_list[0]
        
    # Should not fall through
    raise RuntimeError("No live nodes remaining end of structured coalescent simulation!")


def total_nodes_remaining(nodes: List[List[SeedbankNode]]) -> int:
    """
    Calculate total number of live nodes remaining.

    Returns:
        Number of live nodes remaining.
    """
    return len(nodes[0]) + len(nodes[1])

def update_propensities(migration_prop: List[float], coalesce_prop: List[float], live_nodes: List[List[SeedbankNode]], params) -> float:
    """
    Obtain propensities (instantaneous reaction rates) for coalescence and
    migration events.

    Returns:
        Total reaction propensities.
    """
    total_prop = 0

    N_a = params["N"]
    k_a = len(live_nodes[1])
    k_d = len(live_nodes[0])
    m_ad = params["c"]
    m_da = params["c"] * params["K"]

    coalesce_prop[0] = k_a * (k_a - 1) / (2.0 * N_a)
    total_prop += coalesce_prop[0]

    migration_prop[1] = k_a * m_ad
    migration_prop[0] = k_d * m_da
    total_prop += migration_prop[1]
    total_prop += migration_prop[0]
    return total_prop

# def transition_model(string, c, K, N): # not implemented
#     if string == "getPopSize(1)":
#         return 1
#     elif string == "getBackwardRate(1, 0)":
#         return 2
#     elif string == "getBackwardRate(0, 1)":
#         return 3

def get_next_event(migration_prop: List[float], coalesce_prop: List[float], total_prop: float, t: float) -> SBEvent:
    """
    Obtrain type and location of next reaction.

    Parameters:
        migration_prop: Current migration propensities.
        coalesce_prop: Current coalescence propensities.
        t: Current time

    Returns:
        Event object describing next event.
    """
    # Get time of next event
    if total_prop > 0:
        t += np.random.exponential(1/total_prop)
    else:
        return NullEvent()
    
    # Select event type
    U = np.random.random() * total_prop

    if U < coalesce_prop[0]:
        return CoalescenceEvent(1, t)
    else:
        U -= coalesce_prop[0]
    
    if U < migration_prop[0]:
        return MigrationEvent(0, 1, t)
    else:
        U -= migration_prop[0]
    
    if U < migration_prop[1]:
        return MigrationEvent(1, 0, t)
    else:
        U -= migration_prop[1]

    # Should not fall through.
    raise RuntimeError("Structured coalescence event selection error.")


def update_tree(live_nodes: List[List[SeedbankNode]], event: SBEvent, next_node_nr: int) -> int:
    """
    Update tree with result of latest event.

    Parameters:
        next_node_nr: Integer identifier of last node added to tree.

    Returns:
        Updated next_node_nr.
    """
    if isinstance(event, CoalescenceEvent):
        # Randomly select node pair from active nodes
        daughter = select_random_node(live_nodes[1])
        son = select_random_sibling(live_nodes[1], daughter)

        # Create new parent node with appropriate ID and time
        parent = SeedbankNode()
        parent.nr = next_node_nr
        parent.ID = str(next_node_nr)
        parent.height = event.time
        next_node_nr += 1

        # Connect new parent to children
        parent.left = daughter
        parent.right = son
        son.parent = parent
        daughter.parent = parent

        # Ensure new parent is set to correct color
        parent.node_type = event.from_type

        # Update active nodes
        live_nodes[1].remove(son)
        idx = live_nodes[1].index(daughter)
        live_nodes[1][idx] = parent
    else: # Migration event ... what happens if null event?
        # Randomly select node with chose color
        migrator = select_random_node(live_nodes[event.from_type])

        # Record color change in change lists
        migrator.add_change(event.to_type, event.time)

        # Update active nodes
        live_nodes[event.from_type].remove(migrator)
        live_nodes[event.to_type].append(migrator)
    
    return next_node_nr
        

def select_random_node(node_list: List[SeedbankNode]) -> SeedbankNode:
    """
    Use beast RNG (here we actually use numpy RNG)to select random node from list.

    Returns:
        A randomly selected node.
    """
    indices = np.arange(len(node_list))
    random_index = np.random.choice(indices)
    return node_list[random_index]

def select_random_sibling(node_list: List[SeedbankNode], node: SeedbankNode) -> SeedbankNode:
    """
    Returns random node from list. excluding given node.

    Returns:
        A randomly selected node.
    """
    n = np.random.randint(len(node_list) - 1)
    idx_to_avoid = node_list.index(node)
    if n >= idx_to_avoid:
        n += 1
    return node_list[n]

def newick(root: SeedbankNode):
    """
    Print the newick string of Node root in the standard output.

    Parameter:
        root (SeedbankNode): Root of the tree to be converted to a newick string.
    """
    if root == None:
        return "ERROR"
    
    mig_str = ''
    last_height = root.height
    for i in range(root.n_type_changes):
        type_ = 'active' if root.change_types[i] == 1 else 'dormant'
        mig_str += f'{root.change_times[i]-last_height})[&type={type_}]:'
        last_height = root.change_times[i]

    newick_str = '(' * root.n_type_changes

    if root.left != None:
        newick_str += '(' + newick(root.left) + ',' + newick(root.right) + ')'

    type_ = 'active' if root.node_type == 1 else 'dormant'
    newick_str += f'{root.ID}[&type={type_}]:'
    newick_str += mig_str
    newick_str += "0.0;" if not root.parent else str(root.parent.height - last_height)

    return newick_str

def main():
    if len(sys.argv) != 2:
        print("Call", sys.argv[0], "<path to config file>")
        sys.exit(1)
    
    with open(sys.argv[1]) as f:
        config_data = json.load(f)

    n_leaves = config_data["n_leaves"]
    leaf_names = config_data["leaf_names"]
    leaf_types = config_data["leaf_types"]
    leaf_times = config_data["leaf_times"]
    params = {}
    params['c'] = config_data['c']
    params['K'] = config_data['K']
    params['N'] = config_data['N']

    root = simulate_tree(n_leaves, leaf_names, leaf_types, leaf_times, params)
    string = newick(root)

    with open('seedbank_tree_initialiser_output.txt', 'w') as output:
        output.write(string)


if __name__ == '__main__':
    main()