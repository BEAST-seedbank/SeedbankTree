import sys, os, json
import time
import numpy as np

##########################################################################
# To run this file:
# python3 simulate_tree.py <config.json>
# python3 simulate_tree.py simulate_tree.json
#
# Creates output file:
# simulate_tree_output.txt
#

def newick_serial(id, types, time_vector, anc, dec_1, dec_2):
    type = "active" if types[id] == 1 else "dormant"

    length = 0.0 if anc[id] == -1 else time_vector[anc[id]] - time_vector[id]

    newickStr = ""
    if dec_1[id] != -1:
        newickStr += "(" + newick_serial(dec_1[id], types, time_vector, anc, dec_1, dec_2)
        if dec_2[id] != -1:
            newickStr += "," + newick_serial(dec_2[id], types, time_vector, anc, dec_1, dec_2)
        newickStr += ")"
    newickStr += str(id) + "[&type=\"" + type + "\"]:" + str(length)

    return newickStr

def simulate_tree_serial(n_leaves, leaf_types, leaf_times, transition_rate, relative_seedbank_size, theta, log=None, seed=None):
    if not seed:
        seed = int(time.time() * os.getpid()) % (2 ** 32 - 1)
    np.random.seed(seed)

    # leaf nodes already have spots reserved
    next_node = n_leaves
    anc = [-1] * next_node
    dec_1 = [-1] * next_node
    dec_2 = [-1] * next_node
    types = leaf_types
    time_vector = leaf_times
    sim_time = 0.0
    live_active_lineages = []
    live_dormant_lineages = []

    # first leaf node included (at time 0)
    if types[0] == 1:
        live_active_lineages.append(0)
    else:
        live_dormant_lineages.append(0)
    leaf_i = 1

    if log == []:
        log.append(f"active dormant | c_rate m_rate_1 m_rate_2 total_rate | sim_time coin")
    while leaf_i < n_leaves or len(live_active_lineages) + len(live_dormant_lineages) > 1:
        c_rate_1 = len(live_active_lineages) * (len(live_active_lineages) - 1) / (2 * theta)
        c_rate_2 = 0.0
        m_rate_1 = len(live_active_lineages) * transition_rate
        m_rate_2 = len(live_dormant_lineages) * transition_rate * relative_seedbank_size

        total_rate = c_rate_1 + c_rate_2 + m_rate_1 + m_rate_2
        time_interval = np.random.exponential(1.0 / total_rate)

        # sample
        if leaf_i < n_leaves and sim_time + time_interval >= leaf_times[leaf_i]:
            sim_time = leaf_times[leaf_i]
            if log:
                log.append(f"{len(live_active_lineages)} {len(live_dormant_lineages)} | {c_rate_1} {m_rate_1} {m_rate_2} {total_rate} | {sim_time} SAMPLE")

            if leaf_types[leaf_i] == 1:
                live_active_lineages.append(leaf_i)
            else:
                live_dormant_lineages.append(leaf_i)
            leaf_i += 1
            continue

        sim_time += time_interval
        coin = np.random.uniform()
        if log:
            log.append(f"{len(live_active_lineages)} {len(live_dormant_lineages)} | {c_rate_1} {m_rate_1} {m_rate_2} {total_rate} | {sim_time} {coin}")
        
        if coin < c_rate_1 / total_rate:
            # coalescent event
            child_1 = np.random.randint(len(live_active_lineages))
            child_2 = np.random.randint(len(live_active_lineages))
            while child_2 == child_1:
                child_2 = np.random.randint(len(live_active_lineages))
            time_vector.append(sim_time)
            anc.append(-1)
            types.append(1)
            anc[live_active_lineages[child_1]] = next_node
            anc[live_active_lineages[child_2]] = next_node
            dec_1.append(-1)
            dec_2.append(-1)
            dec_1[next_node] = live_active_lineages[child_1]
            dec_2[next_node] = live_active_lineages[child_2]
            live_active_lineages[min(child_1, child_2)] = next_node
            next_node += 1
            live_active_lineages.pop(max(child_1, child_2))
        elif coin < (total_rate - m_rate_2) / total_rate:
            # active to dormant event
            child_1 = np.random.randint(len(live_active_lineages))
            anc.append(-1)
            types.append(0)
            time_vector.append(sim_time)
            anc[live_active_lineages[child_1]] = next_node
            dec_1.append(-1)
            dec_2.append(-1)
            dec_1[next_node] = live_active_lineages[child_1]
            live_dormant_lineages.append(next_node)
            live_active_lineages.pop(child_1)
            next_node += 1
        else:
            # dormant to active event
            child_1 = np.random.randint(len(live_dormant_lineages))
            anc.append(-1)
            types.append(1)
            time_vector.append(sim_time)
            anc[live_dormant_lineages[child_1]] = next_node
            dec_1.append(-1)
            dec_2.append(-1)
            dec_1[next_node] = live_dormant_lineages[child_1]
            live_active_lineages.append(next_node)
            live_dormant_lineages.pop(child_1)
            next_node += 1

    if log:
        return next_node-1, types, time_vector, anc, dec_1, dec_2, log
    return next_node-1, types, time_vector, anc, dec_1, dec_2

def main():
    if len(sys.argv) != 2:
        print("Call", sys.argv[0], "<path to config file>")
        sys.exit(1)

    with open(sys.argv[1]) as f:
        config_data = json.load(f)

    sample_size_1 = config_data["sample_size_1"]
    sample_size_2 = config_data["sample_size_2"]
    relative_seedbank_size = config_data["relative_seedbank_size"]
    transition_rate = config_data["transition_rate"]
    theta = config_data["theta"]
    log = ["sample_size_1 " + str(sample_size_1), "sample_size_2 " + str(sample_size_2), "relative_seedbank_size " + str(relative_seedbank_size), "transition_rate " + str(transition_rate), "theta " + str(theta)]

    root, island, time_vector, anc, dec_1, dec_2, log = simulate_tree_serial(sample_size_1, sample_size_2, transition_rate, relative_seedbank_size, theta, log)

    with open("simulate_tree_output.txt", "w") as file:
        newick_str = newick_serial(root, island, time_vector, anc, dec_1, dec_2)

        file.write(newick_str + ";\n")
    
    with open("simulate_tree_output_2.txt", "w") as file:
        # update dec_1, dec_2, anc arrays to skip transitions
        stack = [root]
        
        while stack:
            node = stack.pop()

            if (dec_1[node] == -1 and dec_2[node] == -1): #leaf node
                continue

            # update left
            left = dec_1[node]
            while (dec_1[left] != -1 and dec_2[left] == -1):
                left = dec_1[left]
            dec_1[node] = left
            anc[left] = node
            stack.append(left)

            # update right
            right = dec_2[node]
            while (dec_1[right] != -1 and dec_2[right] == -1):
                right = dec_1[right]
            dec_2[node] = right
            anc[right] = node
            stack.append(right)

        newick_str = newick_serial(root, island, time_vector, anc, dec_1, dec_2)
        newick_str = newick_str.replace("dormant", "active")
        file.write(newick_str + ";\n")


    with open("simulate_tree_log.txt", "w") as file:
        file.write("\n".join(log))

if __name__ == "__main__":
    main()