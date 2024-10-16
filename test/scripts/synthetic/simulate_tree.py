import sys, time, os, json
import numpy as np

##########################################################################
# To run this file:
# python3 simulate_tree.py <config.json>
# python3 simulate_tree.py simulate_tree.json
#
# Creates output file:
# simulate_tree_output.txt
#

def newick(id, island, time_vector, anc, dec_1, dec_2):
    type = "active" if island[id] == 0 else "dormant"

    length = 0.0 if anc[id] == -1 else time_vector[anc[id]] - time_vector[id]

    newickStr = ""
    if dec_1[id] != -1:
        newickStr += "(" + newick(dec_1[id], island, time_vector, anc, dec_1, dec_2)
        if dec_2[id] != -1:
            newickStr += "," + newick(dec_2[id], island, time_vector, anc, dec_1, dec_2)
        newickStr += ")"
    newickStr += str(id) + "[&type=\"" + type + "\"]:" + str(length)

    return newickStr

# transition-less newick-string
def newick2(id, island, time_vector, anc, dec_1, dec_2):
    type = "active" if island[id] == 0 else "dormant"

    length = 0.0 if anc[id] == -1 else time_vector[anc[id]] - time_vector[id]

    newickStr = ""
    if dec_1[id] != -1:
        newickStr += "(" + newick(dec_1[id], island, time_vector, anc, dec_1, dec_2)
        if dec_2[id] != -1:
            newickStr += "," + newick(dec_2[id], island, time_vector, anc, dec_1, dec_2)
        newickStr += ")"
    newickStr += str(id) + "[&type=\"" + type + "\"]:" + str(length)

    return newickStr

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
    log = ["sample_size_1 " + str(sample_size_1), "sample_size_2 " + str(sample_size_2), "relative_seedbank_size " + str(relative_seedbank_size), "transition_rate " + str(transition_rate), ""]

    np.random.seed(int(time.time() * os.getpid()) % (2 ** 32 - 1))
    next_parent = sample_size_1 + sample_size_2
    anc = [-1] * next_parent
    dec_1 = [-1] * next_parent
    dec_2 = [-1] * next_parent
    island = [0] * next_parent
    time_vector = [0.0] * next_parent
    sim_time = 0.0
    active_1 = list(range(sample_size_1))
    active_2 = list(range(sample_size_1, sample_size_1 + sample_size_2))

    for i in active_2:
        island[i] = 1

    log.append(f"active dormant | c_rate m_rate_1 m_rate_2 total_rate | sim_time coin")
    while len(active_1) + len(active_2) > 1:
        # elif model == 2:
        c_rate_1 = len(active_1) * (len(active_1) - 1) / 2
        c_rate_2 = 0.0
        m_rate_1 = len(active_1) * transition_rate
        m_rate_2 = len(active_2) * transition_rate * relative_seedbank_size

        total_rate = c_rate_1 + c_rate_2 + m_rate_1 + m_rate_2
        sim_time += np.random.exponential(1.0 / total_rate)
        coin = np.random.uniform()
        log.append(f"{len(active_1)} {len(active_2)} | {c_rate_1} {m_rate_1} {m_rate_2} {total_rate} | {sim_time} {coin}")
        
        if coin < c_rate_1 / total_rate:
            # coalescent event
            child_1 = np.random.randint(len(active_1))
            child_2 = np.random.randint(len(active_1))
            while child_2 == child_1:
                child_2 = np.random.randint(len(active_1))
            time_vector.append(sim_time)
            anc.append(-1)
            island.append(0)
            anc[active_1[child_1]] = next_parent
            anc[active_1[child_2]] = next_parent
            dec_1.append(-1)
            dec_2.append(-1)
            dec_1[next_parent] = active_1[child_1]
            dec_2[next_parent] = active_1[child_2]
            active_1[min(child_1, child_2)] = next_parent
            next_parent += 1
            active_1.pop(max(child_1, child_2))
        elif coin < (total_rate - m_rate_2) / total_rate:
            # active to dormant event
            child_1 = np.random.randint(len(active_1))
            anc.append(-1)
            island.append(1)
            time_vector.append(sim_time)
            anc[active_1[child_1]] = next_parent
            dec_1.append(-1)
            dec_2.append(-1)
            dec_1[next_parent] = active_1[child_1]
            active_2.append(next_parent)
            active_1.pop(child_1)
            next_parent += 1
        else:
            # dormant to active event
            child_1 = np.random.randint(len(active_2))
            anc.append(-1)
            island.append(0)
            time_vector.append(sim_time)
            anc[active_2[child_1]] = next_parent
            dec_1.append(-1)
            dec_2.append(-1)
            dec_1[next_parent] = active_2[child_1]
            active_1.append(next_parent)
            active_2.pop(child_1)
            next_parent += 1

    with open("simulate_tree_output.txt", "w") as file:
        newick_str = newick(next_parent - 1, island, time_vector, anc, dec_1, dec_2)
        file.write(newick_str + ";\n")
    
    with open("simulate_tree_output_2.txt", "w") as file:
        # update dec_1, dec_2, anc arrays to skip transitions
        root = next_parent-1
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

        newick_str = newick(next_parent - 1, island, time_vector, anc, dec_1, dec_2)
        newick_str = newick_str.replace("dormant", "active")
        file.write(newick_str + ";\n")


    with open("simulate_tree_log.txt", "w") as file:
        file.write("\n".join(log))

if __name__ == "__main__":
    main()