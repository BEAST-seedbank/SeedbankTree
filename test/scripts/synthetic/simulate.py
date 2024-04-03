import sys, time, os, json
import numpy as np

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

def main():
    if len(sys.argv) != 2:
        print("Call", sys.argv[0], "<path to config file>")
        sys.exit(1)

    with open(sys.argv[1]) as f:
        config_data = json.load(f)

    weak_seed_bank_mean_delay = config_data["weak_seed_bank_mean_delay"]
    sample_size_1 = config_data["sample_size_1"]
    sample_size_2 = config_data["sample_size_2"]
    model = config_data["model"]
    island_2_relative_size = config_data["island_2_relative_size"]
    migration_rate = config_data["migration_rate"]

    beta = 1.0 / weak_seed_bank_mean_delay
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

    while len(active_1) + len(active_2) > 1:
        if model == 0:
            c_rate_1 = len(active_1) * (len(active_1) - 1) * beta * beta
            c_rate_2 = 0.0
            m_rate_1 = 0.0
            m_rate_2 = 0.0
        elif model == 1:
            c_rate_1 = len(active_1) * (len(active_1) - 1)
            c_rate_2 = len(active_2) * (len(active_2) - 1) * island_2_relative_size
            m_rate_1 = len(active_1) * migration_rate
            m_rate_2 = len(active_2) * migration_rate * island_2_relative_size
        elif model == 2:
            c_rate_1 = len(active_1) * (len(active_1) - 1) / 2
            c_rate_2 = 0.0
            m_rate_1 = len(active_1) * migration_rate
            m_rate_2 = len(active_2) * migration_rate * island_2_relative_size
        else:
            print("Unrecognized model specification.")
            print("See dev.cfg for implemented models.")
            sys.exit(1)

        total_rate = c_rate_1 + c_rate_2 + m_rate_1 + m_rate_2
        sim_time += np.random.exponential(1.0 / total_rate)
        coin = np.random.uniform()
        
        if coin < c_rate_1 / total_rate:
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
        elif coin < (c_rate_1 + c_rate_2) / total_rate:
            child_1 = np.random.randint(len(active_2))
            child_2 = np.random.randint(len(active_2))
            while child_2 == child_1:
                child_2 = np.random.randint(len(active_2))
            time_vector.append(sim_time)
            anc.append(-1)
            island.append(1)
            anc[active_2[child_1]] = next_parent
            anc[active_2[child_2]] = next_parent
            active_2[min(child_1, child_2)] = next_parent
            next_parent += 1
            active_2.pop(max(child_1, child_2))
        elif coin < (total_rate - m_rate_2) / total_rate:
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

    for i in range(next_parent):
        print(i, "|", dec_1[i], "|", dec_2[i])

    print()
    print(newick(next_parent - 1, island, time_vector, anc, dec_1, dec_2))

if __name__ == "__main__":
    main()