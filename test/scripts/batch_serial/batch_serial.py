import numpy as np
import json
import os, shutil
from Bio import Phylo
from treetime.seqgen import SeqGen

from simulate_tree import simulate_tree, newick
from simulate_mutation import conf_to_names, make_model, sanity_check

def main():
    # Get batch configuration
    with open("./batch_serial_config.json") as f:
        batch_serial_config_data = json.load(f)
        batch_size = batch_serial_config_data["batch_size"]
        folder = batch_serial_config_data["folder"]
        shell = batch_serial_config_data["shell"]

    # Get serial configuration
    with open("./simulate_tree_serial_config.json") as f:
        serial_config_data = json.load(f)
        sample_sizes_1 = serial_config_data["sample_sizes_1"]
        sample_sizes_0 = serial_config_data["sample_sizes_0"]
        time_ranges = serial_config_data["time_ranges"]
        c = serial_config_data["c"]
        K = serial_config_data["K"]
        theta = serial_config_data["theta"]

        assert len(sample_sizes_1) == len(sample_sizes_0)
        assert len(sample_sizes_1) == len(time_ranges)
        assert len(sample_sizes_1) != 0

    for i in range(1, batch_size+1):
        print(i)
        try:
            shutil.rmtree(f"./{folder}/{i}")
        except Exception as e:
            pass

        os.mkdir(f"./{folder}/{i}")
        shutil.copyfile(f"./{shell}", f"./{folder}/{i}/{i}.xml")

        # Randomly sample times
        times = np.empty(0)
        types = np.empty(0)
        for j, [lo, hi] in enumerate(time_ranges):
            sim = lo + np.random.random(sample_sizes_1[j] + sample_sizes_0[j]) * (hi-lo)
            times = np.append(times, sim)

            types = np.append(types, np.ones(sample_sizes_1[j]))
            types = np.append(types, np.zeros(sample_sizes_0[j]))

        # Simulate serial tree
        n_leaves = len(times)
        leaf_types = list(map(int, types))
        leaf_times = sorted(list(times - min(times)))
        root, island, time_vector, anc, dec_1, dec_2 = simulate_tree(n_leaves, leaf_types, leaf_times, c, K, theta)
        newick_str = newick(root, island, time_vector, anc, dec_1, dec_2)
        with open(f"./{folder}/{i}/simulate_serial_tree_output.txt", "w") as file:
            file.write(newick_str + "\n")
        
        # Simulate mutations
        with open("./simulate_mutation.json") as config_file:
            mutation_config_data = json.load(config_file)
            seq_len = mutation_config_data["sequence_length"]
        
        with open(f"./{folder}/{i}/simulate_serial_tree_output.txt") as tree_file:
            tree = Phylo.read(tree_file, 'newick', rooted=True)
        conf_to_names(tree)

        gtr, gtr_d = make_model(mutation_config_data)
        sq = SeqGen(seq_len, tree=tree, gtr=gtr, gtr_d = gtr_d)
        sq.evolve_sb()

        sanity_check(tree, sq, mutation_config_data["alpha"], f"./{folder}/{i}/simulate_mutation_output.txt")  

        seq_list = [f"<sequence taxon='{n.name}' value='{''.join(n.ancestral_sequence)}'/>\n\t\t" for n in tree.get_terminals()]
        types_list = [f"{i}={"active" if t==1 else "dormant"}" for i, t in enumerate(leaf_types[:n_leaves])]
        times_list = [f"{i}={t}" for i, t in enumerate(leaf_times[:n_leaves])]

        with open(f"./{folder}/{i}/{i}.xml") as f:
            newText = f.read()
            newText = newText.replace('$$$REPLACE_SEQUENCES$$$', ''.join(seq_list))
            newText = newText.replace('$$$REPLACE_TYPES$$$', ', '.join(types_list))
            newText = newText.replace('$$$REPLACE_TIMES$$$', ', '.join(times_list))
            newText = newText.replace('$$$REPLACE_NEWICK$$$', newick_str.replace("&", "&amp;")+";")

        with open(f"./{folder}/{i}/{i}.xml", "w") as f:
            f.write(newText)
    

if __name__ == "__main__":
    main()
