import json, os, shutil
from Bio import Phylo
from treetime.seqgen import SeqGen
from treetime import GTR

from simulate_tree import simulate_tree, newick
from simulate_mutation import conf_to_names, make_model, sanity_check

def main():
    with open("./batch_config.json") as f:
        batch_serial_config_data = json.load(f)
        batch_size = batch_serial_config_data["batch_size"]
        folder = batch_serial_config_data["folder"]
        shell = batch_serial_config_data["shell"]

    with open('./simulate_tree.json') as f:
            tree_config_data = json.load(f)

            sample_size_1 = tree_config_data["sample_size_1"]
            sample_size_2 = tree_config_data["sample_size_2"]
            relative_seedbank_size = tree_config_data["relative_seedbank_size"]
            transition_rate = tree_config_data["transition_rate"]
            theta = tree_config_data["theta"]
    
    with open("./simulate_mutation.json") as f:
            mutation_config_data = json.load(f)
            seq_len = mutation_config_data["sequence_length"]

    types_list = [f"{i}={"active"}" for i in range(sample_size_1)] + [f"{i}={"dormant"}" for i in range(sample_size_1, sample_size_1 + sample_size_2)]
    times_list = [f"{i}=0.0" for i in range(sample_size_1 + sample_size_2)]

    for i in range(1, batch_size+1):
        print(i)
        try:
            shutil.rmtree(f"./{folder}/{i}")
        except Exception as e:
            pass

        os.mkdir(f"./{folder}/{i}")
        shutil.copyfile(f"./{shell}", f"./{folder}/{i}/{i}.xml")

        root, island, time_vector, anc, dec_1, dec_2, log = simulate_tree(sample_size_1, sample_size_2, transition_rate, relative_seedbank_size, theta)
        with open(f"./{folder}/{i}/simulate_tree_output.txt", "w") as file:
            newick_str = newick(root, island, time_vector, anc, dec_1, dec_2)
            file.write(newick_str + ";\n")

        ##########################################
        # mutations

        with open(f"./{folder}/{i}/simulate_tree_output.txt") as tree_file:
            tree = Phylo.read(tree_file, 'newick', rooted=True)
        conf_to_names(tree)

        gtr, gtr_d = make_model(mutation_config_data)
        sq = SeqGen(seq_len, tree=tree, gtr=gtr, gtr_d = gtr_d)
        sq.evolve_sb()
        
        sanity_check(tree, sq, mutation_config_data["alpha"], f"./{folder}/{i}/simulate_mutation_output.txt")  
        
        seq_list = [f"<sequence taxon='{n.name}' value='{''.join(n.ancestral_sequence)}'/>\n\t\t" for n in tree.get_terminals()]

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