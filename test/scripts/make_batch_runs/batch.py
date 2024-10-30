import sys, json, os, shutil
from Bio import Phylo
from io import StringIO
import numpy as np
from treetime.seqgen import SeqGen
from treetime import GTR

from simulate_tree import simulate_tree, newick
from simulate_mutation import conf_to_names, make_model, sanity_check

def main():
    i = 1
    failed = 0
    while i < 51:
        print(i)
        try:
            shutil.rmtree(f"./batch_10/{i}")
        except Exception as e:
            pass

        os.mkdir(f"./batch_10/{i}")
        shutil.copyfile("./shell3.xml", f"./batch_10/{i}/{i}.xml")

        with open('./simulate_tree.json') as f:
            tree_config_data = json.load(f)

        sample_size_1 = tree_config_data["sample_size_1"]
        sample_size_2 = tree_config_data["sample_size_2"]
        relative_seedbank_size = tree_config_data["relative_seedbank_size"]
        transition_rate = tree_config_data["transition_rate"]
        theta = tree_config_data["theta"]

        root, island, time_vector, anc, dec_1, dec_2, log = simulate_tree(sample_size_1, sample_size_2, transition_rate, relative_seedbank_size, theta)
        with open(f"./batch_10/{i}/simulate_tree_output.txt", "w") as file:
            newick_str = newick(root, island, time_vector, anc, dec_1, dec_2)
            file.write(newick_str + ";\n")

        ##########################################
        ##########################################
        # mutations
        ##########################################

        with open("./simulate_mutation.json") as config_file:
            mutation_config_data = json.load(config_file)
            seq_len = mutation_config_data["sequence_length"]

        with open(f"./batch_10/{i}/simulate_tree_output.txt") as tree_file:
            tree = Phylo.read(tree_file, 'newick', rooted=True)
        conf_to_names(tree)

        #####################
        # dormant_branch_len = 0
        # for n in tree.find_clades():
        #     if (n == tree.root):
        #         continue
        #     if ("dormant" in n.comment):
        #         dormant_branch_len += n.branch_length
        
        # if not (0.25 <= dormant_branch_len/tree.total_branch_length() <= 0.75):
        #     failed += 1
        #     continue
        #####################

        gtr, gtr_d = make_model(mutation_config_data)
        sq = SeqGen(seq_len, tree=tree, gtr=gtr, gtr_d = gtr_d)
        sq.evolve_sb()

        aln = sq.get_aln(False)
        aln_all = sq.get_aln(True)
        depths = tree.depths()
        root_height = max(depths.values())
 
        seq_list = []
        for n in tree.get_terminals():
            seq_list.append(f"<sequence taxon='{n.name}' value='{''.join(n.ancestral_sequence)}'/>\n\t\t")

        with open(f"./batch_10/{i}/{i}.xml") as f:
            newText = f.read()
            newText = newText.replace('$$$REPLACE_TARGET$$$', ''.join(seq_list))
            newText = newText.replace('$$$REPLACE_TARGET_2$$$', newick_str.replace("&", "&amp;")+";")

        with open(f"./batch_10/{i}/{i}.xml", "w") as f:
            f.write(newText)
        
        sanity_check(tree, sq, mutation_config_data["alpha"], f"./batch_10/{i}/simulate_mutation_output.txt")  
        
        i += 1
    
    print(f"NUMBER FAILED TREES {failed}")

if __name__ == "__main__":
    main() 