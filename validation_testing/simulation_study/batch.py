import numpy as np
import json
import os, shutil
from Bio import Phylo
from treetime.seqgen import SeqGen

from simulate_tree_serial import simulate_tree_serial, newick_serial
from simulate_tree_isochronous import simulate_tree_isochronous, newick_isochronous
from simulate_mutation import conf_to_names, make_model, sanity_check
from batch_options import *

def main():
    # Get batch configuration
    with open("./batch_config.json") as f:
        batch_config_data = json.load(f)
        batch_size = batch_config_data["batch_size"]
        rounds = batch_config_data["rounds"]
        folder = batch_config_data["folder"]
        shell = batch_config_data["shell"]

    # Get tree configuration
    with open("./simulate_tree_config.json") as f:
        tree_config_data = json.load(f)
        isochronous = tree_config_data["isochronous"]
        sample_sizes_1 = tree_config_data["sample_sizes_1"]
        sample_sizes_0 = tree_config_data["sample_sizes_0"]
        time_ranges = tree_config_data["time_ranges"]
        c = tree_config_data["c"]
        K = tree_config_data["K"]
        theta = tree_config_data["theta"]
        c_constant = tree_config_data["c_constant"]
        K_constant = tree_config_data["K_constant"]
        theta_constant = tree_config_data["theta_constant"]
        tree_constant = tree_config_data["tree_constant"]

        if not isochronous:
            assert len(sample_sizes_1) == len(sample_sizes_0)
            assert len(sample_sizes_1) == len(time_ranges)
            assert len(sample_sizes_1) != 0
    
    with open("./simulate_mutation.json") as config_file:
        mutation_config_data = json.load(config_file)
        seq_len = mutation_config_data["sequence_length"]
        alpha = mutation_config_data["alpha"]
        mu = mutation_config_data["mu"]
        kappa = mutation_config_data["kappa"]
        alpha_constant = mutation_config_data["alpha_constant"]
        mu_constant = mutation_config_data["mu_constant"]

    # round looping
    for i in range(1, rounds+1):
        try: 
            shutil.rmtree(f"./{folder}/{i}/")
        except Exception as e:
            pass

        os.mkdir(f"./{folder}/{i}/")

        # batch looping
        for j in range(1, batch_size+1):
            print(i, j)
            try:
                shutil.rmtree(f"./{folder}/{i}/{j}")
            except Exception as e:
                pass

            os.mkdir(f"./{folder}/{i}/{j}")
            shutil.copyfile(f"./{shell}", f"./{folder}/{i}/{j}/{j}.xml")

            #####################
            # Simulate tree
            if isochronous:
                sample_size_1 = sample_sizes_1[0]
                sample_size_2 = sample_sizes_0[0]
                types_list = [f"{i}={'active'}" for i in range(sample_size_1)] + [f"{i}={'dormant'}" for i in range(sample_size_1, sample_size_1 + sample_size_2)]
                times_list = [f"{i}=0.0" for i in range(sample_size_1 + sample_size_2)]

                root, island, time_vector, anc, dec_1, dec_2, log = simulate_tree_isochronous(sample_size_1, sample_size_2, c, K, theta)
                newick_str = newick_isochronous(root, island, time_vector, anc, dec_1, dec_2) + ";"

            else: # Randomly sample times
                times = np.empty(0)
                types = np.empty(0)
                for k, [lo, hi] in enumerate(time_ranges):
                    sim = lo + np.random.random(sample_sizes_1[k] + sample_sizes_0[k]) * (hi-lo)
                    times = np.append(times, sim)
                    
                    new_types = np.concatenate((np.ones(sample_sizes_1[k]), np.zeros(sample_sizes_0[k])))
                    np.random.shuffle(new_types)
                    types = np.append(types, new_types)

                n_leaves = len(times)
                leaf_types = list(map(int, types))
                leaf_times = sorted(list(times - min(times)))
                root, island, time_vector, anc, dec_1, dec_2 = simulate_tree_serial(n_leaves, leaf_types, leaf_times, c, K, theta)
                newick_str = newick_serial(root, island, time_vector, anc, dec_1, dec_2) + ";"

            with open(f"./{folder}/{i}/{j}/simulate_tree_output.txt", "w") as file:
                file.write(newick_str + "\n")
            
            ########################
            # Simulate mutations
            
            with open(f"./{folder}/{i}/{j}/simulate_tree_output.txt") as tree_file:
                tree = Phylo.read(tree_file, 'newick', rooted=True)
            conf_to_names(tree)

            gtr, gtr_d = make_model(mutation_config_data)
            sq = SeqGen(seq_len, tree=tree, gtr=gtr, gtr_d = gtr_d)
            sq.evolve_sb()
            
            seq_list = [f"<sequence taxon='{n.name}' value='{''.join(n.ancestral_sequence)}'/>\n\t\t" for n in tree.get_terminals()]
            if not isochronous:
                types_list = [f"{i}={'active' if t==1 else 'dormant'}" for i, t in enumerate(leaf_types[:n_leaves])]
                times_list = [f"{i}={t}" for i, t in enumerate(leaf_times[:n_leaves])]

            with open(f"./{folder}/{i}/{j}/{j}.xml") as f:
                new_text = f.read()
                new_text = new_text.replace('$$$REPLACE_SEQUENCES$$$', ''.join(seq_list))
                new_text = new_text.replace('$$$REPLACE_TYPES$$$', ', '.join(types_list))
                new_text = new_text.replace('$$$REPLACE_TIMES$$$', ', '.join(times_list))

                new_text = new_text.replace('$$$C_VALUE$$$', str(c))
                new_text = new_text.replace('$$$K_VALUE$$$', str(K))
                new_text = new_text.replace('$$$THETA_VALUE$$$', str(theta))
                new_text = new_text.replace('$$$ALPHA_VALUE$$$', str(alpha))
                new_text = new_text.replace('$$$MU_VALUE$$$', str(mu))
                new_text = new_text.replace('$$$KAPPA_VALUE$$$', str(kappa))

                new_text = new_text.replace('$$$C_PRIOR$$$', c_priors[c])
                new_text = new_text.replace('$$$K_PRIOR$$$', K_priors[K])
                new_text = new_text.replace('$$$THETA_PRIOR$$$', theta_priors[theta])
                new_text = new_text.replace('$$$ALPHA_PRIOR$$$', alpha_priors[alpha])
                new_text = new_text.replace('$$$MU_PRIOR$$$', mu_priors[mu])
                new_text = new_text.replace('$$$KAPPA_PRIOR$$$', kappa_priors[kappa])

                new_text = new_text.replace('$$$INITIAL_TREE$$$', 
                                        inital_tree_constant[tree_constant]
                                        .replace('$$$REPLACE_NEWICK$$$', newick_str.replace("&", "&amp;")))

                parameter_scalers_list = [scalers[param] if not constant else "<!-- -->" for (param, constant) in zip(["c", "K", "theta", "alpha", "mu"], [c_constant, K_constant, theta_constant, alpha_constant, mu_constant])]
                new_text = new_text.replace('$$$PARAMETER_SCALERS$$$', '\n\t\t'.join(parameter_scalers_list))

                new_text = new_text.replace('$$$TREE_OPERATORS$$$', tree_operators if not tree_constant else "")

            with open(f"./{folder}/{i}/{j}/{j}.xml", "w") as f:
                f.write(new_text)
        

if __name__ == "__main__":
    main()
