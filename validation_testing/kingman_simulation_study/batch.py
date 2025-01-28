import json
import os, shutil
import re
from Bio import Phylo

def main():
    # Get batch configuration
    with open("./batch_config.json") as f:
        batch_config_data = json.load(f)
        batch_size = batch_config_data["batch_size"]
        rounds = batch_config_data["rounds"]
        folder = batch_config_data["folder"]
        original_batch_directory = batch_config_data["original_batch_directory"]
        shell = batch_config_data["shell"]
        use_effective_mu = batch_config_data["use_effective_mu"]
        alpha = float(batch_config_data["alpha"])
        mu = float(batch_config_data["mu"])

        skip_runs = batch_config_data["skip_runs"]

    # round looping
    for i in range(1, rounds+1):
        try: 
            shutil.rmtree(f"./{folder}/{i}/")
        except Exception as e:
            pass

        os.mkdir(f"./{folder}/{i}/")

        # batch looping
        for j in range(1, batch_size+1):
            if skip_runs and j in skip_runs[i-1]:
                print("skipped")
                continue

            print(i, j)
            try:
                shutil.rmtree(f"./{folder}/{i}/{j}")
            except Exception as e:
                pass

            os.mkdir(f"./{folder}/{i}/{j}")
            shutil.copyfile(f"./{shell}", f"./{folder}/{i}/{j}/{j}.xml")

            #####################
            # Retrieve sequences and times

            seq_line_list = []
            with open(f"./{original_batch_directory}/{i}/{j}/{j}.xml") as f:
                for k, line in enumerate(f):
                    if 3 <= k <= 27:
                        assert "sequence" in line
                        seq_line_list.append(line)
                    elif k == 37:
                        time_values = line
                    elif k > 37:
                        break

            with open(f"./{folder}/{i}/{j}/{j}.xml") as f:
                new_text = f.read()
                new_text = new_text.replace('$$$REPLACE_SEQUENCES$$$', ''.join(seq_line_list))
                new_text = new_text.replace('$$$REPLACE_TIME_VALUES$$$', time_values)

                with open(f"./{original_batch_directory}/{i}/{j}/simulate_tree_output.txt") as tree_file:
                    tree = Phylo.read(tree_file, 'newick', rooted=True)
                
                with open(f"./{original_batch_directory}/{i}/{j}/simulate_tree_output.txt") as tree_file:
                    newick = tree_file.read().rstrip()

                # Calculate effective mu value
                if use_effective_mu:
                    active_branch_len, dormant_branch_len = 0, 0
                    for n in tree.find_clades():
                        if ("active" in n.comment):
                            active_branch_len += n.branch_length
                        elif ("dormant" in n.comment):
                            dormant_branch_len += n.branch_length
                    total_branch_len = active_branch_len + dormant_branch_len
                    effective_mu = (mu * active_branch_len + mu * alpha * dormant_branch_len) / total_branch_len
                    new_text = new_text.replace('$$$REPLACE_MU$$$', str(effective_mu))
                else:
                    new_text = new_text.replace('$$$REPLACE_MU$$$', str(mu))

                # Replace tree (if $$$REPLACE_NEWICK$$$ is in text)
                newick = re.sub(r'\[.*?\]', '', newick)
                newick = re.sub(r'&', '&amp;', newick)
                new_text = new_text.replace("$$$REPLACE_NEWICK$$$", newick)

            with open(f"./{folder}/{i}/{j}/{j}.xml", "w") as f:
                f.write(new_text)
        

if __name__ == "__main__":
    main()
