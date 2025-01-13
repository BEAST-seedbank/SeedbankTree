import json
import os, shutil

def main():
    # Get batch configuration
    with open("./batch_config.json") as f:
        batch_config_data = json.load(f)
        batch_size = batch_config_data["batch_size"]
        rounds = batch_config_data["rounds"]
        folder = batch_config_data["folder"]
        original_batch_directory = batch_config_data["original_batch_directory"]
        shell = batch_config_data["shell"]

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

            with open(f"./{folder}/{i}/{j}/{j}.xml", "w") as f:
                f.write(new_text)
        

if __name__ == "__main__":
    main()
