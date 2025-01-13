#############################################
# Creates a config file for serial samples
# 
# python3 create_serial_config.py
#############################################

import numpy as np
import json

def main():
    active_size = 50
    n_ranges = 5
    ranges = np.array([0, 1, 1, 2, 2, 3, 3, 4, 4, 5]).reshape(n_ranges, 2)
    sizes = np.array([10, 10, 10, 10, 10])
    times = np.empty(0)
    for i in range(n_ranges):
        lo, hi = ranges[i]
        sim = lo + np.random.random(sizes[i]) * (hi-lo)
        times = np.append(times, sim)
    
    dormant_size = 0
    n_ranges = 0
    ranges = np.array([]).reshape(n_ranges, 2)
    sizes = np.array([])
    for i in range(n_ranges):
        lo, hi = ranges[i]
        sim = lo + np.random.random(sizes[i]) * (hi-lo)
        times.append(sim)

    d = {
        "n_leaves": active_size + dormant_size,
        "leaf_names": list(map(str, list(range(active_size + dormant_size)))),
        "leaf_types": active_size * [1] + dormant_size * [0],
        "leaf_times": times.tolist(),
        "c": 0.5,
        "K": 1,
        "N": 1
    }
    
    with open("serial_config.json", "w") as out:
        json.dump(d, out, indent=4)

if __name__ == "__main__":
    main()
