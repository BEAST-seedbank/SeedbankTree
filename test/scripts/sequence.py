import msprime
import numpy as np

def main():
  ts = msprime.sim_ancestry(samples=10, ploidy=1, sequence_length=100, random_seed=8888)

  print(ts)

  kappa = 0.8
  freqs = [0.22, 0.26, 0.21, 0.31]
  hky=msprime.HKY(kappa, freqs)

  mts = msprime.sim_mutations(ts, rate=0.05, random_seed=8888, model=hky)

  print(mts.draw_text())
  print(mts)

  num_muts = np.zeros(mts.num_sites, dtype=int)
  for site in mts.sites():
      num_muts[site.id] = len(site.mutations)  # site.mutations is a list of mutations at the site

  # Print out some info about mutations per site
  for nmuts, count in enumerate(np.bincount(num_muts)):
      info = f"{count} sites"
      if nmuts > 1:
          info += f", with IDs {np.where(num_muts==nmuts)[0]},"
      print(info, f"have {nmuts} mutation" + ("s" if nmuts != 1 else ""))

  print("Genotypes")
  for v in mts.variants():
    print(f"Site {v.site.id}: {v.genotypes}")
    if v.site.id >= 10:  # only print up to site ID 4
        print("...")
        break
    
  print("Genotypes")

  samp_ids = mts.samples()
  print("          ID of individual: ", " ".join([f"{mts.node(s).individual:3}" for s in samp_ids]))
  print("       ID of (sample) node: ", " ".join([f"{s:3}" for s in samp_ids]))
  for v in mts.variants():
      site = v.site
      alleles = np.array(v.alleles)
      print(f"Site {site.id} (ancestral state '{site.ancestral_state}')",  alleles[v.genotypes])
      if site.id >= 4:  # only print up to site ID 4
          print("...")
          break
      

  seqs = [[*seq] for seq in mts.alignments()]
  np.random.seed(1)
  for j in range(100):
     if seqs[0][j] == "N":
        draw = np.random.uniform(0, 1)
        if draw < 0.22:
           nuc = "A"
        elif draw < 0.48:
           nuc = "C"
        elif draw < 0.69:
           nuc = "G"
        else:
           nuc = "T"

        for i in range(10):
           seqs[i][j] = nuc
  
  for seq in seqs:
    print("".join(seq))


if __name__ == "__main__":
  main()