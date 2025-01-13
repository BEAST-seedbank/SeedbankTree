import sys, json
from Bio import Phylo
from io import StringIO
import numpy as np
from treetime.seqgen import SeqGen
from treetime import GTR

##########################################################################
# To run this file:
# python3 simulate_mutation.py <config.json> <tree.txt>
# python3 simulate_mutation.py simulate_mutation.json simulate_tree_output.txt
# python3 simulate_mutation.py simulate_mutation.json serial_beast/seedbank_tree_initialiser_output.txt
#
# Creates output file:
# simulate_mutations_output.txt
# simulate_mutations_xml_output.txt

def conf_to_names (tree):
  clades = [tree.clade]
  while (clades):
    nxt = []
    for c in clades:
      if c.confidence and not c.name:
        c.name=str(c.confidence)
        c.confidence=None
        if c.clades:
          nxt += c.clades
    clades = nxt
  
def make_model(config_data):
  mu = config_data["mu"]
  alpha = config_data["alpha"]
   
  match config_data["model"]:
    case "jc":
      gtr = GTR.standard(model='jc', mu=mu, alphabet="nuc_nogap")
      gtr_d = GTR.standard(model='jc', mu=mu*alpha, alphabet="nuc_nogap")
    case "hky":
      pi = config_data["pi"]
      kappa = config_data["kappa"]
      gtr = GTR.standard(model='hky', mu=mu, pi=np.array(pi), kappa=kappa)
      gtr_d = GTR.standard(model='hky', mu=mu*alpha, pi=np.array(pi), kappa=kappa)
    case "gtr":
        pi = config_data["pi"]
        a,b,c,d,e,f = config_data["W_rates"]
        W = np.array([[0,a,b,c],[a,0,d,e],[b,d,0,f],[c,e,f,0]])
        gtr = GTR(alphabet="nuc_nogap")
        gtr.assign_rates(mu=mu, pi=pi, W=W)

        gtr_d = GTR(alphabet="nuc_nogap")
        gtr_d.assign_rates(mu=mu*alpha, pi=pi, W=W)
    case _:
      raise Exception("Unknown model")

  return gtr, gtr_d

def sanity_check(tree, sq, alpha):
  parent = get_parents(tree)

  active_branch_len = 0
  dormant_branch_len = 0
  active_mutations = 0
  dormant_mutations = 0
  for n in tree.find_clades():
    if (n == tree.root):
      active_branch_len += n.branch_length
      continue
    
    if ("active" in n.comment):
      active_branch_len += n.branch_length
      d = diff(n.ancestral_sequence, parent[n].ancestral_sequence)
      active_mutations += d
      print(f"node {n.name} | length {n.branch_length} | type: active | mutations: {d}")
    elif ("dormant" in n.comment):
      dormant_branch_len += n.branch_length
      d = diff(n.ancestral_sequence, parent[n].ancestral_sequence)
      dormant_mutations += d
      print(f"node {n.name} | length {n.branch_length} | type: dormant | mutations: {d}")
    else:
       print(n.comment)
       raise Exception()
  
  try:
    assert(np.abs(active_branch_len+dormant_branch_len - tree.total_branch_length()) < 1e-5)
  except:
    print(active_branch_len, dormant_branch_len)
    print(active_branch_len+dormant_branch_len)
    print(np.abs(active_branch_len+dormant_branch_len - tree.total_branch_length()))

  print()
  print(f"Summary statistics")
  print(f"active_branch_len: {active_branch_len}, active_mutations: {active_mutations}")
  print(f"dormant_branch_len: {dormant_branch_len}, dormant_mutations: {dormant_mutations}")
  print(f"alpha {alpha}, length {dormant_branch_len/active_branch_len}, mutations {dormant_mutations/active_mutations}")

  with open("simulate_mutation_output.txt", "a") as file:
      file.write("\n")
      file.write(f"Summary statistics\n")
      file.write(f"active_branch_len: {active_branch_len}, active_mutations: {active_mutations}\n")
      file.write(f"dormant_branch_len: {dormant_branch_len}, dormant_mutations: {dormant_mutations}\n")
      file.write(f"alpha {alpha}, length {dormant_branch_len/active_branch_len}, mutations {dormant_mutations/active_mutations}\n")

  return

def diff(seq1, seq2):
   return sum(1 for a, b in zip(seq1, seq2) if a != b)

def get_parents(tree):
    parents = {}
    for clade in tree.find_clades(order="level"):
        for child in clade:
            parents[child] = clade
    return parents



def main():
  if len(sys.argv) != 3:
        print("Call", sys.argv[0], "<path to config file> <path to tree file>")
        sys.exit(1)

  with open(sys.argv[1]) as config_file:
    config_data = json.load(config_file)
    seq_len = config_data["sequence_length"]

  with open(sys.argv[2]) as tree_file:
    tree = Phylo.read(tree_file, 'newick', rooted=True)
  conf_to_names(tree)

  gtr, gtr_d = make_model(config_data)
  sq = SeqGen(seq_len, tree=tree, gtr=gtr, gtr_d = gtr_d)
  sq.evolve_sb()

  aln = sq.get_aln(False)
  aln_all = sq.get_aln(True)
  depths = tree.depths()
  root_height = max(depths.values())

  with open("simulate_mutation_output.txt", "w") as file:
      file.write(f"sequence_length {seq_len}\n")
      file.write(f"number_of_leaves {tree.count_terminals()}\n")
      file.write(f"number_of_nodes {tree.count_terminals() + len((list(tree.get_nonterminals())))}\n")

      file.write("leaves\n")
      file.write("node_id node_height node_seq\n")
      for seqrecord in aln._records:
          file.write(f"{seqrecord.id} {seqrecord.seq}\n")

      file.write("\nall\n")
      file.write("node_id node_seq\n")
      for seqrecord in aln_all._records:
          file.write(f"{seqrecord.id} {seqrecord.seq}\n")

      file.write("\n\n")
      for n in tree.get_terminals():
         file.write(f"{n.name} {root_height - depths[n]}\n")
  
  with open("simulate_mutation_xml_output.txt", "w") as file:
     type_ = []
     time_ = []
     for n in tree.get_terminals():
        file.write(f"<sequence taxon='sample{n.name}' value='{''.join(n.ancestral_sequence)}'/>\n\t\t")
        type_.append(f"sample{n.name}=active" if "active" in n.comment else f"sample{n.name}=dormant")
        time_.append(f"sample{n.name}={root_height - depths[n]}")
     file.write("\n")
     file.write(f"<typeTraitSet spec='TraitSet' id='typeTraitSet' traitname='type' \n\t\tvalue='{', '.join(type_)}'>")
     file.write("\n")
     file.write(f"<timeTraitSet spec='TraitSet' id='timeTraitSet' traitname='date-backward' \n\t\tvalue='{', '.join(time_)}'>")
  
  sanity_check(tree, sq, config_data["alpha"])

if __name__ == "__main__":
    main() 