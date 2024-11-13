import json
import re
import pickle
import numpy as np
import msprime
import hopkins
import argparse

parser = argparse.ArgumentParser(
                    prog='msprime_simulator')
parser.add_argument('--output', type=str, default="output")
parser.add_argument('--Ne', type=int, default=1)
parser.add_argument('--length', type=int, default=1000)
parser.add_argument('--track_length', type=int, default=1)
parser.add_argument('--nsample', type=int, default=500)
parser.add_argument('--mu', type=float, default=0.0000006)
parser.add_argument('--r_m', type=float, default=0.1)
parser.add_argument('--model', type=str, default="beta1.0001")
parser.add_argument('--continuous_genome', action='store_true')

def save_metadata(params):
    params_dict = vars(params)
    with open(args.output + ".json", 'w') as metadata_file:
        json.dump(params_dict, metadata_file, indent=4)

args = parser.parse_args()

l = args.length  # number of genes
t = args.track_length  # tract length
nsample = args.nsample  # the number of genomes sampled
mu = args.mu  # mutation rate
r_m = args.r_m

r = r_m * mu

match = re.match(r"^beta([0-9]*\.[0-9]+)$", args.model)
if match:
    a = float(match.group(1))
    model = msprime.BetaCoalescent(alpha=a)
elif  args.model == "kingman":
    model = None
else:
    raise ValueError(f"Invalid model argument: {input_str}")

ts = msprime.sim_ancestry(nsample,
                          model=model,
                          population_size=args.Ne,
                          ploidy=1,
                          sequence_length=l,
                          gene_conversion_rate=r,
                          gene_conversion_tract_length=t,
                          discrete_genome=not args.continuous_genome,
                          )

mts = msprime.sim_mutations(ts, rate=mu)

pairs = [(i, j) for i in range(nsample) for j in range(i)]
distance_list = mts.diversity(pairs, mode='site')
distance_matrix = np.zeros((nsample,nsample))

for pair_index, distance in enumerate(distance_list):
    (i,j) = pairs[pair_index]
    distance_matrix[i, j] = distance
    distance_matrix[j, i] = distance

mts.dump(args.output)
with open(args.output + "_ds", 'wb') as file:
    pickle.dump(distance_list, file)
    pickle.dump(distance_matrix, file)

save_metadata(args)
