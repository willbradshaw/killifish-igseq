#!/usr/bin/env python

import numpy as np, pandas as pd, argparse, os, sys
print(sys.version)
from pygor.models.genmodel import GenModel


#####################
## Parse arguments ##
#####################

parser = argparse.ArgumentParser()
parser.add_argument("indir")
parser.add_argument("outdir")
parser.add_argument("id")
args = parser.parse_args()

marginals_path = os.path.join(os.getcwd(), args.indir, "final_marginals.txt")
parms_path = os.path.join(os.getcwd(), args.indir, "final_parms.txt")

##################
## Define model ##
##################

g = GenModel(parms_path, marginals_path)

def iterate_event(nickname, attr):
    """Iterate over the realizations in an event for a specific attribute."""
    R = g.get_event(nickname, True).realizations
    return np.array([getattr(r, attr) for r in R])

def event_order(nickname):
    """Get the sort order for realizations in an event."""
    return np.argsort(iterate_event(nickname, "index"))

def event_names(nickname):
    """Get the sorted names of realizations in an event."""
    return iterate_event(nickname, "name")[event_order(nickname)]

def event_values(nickname):
    """Get the sorted values of realizations in an event."""
    return iterate_event(nickname, "value")[event_order(nickname)]

def event_marginals(nickname):
    """Get the marginal probabilities of realizations in an event."""
    return g.marginals[0][nickname]

def make_distr_df(dct):
    """Make a distribution dataframe from an input dictionary."""
    dct.update({"id":args.id})
    return pd.DataFrame(dct)

distributions = {}

#################################
## Get V-deletion distribution ##
#################################

# Get embedded distribution
v3_values = event_values("v_3_del")
v3v_marginals = event_marginals("v_3_del")
# Normalise by V-frequency
v_marginals = event_marginals("v_choice")
v3_marginals = np.sum(v3v_marginals.T * v_marginals, 1)
# Save distribution
distributions["V3"] = make_distr_df({"n":v3_values, "p":v3_marginals, 
    "event":"v_3_del"})

#################################
## Get J-deletion distribution ##
#################################

# Get embedded distribution
j5_values = event_values("j_5_del")
j5j_marginals = event_marginals("j_5_del")
# Normalise by J-frequency
vj_marginals = event_marginals("j_choice")
j_marginals = np.sum(vj_marginals.T * v_marginals.astype(float), 1)
j5_marginals = np.sum(j5j_marginals.T * j_marginals, 1)
# Save distribution
distributions["J5"] = make_distr_df({"n":j5_values, "p":j5_marginals,
    "event":"j_5_del"})

###########################################
## Get D-deletion distribution (5-prime) ##
###########################################

# Get embedded distribution
d5_values = event_values("d_5_del")
d5d_marginals = event_marginals("d_5_del")
# Normalise by D-frequency
vdj_marginals = event_marginals("d_gene")
d_marginals = np.einsum("ijk -> i",
    vdj_marginals.T * v_marginals * j_marginals[:, np.newaxis])
d5_marginals = np.sum(d5d_marginals.T * d_marginals, 1)
# Save distribution
distributions["D5"] = make_distr_df({"n":d5_values, "p":d5_marginals,
    "event":"d_5_del"})

###########################################
## Get D-deletion distribution (5-prime) ##
###########################################

# Get embedded distribution
d3_values = event_values("d_3_del")
d3d_marginals = event_marginals("d_3_del")
# Normalise by D-frequency and 5' deletions
d3_marginals = np.einsum("ijk -> i",
        d3d_marginals.T * d5d_marginals.T * d_marginals)

# Save distribution
distributions["D3"] = make_distr_df({"n":d3_values, "p":d3_marginals,
    "event":"d_3_del"})

#################################
## Write distributions to file ##
#################################

output = pd.DataFrame()
for k in distributions.keys():
    output = output.append(distributions[k])
    print len(distributions[k]), len(output)

assert len(output) > 0

outpath = os.path.join(args.outdir, "deletes.tsv")
output.to_csv(outpath, index=False, sep = "\t")
