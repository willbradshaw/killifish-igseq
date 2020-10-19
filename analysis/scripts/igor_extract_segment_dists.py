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
print marginals_path
print parms_path

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

def event_marginals(nickname):
    """Get the marginal probabilities of realizations in an event."""
    return g.marginals[0][nickname]

def make_distr_df(dct):
    """Make a distribution dataframe from an input dictionary."""
    dct.update({"id":args.id})
    return pd.DataFrame(dct)

########################
## Get V distribution ##
########################

distributions = {}

v_names = event_names("v_choice")
v_marginals = event_marginals("v_choice")

distributions["V"] = make_distr_df({"v":v_names, "p":v_marginals})

################################
## Get VJ and J distributions ##
################################

j_names = event_names("j_choice")
vj_conditionals = event_marginals("j_choice") # J marginals *given* V choice

# Get marginal J probabilities
vj_marginals = vj_conditionals.T * v_marginals.astype(float)
j_marginals = np.sum(vj_marginals, 1)
distributions["J"] = make_distr_df({"j":j_names, "p":j_marginals})

# Get joint V/J probabilities
distributions["VJ"] = make_distr_df({"v": np.tile(v_names, len(j_names)),
    "j": np.repeat(j_names, len(v_names)),
    "p": vj_marginals.T.flatten()
    })

#################################
## Get VDJ and D distributions ##
#################################

d_names = event_names("d_gene")
vdj_conditionals = event_marginals("d_gene")

# Get marginal D probabilities
vdj_marginals = vdj_conditionals.T * vj_marginals
d_marginals = np.einsum("ijk -> i", vdj_marginals)
distributions["D"] = make_distr_df({"d":d_names, "p":d_marginals})

# Get joint V/D/J probabilities
distributions["VDJ"] = make_distr_df({
    "v": np.tile(v_names, len(d_names)*len(j_names)),
    "d": np.repeat(d_names, len(v_names)*len(j_names)),
    "j": np.tile(np.repeat(j_names, len(v_names)), len(d_names)),
    "p": vdj_marginals.flatten()
    })

#################################
## Write distributions to file ##
#################################

for k in distributions.keys():
    outpath = os.path.join(args.outdir, "{}.tsv".format(k))
    distributions[k].to_csv(outpath, index=False, sep="\t")
