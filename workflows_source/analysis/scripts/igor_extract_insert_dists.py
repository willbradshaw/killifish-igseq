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

################################
## Get VD insert distribution ##
################################

vd_values = event_values("vd_ins")
vd_marginals = event_marginals("vd_ins")
distributions["VD"] = make_distr_df({"n":vd_values, "p":vd_marginals,
    "event":"vd_ins"})

################################
## Get DJ insert distribution ##
################################

dj_values = event_values("dj_ins")
dj_marginals = event_marginals("dj_ins")
distributions["DJ"] = make_distr_df({"n":dj_values, "p":dj_marginals,
    "event":"dj_ins"})

#################################
## Write distributions to file ##
#################################

output = pd.DataFrame()
for k in distributions.keys():
    output = output.append(distributions[k])
    print len(distributions[k]), len(output)

assert len(output) > 0

outpath = os.path.join(args.outdir, "inserts.tsv")
output.to_csv(outpath, index=False, sep = "\t")
