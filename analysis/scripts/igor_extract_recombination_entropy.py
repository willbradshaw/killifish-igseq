#!/usr/bin/env python

import numpy as np, pandas as pd, argparse, os, sys
print(sys.version)
from pygor.models.genmodel import GenModel
from pygor.models.entropy import compute_event_entropy
from pygor.models.entropy import compute_model_entropy


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

# TODO: Delete unneeded aux functions
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

###############################
## Get overall model entropy ##
###############################

print "Extracting overall model entropy..."
h_model = {"event":["model"], "h": [compute_model_entropy(g)], "id":[args.id]}
h_model = pd.DataFrame(h_model)
print "done."
print "Model entropy:", h_model["h"]

###################################
## Get different event entropies ##
###################################

print "Extracting separate event entropies..."
events = [e.nickname for e in g.events]
print "   ", "Events:", events
entropies = [compute_event_entropy(g.get_event(e,True).name, g) for e in events]
print "   ", "Event entropies:", entropies
h_events = {"event":events, "h":entropies, "id":args.id}
h_events = pd.DataFrame(h_events)
print "done."

#######################################
## Combine model and event entropies ##
#######################################

print "Combining event and model entropies..."
h_all = pd.concat([h_model, h_events])
print "done."

#############################
## Write entropies to file ##
#############################

outpath = os.path.join(args.outdir, "entropies.tsv")
h_all.to_csv(outpath, index=False, sep = "\t")
