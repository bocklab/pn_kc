#
# This is a modification of the previous startup.py to work as an import.
#

import json
import os
import sys

pnkc_root = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
data_dir = os.path.join(pnkc_root, "data")

# Add the local_path to make mushroom_2to3 visible
sys.path.insert(1, os.path.join(sys.path[0], pnkc_root))

import mushroom_2to3.connect as cc
import mushroom_2to3.analysis_routine as ar

from mushroom_2to3.shuffle import *
from mushroom_2to3.build_connectivity import *
from mushroom_2to3.detect_community import *


def load_json(path):
    with open(path) as outfile:
        r = json.load(outfile)
    return r

pn_skids = load_json(os.path.join(data_dir,  "skids/pn"))
rd = load_json(os.path.join(data_dir,  "skids/RandomDraw"))
t1p = load_json(os.path.join(data_dir, "skids/t1p"))
bundle = load_json(os.path.join(data_dir, "skids/bundle"))

try:
    with open(os.path.join(data_dir, "node_to_segment_cache.pkl"), "rb") as f:
        node_to_segment_cache = pickle.load(f)
except:
    print("Warning: Could not read existing node_to_segment_cache.")
    node_to_segment_cache = {}

ana_all_rd = None
ana_rd = None

def get_ana_rd():
    global ana_rd
    if ana_rd is None:
        ana_rd = ar.Analysis.init_connectivity(data_dir, pn_skids, rd, 'pn_rd_kc', node_to_segment_cache)
    return ana_rd

def get_ana_all_rd():
    global ana_all_rd
    if ana_all_rd is None:        
        ana_all_rd = ar.Analysis.init_connectivity(data_dir, pn_skids, rd + t1p, 'pn_all_kc', node_to_segment_cache)
    return ana_all_rd


def save_cache(filename="node_to_segment_cache.pkl"):
    if len(node_to_segment_cache) > 0:
        with open(os.path.join(data_dir, "node_to_segment_cache.pkl"), "wb") as f:
            pickle.dump(node_to_segment_cache, f)


# exec(open(local_path + "/connectivity/load_pn_tbl.py").read())
