#
# This is a combination of several of files to have this code work on
# a generic python instance (e.g. binder.)
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


# exec(open(local_path + "/connectivity/load_pn_metadata_v2.py").read())

pn_skids = load_json(os.path.join(data_dir,  "skids/pn"))
rd = load_json(os.path.join(data_dir,  "skids/RandomDraw"))
t1p = load_json(os.path.join(data_dir, "skids/t1p"))
bundle = load_json(os.path.join(data_dir, "skids/bundle"))

ana_all_rd = ar.Analysis.init_connectivity(data_dir, pn_skids, rd + t1p, 'pn_all_kc')
ana_rd = ar.Analysis.init_connectivity(data_dir, pn_skids, rd, 'pn_rd_kc')

# exec(open(local_path + "/connectivity/load_pn_tbl.py").read())
