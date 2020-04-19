import sys
import os
sys.path.append('/Users/zhengz11/myscripts/git_clone/pn_kc/')

import json
import mushroom_2to3.connect_path as cp
import mushroom_2to3.analysis_routine as ar


# credential, to delete when push to remote
sys.path.append('/Users/zhengz11/myscripts/mushroom_v9/credential/')
from fafb_tokens import token
fafb_c = cp.fafb_connection(token)


def save_json(js, file_name):
    with open(file_name, 'w+') as file:
        json.dump(js, file)
    file.close()

def load_json(path):
    with open(path) as outfile:
        r = json.load(outfile)
    return r

path = "/Users/zhengz11/myscripts/git_clone/pn_kc/"

##----------------------------------------
import mushroom_2to3.connect as cc
pn_skids = cc.get_skids_from_annos(
    fafb_c, [['right_calyx_PN'], ['has_bouton']], ["multiglomerular PN"])

# get KC skeleton ids from CATMAID
# rd is random draw manually traced KCs
rd = cc.get_skids_from_annos(fafb_c,
                             [['Random Draw 1 KC', 'Random Draw 2 KC'], ['Complete']],
                             ['KCaBp', 'KCyd'])

t1p = cc.get_skids_from_annos(fafb_c,
                             [['T1+ Complete']])


bundle = cc.get_skids_from_annos(fafb_c,
    [['Bundle 1 Seed', 'Different Tracing Protocol in Bundle 1'], ['Complete']], ['KCaBp', 'KCyd'])

save_path = path + "data/skids/"

if not os.path.exists(save_path):
    os.makedirs(save_path)

save_json(pn_skids, save_path + "PN")
save_json(rd, save_path + "RandomDraw")
save_json(t1p, save_path + "t1p")
save_json(bundle, save_path + "bundle")


# download the different KC subtypes
y = cc.get_skids_from_annos(fafb_c,[['KCy']])

prime = cc.get_skids_from_annos(fafb_c,
                             [["KCa'B'", "KCa'B'ap", "KCa'B'm", "KCa'B'x"]])

ab = cc.get_skids_from_annos(fafb_c,
                             [['KCaBc', 'KCaBs', 'KCaBx']])

save_json(y, save_path + "kcy")
save_json(prime, save_path + "kcprime")
save_json(ab, save_path + "kcab")

# pn_skids = load_json(save_path + "PN")
# rd = load_json(save_path + "RandomDraw")
# bundle = load_json(save_path + "bundle")
# load_json(save_path + "t1p.txt")
# load_json(save_path + "bundle.txt")

all_skids = pn_skids + rd + t1p + bundle

for i in all_skids:
    cp.save_compact_sk(fafb_c, i, path)

cp.save_annotations_for_skeleton(fafb_c, all_skids, path)
cp.save_neurons_names(fafb_c, all_skids, path)
cp.save_root_node(fafb_c, all_skids, path)
cp.save_annotated_annotations(fafb_c, 'glom_class', 'id', path)
cp.save_annotated_annotations(fafb_c, 'kc_class', 'id', path)
cp.save_pre_post_info(fafb_c, pn_skids, rd + t1p, path, 'pn_all_kc')

# 200329
cp.save_pre_post_info(fafb_c, pn_skids, rd + t1p, path, 'pn_rd_kc')

## testing analysis
ana_all_rd = ar.Analysis.init_connectivity(path, pn_skids, rd + t1p, 'pn_all_kc')

ana_rd = ar.Analysis.init_connectivity(path, pn_skids, rd + t1p, 'pn_rd_kc')


# save pickle python object for the code in PN bouton clustering

import pandas as pd
import numpy as np
from pymaid import rmaid
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import pymaid
import mushroom_2to3.neurogenesis as neurogenesis
# %run startup_py3.py
# %run load_pn_metadata_v1.py

# pn_skids = cc.get_skids_from_annos(fafb_c, [['right_calyx_PN'], ['has_bouton']], ["multiglomerular PN"])

pns_ms = neurogenesis.init_from_skid_list(fafb_c, pn_skids)

import pickle
path = local_path + "data/pn_bouton_clusters/"

with open(path + "pns_ms.pkl", 'wb') as f:
    pickle.dump(pns_ms, f, -1)
    
nl = [pymaid.get_neuron([str(i) for i in j]) for j in [pn_skids[:40], pn_skids[40:80], pn_skids[80:]]]
pns_pm = nl[0] + nl[1] + nl [2]

with open(path + "pns_pm.pkl", 'wb') as f:
    pickle.dump(pns_pm, f, -1)

ca = pymaid.get_volume('MB_CA_R')

with open(path + "ca.pkl", 'wb') as f:
    pickle.dump(ca, f, -1)

df = pd.read_excel(local_path + 'data/180613-pn_subtypes.xlsx')


# for loading the pickle pymaid neuronlist data
path = local_path + "data/pn_bouton_clusters/"
with open(path + "pns_pm.pkl", 'rb') as f:
    pns_pm = pickle.load(f)

    
with open(path + "pns_ms.pkl", 'rb') as f:
    pns_ms = pickle.load(f)
    
with open(path + "ca.pkl", 'rb') as f:
    ca = pickle.load(f)
    
    
# trying
import pickle
path = local_path + "data/pn_bouton_clusters/"
for i in pns_pm:
    with open(path + "pns_pm/" + "{}.pkl".format(i.skeleton_id), 'wb') as f:
        pickle.dump(i, f, -1)