
# start a connection to fafb_c to download data
import sys
import mushroom_2to3.connect_path as cp
import mushroom_2to3.analysis_routine as ar
sys.path.append('/Users/zhengz11/myscripts/git_clone/pn_kc/')
# import mushroom_2to3.connect as cc

# credential, to delete when push to remote
sys.path.append('/Users/zhengz11/myscripts/mushroom_v9/credential/')
from fafb_tokens import token
fafb_c = cp.fafb_connection(token)

import mushroom_2to3.connect_path as cp
import json

path = "/Users/zhengz11/myscripts/git_clone/pn_kc/test/"

def save_json(js, file_name):
    with open(file_name, 'w+') as file:
        json.dump(js, file)
    file.close()

def load_json(path):
    with open(path) as outfile:
        r = json.load(outfile)
    return r

# save_json(pn_skids, save_path + "PN")
# save_json(rd, save_path + "RandomDraw")
# save_json(t1p, save_path + "t1p")
# save_json(bundle, save_path + "bundle")

save_path = "/Users/zhengz11/myscripts/git_clone/pn_kc/data/skids/"
pn_skids = load_json(save_path + "pn")
rd = load_json(save_path + "RandomDraw")


t_pn_skids = pn_skids
t_kc_skids = rd[:30]

t_skids = t_pn_skids + t_kc_skids


for i in t_skids:
    cp.save_compact_sk(fafb_c, i, path)

cp.save_annotations_for_skeleton(fafb_c, t_skids, path)
cp.save_neurons_names(fafb_c, t_skids, path)
cp.save_root_node(fafb_c, t_skids, path)
cp.save_annotated_annotations(fafb_c, 'glom_class', 'id', path)
cp.save_annotated_annotations(fafb_c, 'kc_class', 'id', path)
cp.save_pre_post_info(fafb_c, t_pn_skids, t_kc_skids, path, 'testing')


# testing

def load_json(path):
    with open(path) as outfile:
        r = json.load(outfile)
    return r

ana_t = ar.Analysis.init_connectivity(path, t_pn_skids, t_kc_skids, 'testing')







# test the get function with local paths
import mushroom_2to3.connect_path as cp
import mushroom_2to3.analysis_routine as ar
path = "/Users/zhengz11/myscripts/git_clone/pn_kc/test/"
t1 = cp.get_root_node(path, t_pn_skids[0])


# save binary connectivity matrix between PNs and KCs
pk_conn = ana.conn_data['pn_kc_in_claw_unit'].conn['1s'].copy()
# np.save(local_path + "data/200514_pn_kc_bi_conn.npy", pk_conn)