
import sys
sys.path.append('/Users/zhengz11/myscripts/git_clone/pn_kc/')
import mushroom_2to3.connect as cc

# credential, to delete when push to remote
sys.path.append('/Users/zhengz11/myscripts/mushroom_v9/credential/')
from fafb_tokens import token
fafb_c = cc.fafb_connection(token)

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


import json
save_path = "/Users/zhengz11/myscripts/git_clone/pn_kc/data/skids/"

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
# load_json(save_path + "t1p.txt")
# load_json(save_path + "bundle.txt")
t_pn_skids = pn_skids
t_kc_skids = rd[:30]

t_skids = t_pn_skids + t_kc_skids








path =

for all skids
save_compact_sk(connection, skid, path)


# run it for
'glom_class'
'kc_class'
save_annotated_annotations



# need to call for specific pairing of pre_skids, and post_skids
save_pre_post_info(connection, pre_skids, post_skids):
