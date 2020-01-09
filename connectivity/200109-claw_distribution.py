

# 200109 copy from /Users/zhengz11/myscripts/bocklab_git/bocklab/zhihao/mushroom_py/v10/191211-update_claw_distribution.py
# 200109 this will generate supplementary figure 1 (claw distribution)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

%run startup_py3.py
%run load_pn_metadata_v2.py
%run medium_term_functions.py

save_path = "/Users/zhengz11/myscripts/data_results/191112-covariance_matrices/191202_updated/"

fafb_c = cc.fafb_connection(token)

pn_skids = cc.get_skids_from_annos(
    fafb_c, [['right_calyx_PN'], ['has_bouton']], ["multiglomerular PN"])

rd = cc.get_skids_from_annos(fafb_c,
                             [['Random Draw 1 KC', 'Random Draw 2 KC'], ['Complete']],
                             ['KCaBp', 'KCyd'])

t1p = cc.get_skids_from_annos(fafb_c,
                             [['T1+ Complete']])

# ana_rd = ar.Analysis.init_connectivity(fafb_c, pn_skids, rd)
ana_all_rd = ar.Analysis.init_connectivity(fafb_c, pn_skids, rd + t1p)

ana = ana_all_rd

claw_num = [len(i.segments.ids) for i in list(ana.row_neurons.values())]

ax = plt.subplot(1, 1, 1)
sns.distplot(claw_num, bins=13, kde=False)
# save_path = '/Users/zhengz11/myscripts/data_results/191209-modeling_claw_null_model/'
# plt.savefig(save_path + "191211-claw_num_dist_AllKCs.png", dpi=600, bbox_inches='tight')
