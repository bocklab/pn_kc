

# 200109 copy from /Users/zhengz11/myscripts/bocklab_git/bocklab/zhihao/mushroom_py/v10/191211-update_claw_distribution.py
# 200109 this will generate supplementary figure 2 (claw distribution), 200326, PNKC2019_v9_fig_200313DB-ZZfixedSuppl6B.pptx

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ana = ana_all_rd

claw_num = [len(i.segments.ids) for i in list(ana.row_neurons.values())]

ax = plt.subplot(1, 1, 1)
sns.distplot(claw_num, bins=13, kde=False)
# save_path = '/Users/zhengz11/myscripts/data_results/191209-modeling_claw_null_model/'
# plt.savefig(save_path + "191211-claw_num_dist_AllKCs.png", dpi=600, bbox_inches='tight')
