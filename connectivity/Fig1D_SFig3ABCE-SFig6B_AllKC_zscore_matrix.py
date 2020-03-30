# Figre 1D (200326, PNKC2019_v9_fig_200313DB-ZZfixedSuppl6B.pptx), z-score matrix with PNs and all KCs (random draw and t1p), random bouton null model

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# %run startup_py3.py
# %run load_pn_metadata_v2.py
# %run medium_term_functions.py


##-------------------------------------------------
ana = ana_all_rd

conn_data = ana.conn_data['glom_kc_in_claw_unit']
ob_conn, glom_prob, glom_idx_ids = get_conn_prob_idx(conn_data)

num_exp = 1000

stat = [get_raw_inputs(shuffle_glom_kc_w_prob(ob_conn, glom_prob)) for i in range(num_exp)]
stat = np.array(stat)
sd = np.nanstd(stat, axis=0)
avg = np.nanmean(stat, axis=0)

ob_ci = get_raw_inputs(ob_conn)
comm_zscore = np.divide(np.subtract(ob_ci, avg), sd)

cm_zs = PairMatrix('', comm_zscore, glom_idx_ids)

# reorder_idx = km_cluster(cm_zs.conn)

# this is the order used in final figures 'ClusterOrder0707'
reorder_idx = reorder(ClusterOrder0707, glom_idx_ids)
t1_zs = cm_zs.reorder(reorder_idx, return_new=True)

# plotting z score matrix
fig, ax1 = plt.subplots()
t1 = t1_zs;
gloms = df_lookup('glom_id',t1.col_ids,'short_glom_name',glom_btn_table)
sns.heatmap(t1.conn, xticklabels=gloms, yticklabels=gloms, ax=ax1, vmin=-8.53, vmax=8.53, cmap="RdBu_r")

ax1.tick_params(bottom=False,labeltop=True, top=True, labelbottom=False)
ax1.tick_params(axis='x',labelrotation=90)

# run and get tbl in 191029-bouton-KC-representations_per_PN.py
col_list = t1.col_ids
col_colors = df_lookup('short_glom_name', gloms, 'color', tbl)

for x in [ax1.get_xticklabels(), ax1.get_yticklabels()]:
    for idx, tick in enumerate(x):
        tick.set_color(col_colors[idx])
        if col_list[idx] in comm_ids:
            tick.set_weight("extra bold")
#            tick.set_bbox(dict(ec='green', fc=None, alpha=0.05))

ax1.set_aspect("equal")
fig.set_size_inches(16,12)
plt.show()

# save_path = '/Users/zhengz11/myscripts/data_results/191119-rerun_main_figs/'
# fig.savefig(save_path + "191128-allKCs_CondInputMatrix_LabelColored.png", bbox_inches='tight')


# SFig3A, B, C
# the following is copied from 200206-main_matrices_zscore_distribution.py
# 200219, distribution of z scores for the community (first 10), wider community (first 19), and remaining entries
save_path = "/Users/zhengz11/myscripts/data_results/200206-zscore_distribution/"

t1 = t1_zs.conn[:10,:10]
t1_f = t1.flatten()
# mean 5.7 std 2.9

t2 = t1_zs.conn[:19,:19]
t2_idx = np.ones(t2.shape)
t2_idx[:10,:10] = 0
t2_f = t2[np.where(t2_idx)]
# mean ± std 2.3 ± 1.9

t3 = t1_zs.conn.copy()
t3_idx = np.ones(t3.shape)
t3_idx[:19,:19] = 0
t3_f = t3[np.where(t3_idx)]
# mean ± std -0.5 ± 1.5

bin_set = np.arange(-5,16)
xtick_labels = np.arange(-5,16,2)
fig, ax = plt.subplots()
sns.distplot(t1_f, kde=False, bins=bin_set, ax=ax)
ax.set_xticks(xtick_labels)
ax.set_yticks(np.arange(0,22,4))
fig.set_size_inches(12,8)
# fig.savefig(save_path + '200219-zscore_distr_commmunity.png', bbox_inches='tight')

fig, ax = plt.subplots()
sns.distplot(t2_f, kde=False, bins=bin_set, ax=ax)
ax.set_xticks(xtick_labels)
ax.set_yticks(np.arange(0,70,10))
fig.set_size_inches(12,8)
# fig.savefig(save_path + '200219-zscore_distr_bigger_commmunity.png', bbox_inches='tight')

fig, ax = plt.subplots()
sns.distplot(t3_f, kde=False, bins=bin_set, ax=ax)
ax.set_xticks(xtick_labels)
ax.set_yticks(np.arange(0,800,100))
fig.set_size_inches(12,8)
# fig.savefig(save_path + '200219-zscore_distr_non_commmunity.png', bbox_inches='tight')


import scipy as sc
sc.stats.ks_2samp(t1_f, t2_f)
# Ks_2sampResult(statistic=0.5960919540229885, pvalue=2.220446049250313e-16)


sc.stats.ks_2samp(t1_f, t3_f)
# Ks_2sampResult(statistic=0.8619373776908024, pvalue=2.3314683517128287e-15)

sc.stats.ks_2samp(t2_f, t3_f)
# Ks_2sampResult(statistic=0.624372614736337, pvalue=1.3322676295501878e-15)



# SFig3E
# copied from 191127-rerun_main_matrices.py
##--------------------------------------------------
# re-run the control matrix
sfl_conn = shuffle_glom_kc_w_prob(ob_conn, glom_prob)
stat = [get_raw_inputs(shuffle_glom_kc_w_prob(sfl_conn, glom_prob)) for i in range(num_exp)]
stat = np.array(stat)
sd = np.nanstd(stat, axis=0)
avg = np.nanmean(stat, axis=0)

ob_ci = get_raw_inputs(sfl_conn)
comm_zscore = np.divide(np.subtract(ob_ci, avg), sd)

cm_zs = PairMatrix('', comm_zscore, glom_idx_ids)

reorder_idx = km_cluster(cm_zs.conn)
t1_zs = cm_zs.reorder(reorder_idx, return_new=True)

# plotting z score matrix
fig, ax1 = plt.subplots()
t1 = t1_zs;
gloms = df_lookup('glom_id',t1.col_ids,'short_glom_name',glom_btn_table)
sns.heatmap(t1.conn, xticklabels=gloms, yticklabels=gloms, ax=ax1, vmin=-8.53, vmax=8.53, cmap="RdBu_r")

ax1.tick_params(bottom=False,labeltop=True, top=True, labelbottom=False)
ax1.tick_params(axis='x',labelrotation=90)

col_list = t1.col_ids
col_colors = df_lookup('short_glom_name', gloms, 'color', tbl)

for x in [ax1.get_xticklabels(), ax1.get_yticklabels()]:
    for idx, tick in enumerate(x):
        tick.set_color(col_colors[idx])
        if col_list[idx] in comm_ids:
            tick.set_weight("extra bold")

ax1.set_aspect("equal")
fig.set_size_inches(16,12)
plt.show()
# fig.savefig(save_path + "191128-allKCs_CondInputMatrix_LabelColored_NullModel.png", bbox_inches='tight')



# SFig6B (200326, PNKC2019_v9_fig_200313DB-ZZfixedSuppl6B.pptx)
# copied from 200206-main_matrices_zscore_distribution.py
# (1) The following is copied from above
ana = ana_all_rd
conn_data = ana.conn_data['glom_kc_in_claw_unit']
ob_conn, glom_prob, glom_idx_ids = get_conn_prob_idx(conn_data)
num_exp = 1000
stat = [get_raw_inputs(shuffle_glom_kc_w_prob(ob_conn, glom_prob)) for i in range(num_exp)]
stat = np.array(stat)
sd = np.nanstd(stat, axis=0)
avg = np.nanmean(stat, axis=0)
ob_ci = get_raw_inputs(ob_conn)
comm_zscore = np.divide(np.subtract(ob_ci, avg), sd)
cm_zs = PairMatrix('', comm_zscore, glom_idx_ids)
zs_t1 = cm_zs.conn.flatten().copy()

#(2) The following is basically copied from Fig2D_SFig6A_ob_vs_RandomClawModel.py
ana = ana_all_rd
conn_data = ana.conn_data['glom_kc_in_claw_unit']
ob_conn, glom_prob, glom_idx_ids = get_conn_prob_idx(conn_data)
stat = [get_raw_inputs(i) for i in shuffle_glom_kc_iterate(ob_conn, 1000)]
stat = np.array(stat)
sd = np.nanstd(stat, axis=0)
avg = np.nanmean(stat, axis=0)
ob_ci = get_raw_inputs(ob_conn)
comm_zscore = np.divide(np.subtract(ob_ci, avg), sd)
# clustering
claw_cm_zs = PairMatrix('', comm_zscore, glom_idx_ids)
zs_t2 = claw_cm_zs.conn.flatten().copy()

# save_path = "/Users/zhengz11/myscripts/data_results/200206-zscore_distribution/200314-zscore_compilation/"

fig, ax = plt.subplots()
bin_set = np.arange(-5,16,0.5)
xtick_labels = np.arange(-5,16,1)
ax.hist([zs_t1, zs_t2], bins=bin_set, label=['observed vs. random bouton', 'observed vs. random claw'])
ax.legend(loc='upper right')
ax.set_xticks(xtick_labels)
ax.set_xlabel('Z scores')
fig.set_size_inches(12,8)
# fig.savefig(save_path + '200314-zscore_hist_ObvsRandBtn_ObvsRandClaw.png', bbox_inches='tight')

# Ob vs. random bouton, mean -0.044, std 2.11
# Ob vs. random claw, mean -0.058, std 1.47
sc.stats.ks_2samp(zs_t1, zs_t2)
# P < 1e-10
