

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

%run startup_py3.py
%run load_pn_metadata_v2.py
%run medium_term_functions.py

# 200109 copy from:
# /Users/zhengz11/myscripts/bocklab_git/bocklab/zhihao/mushroom_py/v10/191127-rerun_main_matrices.py
# 200109 currently generates:
# Fig 1D,E (conditional input matrix bouton-based null model) and
# Fig 2F,G (conditional input matrix claw-based null model)

##-------------------------------------------------
# 191127 re-run the matrix with glomerular names colored according to behavioral Significance
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
save_path = '/Users/zhengz11/myscripts/data_results/191119-rerun_main_figs/'
fig.savefig(save_path + "191128-allKCs_CondInputMatrix_LabelColored.png", bbox_inches='tight')


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
fig.savefig(save_path + "191128-allKCs_CondInputMatrix_LabelColored_NullModel.png", bbox_inches='tight')


##-----------------------------------------------------
# Preferential convergence: use number of claws as null models
input_prob = conn_data.conn['1s'].sum(0)
input_prob = input_prob / sum(input_prob)

stat = [get_raw_inputs(shuffle_glom_kc_w_prob(ob_conn, input_prob)) for i in range(num_exp)]

stat = np.array(stat)
sd = np.nanstd(stat, axis=0)
avg = np.nanmean(stat, axis=0)

ob_ci = get_raw_inputs(ob_conn)

comm_zscore = np.divide(np.subtract(ob_ci, avg), sd)

# clustering
cm_zs = PairMatrix('', comm_zscore, glom_idx_ids)

reorder_idx = km_cluster(cm_zs.conn)
# reorder_idx = reorder(NewClusterOrder, glom_idx_ids)
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
fig.savefig(save_path + "191128-AllKCs_ClawNullModel_LabelColored.png", bbox_inches='tight')


##------------------------------------------
# null model for the Preferential convergence result
null_conn = shuffle_glom_kc_w_prob(ob_conn, input_prob)
ob_ci = get_raw_inputs(null_conn)
comm_zscore = np.divide(np.subtract(ob_ci, avg), sd)

# clustering
cm_zs = PairMatrix('', comm_zscore, glom_idx_ids)

reorder_idx = km_cluster(cm_zs.conn)
# reorder_idx = reorder(NewClusterOrder, glom_idx_ids)
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

fig.savefig(save_path + "191128-AllKCs_ClawNullModel_NullModel_LabelColored.png", bbox_inches='tight')
