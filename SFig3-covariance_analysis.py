
# this will generate suppl fig 4 A - D (4 matrices of covariance analysis)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
# need ana_rd and ana_all_rd from analysis.py

ana = ana_all_rd
conn_data = ana.conn_data['glom_kc_in_claw_unit']
ob_conn, glom_prob, glom_idx_ids = get_conn_prob_idx(conn_data)
num_exp = 1000
t11 = [np.cov(shuffle_glom_kc_w_prob(ob_conn, glom_prob),rowvar=False) for i in range(num_exp)]
t12 = np.stack(t11,2)

ob_cov = np.cov(ob_conn,rowvar=False)
idx = np.triu_indices(len(ob_cov))
results = np.zeros(ob_cov.shape)
for i,j1 in enumerate(idx[0]):
    j2 = idx[1][i]
    t1 = t12[j1,j2,:]
    j3 = ob_cov[j1,j2]
    p = sum(t1>j3)/num_exp
    results[j1,j2]=p
    results[j2,j1]=p

pvs = results.copy()
pvs[pvs>=0.05]=1
pvs = 1-pvs
pvs[np.where(pvs)]=1
# re-cluster the P binary matrix
cm_cv = PairMatrix('', pvs.copy(), glom_idx_ids)

# reorder_idx = km_cluster(pvs)
reorder_idx = reorder(covar_order_ids, glom_idx_ids)

t1_cv = cm_cv.reorder(reorder_idx, return_new=True)

fig, ax1 = plt.subplots()
t1 = t1_cv;
gloms = df_lookup('glom_id',t1.col_ids,'short_glom_name',glom_btn_table)
sns.heatmap(t1.conn, xticklabels=gloms, yticklabels=gloms, ax=ax1, cmap=cm.get_cmap('viridis', 2))

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
fig.set_size_inches(18,12)
plt.show()
fig.savefig(save_path + "191202-allKC_cov_BinaryPval_ClusterP_over-rep_wGlomAnno.png", bbox_inches='tight')


# use original cluster order
cm_cv = PairMatrix('', pvs.copy(), glom_idx_ids)
# reorder_idx = km_cluster(cm_zs.conn)
reorder_idx = reorder(ClusterOrder0707, glom_idx_ids)

t1_cv = cm_cv.reorder(reorder_idx, return_new=True)

fig, ax1 = plt.subplots()
t1 = t1_cv;
gloms = df_lookup('glom_id',t1.col_ids,'short_glom_name',glom_btn_table)

sns.heatmap(t1.conn, xticklabels=gloms, yticklabels=gloms, ax=ax1, cmap=cm.get_cmap('viridis', 2))

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
fig.set_size_inches(18,12)
plt.show()
fig.savefig(save_path + "191202-allKC_cov_BinaryPval_OriginalOrder_over-rep_wGlomAnno.png", bbox_inches='tight')


# calculate covariance matrices with synapse numbers
def shuffle_glom_kc_w_prob_syn(gk_conn, col_prob, syn_d):
    '''Given a glomerulus-KC connectivity matrix, shuffle the connection
    while maintain the numbers CLAWS ONLY and return the shuffled matrix.
    Note that ndividual claw connection is not identifiable but as numbers in the glom-KC cell
    e.g. 2 in a cell means the KC and glom connects with 2 claws
    This one with probability of choice for each glomerulus (eacg column)

    191112 add synapse distribution to generate a shuffled connectivity matrix'''
    sfl_conn = np.zeros(gk_conn.shape)
    num_col = sfl_conn.shape[1]
    for i in range(sfl_conn.shape[0]):
        t1 = np.random.choice(int(num_col), size=int(sum(gk_conn[i,:])), p=col_prob)
        for j in t1:
            sfl_conn[i, j] += np.random.choice(syn_d)
    return sfl_conn


# get the distribution of all synapses from random draw KCs
syn_conn = ana_rd.conn_data['bouton_claw'].conn['5s']
syn_d = syn_conn[np.where(syn_conn)]

conn_data = ana.conn_data['glom_kc_in_claw_unit']
ob_conn, glom_prob, glom_idx_ids = get_conn_prob_idx(conn_data)

num_exp = 1000

t11 = [np.cov(shuffle_glom_kc_w_prob_syn(ob_conn, glom_prob, syn_d),rowvar=False) for i in range(num_exp)]
t12 = np.stack(t11,2)

# get the observed covariance matrix of synapse level connections
conn = ana.conn_data['glom_kc_contracted'].conn['5s']
ob_cov = np.cov(conn, rowvar=False)

idx = np.triu_indices(len(ob_cov))
results = np.zeros(ob_cov.shape)
for i,j1 in enumerate(idx[0]):
    j2 = idx[1][i]
    t1 = t12[j1,j2,:]
    j3 = ob_cov[j1,j2]
    p = sum(t1>j3)/num_exp
    results[j1,j2]=p
    results[j2,j1]=p

pvs = results.copy()
pvs[pvs>=0.05]=1
pvs = 1-pvs
pvs[np.where(pvs)]=1

cm_cv = PairMatrix('', pvs.copy(), glom_idx_ids)
# reorder_idx = km_cluster(1-pvs)
reorder_idx = reorder(covar_order_syn_ids, glom_idx_ids)

t1_cv = cm_cv.reorder(reorder_idx, return_new=True)

fig, ax1 = plt.subplots()
t1 = t1_cv;
gloms = df_lookup('glom_id',t1.col_ids,'short_glom_name',glom_btn_table)

sns.heatmap(t1.conn, xticklabels=gloms, yticklabels=gloms, ax=ax1, cmap=cm.get_cmap('viridis', 2))

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
fig.set_size_inches(18,12)
plt.show()
fig.savefig(save_path + "191202-allKC_cov_BinaryPval_ClusterP_SynDistr_over-rep_wGlomAnno.png", bbox_inches='tight')



# same as above but use original clustering order
cm_cv = PairMatrix('', pvs.copy(), glom_idx_ids)
# reorder_idx = km_cluster(1-results)
reorder_idx = reorder(ClusterOrder0707, glom_idx_ids)
t1_cv = cm_cv.reorder(reorder_idx, return_new=True)
fig, ax1 = plt.subplots()
t1 = t1_cv;
gloms = df_lookup('glom_id',t1.col_ids,'short_glom_name',glom_btn_table)
sns.heatmap(t1.conn, xticklabels=gloms, yticklabels=gloms, ax=ax1, cmap=cm.get_cmap('viridis', 2))

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
fig.set_size_inches(18,12)
plt.show()
fig.savefig(save_path + "191202-allKC_cov_BinaryPval_OriginalOrder_SynDistr_over-rep_wGlomAnno.png", bbox_inches='tight')
# fig.savefig(save_path + "191112-allKC_cov_pval_OriginCluster_SynDistr_over-representation.png", bbox_inches='tight')

# below are old comments
#------------------------------------------------------------------
#copy from /Users/zhengz11/myscripts/bocklab_git/bocklab/zhihao/mushroom_py/v10/191202-replot_covariance_matrix.py
# change date to 200110

# save_path = "/Users/zhengz11/myscripts/data_results/191112-covariance_matrices/191202_updated/"
