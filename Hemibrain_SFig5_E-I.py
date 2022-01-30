
date = "200911"

csv_tbl = pd.read_csv(wd + "tables/200628_pn_kc_conn", header=0, index_col=0, sep=" ")

exclude = []
for i in ['VP1l', 'VP1m', 'VP4', 'VP2', 'VP3']:
    for idx,j in enumerate(csv_tbl.index.values):
        if i in j:
            exclude.append(idx)

hemi_tbl = csv_tbl.drop(csv_tbl.index[exclude])

hemi_conn = np.transpose(np.array(hemi_tbl))

hemi_pn_names = hemi_tbl.index.values
hemi_short_glom_names = [i.split(" ")[0] for i in hemi_tbl.index.values]

hemi_id_map = pd.DataFrame({'PN': hemi_pn_names, 'glom': hemi_short_glom_names})

hemi_pn_kc = ConnectivityMatrix('hemi_pn_kc_conn', hemi_conn, hemi_pn_names, hemi_tbl.columns.values)

hemi_glom_kc = hemi_pn_kc.combine_new(hemi_id_map, output_name="hemi_glom_kc_conn", axis=1)


conn = hemi_pn_kc.conn['1s'].copy()
conn[conn<3]=0
conn[conn>0]=1
hemi_pn_kc_cu = ConnectivityMatrix('hemi_pn_kc_conn_in_claw_units', conn, hemi_pn_names, hemi_tbl.columns.values)
hemi_glom_kc_cu = hemi_pn_kc_cu.combine_new(hemi_id_map, output_name="hemi_glom_kc_conn_in_claw_units", axis=1)

glom_ids_hemi = df_lookup(
        'short_glom_name', hemi_glom_kc_cu.col_ids, 'glom_id', glom_btn_table)

    # get glom_prob
glom_prob_hemi = df_lookup('short_glom_name', hemi_glom_kc_cu.col_ids,
                          'norm_bouton_count', glom_btn_table)




ob_conn = hemi_glom_kc_cu.conn['1s'].copy()

glom_prob = glom_prob_hemi / np.sum(glom_prob_hemi)
glom_idx_ids = glom_ids_hemi


# random bouton model
stat = [get_raw_inputs(shuffle_glom_kc_w_prob(ob_conn, glom_prob)) for i in range(10000)]
stat = np.array(stat)
sd = np.nanstd(stat, axis=0)
avg = np.nanmean(stat, axis=0)
ob_ci = get_raw_inputs(ob_conn)
comm_zscore = np.divide(np.subtract(ob_ci, avg), sd)
cm_zs = PairMatrix('', comm_zscore, glom_idx_ids)
conn = cm_zs.conn.copy()
gloms = df_lookup('glom_id',cm_zs.col_ids,'short_glom_name',glom_btn_table)


zs_t1 = comm_zscore.flatten().copy()

#--------------------------------------------------------
# random bouton model
# use original cluster used by PNKC paper
reorder_idx = reorder(ClusterOrder0707, glom_idx_ids)
t1_zs = cm_zs.reorder(reorder_idx, return_new=True)

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

# fig.savefig(save_path + date + "-HemibrainConn_RandomBoutonModel_FAFBCluster_10kRandomizations.png", bbox_inches='tight')



#--------------------------------------------------------
# re-cluster k-means with the random bouton model

num_clusters = 4
fi = KMeans(n_clusters=num_clusters, n_init=100).fit(cm_zs.conn.copy())

glom_seq = df_lookup('glom_id',cm_zs.col_ids,'short_glom_name',glom_btn_table)

reorder_idx = np.argsort(fi.labels_)
t1_zs = cm_zs.reorder(reorder_idx, return_new=True)

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
# fig.savefig(save_path + date + "-HemibrainConn_RandomBoutonModel_ReCluster_OptimalNumCluster{}_10kRd.png".format(num_clusters), bbox_inches='tight')


##-----------------------------------------------
# random claw null model
stat = [get_raw_inputs(i) for i in shuffle_glom_kc_iterate(ob_conn, 10000)]
stat = np.array(stat)
sd = np.nanstd(stat, axis=0)
avg = np.nanmean(stat, axis=0)
ob_ci = get_raw_inputs(ob_conn)
comm_zscore = np.divide(np.subtract(ob_ci, avg), sd)
cm_zs = PairMatrix('', comm_zscore, glom_idx_ids)
conn = cm_zs.conn.copy()
gloms = df_lookup('glom_id',cm_zs.col_ids,'short_glom_name',glom_btn_table)


t1 = df_lookup('short_glom_name', gloms,'glom_id',glom_btn_table)


# show the conditional input matrix with original cluster order
cm_zs = PairMatrix('', comm_zscore, t1.copy())

zs_t2 = comm_zscore.flatten().copy()

conn = cm_zs.conn.copy()

reorder_idx = reorder(ClusterOrder0707, glom_idx_ids)
t1_zs = cm_zs.reorder(reorder_idx, return_new=True)



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
# fig.savefig(save_path + date + "-HemibrainConn_RandomClawModel_FAFBClusterOrder_10kRd.png", bbox_inches='tight')




# random claw model
# kmeans recluster with optimal number of cluster 4
num_clusters = 4

fi = KMeans(n_clusters=num_clusters, n_init=100).fit(cm_zs.conn.copy())

glom_seq = df_lookup('glom_id',cm_zs.col_ids,'short_glom_name',glom_btn_table)

# save cluster labels
pd.DataFrame({'glom':glom_seq, 'cluster_labels': fi.labels_}).to_csv(save_path + date + "-HemibrainConn_RandomClawModel_10kRd_cluster_labels.csv")

reorder_idx = np.argsort(fi.labels_)

t1_zs = cm_zs.reorder(reorder_idx, return_new=True)

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
# fig.savefig(save_path + date + "-HemibrainConn_RandomClawModel_ReCluster_OptimalNumCluster{}_10kRd.png".format(num_clusters), bbox_inches='tight')


# 210524
# compare distributions of random bouton model vs. random claw model

zs_t1 = comm_zscore.flatten().copy()

#(2)
zs_t2 = comm_zscore.flatten().copy()


# adjust fig size and font size
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 13})

fig, ax = plt.subplots()
bin_set = np.arange(-10,14,0.5)
xtick_labels = np.arange(-10,14,1)
ax.hist([zs_t1, zs_t2], bins=bin_set, label=['observed vs. random bouton', 'observed vs. random claw'])
ax.legend(loc='upper right')
ax.set_xticks(xtick_labels)
ax.set_xlabel('Z scores',fontdict={'fontsize': 13})
plt.xlim(-10,13)
fname = '210524-zscore_hist_hemibrain_ObvsRandBtn_ObvsRandClaw_adj.png'
# fname = '210621-zscore_hist_hemibrain_ObvsRandBtn_ObvsRandClaw_adj.png'
ax.set_title(fname, fontdict={'fontsize': 10}, pad=15)
fig.set_size_inches(12,10)

fig.savefig(save_path + fname, bbox_inches='tight')

'''
zs_t1 mean: 0.2243
zs_t1 std: 2.7555
zs_t2 mean: 0.1776
zs_t2 std: 1.7695
'''
