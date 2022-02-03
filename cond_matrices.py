
##----------------------------
# Fig3B
# Ob vs. Random bouton
ana = ana_all_rd

conn_data = ana.conn_data['glom_kc_in_claw_unit']
ob_conn, glom_prob, glom_idx_ids = get_conn_prob_idx(conn_data)

num_exp = 10000

stat = [get_raw_inputs(shuffle_glom_kc_w_prob(ob_conn, glom_prob)) for i in range(num_exp)]
stat = np.array(stat)
sd = np.nanstd(stat, axis=0)
avg = np.nanmean(stat, axis=0)

ob_ci = get_raw_inputs(ob_conn)
comm_zscore = np.divide(np.subtract(ob_ci, avg), sd)

cm_zs = PairMatrix('', comm_zscore, glom_idx_ids)

# reorder_idx = km_cluster(cm_zs.conn)
reorder_idx = reorder(ClusterOrder0707, glom_idx_ids)
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
#            tick.set_bbox(dict(ec='green', fc=None, alpha=0.05))

ax1.set_aspect("equal")
fig.set_size_inches(16,12)
plt.show()


##----------------------------
# Fig4B
# observed vs. random claw
conn_data = ana.conn_data['glom_kc_in_claw_unit']

ob_conn, glom_prob, glom_idx_ids = get_conn_prob_idx(conn_data)

stat = [get_raw_inputs(i) for i in shuffle_glom_kc_iterate(ob_conn, 10000)]

stat = np.array(stat)
sd = np.nanstd(stat, axis=0)
avg = np.nanmean(stat, axis=0)

ob_ci = get_raw_inputs(ob_conn)

comm_zscore = np.divide(np.subtract(ob_ci, avg), sd)

cm_zs = PairMatrix('', comm_zscore, glom_idx_ids)

gloms = df_lookup('glom_id', cm_zs.col_ids,'short_glom_name',glom_btn_table)
# np.savetxt(save_path + date + "-PNKC_10kRD_comm_zscore_ObvsRandomClaw.csv", comm_zscore, '%1.3f', delimiter=',')
# np.savetxt(save_path + date + "-PNKC_10kRD_comm_zscore_RandomClaw_ids.csv", gloms, '%s', delimiter=',')

t1 = df_lookup('short_glom_name', gloms,'glom_id',glom_btn_table)

# clustering
cm_zs = PairMatrix('', comm_zscore, glom_idx_ids)

k = 4

num_clusters = k
fi = KMeans(n_clusters=num_clusters, n_init=1000).fit(cm_zs.conn.copy())
# order k-means group based on their average z-scores
avg_list = []
for i in range(num_clusters):
    ti = np.nonzero(fi.labels_ == i)[0]
    avg_list.append(np.mean(cm_zs.conn[ti,:][:,ti]))
reorder_idx = np.concatenate([np.nonzero(fi.labels_==i)[0] for i in np.flip(np.argsort(avg_list))])

# reorder_idx = reorder(ClusterOrder0707, glom_idx_ids)
t1_zs = cm_zs.reorder(reorder_idx, return_new=True)

fname = '210922-ob_vs_randomclaw_recluster_k{}.png'.format(num_clusters)

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
ax1.set_title(fname, pad = 50)
fig.set_size_inches(16,12)
plt.show()
# fig.savefig(save_path + fname, bbox_inches='tight')


##----------------------------
# Fig5F observed vs. local random model
n = 5
num_exp = 10000
ana = ana_all_rd
geom_conn = ana.geo_data['centroids'].conn

row_ids = segid_to_skid(ana.geo_data['centroids'].row_ids)
col_ids = segid_to_skid(ana.geo_data['centroids'].col_ids)
gloms = [ana.pn_mapping.skids_to_types(i)[0] for i in col_ids]

conn_data = ana.conn_data['glom_kc_in_claw_unit']
ob_conn, glom_prob, glom_idx_ids = get_conn_prob_idx(conn_data)
sfl_cols = conn_data.col_ids

t21 = ana.geo_data['centroids'].row_ids
t22 = ana.conn_data['bouton_claw'].row_ids
t21 == t22
conn_t = ana.conn_data['bouton_claw'].conn['5s']

select = np.where(np.sum(conn_t,1))[0]
row_ids = [row_ids[i] for i in select]
geom_conn = geom_conn[select,:]
conn_t = conn_t[select,:].copy()

ip_gloms = [sfl_cols.index(ana.pn_mapping.skids_to_types(col_ids[np.argmax(conn_t[i,:])])) for i in range(conn_t.shape[0])]

# importantly sfl_cols is from conn_data, and the glom_ids for the main matrix is from conn_data, so col_ids orders are the same, no need to re-order

t1 = []
for i,v in enumerate(row_ids):
    geom_sort = np.argsort(geom_conn[i,:])
    to_pick = [sfl_cols.index(gloms[j]) for j in geom_sort[:n]]
    t1.append(to_pick)
glom_to_pick = np.stack(t1)
rows = np.array(row_ids)



def shuffle_glom_kc_w_nearest5(gk_conn, row_ids, rows, glom_to_pick, n=5):
    # row_ids - skids for KCs in rows of gk_conn
    # rows - skids for each claw
    # glom_to_pick the 5 nearest boutons belong to what gloms for each claw
    # The index in glom_to_pick is column index for gk_conn
    sfl_conn = np.zeros(gk_conn.shape)
    for i,v in enumerate(row_ids):
        for j in np.where(rows==v)[0]:
            sfl_conn[i, glom_to_pick[j,np.random.randint(n)]] += 1
    return sfl_conn



sfl_conn_t = shuffle_glom_kc_w_nearest5(ob_conn, conn_data.row_ids, rows, glom_to_pick)

row_ids = conn_data.row_ids.copy()
stat = [get_raw_inputs(shuffle_glom_kc_w_nearest5(ob_conn, row_ids, rows, glom_to_pick)) for i in range(num_exp)]
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
