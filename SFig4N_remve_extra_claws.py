

ana = ana_all_rd

conn_data = ana.conn_data['glom_kc_in_claw_unit']
ob_conn, glom_prob, glom_idx_ids = get_conn_prob_idx(conn_data)

row_ids = ana.conn_data['glom_kc_in_claw_unit'].row_ids

all_claws = [len(ana.row_neurons[i].segments.ids) for i in row_ids]

conn_fill = np.concatenate([ob_conn,np.expand_dims(np.zeros(ob_conn.shape[0]),1)],1)
for i,v in enumerate(np.sum(ob_conn,1)):
    if v < all_claws[i]:
        conn_fill[i,-1] = all_claws[i] - v

idx = find_elements(glom_idx_ids, comm_ids)
comm_conn = conn_fill[:,idx]

conn_c2 = conn_fill[np.where(np.sum(comm_conn,1)>=2)[0],:]
conn_c1 = conn_fill[np.where(np.sum(comm_conn,1)<2)[0],:]

c2_avg = np.mean(np.sum(conn_c2,1))
c1_avg = np.mean(np.sum(conn_c1,1))

conn_trim = conn_c2.copy()
while np.mean(np.sum(conn_trim,1)) >  c1_avg:
    nz = np.where(conn_trim)
    t1 = np.random.randint(len(nz[0])-1)
    conn_trim[nz[0][t1],nz[1][t1]] -= 1

conn_rt = np.concatenate([conn_trim,conn_c1])[:,:-1]

conn = conn_rt.copy()

stat = [get_raw_inputs(shuffle_glom_kc_w_prob(conn, glom_prob)) for i in range(num_exp)]
stat = np.array(stat)

sd = np.nanstd(stat, axis=0)
avg = np.nanmean(stat, axis=0)

# avg = np.load(save_path + "200712-10k_allKCs_avg.npy")
# sd = np.load(save_path + "200712-10k_allKCs_sd.npy")

ob_ci = get_raw_inputs(conn)
comm_zscore = np.divide(np.subtract(ob_ci, avg), sd)

cm_zs = PairMatrix('', comm_zscore, glom_idx_ids)

num_clusters = 4
fi = KMeans(n_clusters=num_clusters, n_init=100).fit(cm_zs.conn.copy())
# order k-means group based on their average z-scores
avg_list = []
for i in range(num_clusters):
    ti = np.nonzero(fi.labels_ == i)[0]
    avg_list.append(np.mean(cm_zs.conn[ti,:][:,ti]))
reorder_idx = np.concatenate([np.nonzero(fi.labels_==i)[0] for i in np.flip(np.argsort(avg_list))])

# reorder_idx = reorder(ClusterOrder0707, glom_idx_ids)

t1_zs = cm_zs.reorder(reorder_idx, return_new=True)

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

fname = '211024-' + 'Remove_2commKCs_claws.png'
ax1.set_aspect("equal")
ax1.set_title(fname, pad = 50)
fig.set_size_inches(16,12)
plt.show()

# save_path = "/Users/zhengz11/myscripts/data_results/211022-KC_claws/"
fig.savefig(save_path + fname, bbox_inches='tight')
