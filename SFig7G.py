
# codes of various models are copied from:
# 210218-RandomGlomBoutonClaws_LocalRandom_rerun_zscores.py

num_exp = 10000
ana = ana_all_rd
conn_data = ana.conn_data['glom_kc_in_claw_unit']
ob_conn, glom_prob, glom_idx_ids = get_conn_prob_idx(conn_data)
comm_gloms = df_lookup("glom_id", comm_ids, "glom_anno_id", glom_btn_table)
idx = find_elements(glom_idx_ids,comm_ids)

def get_fraction(conn, comm_idx):
    comm_conn = conn[:,comm_idx]
    target_kcs=np.where(np.sum(comm_conn,1)>=2)[0]
    return np.sum(conn[target_kcs,:][:,comm_idx])/np.sum(conn[target_kcs,:])

observed = get_fraction(ob_conn, idx)
# 0.554254219772649


## random glom
#--------------------------------------------------------
RandomGlom = []
n = len(glom_prob)
prob = [1/n]*n
for i in range(num_exp):
    conn = shuffle_glom_kc_w_prob(ob_conn, prob)
    RandomGlom.append(get_fraction(conn, idx))



## random bouton
#--------------------------------------------------------
RandomBouton = []
for i in range(num_exp):
    conn = shuffle_glom_kc_w_prob(ob_conn, glom_prob)
    RandomBouton.append(get_fraction(conn, idx))


## random claws
#--------------------------------------------------------
RandomClaw = []
for i in range(num_exp):
    conn = shuffle_glom_kc_iterate(ob_conn, 1)[0]
    RandomClaw.append(get_fraction(conn, idx))



## local random
#--------------------------------------------------------
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  this part is copied from 200212-Ob_vs_LocalBtnRdModel.py
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
n=5
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
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

LocalRandom = []
for i in range(num_exp):
    conn = shuffle_glom_kc_w_nearest5(ob_conn, conn_data.row_ids, rows, glom_to_pick)
    LocalRandom.append(get_fraction(conn, idx))

'''
np.save(save_path + 'RandomGlom_10k.npy', np.array(RandomGlom))
np.save(save_path + 'RandomBouton_10k.npy', np.array(RandomBouton))
np.save(save_path + 'RandomClaw_10k.npy', np.array(RandomClaw))
np.save(save_path + 'LocalRandom_10k.npy', np.array(LocalRandom))
'''


means = [observed]
sds = [0]
for i in [RandomGlom, RandomBouton, RandomClaw, LocalRandom]:
    means.append(np.mean(i))
    sds.append(np.std(i))

labels = ["observed", "RandomGlom", "RandomBouton", "RandomClaw", "LocalRandom"]

fname = "211031-2orMoreCommKC_ClawFraction.png"

fig, ax = plt.subplots()
ax.bar(labels, means, yerr=sds, align='center', color=sns.color_palette("tab10")[0], ecolor='black', capsize=6)
ax.set_ylabel('Fraction of core community claws', fontsize=16)
ax.set_xticklabels(labels, rotation = 45, ha="right")
ax.set_title(fname, pad=18)
ax.set_ylim(0.4, 0.5)
ax.tick_params(axis="both", labelsize=16)
plt.yticks(np.arange(0.4,0.6,0.04),[i/100 for i in range(40,60,4)])
# Save the figure and show
plt.show()
fig.set_size_inches(8,6)
# fig.savefig(save_path + fname + ".png", bbox_inches="tight")
