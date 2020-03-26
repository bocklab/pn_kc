
# copy from 191028-local_randomization_conditional_input_matrix.py
# 191119 number of community claws with all randomizations, local randomizations, and observed
# 200211---------------begin (local random bouton set)


num_exp = 1000
ana = ana_all_rd
geom_conn = ana.geo_data['centroids'].conn

row_ids = segid_to_skid(ana.geo_data['centroids'].row_ids)
col_ids = segid_to_skid(ana.geo_data['centroids'].col_ids)
gloms = [ana.pn_mapping.skids_to_types(i)[0] for i in col_ids]

conn_data = ana.conn_data['glom_kc_in_claw_unit']
ob_conn, glom_prob, glom_idx_ids = get_conn_prob_idx(conn_data)
sfl_cols = conn_data.col_ids


# check row_ids from geo_data are the same as row_ids from bouton_claw
# and trim it down to claws that receive inputs from unilgomerular PNs
t21 = ana.geo_data['centroids'].row_ids
t22 = ana.conn_data['bouton_claw'].row_ids
t21 == t22
conn_t = ana.conn_data['bouton_claw'].conn['5s']

select = np.where(np.sum(conn_t,1))[0]
row_ids = [row_ids[i] for i in select]
geom_conn = geom_conn[select,:]
conn_t = conn_t[select,:].copy()

ip_gloms = [sfl_cols.index(ana.pn_mapping.skids_to_types(col_ids[np.argmax(conn_t[i,:])])) for i in range(conn_t.shape[0])]


# glom_to_pick is an array; each row represents a claw, ordered as in the KC skeleton id (row_ids above); each row has 5 glom indices that represent 5 boutons closest to the claw
n = 5
including_connected = True

t1 = []

if including_connected:
    for i,v in enumerate(row_ids):
        geom_sort = np.argsort(geom_conn[i,:])
        to_pick = [sfl_cols.index(gloms[j]) for j in geom_sort[:n]]
        t1.append(to_pick)
else:
    for i,v in enumerate(row_ids):
        geom_sort = np.argsort(geom_conn[i,:])
        to_pick = [sfl_cols.index(gloms[j]) for j in geom_sort[:n+1]]
        t2 = ip_gloms[i]
        if t2 in to_pick:
            to_pick = [to_pick[l] for l in range(n+1) if l != to_pick.index(t2)]
        else:
            to_pick = to_pick[:n]
        t1.append(to_pick)

glom_to_pick = np.stack(t1)
rows = np.array(row_ids)
# 200211---------------end (local random bouton set)



comm_ids = [22, 44, 32, 42, 34, 38, 30, 23, 47, 4]
comm_anno_ids = df_lookup('glom_id',comm_ids,'glom_anno_id',glom_btn_table)
ana_obj = ana.conn_data['glom_kc_in_claw_unit']

comm_idx = find_elements(sfl_cols,comm_anno_ids)

local_claw_set = []
stat = []
for t11 in range(num_exp):
    sfl_conn = np.zeros(ob_conn.shape)
    for i,v in enumerate(conn_data.row_ids):
        for j in np.where(rows==v)[0]:
            sfl_conn[i, glom_to_pick[j,np.random.randint(n)]] += 1
    local_claw_set.append(np.sum(sfl_conn[:,comm_idx]))
    print(t11)


# global_claw_set
conn_data = ana.conn_data['glom_kc_in_claw_unit']
ob_conn, glom_prob, glom_idx_ids = get_conn_prob_idx(conn_data)

comm_prob = sum([glom_prob[i] for i in find_elements(glom_idx_ids, comm_ids)])

claw_num = np.sum(ob_conn,1)

global_claw_set = []
sfl_results = []
for i in range(num_exp):
    sfl_conn = get_random_conn(claw_num.shape[0], claw_num, prob=[comm_prob, 1-comm_prob])
    global_claw_set.append(np.sum(sfl_conn,0)[0])

# ob
ob_claw_set = np.sum(ob_conn[:,comm_idx])


# 200326, copy from 200212-Ob_vs_LocalBtnRdModel.py
# Fig3B (200326, PNKC2019_v9_fig_200313DB-ZZfixedSuppl6B.pptx)
fig, ax = plt.subplots()
cols = ['#1f77b4', '#2ca02c']
for i, x in enumerate([global_claw_set, local_claw_set]):
    ax.hist(x, bins=range(1300,2000,15), color = cols[i], edgecolor='black', linewidth=0.8)
ax.plot(ob_claw_set, 4, 'ro', ms=8)

fig.set_size_inches([12,6])
ax.legend(['observed','random bouton model','local random claw model'])
# save_path = '/Users/zhengz11/myscripts/data_results/200226-rerun_comm_claw_hist/'
# fig.savefig(save_path + '200226-comm_claws_hist_LocalvsGlobalvsObserved_dpi600.png', bbox_inches='tight', dpi=600)



# number of community claws in random bouton model vs. observed
# Fig2B (200326, PNKC2019_v9_fig_200313DB-ZZfixedSuppl6B.pptx)
fig, ax = plt.subplots()
ax.hist(global_claw_set, bins=range(1300,2000,15), color = '#1f77b4', edgecolor='black', linewidth=0.8)
ax.plot(ob_claw_set, 4, 'ro', ms=8)
# ax.set_ylim([0,300])

fig.set_size_inches([12,6])
fig.savefig(save_path + '200226-comm_claws_hist_GlobalvsObserved_dpi600.png', bbox_inches='tight', dpi=600)
