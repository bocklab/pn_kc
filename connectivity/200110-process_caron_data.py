
# todo: need to add process_caron_data_v2

# 200110 copy from /Users/zhengz11/myscripts/bocklab_git/bocklab/zhihao/mushroom_py/v10/191202-process_caron_data.py
# this generates Suppl Fig4 A-D (re-analysis of Caron et al. result)
# change date to 200110

save_path = "/Users/zhengz11/myscripts/data_results/191202-redo_caron_comparison/"

from process_caron_data_v2 import *

num_exp = 1000

# make all gloms use this universal ids
glom_id_table = pd.read_excel(
    "/Users/zhengz11/myscripts/data_results/171012-1D_olfactory_space/171012-glom_index_list.xlsx")

def get_zscores(conn, prob, num_exp=1000):
    stat = [get_raw_inputs(shuffle_glom_kc_w_prob(conn, prob)) for i in range(num_exp)]
    stat = np.array(stat)
    sd = np.nanstd(stat, axis=0)
    avg = np.nanmean(stat, axis=0)
    ob_ci = get_raw_inputs(conn)
    comm_zscore = np.divide(np.subtract(ob_ci, avg), sd)
    return comm_zscore

tbl_wc = tbl.copy()
ana_ca = CaronAnalysis.init_processing()
caron_obj = ana_ca.conn_data['caron_glom_kc']
ca_conn = caron_obj.conn['1s'].copy()
caron_obj.col_ids = df_lookup(
    'glom', ana_ca.glom_tbl['glom_class'], 'id', glom_id_table)

# Combine these to a single column in the connectivity matrix and left it out from displaying and clustering
# ['DL6', 'VC3', 'VM6', 'cold', 'heat', 'other']

##--------------------------------------------------
# get the bouton probability for ca_conn
glom_list = glom_btn_table['glom_id'].tolist()
caron_prob = []
no_btn_list = []
for i in caron_obj.col_ids:
    if i not in glom_list:
        print(i)
        caron_prob.append(0)
        no_btn_list.append(i)
    else:
        j = df_lookup('glom_id', [i], 'norm_bouton_count', glom_btn_table)
        caron_prob.append(j[0])

# set the bouton probability of VC3 as the sum of VC3l and VC3m
vc3 = glom_btn_table.query(
    'short_glom_name=="VC3m" or short_glom_name=="VC3l"').norm_bouton_count.sum()

# must have a simpler way to do this?
# get VC3 id, get VC3 index from the id, then get
caron_prob[caron_obj.col_ids.index(
    glom_id_table.loc[glom_id_table.glom == 'VC3', 'id'].tolist()[0])] = vc3

# divide the remaining probability to the zero gloms
nnz = np.where(caron_prob)[0]
zeros = np.where(np.array(caron_prob)==0)[0]
ca_conn = np.concatenate([ca_conn[:,nnz],
                            np.sum(ca_conn[:,zeros],1).reshape(ca_conn.shape[0],1)],1)
caron_prob = [v for i,v in enumerate(caron_prob) if i in nnz]
caron_prob.append(1 - sum(caron_prob))
caron_col_ids = caron_obj.col_ids.copy()
caron_col_ids = [v for i,v in enumerate(caron_col_ids) if i in nnz]
caron_col_ids.append(100)

comm_zscore = get_zscores(ca_conn, caron_prob, num_exp=1000)

# ca_conn, caron_prob, caron_col_ids
# remove the final extra column

cm_zs = PairMatrix('', comm_zscore[:-1,:-1].copy(), caron_col_ids[:-1].copy())

reorder_idx = km_cluster(cm_zs.conn)
t1_zs = cm_zs.reorder(reorder_idx, return_new=True)

t1 = t1_zs
fig, ax1 = plt.subplots()
col_ids = t1.col_ids.copy()
t55 = col_ids.index(55)

gloms = df_lookup('id', t1.col_ids, 'glom', glom_id_table)
sns.heatmap(t1.conn, xticklabels=gloms, yticklabels=gloms, ax=ax1, vmin=-8.53, vmax=8.53, cmap="RdBu_r")

ax1.tick_params(bottom=False,labeltop=True, top=True, labelbottom=False)
ax1.tick_params(axis='x',labelrotation=90)

# run and get tbl in 191029-bouton-KC-representations_per_PN.py
col_list = t1.col_ids
col_colors = []
existing_gloms = pd.unique(tbl.short_glom_name)
for glom in gloms:
    if glom in existing_gloms:
        col_colors.append(df_lookup('short_glom_name', [glom], 'color', tbl)[0])
    else:
        col_colors.append('black')

for x in [ax1.get_xticklabels(), ax1.get_yticklabels()]:
    for idx, tick in enumerate(x):
        tick.set_color(col_colors[idx])
        if col_list[idx] in comm_ids:
            tick.set_weight("extra bold")

ax1.set_aspect("equal")
fig.set_size_inches(16,12)
plt.show()

fname = '191204-Analyze_caron_conn_BtnNullModel'
fig.savefig(save_path + fname + '.png', bbox_inches='tight')
pd.DataFrame({'glom_seq_names': gloms, 'glom_seq_ids': t1.col_ids}).to_csv(save_path + fname + '_GlomSeq'+ '.csv')


##-----------------------------------------------
# run Caron et al. connectivity with claw null models
t1 = np.sum(ca_conn,0)
ob_prob = t1/np.sum(t1)
comm_zscore = get_zscores(ca_conn, ob_prob, num_exp=1000)

# remove the final extra column

cm_zs = PairMatrix('', comm_zscore[:-1,:-1].copy(), caron_col_ids[:-1].copy())

reorder_idx = km_cluster(cm_zs.conn)
t1_zs = cm_zs.reorder(reorder_idx, return_new=True)

t1 = t1_zs
fig, ax1 = plt.subplots()
col_ids = t1.col_ids.copy()
t55 = col_ids.index(55)

gloms = df_lookup('id', t1.col_ids, 'glom', glom_id_table)
sns.heatmap(t1.conn, xticklabels=gloms, yticklabels=gloms, ax=ax1, vmin=-8.53, vmax=8.53, cmap="RdBu_r")

ax1.tick_params(bottom=False,labeltop=True, top=True, labelbottom=False)
ax1.tick_params(axis='x',labelrotation=90)

# run and get tbl in 191029-bouton-KC-representations_per_PN.py
col_list = t1.col_ids
col_colors = []
existing_gloms = pd.unique(tbl.short_glom_name)
for glom in gloms:
    if glom in existing_gloms:
        col_colors.append(df_lookup('short_glom_name', [glom], 'color', tbl)[0])
    else:
        col_colors.append('black')

for x in [ax1.get_xticklabels(), ax1.get_yticklabels()]:
    for idx, tick in enumerate(x):
        tick.set_color(col_colors[idx])
        if col_list[idx] in comm_ids:
            tick.set_weight("extra bold")

ax1.set_aspect("equal")
fig.set_size_inches(16,12)
plt.show()

fname = '191204-Analyze_caron_conn_ClawNullModel'
fig.savefig(save_path + fname + '.png', bbox_inches='tight')
pd.DataFrame({'glom_seq_names': gloms, 'glom_seq_ids': t1.col_ids}).to_csv(save_path + fname + '_GlomSeq'+ '.csv')



# caron_obj.col_ids
# ca_conn, ca_prob

##-----------------------------------------------
# what is the claw distribution of the subsampled conn from our observation
# sbs_conn
# sbs_prob
# sbs_col_ids
# ob_sbs_claws

ana = ana_all_rd
ana_ca = CaronAnalysis.init_processing()
caron_obj = ana_ca.conn_data['caron_glom_kc']
ca_conn = caron_obj.conn['1s'].copy()
t1 = ca_conn.sum(0)

ca_tbl = ana_ca.suppl_tbl.copy()

# in ana_all_rd, How many KCs have known subtypes? What are their subtypes?
y = cc.get_skids_from_annos(fafb_c,[['KCy']])

prime = cc.get_skids_from_annos(fafb_c,
                             [["KCa'B'", "KCa'B'ap", "KCa'B'm", "KCa'B'x"]])

ab = cc.get_skids_from_annos(fafb_c,
                             [['KCaBc', 'KCaBs', 'KCaBx']])

kcs = ana.conn_data['glom_kc_in_claw_unit'].row_ids
conn_obj = ana.conn_data['glom_kc_in_claw_unit']
kci = [i for i,v in enumerate(kcs) if v in y + prime + ab]

kci_row_ids = [kcs[i] for i in kci]
kci_conn = conn_obj.conn['1s'][kci,:]
kci_claws = np.sum(kci_conn,1)
kci_types = []
for i in kci_row_ids:
    if i in y:
        m = 'y'
    elif i in ab:
        m = 'ab'
    else:
        m = 'prime'
    kci_types.append(m)

kci_tbl = pd.DataFrame({'kci_row_ids':kci_row_ids,
                        'kci_claws': kci_claws,
                        'kci_types': kci_types})
kci_tbl['index']=kci_tbl.index
type_dict = dict(zip(["alpha'/beta'", 'alpha/beta', 'gamma'],
["prime", 'ab', 'y']))

index_list = []
claw_list = []
for t in pd.unique(ca_tbl.KC_subclass):
    t1 = type_dict[t]
    claws = ca_tbl.query('KC_subclass == @t')['claws_filled'].tolist()
    sub_tbl = kci_tbl.query('kci_types == @t1')
    perm = np.random.permutation(sub_tbl.index.values)
    t6 = []
    t8 = []
    ct = 0
    kci_claws = sub_tbl['kci_claws']
    kci_index = sub_tbl.index
    while len(t6)<len(claws):
        if kci_claws.iloc[ct] > claws[len(t6)]:
            t6.append(claws[len(t6)])
            t8.append(kci_index[ct])
        print(ct)
        ct += 1
    claw_list.extend(t6)
    index_list.extend(t8)

picked_tbl = pd.DataFrame({'index': index_list, 'caron_claws': claw_list}).merge(kci_tbl, on='index', how='left').assign(delta_claws=lambda x: x.kci_claws - x.caron_claws)


rows = [conn_obj.row_ids.index(i) for i in picked_tbl.kci_row_ids]

sbs_conn = conn_obj.conn['1s'][rows,:]

for j in range(sbs_conn.shape[0]):
    t5 = []
    for i in np.where(sbs_conn[j,:])[0]:
        t5.extend(int(sbs_conn[j,i])*[i])
    t1 = np.random.choice(t5, size=int(picked_tbl['delta_claws'].iloc[j]), replace=False)
    for l in t1:
        sbs_conn[j,l] -= 1

sbs_col_ids = conn_obj.col_ids
ob_sbs_claws = np.sum(sbs_conn,0)
sbs_prob = df_lookup('glom_anno_id', sbs_col_ids, 'norm_bouton_count', glom_btn_table)

comm_zscore = get_zscores(sbs_conn, sbs_prob)
cm_zs = PairMatrix('', comm_zscore, sbs_col_ids)

reorder_idx = km_cluster(cm_zs.conn)
t1_zs = cm_zs.reorder(reorder_idx, return_new=True)

t1 = t1_zs;

fig, ax1 = plt.subplots()

# gloms is for annotating the axes
gloms = df_lookup('glom_anno_id',t1.col_ids,'short_glom_name',glom_btn_table)

# col_list is for knowing which one is community PNs
col_list = df_lookup('glom_anno_id', t1.col_ids,'glom_id',glom_btn_table)

# for assigning colors to the labels in each axis
col_colors = df_lookup('short_glom_name', gloms, 'color', tbl)

sns.heatmap(t1.conn, xticklabels=gloms, yticklabels=gloms, ax=ax1, vmin=-8.53, vmax=8.53, cmap="RdBu_r")
ax1.tick_params(bottom=False,labeltop=True, top=True, labelbottom=False)
ax1.tick_params(axis='x',labelrotation=90)

for x in [ax1.get_xticklabels(), ax1.get_yticklabels()]:
    for idx, tick in enumerate(x):
        tick.set_color(col_colors[idx])
        if col_list[idx] in comm_ids:
            tick.set_weight("extra bold")

ax1.set_aspect("equal")
fig.set_size_inches(16,12)
plt.show()

fname = '191204-Analyze_subsampled_conn_BtnNullModel'
fig.savefig(save_path + fname + '.png', bbox_inches='tight')
pd.DataFrame({'glom_seq_names': gloms, 'glom_seq_ids': t1.col_ids}).to_csv(save_path + fname + '_GlomSeq'+ '.csv')



##---------------------------------------------------
# analyze subsampled connectivity with claw null models
t1 = np.sum(sbs_conn,0)
ob_prob = t1 / np.sum(t1)
comm_zscore = get_zscores(sbs_conn, ob_prob)
cm_zs = PairMatrix('', comm_zscore, sbs_col_ids)

reorder_idx = km_cluster(cm_zs.conn)
t1_zs = cm_zs.reorder(reorder_idx, return_new=True)

t1 = t1_zs;

fig, ax1 = plt.subplots()

# gloms is for annotating the axes
gloms = df_lookup('glom_anno_id',t1.col_ids,'short_glom_name',glom_btn_table)

# col_list is for knowing which one is community PNs
col_list = df_lookup('glom_anno_id', t1.col_ids,'glom_id',glom_btn_table)

# for assigning colors to the labels in each axis
col_colors = df_lookup('short_glom_name', gloms, 'color', tbl)

sns.heatmap(t1.conn, xticklabels=gloms, yticklabels=gloms, ax=ax1, vmin=-8.53, vmax=8.53, cmap="RdBu_r")
ax1.tick_params(bottom=False,labeltop=True, top=True, labelbottom=False)
ax1.tick_params(axis='x',labelrotation=90)

for x in [ax1.get_xticklabels(), ax1.get_yticklabels()]:
    for idx, tick in enumerate(x):
        tick.set_color(col_colors[idx])
        if col_list[idx] in comm_ids:
            tick.set_weight("extra bold")

ax1.set_aspect("equal")
fig.set_size_inches(16,12)
plt.show()

fname = '191204-Analyze_subsampled_conn_ClawNullModel'
fig.savefig(save_path + fname + '.png', bbox_inches='tight')
pd.DataFrame({'glom_seq_names': gloms, 'glom_seq_ids': t1.col_ids}).to_csv(save_path + fname + '_GlomSeq'+ '.csv')

pd.DataFrame(sbs_conn).to_csv(save_path + fname + 'sbs_conn.csv')
pd.DataFrame(sbs_col_ids).to_csv(save_path + fname + 'sbs_conn_col_ids.csv')



##---------------------------------------------
# ca_conn, caron_prob, caron_col_ids, ca_ob
stat = [np.sum(shuffle_glom_kc_w_prob(ca_conn, caron_prob),0) for i in range(num_exp)]
stat = np.array(stat)
sd = np.nanstd(stat, axis=0)
avg = np.nanmean(stat, axis=0)

ca_ob = np.sum(ca_conn,0)

glom_seq = pd.read_csv(save_path + '191204-Analyze_caron_conn_BtnNullModel_GlomSeq.csv').glom_seq_names.tolist()

gloms = df_lookup('id', caron_col_ids[:-1], 'glom', glom_id_table)

ro_index = [gloms.index(i) for i in glom_seq]
gloms_ro = [gloms[j] for j in ro_index]
ca_ob_ro = ca_ob[ro_index]
avg_ro = avg[ro_index]
sd_ro = sd[ro_index]
# glom_seq

fig, ax = plt.subplots()
ind = np.arange(len(glom_seq))
width = 0.4

r1 = ax.bar(ind + 0.2, avg_ro, width, color=sns.xkcd_palette(['windows blue']), align='center')

r2 = ax.bar(ind + 0.6, ca_ob_ro, width, color=sns.xkcd_palette(['amber']), align='center')

r3 = ax.errorbar(ind + 0.2, avg_ro, yerr=sd_ro, fmt = 'none', color=sns.xkcd_palette(['windows blue']))

ax.legend((r1[0], r2[0]), ['bouton null model', 'observed in Caron et al.'], fontsize=12, loc='upper right')

plt.xticks(ind + 0.4, glom_seq)
ax.tick_params(axis='both', which='major', labelsize=11)
plt.xlim(-0.2,49)
fig.set_size_inches(30,8)

fname = '191204-caron_conn_sampling_tick'
fig.savefig(save_path + fname + '.png', bbox_inches='tight')
