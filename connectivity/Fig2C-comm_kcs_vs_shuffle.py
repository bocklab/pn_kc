
#
# copy from  part of 191117-comm_claws-kcs_vs_shuffle.py
# Fig 2C (200326, PNKC2019_v9_fig_200313DB-ZZfixedSuppl6B.pptx)
# save_path = '/Users/zhengz11/myscripts/data_results/191117-community_expansion/'
ana = ana_all_rd

num_exp = 1000

comm_anno_ids = df_lookup('glom_id',comm_ids,'glom_anno_id',glom_btn_table)

ana_obj = ana.conn_data['glom_kc_in_claw_unit']

all_ids = ana_obj.col_ids

comm_idx = find_elements(all_ids,comm_anno_ids)

reorder_glom = comm_idx + [i for i in range(len(all_ids)) if i not in comm_idx]

ob_comm = ana_obj.conn['1s'][:,comm_idx]
ob_claw_set = np.sum(ob_comm)

ob_conn, glom_prob, glom_idx_ids = get_conn_prob_idx(ana_obj)

ro_prob = [glom_prob[i] for i in reorder_glom]

ro_conn = ob_conn[:,reorder_glom]
ob_claw = np.sum(ro_conn[:,:10])
ob_kc = np.count_nonzero(np.sum(ro_conn[:,:10],1))

sfl_claw_set = []
sfl_kc_set = []
for i in range(num_exp):
    sfl_conn_g = shuffle_glom_kc_w_prob(ro_conn, ro_prob)[:,:10]
    sfl_claw_set.append(np.sum(sfl_conn_g))
    sfl_kc_set.append(np.count_nonzero(np.sum(sfl_conn_g,1)))

plt.style.use('seaborn-deep')
fig, ax1 = plt.subplots()
clawhist = ax1.hist(sfl_claw_set)
ax1.plot(ob_claw, 4, 'ro', ms=12)
fig.set_size_inches([10,6])
# plt.savefig(save_path + '191117-num_claws_for_comm_vs_noncomm_allKCs.png')
# mean - 1412.5
# std - 34.0
# observed - 1901
# 14.3 std
# (ob_claw - np.mean(sfl_claw_set))/np.std(sfl_claw_set)

fig = plt.figure()
plt.hist(sfl_kc_set)
plt.plot(ob_kc, 4, 'ro', ms=12)
fig.set_size_inches([10,6])
# plt.savefig(save_path + '191117-num_kcs_for_comm_vs_noncomm_allKCs.png')
# mean - 903.4
# std - 17.3
# observed - 844
# 3.4 std

# 200226 re-plot the above figure with the same data
fig, ax = plt.subplots()
ax.hist(sfl_kc_set, bins=range(840,970,8), edgecolor='black')
ax.plot(ob_kc, 4, 'ro', ms=12)
fig.set_size_inches([12,6])
# save_path = '/Users/zhengz11/myscripts/data_results/200226-rerun_comm_claw_hist/'
# fig.savefig(save_path + '200226-num_kcs_for_comm_vs_noncomm_allKCs_morebins.png', bbox_inches='tight', dpi=600)
