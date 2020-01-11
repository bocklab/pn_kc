
# copy from /Users/zhengz11/myscripts/bocklab_git/bocklab/zhihao/mushroom_py/v10/191117-comm_claws-kcs_vs_shuffle.py

# 200110 change date in name to 200110
# this will generate Fig 2A and Fig 2E (community input frequency between observed and claw-based null model)
# this will generate Fig 2B,C
# 191117-num_claws_for_comm_vs_noncomm_allKCs.png
# 191117-num_kcs_for_comm_vs_noncomm_allKCs.png

save_path = '/Users/zhengz11/myscripts/data_results/191117-community_expansion/'
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
plt.savefig(save_path + '191117-num_claws_for_comm_vs_noncomm_allKCs.png')
# mean - 1412.5
# std - 34.0
# observed - 1901
# 14.3 std
# (ob_claw - np.mean(sfl_claw_set))/np.std(sfl_claw_set)

fig = plt.figure()
plt.hist(sfl_kc_set)
plt.plot(ob_kc, 4, 'ro', ms=12)
fig.set_size_inches([10,6])
plt.savefig(save_path + '191117-num_kcs_for_comm_vs_noncomm_allKCs.png')
# mean - 903.4
# std - 17.3
# observed - 844
# 3.4 std


##--------------
# compare distribution of synapses per claw
# between community claws vs non-community claws
ana_obj = ana.conn_data['glom_claw_contracted']
all_ids = ana_obj.col_ids
omm_idx = find_elements(all_ids,comm_anno_ids)
reorder_glom = comm_idx + [i for i in range(len(all_ids)) if i not in comm_idx]
ro_conn = ana_obj.conn['1s'][:,reorder_glom]

comm_conn = ro_conn[:,:10].copy()
comm_set = comm_conn[np.where(comm_conn)]

noncomm_conn = ro_conn[:,10:].copy()
noncomm_set = noncomm_conn[np.where(noncomm_conn)]


##---------------------
# run to get tbl from 191029-bouton-KC-representations_per_PN
num_exp = 1000

total_claws = int(sum(tbl.num_claws))
num_pn = tbl.shape[0]
btn_prob = tbl.norm_num_boutons

sfl_set = np.zeros([num_pn,num_exp])

for i in range(num_exp):
    t1 = np.random.choice(num_pn,total_claws,p=btn_prob)
    for j in t1:
        sfl_set[j,i] += 1

tbl['rd_btn_mean'] = np.nanmean(sfl_set,1)
tbl['rd_btn_sd'] = np.nanstd(sfl_set,1)

plot_tbl = tbl.copy()
fig, ax = plt.subplots()
ind = np.arange(plot_tbl.shape[0])
width = 0.4

r1 = ax.bar(ind + 0.2, plot_tbl['rd_btn_mean'], width, color=sns.xkcd_palette(['windows blue']), align='center')
r2 = ax.bar(ind + 0.6, plot_tbl['num_claws'], width, color=sns.xkcd_palette(['amber']), align='center')
r3 = ax.errorbar(ind + 0.2, plot_tbl['rd_btn_mean'], yerr=plot_tbl['rd_btn_sd'], fmt = 'none', color=sns.xkcd_palette(['windows blue']))

ax.legend((r1[0], r2[0]), ['null model', 'observed'], fontsize=12, loc='upper right')

plt.xticks(ind + 0.5, plot_tbl.short_glom_name, rotation='vertical')

ticks = ax.get_xticklabels()

for i,tick in enumerate(ticks):
        tick.set_color(plot_tbl['color'].iloc[i])

ax.set_ylabel("claws")
ax.set_xlim([-1,114])
fig.set_size_inches(30,8)
plt.show()
fig.savefig(save_path + '191117-allKCs_claws_ObvsBtnNull.png', bbox_inches='tight')


# 191118 to do!!!!
# more synapses per KC between community connections and non-community connections?
ana = ana_all_rd
ana_obj = ana.conn_data['glom_claw_contracted']
all_ids = ana_obj.col_ids
omm_idx = find_elements(all_ids,comm_anno_ids)
reorder_glom = comm_idx + [i for i in range(len(all_ids)) if i not in comm_idx]
ro_conn = ana_obj.conn['1s'][:,reorder_glom]


#191119
# randomize connectivity based on the number of claws each PN has
# More KCs with multiple claws receiving inputs from community PNs?
# some minor modifications of 190428-random_draw_community_effect_size.py

ana = ana_all_rd
conn_data = ana.conn_data['glom_kc_in_claw_unit']
ob_conn, glom_prob, glom_idx_ids = get_conn_prob_idx(conn_data)

reordered_glom = comm_ids + [i for i in glom_idx_ids if i not in comm_ids]
reorder_idx = reorder(reordered_glom, glom_idx_ids)
ro_conn = ob_conn[:, reorder_idx]

comm_prob = np.sum(ro_conn[:,:10])/np.sum(ro_conn)

claw_num = np.sum(ob_conn,1)

def get_random_conn(num_row, claw_set, prob):
    # given num_row, claw_set,
    num_col = 2
    num_row = int(num_row)
    conn = np.zeros((num_row,num_col))
    for i in range(int(num_row)):
        for j in np.random.choice([0,1],size=int(np.random.choice(claw_set)),p=prob).astype(int):
            conn[i,j] += 1
    return  conn

prob = comm_prob
sfl_results = []
for i in range(1000):
    sfl_conn = get_random_conn(claw_num.shape[0], claw_num, prob=[prob, 1-prob])
    sfl_results.append(pd.Series(sfl_conn[:,0]).value_counts(sort=False))

summary = pd.DataFrame(sfl_results)
s_mean = summary.mean(0)
s_std = summary.std(0)

# observed community frequency
target_freq = ro_conn[:,:10].sum(axis=1)
ob = pd.Series(target_freq).value_counts()


fig, ax = plt.subplots(figsize=(12, 8))
ind = ob.index
width = 0.4

r1 = ax.bar(ind + 0.4, ob, width, color='r', align='center')
r2 = ax.bar(ind+ 0.4 + width, s_mean, width, color='b', align='center', yerr=s_std, error_kw=dict(ecolor='black', lw=3, capsize=5, capthick=2))
# ax.legend((r1[0], r2[0]), ('random draw KCs', 'input-based null model'), fontsize=15)
ax.legend((r1[0], r2[0]), ('observed KCs', 'claw null model'), fontsize=15)
ax.tick_params(axis='both', labelsize=25)
ax.set_xlim(0,max(ind)+1)
ax.set_ylim(top=550)
plt.xticks(ind + 0.2 + width, [int(i) for i in ind])
plt.show()
save_path = '/Users/zhengz11/myscripts/data_results/191119-rerun_main_figs/'
# fig.savefig(save_path + "191119-CommEffectSize_RandomDrawKCs_vs_InputNullModel.png", bbox_inches='tight')
fig.savefig(save_path + "191128-CommEffectSize_RandomDrawKCs_vs_ClawNullModel.png", bbox_inches='tight')

# chi square test p < 1e-10
from scipy.stats import chisquare
t21 = chisquare(ob[:5], s_mean[:5])
# Power_divergenceResult(statistic=401.05640404183276, pvalue=1.7562241520029914e-84)

##-----------------------------------------
# note for the effect size comparison with bouton null model,
# please refer to 190428-random_draw_community_effect_size.py
