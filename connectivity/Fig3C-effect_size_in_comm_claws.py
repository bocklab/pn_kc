
# Fig3C, 200326 PNKC2019_v9_fig_200313DB-ZZfixedSuppl6B.pptx
from mushroom_2to3 import detect_community as dc
ana = ana_all_rd

comm_ids = [22, 44, 32, 42, 34, 38, 30, 23, 47, 4]

conn_data = ana.conn_data['glom_kc_in_claw_unit']
ob_conn, glom_prob, glom_idx_ids = get_conn_prob_idx(conn_data)

reordered_glom = comm_ids + [i for i in glom_idx_ids if i not in comm_ids]
reorder_idx = dc.reorder(reordered_glom, glom_idx_ids)
ro_conn = ob_conn[:, reorder_idx]

btn_comm_prob = [glom_prob[i] for i in find_elements(comm_ids, glom_idx_ids)]
claw_comm_prob = np.sum(ro_conn[:,:10])/np.sum(ro_conn)

claw_num = np.sum(ob_conn,1)

def get_random_conn(num_row, claw_set, prob):
    # given num_row, claw_set,
    # prob=sum(comm_prob)
    num_col = 2
    num_row = int(num_row)
    conn = np.zeros((num_row,num_col))
    for i in range(int(num_row)):
        for j in np.random.choice([0,1],size=int(np.random.choice(claw_set)),p=prob).astype(int):
            conn[i,j] += 1
    return  conn


prob = sum(btn_comm_prob)
sfl_results = []
for i in range(1000):
    sfl_conn = get_random_conn(claw_num.shape[0], claw_num, prob=[prob, 1-prob])
    sfl_results.append(pd.Series(sfl_conn[:,0]).value_counts(sort=False))
summary = pd.DataFrame(sfl_results)
btn_s_mean = summary.mean(0).copy()
btn_s_std = summary.std(0).copy()


prob = np.sum(claw_comm_prob)
sfl_results = []
for i in range(1000):
    sfl_conn = get_random_conn(claw_num.shape[0], claw_num, prob=[prob, 1-prob])
    sfl_results.append(pd.Series(sfl_conn[:,0]).value_counts(sort=False))
summary = pd.DataFrame(sfl_results)
claw_s_mean = summary.mean(0).copy()
claw_s_std = summary.std(0).copy()



# observed community frequency
target_freq = ro_conn[:,:10].sum(axis=1)
ob = pd.Series(target_freq).value_counts()

# chi square test p < 1e-10
from scipy.stats import chisquare
# t21 = chisquare(ob[:5], s_mean[:5])


## run 200211-- begin-to-end (local random bouton set)
# from 191028-local_randomization_conditional_input_matrix.py

comm_ids = [22, 44, 32, 42, 34, 38, 30, 23, 47, 4]
comm_anno_ids = df_lookup('glom_id',comm_ids,'glom_anno_id',glom_btn_table)
ana_obj = ana.conn_data['glom_kc_in_claw_unit']

comm_idx = find_elements(sfl_cols,comm_anno_ids)

sfl_results = []
local_claw_set = []
stat = []
for t11 in range(num_exp):
    sfl_conn = np.zeros(ob_conn.shape)
    for i,v in enumerate(conn_data.row_ids):
        for j in np.where(rows==v)[0]:
            sfl_conn[i, glom_to_pick[j,np.random.randint(n)]] += 1
    sfl_results.append(pd.Series(sfl_conn[:,comm_idx].sum(1)).value_counts(sort=False))
    print(t11)

summary = pd.DataFrame(sfl_results)
local_s_mean = summary.mean(0).copy()
local_s_std = summary.std(0).copy()


# temporarily make up for n=9 scenario
ob_p = ob.append(pd.Series(0, index=[9]))
btn_s_mean_p = btn_s_mean.append(pd.Series(0, index=[9]))
btn_s_std_p = btn_s_std.append(pd.Series(0, index=[9]))
claw_s_mean_p = claw_s_mean.append(pd.Series(0, index=[9]))
claw_s_std_p = claw_s_std.append(pd.Series(0, index=[9]))

'''
# some script to make all with the same nember of bins
all_set = [ob, btn_s_mean, btn_s_std, claw_s_mean, claw_s_std, local_s_mean, local_s_std]
max_size = max([len(i) for i in all_set])
all_set_p = []
for cts_set in all_set:
    if len(cts_set) < max_set:
        cts_set.append(pd.Series(0, index=[max_size-len(cts_set)]))
    all_set_p.append(cts_set)
'''

fig, ax = plt.subplots(figsize=(12, 8))
ind = ob_p.index
width = 0.2

first = ind + width/2;
er_bars = dict(ecolor='black', lw=1.5, capsize=4, capthick=1.5);

r1 = ax.bar(first, ob_p, width, color='r', align='center')
r2 = ax.bar(first + width, btn_s_mean_p, width, color='b', align='center', yerr=btn_s_std_p, error_kw=er_bars)
r3 = ax.bar(first + width*2, claw_s_mean_p, width, color='y', align='center', yerr=claw_s_std_p, error_kw=er_bars)
r4 = ax.bar(first + width*3, local_s_mean, width, color='c', align='center', yerr=local_s_std, error_kw=er_bars)
ax.legend((r1[0], r2[0], r3[0], r4[0]), ('observed', 'bouton random model', 'claw random model', 'local bouton random model'), fontsize=15)
ax.tick_params(axis='both', labelsize=20)
ax.set_xlim(0,max(ind)+1)
ax.set_ylim(top=550)
plt.xticks(ind + width*2, [int(i) for i in ind])

# fig.savefig(save_path + "200212-CommEffectSize_BtnClawLocalbtnNullModel.png", bbox_inches='tight')

from scipy.stats import chisquare
t21 = chisquare(ob[:9], local_s_mean[:9])



# Chi-squared statistical test
# 200212
from scipy.stats import chisquare
t21 = chisquare(ob, local_s_mean)
t21
# Power_divergenceResult(statistic=23.845043068217716, pvalue=0.004553135508363728)

t19 = chisquare(ob[:8], btn_s_mean[:8])
# Power_divergenceResult(statistic=972.2573265664316, pvalue=1.4514074337354426e-204)

t19 = chisquare(ob[:8], claw_s_mean[:8])
# Power_divergenceResult(statistic=444.9229345486355, pvalue=4.5257981764215116e-91)

# below are old comments
#--------------------------------------------------
# copy from 200210-effect_size_vs_BtnClawRandomModel_newFig2E.py
# Fig3C, 200326 PNKC2019_v9_fig_200313DB-ZZfixedSuppl6B.pptx
# copy and modify from 190428-random_draw_community_effect_size.py
# and 191117-comm_claws-kcs_vs_shuffle.py
