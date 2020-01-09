
# 200109 copy from /Users/zhengz11/myscripts/bocklab_git/bocklab/zhihao/mushroom_py/v10/191201-inter_bouton_dist.py
# rename to 200109-inter_bouton_dist
# this should generate distance between boutons (currently Fig 4G)
import scipy as sc

save_path = '/Users/zhengz11/myscripts/data_results/191201-inter_bouton_dist/'

comm_ids = [22, 44, 32, 42, 34, 38, 30, 23, 47, 4]
comm_anno_ids = df_lookup('glom_id',comm_ids,'glom_anno_id',glom_btn_table)
comm_pns = ana.pn_mapping.types_to_skids(comm_anno_ids)

ana = ana_all_rd

all_pns = ana.col_neurons.values()
btn_ids = []
btn_ctds = []
for pn in all_pns:
    btn_ids.extend(pn.segments.ids)
    btn_ctds.extend(pn.segments.centroids)

nbtn = len(btn_ids)
dist_conn =np.zeros((nbtn,nbtn))
t1 = range(nbtn)
for i,j in [(i,j) for i in t1 for j in t1 if i>j]:
    d = np.linalg.norm(btn_ctds[i] - btn_ctds[j])
    dist_conn[i,j] = d
    dist_conn[j,i] = d

pn_ids = segid_to_skid(btn_ids)

upn = np.unique(pn_ids)
# n_pn = range(len(upn))
pb_dist = np.zeros((len(upn),len(pn_ids)))

for idx,val in enumerate(upn):
    for idx_btn in range(pb_dist.shape[1]):
        pb_dist[idx,idx_btn] = np.min(dist_conn[find_elements(pn_ids,[val]),idx_btn])

# output a dist_conn with boutons as columns (pn_ids, the pn ids of every bouton) and PN as rows (upn, PN id of every PN)

# reorder columns to put the community boutons (first 110 boutons) at first
comm_idx = find_elements(pn_ids, comm_pns)
reorder = comm_idx + [i for i in range(len(pn_ids)) if i not in comm_idx]
pb_dist = pb_dist[:,reorder]
new_pn_ids = [pn_ids[i] for i in reorder]


# reorder rows to put 16 community PNs in the first 16 rows
new_upn = pd.unique(new_pn_ids)
reorder = [upn.tolist().index(i) for i in new_upn]

pb_dist = pb_dist[reorder,:]


rev = pb_dist.copy()

import seaborn as sns
sns.set_style('whitegrid')
t4 = rev.copy()
t5 = t4[:,:110][:16,:]
t5 = t5.flatten()
t5 = t5[np.where(t5)]
plt.tick_params(labelsize=12)
plt.hist(t5, bins=range(0,50000,2500))
plt.xlim(left=0)
# plt.savefig(save_path + '191201-comm_inter-subtype_btn_dist.png')
# save t5 for doing statistical test
CommToComm = t5.copy()


t4 = rev.copy()
t6 = t4[:,:110][16:,:]
plt.tick_params(labelsize=12)
plt.hist(t6.flatten(), bins=range(0,50000,2500))
plt.xlim(left=0)
# plt.savefig(save_path + '191201-CommToNoncomm_inter-subtype_btn_dist.png')
CommToNoncomm = t6.flatten().copy()

t4 = rev.copy()
t6 = t4[:,110:][16:,:]
t6 = t6.flatten()
t6 = t6[np.where(t6)]
plt.tick_params(labelsize=12)
plt.hist(t6.flatten(), bins=range(0,50000,2500))
plt.xlim(left=0)
# plt.savefig(save_path + '191201-NoncommToNoncomm_inter-subtype_btn_dist.png')
NoncommToNoncomm = t6.flatten().copy()

# https://data.library.virginia.edu/the-wilcoxon-rank-sum-test/
CC_to_CN = sc.stats.ks_2samp(CommToComm, CommToNoncomm)
# Ks_2sampResult(statistic=0.50893470790378, pvalue=6.4e-322)

CC_to_NN = sc.stats.ks_2samp(CommToComm, NoncommToNoncomm)
# Ks_2sampResult(statistic=0.3354948216340621, pvalue=5.179808697970495e-155)


##------------------------
# 191201 plot the two comparison distribution together
sns.set_style('white')
colors = sns.color_palette()
#sets up the axis and gets histogram data
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
t1 = ax1.hist([CommToComm, CommToNoncomm])
n, bins, patches = ax1.hist([CommToComm, CommToNoncomm],bs)
ax1.cla() #clear the axis

#plots the histogram data
width = (bins[1] - bins[0]) * 0.4
bins_shifted = bins + width
ax1.bar(bins[:-1], n[0], width, align='edge', color=colors[0])
ax2.bar(bins_shifted[:-1], n[1], width, align='edge', color=colors[1])

#finishes the plot
ax1.set_ylabel("Count", fontsize=16, color=colors[0])
ax2.set_ylabel("Count", fontsize=16, color=colors[1])
ax1.tick_params('y', colors=colors[0])
ax2.tick_params('y', colors=colors[1])
ax1.tick_params(labelsize=16)
ax2.tick_params(labelsize=16)
ax1.tick_params('x', bottom=True)
plt.xlim(left=0)
plt.xticks([0]+[i-250 for i in range(10000,50000,10000)], range(0,50000,10000))
fig.set_size_inches([12,6])
plt.tight_layout()
plt.show()
fig.savefig(save_path + '191202-CommComm_vs_NoncommNoncomm_hist.png', bbox_inches='tight', dpi=600)
