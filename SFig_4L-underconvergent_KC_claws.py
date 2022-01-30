
ana = ana_all_rd

negs = ["VP2", "V", "VA7m", "DA3", "DA2", "DA1", "VP1", "VL1", "DA4m", "DL4"]
neg_gloms = df_lookup("short_glom_name", negs, "glom_anno_id", glom_btn_table)

# neg_gloms = df_lookup("short_glom_name", t3.query("counts_NumCluster4>800").short_glom_name,"glom_anno_id", glom_btn_table)

conn_obj = ana.conn_data['glom_kc_contracted']
neg_kcs = [conn_obj.row_ids[j] for j in np.unique(np.nonzero(conn_obj.conn['5s'][:,find_elements(conn_obj.col_ids,neg_gloms)])[0])]

sc_set = [len(ana.row_neurons[i].segments.ids) for i in neg_kcs]

# np.save(save_path + "211031-UnderconvergentTotalClaws.npy", np.array(sc_set))

noncold_set = [len(ana.row_neurons[i].segments.ids) for i in rd + t1p if i not in neg_kcs]

# np.save(save_path + "211031-NonUnderconvergentTotalClaws.npy", np.array(noncold_set))

# plotting
# claw distribution for KCs downstream of the cold cluster
sns.set_style("white")
claw_num = sc_set
freq = pd.Series(claw_num).value_counts()
t1 = [i for i in range(0, max(freq.index)) if i not in freq.index]
freq = pd.concat([pd.Series([0]*len(t1),index=t1),freq]).sort_index()
freq_cold = freq

ana = ana_all_rd
sk = ana.row_neurons.keys()
claw_num_all = [len(ana.row_neurons[i].segments.ids) for i in sk]


freq = noncold_set
freq = pd.Series(noncold_set).value_counts()
t1 = [i for i in range(0, max(freq.index)) if i not in freq.index]
freq = pd.concat([pd.Series([0]*len(t1),index=t1),freq]).sort_index()
freq_noncold = freq

freq_cold[13]=0
freq_noncold.drop(0, inplace=True)
freq_cold.drop(0, inplace=True)


col = sns.color_palette("tab10")
fig, ax1 = plt.subplots()
width = 0.4
ax2 = ax1.twinx()
r1 = ax1.bar(freq_cold.index - 0.2, freq_cold, width, color=col[0], align='center', linewidth=0)

r2 = ax2.bar(freq_noncold.index + 0.2, freq_noncold, width, color=col[1], align='center', linewidth=0)

plt.xticks(np.arange(1, 14), range(1, 14), fontsize=16)
ax1.tick_params(axis="both", labelsize=16)
ax2.tick_params(axis="y", labelsize=16, colors=col[1])
ax2.set_ylim(0,250)
ax1.set_ylim(0,180)
ax1.tick_params('y', colors=col[0])
ax1.set_xlim(0.4,13.7)
fname = "211031-UnderconvergentClaws_vs_otherClaws.png"

ax1.set_title(fname, pad=18)

def autolabel(rects, ax,xpos='center'):
    """
    Attach a text label above each bar in *rects*, displaying its height.

    *xpos* indicates which side to place the text w.r.t. the center of
    the bar. It can be one of the following {'center', 'right', 'left'}.
    """

    xpos = xpos.lower()  # normalize the case of the parameter
    ha = {'center': 'center', 'right': 'left', 'left': 'right'}
    offset = {'center': 0.5, 'right': 0.57, 'left': 0.43}  # x_txt = x + w*off

    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()*offset[xpos], 1.01*height,
                '{}'.format(height), ha=ha[xpos], va='bottom', fontsize=14)

autolabel(r1, ax1)
autolabel(r2, ax2)
fig.set_size_inches(13,6)

save_path = "/Users/zhengz11/myscripts/data_results/211031-bar_graph_replots/"
# fig.savefig(save_path + fname + ".png", bbox_inches="tight")

from scipy import stats
print(stats.ks_2samp(sc_set, noncold_set))

import numpy as np
for i in  [np.mean(sc_set),
np.std(sc_set),
np.mean(noncold_set),
np.std(noncold_set)]:
    print(i)
