# Suppl Fig 2 B,C from PNKC2019_figs_v10_200406DB.pptx
# This script is copied from PNKC2019_figs_v10_200406DB.pptx

# save_path = "/Users/zhengz11/myscripts/data_results/200331-hemibrain_data/"
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

pn_kc_path = "/Users/zhengz11/myscripts/git_clone/pn_kc/"
# exec(open(pn_kc_path + "/connectivity/analysis.py").read())
# exec(open(pn_kc_path + "/connectivity/startup.py").read())

local_path = "/Users/zhengz11/myscripts/git_clone/pn_kc/"
exec(open(local_path + "connectivity/load_pn_metadata_v2.py").read())

data_path = local_path + "data/skids/"
pn_skids = load_json(data_path + "pn")
rd = load_json(data_path + "RandomDraw")
t1p = load_json(data_path + "t1p")

path = '/Users/zhengz11/myscripts/git_clone/pn_kc/data/'
ana_all_rd = ar.Analysis.init_connectivity(path, pn_skids, rd + t1p, 'pn_all_kc')

exec(open(local_path + "connectivity/load_pn_tbl.py").read())

# this piece is copied from 191022-rerun_figs.py
##--------------------------------------------
conn_data = ana.conn_data['pn_kc']
pn_gloms = df_lookup('glom_anno_id',
[ana.pn_mapping.skids_to_types(i)[0] for i in conn_data.col_ids],'short_glom_name',glom_btn_table)

t11 = np.copy(conn_data.conn['5s'])
t11[np.where(t11)]=1
pn_connections = t11.sum(0)

fafb_pn_tbl = pd.DataFrame({'fafb_glom': pn_gloms, 'fafb_connections': pn_connections})
fafb_glom_tbl = fafb_pn_tbl.groupby(['fafb_glom']).sum().assign(fafb_perc_connections=lambda x: x.fafb_connections/sum(x.fafb_connections))
fafb_glom_tbl['glom'] = fafb_glom_tbl.index
fafb_glom_tbl['fafb_glom'] = fafb_glom_tbl.index

# read the hemibrain downloaded data

# hemi_path = local_path + "/Users/zhengz11/myscripts/data_results/200331-hemibrain_data/"
hemi_glom_tbl = pd.read_excel(local_path + 'data/200331-hemibrain_pn_connections_3syn.xlsx', sheet_name="glom")

glom_merged_tbl =pd.merge(hemi_glom_tbl, fafb_glom_tbl, how='outer', on='glom')[:51].sort_values('glom')

comm_gloms = df_lookup('glom_id', comm_ids, 'short_glom_name', glom_btn_table)

p_tbl = glom_merged_tbl.copy()

fig, ax = plt.subplots()
ind = p_tbl.index
width = 0.4

r1 = ax.bar(ind + 0.2, p_tbl.fafb_perc_connections, width, color=sns.xkcd_palette(['windows blue']), align='center')

r2 = ax.bar(ind + 0.6, p_tbl.glom_perc_conn, width, color=sns.xkcd_palette(['amber']), align='center')

ax.legend((r1[0], r2[0]), ['FAFB (the current study)', 'Hemibrain (Xu et al. 2020)'], fontsize=20, loc='upper right')

plt.xticks(ind + 0.4, p_tbl.glom)

ax.set_ylabel('% connections', fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=15)
plt.xlim(-0.2,50)
fig.set_size_inches(40,10)
# fig.savefig(save_path + "200408-compare_fafb_hemibrain_connections.png",dpi=600, bbox_inches='tight')


# produce one with R2 and p value
from scipy import stats
def r2(x, y):
    return stats.pearsonr(x, y)[0] ** 2
sns.jointplot(x="fafb_perc_connections", y="glom_perc_conn", data=p_tbl, kind="reg", stat_func=r2)
# r-squred value: 0.83

## 200406 scatter plot without colored dots
fig, ax = plt.subplots()
# sns.regplot(x="fafb_perc_connections", y="glom_perc_conn", data=p_tbl, scatter_kws={"s": 20, 'edgecolors':'none', 'facecolors':p_tbl_wcol['color']})
sns.regplot(x="fafb_perc_connections", y="glom_perc_conn", data=p_tbl,
scatter_kws={"s": 20, 'edgecolors':'none'})
ax.set_ylim([0,0.055])
ax.set_xlim([0,0.055])
ax.set_xlabel('FAFB (in the current study)')
ax.set_ylabel('Hemibrain (Xu et al. 2020)')
ax.set_title("Percentage of connections for each PN type")
fig.savefig(hemi_path + "200406-compare_fafb_hemibrain_connections_scatterplot_3syn_woCol.png",dpi=600, bbox_inches='tight')


import scipy

slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(p_tbl_wcol.glom_perc_conn, p_tbl_wcol.fafb_perc_connections)
# r-squared value r_value**2 = 0.83

# https://stats.stackexchange.com/questions/146804/difference-between-statsmodel-ols-and-scikit-linear-regression

# https://stats.stackexchange.com/questions/249892/wildly-different-r2-between-statsmodels-linear-regression-and-sklearn-linear?noredirect=1&lq=1

# https://stackoverflow.com/questions/29664471/why-would-r-squared-decrease-when-i-add-an-exogenous-variable-in-ols-using-pytho

'''
# in case if I want to color the dots according to Fig 1D
col_tbl = tbl.loc[:,['short_glom_name','color']]
col_tbl['glom']=col_tbl.short_glom_name
col_tbl.drop_duplicates(subset='short_glom_name', inplace=True)
p_tbl_wcol = pd.merge(p_tbl, col_tbl, on='glom', how='left')

fig, ax = plt.subplots()
sns.regplot(x="fafb_perc_connections", y="glom_perc_conn", data=p_tbl_wcol, scatter_kws={"s": 20, 'edgecolors':'none', 'facecolors':p_tbl_wcol['color']})
ax.set_ylim([0,0.055])
ax.set_xlim([0,0.055])
ax.set_xlabel('FAFB (in the current study)')
ax.set_ylabel('Hemibrain (Xu et al. 2020)')
ax.set_title("Percentage of connections for each PN type")
# fig.savefig(hemi_path + "200331-compare_fafb_hemibrain_connections_scatterplot_3syn_wCol.png",dpi=600, bbox_inches='tight')


# a version of SFig2C (compare hemibrain and our data with the bars) PNKC2019_figs_v10_200407DB.pptx
fig, ax = plt.subplots()
ind = p_tbl.index
width = 0.4

r1 = ax.bar(ind + 0.2, p_tbl.fafb_perc_connections, width, color=sns.xkcd_palette(['windows blue']), align='center')

r2 = ax.bar(ind + 0.6, p_tbl.glom_perc_conn, width, color=sns.xkcd_palette(['amber']), align='center')

ax.legend((r1[0], r2[0]), ['FAFB (the current study)', 'Hemibrain (Xu et al. 2020)'], fontsize=12, loc='upper right')

plt.xticks(ind + 0.4, p_tbl.glom)

ticks = ax.get_xticklabels()
for i,tick in enumerate(ticks):
        tick.set_color(tbl.query("short_glom_name==@tick.get_text()").color.iloc[0])
#        if tick.get_text() in comm_gloms:
#            tick.set_weight("extra bold")


ax.set_ylabel('% connections', fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=14)
plt.xlim(-0.2,50)
fig.set_size_inches(40,10)
# fig.savefig(hemi_path + "200331-compare_fafb_hemibrain_connections.png",dpi=600, bbox_inches='tight')
'''
