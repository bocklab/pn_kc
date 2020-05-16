
# This will produce the table for Fig 2A

ana = ana_all_rd

# y number of boutons
# y number of downstream KCs
# per PN, also per glomerulus

pn_skids = load_json(local_path + "data/skids/pn")
# pn_skids = cc.get_skids_from_annos(fafb_c, [['right_calyx_PN'], ['has_bouton']], ["multiglomerular PN"])

boutons_per_pn = [len(ana.col_neurons[i].segments.ids) for i in pn_skids]
gloms = df_lookup('glom_anno_id', ana.pn_mapping.skids_to_types(pn_skids),
'short_glom_name', glom_btn_table)

meta_tbl = pn_meta_table.copy()
t1 = meta_tbl.Significance.copy()
t1[pd.isna(t1)] = 'Unknown'
meta_tbl.Significance = t1

tbl = pd.DataFrame({'pn_skids':pn_skids,
'num_boutons': boutons_per_pn,
'short_glom_name': gloms,
'significance': pd.concat([meta_tbl.query('glom==@i').Significance for i in gloms])})


# add connections, namely the number of KCs per PN
conn_data = ana.conn_data['pn_kc_contracted']
pn_gloms = conn_data.col_ids
t11 = np.copy(conn_data.conn['5s'])
t11[np.where(t11)]=1
pn_connections = t11.sum(0)
tbl = tbl.merge(pd.DataFrame({'pn_skids': conn_data.col_ids,'outdegree': pn_connections}), how='outer',on='pn_skids')


# add claws per PNs
t16 = ana.conn_data['pn_claw_contracted'].conn['1s'].copy()
t16[np.where(t16<3)]=0
t16[np.where(t16)]=1
ds_claws = t16.sum(0)
tbl = tbl.merge(pd.DataFrame({'pn_skids': conn_data.col_ids,'num_claws': ds_claws}), how='outer',on='pn_skids')

tbl = tbl.assign(norm_num_boutons=lambda x: x.num_boutons/sum(tbl.num_boutons), norm_num_kcs=lambda x: x.outdegree/sum(tbl.outdegree), norm_num_claws=lambda x: x.num_claws/sum(tbl.num_claws))

tbl.significance = pd.Categorical(tbl.significance, ordered=True, categories=['Food','Aversive','Pheromonal','Egg-laying','Unknown'])
tbl = tbl.sort_values(by=['significance','short_glom_name'])

t1_dict = dict(zip(['Food', 'Aversive', 'Pheromonal', 'Egg-laying', 'Unknown'],['green','red','purple','blue','black']))

tbl['color'] = tbl['significance'].apply(lambda x: t1_dict.get(x))

comm_gloms = df_lookup('glom_id', comm_ids,'short_glom_name', glom_btn_table)
tbl['community'] = [True if i in comm_gloms else False for i in tbl['short_glom_name']]

## to plot
plot_tbl = tbl.copy()

t1_dict = dict(zip(['Food', 'Aversive', 'Pheromonal', 'Egg-laying', 'Unknown'],['green','red','purple','blue','black']))

plot_tbl['color'] = plot_tbl['significance'].apply(lambda x: t1_dict.get(x))


# below are old comments
#------------------------------------------------------------------------
# copy from /Users/zhengz11/myscripts/bocklab_git/bocklab/zhihao/mushroom_py/v10/191029-bouton-KC-representations_per_PN.py
# save_path = '/Users/zhengz11/myscripts/data_results/191028-bouton-KC_representations/'
# 200110 change name to 200110-pn_conn_tbl.py
# 200326 change name to load_pn_tbl.py
