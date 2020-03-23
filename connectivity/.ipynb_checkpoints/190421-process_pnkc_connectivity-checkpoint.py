
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

%run startup_py3.py
# %run load_pn_metadata_v1.py
%run medium_term_functions.py

# f_path = '/Users/zhengz11/myscripts/data_results/190412-varying_structure_movie/'

# get PN skeleton ids from CATMAID
pn_skids = cc.get_skids_from_annos(
    fafb_c, [['right_calyx_PN'], ['has_bouton']], ["multiglomerular PN"])

# get KC skeleton ids from CATMAID
# rd is random draw manually traced KCs
rd = cc.get_skids_from_annos(fafb_c,
                             [['Random Draw 1 KC', 'Random Draw 2 KC'], ['Complete']],
                             ['KCaBp', 'KCyd'])

# an analysis object that includes everything I need for the analysis of the connectivity between PNs and KCs (e.g. boutons, claws, connectivity matrices of different forms, etc.)
# with manually traced KCs (random draw) only
ana_rd = ar.Analysis.init_connectivity(fafb_c, pn_skids, rd)


#-------------------------------------------------
# part 1
# How I analyze the PN-KC connectivity and visualize the cluster of community PNs
ana = ana_rd

# this is a connectivity matrix
conn_data = ana.conn_data['glom_kc_in_claw_unit']

ob_conn, glom_prob, glom_idx_ids = get_conn_prob_idx(conn_data)

stat = [get_raw_inputs(shuffle_glom_kc_w_prob(ob_conn, glom_prob)) for i in range(num_exp)]

stat = np.array(stat)
sd = np.nanstd(stat, axis=0)
avg = np.nanmean(stat, axis=0)

ob_ci = get_raw_inputs(ob_conn)
comm_zscore = np.divide(np.subtract(ob_ci, avg), sd)

# clustering
cm_zs = PairMatrix('', comm_zscore, glom_idx_ids)
reorder_idx = km_cluster(cm_zs.conn)
# reorder_idx = reorder(NewClusterOrder, glom_idx_ids)

t1_zs = cm_zs.reorder(reorder_idx, return_new=True)
t1 = t1_zs
fig, ax1 = plt.subplots()
t1 = sns.heatmap(t1.conn, xticklabels=False, yticklabels=False, square=True, ax=ax1, vmin=-8.53, vmax=8.53, cbar=False, cmap="RdBu_r")
# fig.savefig()


#-------------------------------------------------
# part 2
# save a PN to KC connectivity matrix for the modeling
# the first 16 PNs of the columns are community PNs (given by part 1 results)

# for modeling, using the manually traced KCs only, so they have all synaptic inputts from PNs
ana = ana_rd

# ids of the community glomeruli (look up ids in data/180320-glom_name_id.xlsx)
comm_ids = [22, 44, 32, 42, 34, 38, 30, 23, 47, 4]

# ids of the glomeruli that are "blue" in the k-means clustering results
blue = [27, 35, 31, 7, 8, 28, 46, 18, 2, 0]

# map these glomeruli onto PN ids
comm_anno_ids = df_lookup('glom_id',comm_ids,'glom_anno_id',glom_btn_table)
comm_pns = ana.pn_mapping.types_to_skids(comm_anno_ids)
blue_anno_ids = df_lookup('glom_id',blue,'glom_anno_id',glom_btn_table)
blue_pns = ana.pn_mapping.types_to_skids(blue_anno_ids)


# ids of the PN-KC connectivity matrix (each column is a PN, each row is a KC, ids stored in col_ids and row_ids)
all_sks = ana.conn_data['pn_kc_contracted'].col_ids


# reorder the connectivity matrix so that:
# first 16 PNs are community PNs
# last 22 PNs are blue PNs
comm_idx = find_elements(all_sks,comm_pns)
blue_idx = find_elements(all_sks,blue_pns)
reorder_pn = comm_idx + [i for i in range(len(all_sks)) if i not in comm_idx + blue_idx] + blue_idx

pk_conn = ana.conn_data['pn_kc_contracted'].conn['1s'][:,reorder_pn]

col_names = ana.skid_to_name([ana.conn_data['pn_kc_contracted'].col_ids[i] for i in reorder_pn])

row_names = ana.skid_to_name(ana.conn_data['pn_kc_contracted'].row_ids)

# save the connectivity matrix
# t8 = pd.DataFrame(pk_conn, index=row_names, columns=col_names).to_csv('data/190319-440RDKC_PNfirst16comm_last22blue.csv')
