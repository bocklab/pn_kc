
import pandas as pd

wd = "~/myscripts/git_clone/pn_kc/pn_kc/"
data_path = ""

glom_id_table = pd.read_excel(wd + "tables/171012-glom_index_list.xlsx")
glom_btn_table = pd.read_excel(wd + "tables/180320-glom_name_id.xlsx")
pn_meta_table = pd.read_excel(wd + "tables/191029-processed_pn_metadata.xlsx")
tbl = pd.read_excel(wd + "tables/191029_tbl.xlsx", index_col=0)


def get_conn_prob_idx(conn_obj, glom_data_table=glom_btn_table):
    '''
    Parameters
    ----------
    conn_obj:            a connectivity matrix object
                         e.g. conn_obj = ana.conn_data['glom_kc_in_claw_unit']
    glom_data_table :    a table with meta information (annotation ids, short
                         glomerular names, annotation names, index ids used for matrix columns and rows, glom_prob based on bouton counts, etc.) on PNs
    Returns
    -------
    conn:                the PN (gloms) - KCs connectivity matrix
    glom_prob:           meta info in the PN table
    glom_idx_ids:        the indices I used for showing the conditional matrix
    '''
    glom_anno_ids = conn_obj.col_ids

    glom_idx_ids = df_lookup(
        'glom_anno_id', glom_anno_ids, 'glom_id', glom_btn_table)

    # get glom_prob
    glom_prob = df_lookup('glom_anno_id', glom_anno_ids,
                          'norm_bouton_count', glom_btn_table)

    conn = conn_obj.conn['1s']
    return conn, glom_prob, glom_idx_ids






original_order = [53,10,43,37,36,35,33,31,52,25,
20,15,26,3,6,9,30,22,23,34,
2,32,29,38,42,13,4,44,19,16,
47,8,40,45,46,48,1,5,28,27,
7,18,12,17,0,39,41,24,21,11,
14,49,50,51]

NewClusterOrder = [13,32,30,38,23,22,19,47,16,42,34,4,8,11,1,36,7,51,27,25,24,12,6,33,9,48,49,10,40,2,31,50,44,59,43,0,26,29,3,5,14,15,17,35,53,21,52,20,46,41,39,37,18,28]

# old community glomeruli
# CommunityGlom = [4, 13, 16, 19, 22, 23, 30, 32, 34, 38, 42, 47]

ClusterOrder0707 = [22,44,32,42,34,
38,30,23,47,4,19,
16,29,13,11,21,
50,3,26,20,17,
48,53,15,12,10,
1,9,5,6,36,
25,24,37,51,49,
39,59,43,41,52,
40,14,33,27,35,
31,7,8,28,46,
18,2,0]

comm_ids = [22, 44, 32, 42, 34, 38, 30, 23, 47, 4]

'''
FAFB2018_olfactory_PNs
FAFB2018_KCs
FAFB2018_non-mALT_PNs
FAFB2018_Subesophogeal_PN
FAFB2018_other_PNs
exclude:
"multiglomerular PN"
based on "20170519 PN Comparison With 20180211.xlsx"
'''

nblast_order = ["DL2v","DL2d","DM2","V", "VM3","VA4","DM1","VM2","VA2","DM3","DM4",
"VP3", "VL1", "DC2", "VP1", "DA4m", "VA5", "DA3", "VA3", "DM5", "VM1", "DA2", "VA7m", "DA4l", "VA7l", "VC1", "VC2", "DP1l", "DP1m", "VL2p",
"DA1", "VA1v", "DL4", "VA1d", "D", "DC3", "VM7v", "VC4", "DL1", "DC1", "VA6", "VL2a",
"VC3m", "VM4", "VP2", "VC5", "DL3", "VM5d", "VM5v", "VM7d", "DL5", "DC4", "DM6", "VC3l"]
