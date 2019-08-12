
# %run startup.py first
from load_pn_metadata_v2 import *
from mushroom_2to3.detect_community import *

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
