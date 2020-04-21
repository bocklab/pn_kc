
import pandas as pd
import numpy as np

fpath = local_path + "data/160928-Caron_suppl_table.xlsx"

class CaronAnalysis(object):
    def __init__(self):
        self.conn_data = {}

    @classmethod
    def init_processing(cls):
        c = cls()
        c.suppl_tbl, c.glom_tbl = import_caron_data(fpath)
        c.conn_data, c.kc_mapping = parse_data(c.suppl_tbl, c.glom_tbl)

        co_matrix = c.conn_data['caron_glom_kc'].co_matrix['1s']
        # pd.DataFrame(caron_gk).to_excel('160928_caron_fig3.xlsx')
        co_matrix = c.conn_data['caron_glom_claw'].fill_co_matrix_diagonal(c.kc_mapping, co_matrix, syn='1s')
        return c

def import_caron_data(path=fpath):
    return (pd.read_excel(path, 'combined'), pd.read_excel(path, 'glom_seq'))

def parse_data(suppl_tbl, glom_tbl):
    num_kc = suppl_tbl.shape[0]
    num_glom = glom_tbl.shape[0]

    claws = suppl_tbl.filter(like='glom').values.flatten()
    num_claw = claws[np.invert(np.isnan(claws))].shape[0]

    caron_gk = np.zeros((num_kc, num_glom))
    caron_gc = np.zeros((num_claw,num_glom))

    claw_ids = []
    kc_ids = []
    j = 0
    for row_i in range(num_kc-1):
        col_i = suppl_tbl.filter(like='glom').loc[row_i,]
        col_i = col_i[np.invert(np.isnan(col_i))]
        for i in col_i:
            caron_gk[row_i, int(i-1)] += 1
            caron_gc[j, int(i-1)] += 5
            kc_id = suppl_tbl.loc[row_i,'cell_number']
            claw_ids.append(str(kc_id) + '.' + str(j))
            kc_ids.append(kc_id)
            j += 1

    glom_ids = glom_tbl['glom_number'].tolist()

    conn_data = {}
    conn_data['caron_glom_kc'] = ConnectivityMatrix('caron_glom_kc', caron_gk,
                    glom_ids, suppl_tbl['cell_number'].tolist())

    conn_data['caron_glom_claw'] = ConnectivityMatrix('caron_glom_claw', caron_gc,
                    glom_ids, claw_ids)

    kc_mapping = pd.DataFrame({'segments': claw_ids, 'skids': kc_ids}, columns=['segments','skids'])
    return conn_data, kc_mapping

def freq_identical_pairs(co_m):
# frequency distribution in the diagonal of the glom-glom matrix (Fig 4b)
    diagonal = np.diagonal(co_m)
    return pd.Series(diagonal[diagonal!=0]).value_counts()

def freq_non_identical_pairs(co_m):
    t2 = co_m[np.triu_indices(co_m.shape[0], k=1)]
    t2 = t2[t2>1]
    return pd.Series(t2).value_counts()

def get_caron_fig4a(conn, w=False):
    sum_inputs = sum_conditioned_inputs(conn, w, False, threshold=1)
    sum_inputs[np.diag_indices_from(sum_inputs)] = 0
    for i in np.where(conn>1)[1]:
        sum_inputs[i,i] += 1
    return sum_inputs

def get_caron_FigS2():
    c_suppl_tbl, c_glom = import_caron_data(fpath)
    r = {}
    for t in ['all', 'alpha/beta', "alpha'/beta'", 'gamma']:
        if t is 'all':
            t1 = c_suppl_tbl
        else:
            t1 = c_suppl_tbl.query("KC_subclass==@t")
        t2 = t1[[col for col in t1 if col.startswith('glom')]]
        t2 = np.array(t2).flatten()
        t3 = pd.Series(t2[~np.isnan(t2)]).value_counts()
        zero_idx = [i for i in c_glom['glom_number'] if i not in t3.index]
        glom_names = [c_glom.loc[c_glom.glom_number==i,'glom_class'].to_string(header=False, index=False) for i in t3.index.tolist() + zero_idx]
        r[t] = pd.DataFrame({'inputs': t3.tolist() + [0]*len(zero_idx), 'glom': glom_names})
    return r

# old comments
#---------------------------------------------------------------
# from mushroom_2to3.shuffle import *
# from mushroom_2to3.build_connectivity import *
# from mushroom_2to3.detect_community import *
# likely need to import the mushroom package

# fpath = "/Users/zhengz11/myscripts/data_results/160928-caron_equiv/160928-Caron_suppl_table.xlsx"
# fpath = local_path + "data/160928-Caron_suppl_table.xlsx"
# abKC_path = "/Users/zhengz11/myscripts/data_results/171012-1D_olfactory_space/171025-Caron_suppl_table_abonly.xlsx"
