import time
import pandas as pd
import numpy as np
# import seaborn as sns
from . import connect as cc
from .build_connectivity import *

from igraph import *

import scipy as sp
from scipy.stats import chi2_contingency
from functools import reduce


class GlomKcConnectivityMatrix(ConnectivityMatrix):
    """docstring for ."""

    def __init__(self, name, conn, col_ids=[], row_ids=[]):
        ConnectivityMatrix.__init__(self, name, conn, col_ids, row_ids)
        self.fill_diagonal()

    def fill_diagonal(self, syn='1s'):
        conn = self.conn[syn]
        idx = np.where(conn >= 2)
        for i in range(idx[0].shape[0]):
            row = idx[0][i]
            col = int(idx[1][i])
            n = conn[row, col]
            self.co_matrix[syn][col, col] += n * (n - 1) / 2


def get_co_matrix(conn):
    """
    A much simpler implementation of class GlomKcConnectivityMatrix(ConnectivityMatrix).
    Given conn, produce co_matrix with simple functions
    """
    co_matrix = get_ci_matrix(conn.transpose())
    idx = np.where(conn >= 2)
    for i in range(idx[0].shape[0]):
        row = idx[0][i]
        col = int(idx[1][i])
        n = conn[row, col]
        co_matrix[col, col] += n * (n - 1) / 2
    return co_matrix


def freq_identical_pairs(co_m):
    # frequency distribution in the diagonal of the glom-glom matrix (Fig 4b)
    diagonal = np.diagonal(co_m)
    return pd.Series(diagonal[diagonal != 0]).value_counts(sort=False)


def freq_non_identical_pairs(co_m):
    t2 = co_m[np.triu_indices(co_m.shape[0], k=1)]
    t2 = t2[t2 > 1]
    return pd.Series(t2).value_counts(sort=False)


def to_df(ser, col_i): return pd.DataFrame(ser, columns=[col_i])


def join_freq_table(tbl):
    # join a list of frequency tables to a single DataFrame
    # and calculate mean and std
    r = to_df(tbl[0], 0).join([to_df(e, i)
                               for i, e in enumerate(tbl) if i > 0], how='outer')
    r[pd.isnull(r)] = 0
    ci_lower = []
    ci_upper = []
    for i, row in r.iterrows():
        s = np.sort(row)
        ci_lower.append(s[int(np.floor(r.shape[1] * 0.025)) - 1])
        ci_upper.append(s[int(np.ceil(r.shape[1] * 0.975)) - 1])

    r = pd.DataFrame({'mean': r.mean(axis=1), 'sd': r.std(axis=1),
                      'ci_lower': ci_lower, 'ci_upper': ci_upper},
                     columns=['mean', 'sd', 'ci_lower', 'ci_upper'])
    return r


def shuffle_glom_claw_conn(gc_obj, seg_id_mapping):
    # 1) take a glom_claw conn object (note that each claw only has 1 glom input)
    # 2) permute the connections as in Caron et al. to construct permuted glom_claw conn
    # 3) combine from glom_claw to glom_kc matrix
    # 4) get common output matrix and fill in the diagonal
    gc_conn = gc_obj.conn['5s']
    perm = np.zeros(gc_conn.shape)
    perm[np.nonzero(gc_conn)[0], np.random.permutation(
        np.nonzero(gc_conn)[1])] = 5

    # combine from glom_claw_conn to glom_kc_conn
    gc_perm = ConnectivityMatrix('perm', perm, gc_obj.col_ids, gc_obj.row_ids)
    gk_perm = gc_perm.combine(seg_id_mapping, 'combined', 0, '1s')

    # fill in the diagonal, plot 4b, 4c
    co_matrix = gk_perm.co_matrix['1s']
    co_matrix = gc_perm.fill_co_matrix_diagonal(
        seg_id_mapping, co_matrix, syn='1s')
    return {'glom_kc': gk_perm, 'glom_claw': gc_perm}


def make_freq_table(pm, ob, group):
    # join permuted and observed and assign a group type
    pm = join_freq_table(pm)
    ob = pd.DataFrame(ob, columns=['observed_pairs'])
    pm = pd.concat([pm, pd.DataFrame(pm.index.values, columns=['num_kc'],
                                     index=pm.index.values)], axis=1)
    r = ob.join(pm, how='outer')
    r = r.where(~pd.isnull(r), other=0)
    r['num_kc'] = r.index.values.astype(int)
    r.assign(type=group)
    return r


def reorder_glom_seq(glom_ids):
    # glom_ids - a list of glomerular ids. For example, ana.conn_data['glom_kc_contracted'].col_ids
    # glom_ids = ana.conn_data['glom_kc_contracted'].col_ids
    fafb_c = cc.fafb_connection()
    t1 = [str(i) for i in cc.ids_to_annotations(fafb_c, glom_ids)]
    fpath = '/Users/zhengz11/myscripts/data_results/161021-caron_equiv/161024_FAFB_glom_seq.xlsx'
    glom_seq = pd.read_excel(fpath).fafb.tolist()
    glom_seq.pop(17)
    # glom_seq += [i for i in t1 if i not in glom_seq]

    reorder_idx = [t1.index(i) for i in glom_seq]
    reorder_glom = [t1[i] for i in reorder_idx]
    reorder_glom_id = [glom_ids[i] for i in reorder_idx]
    return reorder_idx, reorder_glom, reorder_glom_id


def observed_vs_shuffle(co_matrix, gc_obj, seg_skid, num_exp=10):
    # gc_obj - a glom_claw connectivity object. For example, ana.conn_data['glom_claw_contracted']
    # seg_skid - a mapping from segments to skids for KCs. For example, ana.kc_mapping.segment_skid

    fafb_ob_f4b = freq_identical_pairs(co_matrix)
    fafb_ob_f4c = freq_non_identical_pairs(co_matrix)

    fafb_pm_f4b = []
    fafb_pm_f4c = []

    for i in range(num_exp):
        # permute and build glom claw conn
        perm = shuffle_glom_claw_conn(gc_obj, seg_skid)
        co_matrix = perm['glom_kc'].co_matrix['1s']

        fafb_pm_f4b.append(freq_identical_pairs(co_matrix))
        fafb_pm_f4c.append(freq_non_identical_pairs(co_matrix))
        print(i)

    f4b = make_freq_table(fafb_pm_f4b, fafb_ob_f4b, 'identical')
    f4c = make_freq_table(fafb_pm_f4c, fafb_ob_f4c, 'non-identical')
    return f4b, f4c


'''
whenever re-run, set diagonal to zeros
co_matrix = ana.conn_data['glom_kc_contracted'].co_matrix['5s']
co_matrix[np.diag_indices(co_matrix.shape[0])]=0

co_matrix = ana.conn_data['glom_kc_contracted'].co_matrix['5s']
co_matrix = ana.conn_data['glom_claw_contracted'].fill_co_matrix_diagonal(ana.kc_mapping.segment_skid, co_matrix, syn='5s')
'''


def group_division(total, num_group):
    r = []
    for i in range(num_group - 1):
        r.extend([i] * int(round(float(total) / num_group, 0)))
    r.extend([num_group - 1] * (total - len(r)))
    return r

def get_gk_conn_list(gk_conn):
    '''Given a glomerulus-KC connectivity matrix, generate a list (idx) of 2 elements.
    idx[0] - represent each output instance of a bouton
    idx[1] - represent each input instance of a claw. Namely, a claw since one claw only receives input from 1 bouton
    '''
    idx_list = [np.nonzero(gk_conn)]
    for i in range(2, int(np.max(gk_conn) + 1)):
        idx = np.where(gk_conn == i)
        idx = [np.repeat(idx[j], i - 1) for j in (0, 1)]
        idx_list.append(idx)
    idx = [np.concatenate([idx_list[i][j]
                           for i in range(len(idx_list))]) for j in (0, 1)]
    return idx

def shuffle_glom_kc(gk_conn):
    '''Given a glomerulus-KC connectivity matrix, shuffle the connection
    while maintain the numbers of boutons and claws and return the shuffled matrix.
    Note that ndividual claw connection is not identifiable but as numbers in the glom-KC cell
    e.g. 2 in a cell means the KC and glom connects with 2 claws'''
    idx = get_gk_conn_list(gk_conn)
    shuffled_conn = np.zeros(gk_conn.shape)
    idx[1] = np.random.permutation(idx[1])
    for i in range(len(idx[0])):
        shuffled_conn[idx[0][i], idx[1][i]] += 1
    return shuffled_conn

def shuffle_glom_kc_iterate(gk_conn, num_exp):
    '''same as shuffle_glom_kc but add num_exp'''
    idx = get_gk_conn_list(gk_conn)
    r = []
    for j in range(num_exp):
        shuffled_conn = np.zeros(gk_conn.shape)
        idx[1] = np.random.permutation(idx[1])
        for i in range(len(idx[0])):
            shuffled_conn[idx[0][i], idx[1][i]] += 1
        r.append(shuffled_conn)
    return r

def shuffle_glom_kc_w_prob(gk_conn, col_prob):
    '''Given a glomerulus-KC connectivity matrix, shuffle the connection
    while maintain the numbers CLAWS ONLY and return the shuffled matrix.
    Note that ndividual claw connection is not identifiable but as numbers in the glom-KC cell
    e.g. 2 in a cell means the KC and glom connects with 2 claws
    This one with probability of choice for each glomerulus (eacg column)'''
    sfl_conn = np.zeros(gk_conn.shape)
    num_col = sfl_conn.shape[1]
    for i in range(sfl_conn.shape[0]):
        t1 = np.random.choice(int(num_col), size=int(sum(gk_conn[i,:])), p=col_prob)
        for j in t1:
            sfl_conn[i, j] += 1
    return sfl_conn

def simulated_vs_shuffle(gk_obj, num_exp=10):
    # gk_obj - a glom_KC connectivity object. For example, ana.conn_data['glom_claw_contracted']
    co_matrix = gk_obj.co_matrix['1s']
    fafb_ob_f4b = freq_identical_pairs(co_matrix)
    fafb_ob_f4c = freq_non_identical_pairs(co_matrix)

    fafb_pm_f4b = []
    fafb_pm_f4c = []

    for i in range(num_exp):
        # permute and build glom claw conn
        perm_obj = GlomKcConnectivityMatrix(
            'shuffled', shuffle_glom_kc(gk_obj.conn['1s']))

        fafb_pm_f4b.append(freq_identical_pairs(perm_obj.co_matrix['1s']))
        fafb_pm_f4c.append(freq_non_identical_pairs(perm_obj.co_matrix['1s']))
#        print(i)

    f4b = make_freq_table(fafb_pm_f4b, fafb_ob_f4b, 'identical')
    f4c = make_freq_table(fafb_pm_f4c, fafb_ob_f4c, 'non-identical')
    return f4b, f4c


def simulated_vs_shuffle_simple(gk_conn):
    # for example, gk_conn=gk_obj.conn['1s']
    # similar to simulated_vs_shuffle but only permute once and therefore trim all sd, ci, mean, etc.
    gk_co=get_co_matrix(gk_conn)
    perm_co = get_co_matrix(shuffle_glom_kc(gk_conn))

    r_f4b = pd.DataFrame({'observed': freq_identical_pairs(
        gk_co), 'permuted': freq_identical_pairs(perm_co)})
    r_f4c = pd.DataFrame({'observed': freq_non_identical_pairs(
        gk_co), 'permuted': freq_non_identical_pairs(perm_co)})
    r_f4b = r_f4b.where(~pd.isnull(r_f4b), other=0)
    r_f4c = r_f4c.where(~pd.isnull(r_f4c), other=0)
    return r_f4b, r_f4c


def get_weighted_transitivity(ci_matrix, conn):
    simi = get_ci_similarity(ci_matrix, conn)
    return get_local_transitivity(simi, 'weighted')


def get_local_transitivity(ci_matrix, graph_type='binary'):
    ci_matrix = ci_matrix.copy()
    if graph_type == 'binary':
        ci_matrix[ci_matrix > 0] = 1
        g = Graph.Adjacency(ci_matrix.tolist(), mode=ADJ_UPPER)
    elif graph_type == 'weighted':
        g = Graph.Weighted_Adjacency(ci_matrix.tolist(), mode=ADJ_UPPER)
    return g.transitivity_local_undirected(mode="zero", weights="weight")


def pick_random_neighbour(conn_row, geom_row, threshold=3000):
    if np.count_nonzero(conn_row) > 0:
        data = np.zeros((geom_row.shape))
#		conn_idx = np.nonzero(conn_row)[0]
        to_sample = np.where(geom_row < threshold)[0]
    if len(to_sample) == 0:
        data = conn_row
    else:
        data[np.random.choice(to_sample, size=1)[0]] = 5
    return data


def pick_partner(conn, geom, func):
    perm_matrix = np.zeros((conn.shape))
    for i in range(conn.shape[0]):
        if np.count_nonzero(conn[i, :]) > 0:
            perm_matrix[i, :] = func(conn[i, :], geom[i, :])
    return perm_matrix


def pick_random_bouton(conn_row, geom_row):
    if np.count_nonzero(conn_row) > 0:
        data = np.zeros((geom_row.shape))
    #	conn_idx = np.nonzero(conn_row)[0]
    #	to_sample = np.setdiff1d(range(len(conn_row)), conn_idx, True)
        data[np.random.choice(list(range(len(conn_row))), size=1)[0]] = 5
    return data

def pick_next_neighbour(conn_row, geom_row):
    if np.count_nonzero(conn_row) > 0:
        data = np.zeros((geom_row.shape))
        geom_sort = np.argsort(geom_row)
        conn_idx = np.nonzero(conn_row)[0]
        to_sample = np.setdiff1d(geom_sort[:5], conn_idx, True)
        data[to_sample[0]]=5
    return data

def pick_random_from_neighbours(conn_row, geom_row, n=5):
    if np.count_nonzero(conn_row) > 0:
        data = np.zeros((geom_row.shape))
        geom_sort = np.argsort(geom_row)
        conn_idx = np.nonzero(conn_row)[0]
        to_sample = np.setdiff1d(geom_sort[:n], conn_idx, True)
        data[np.random.choice(to_sample, 1)]=5
    return data

def detect_structured_matrix(sampled_kc=200, ratio=0.2, nglom=52, ngroup=5, nkc=2000, p=1, num_exp=1000, unfilled_claw=True):
    fpath = "/Users/zhengz11/myscripts/data_results/160928-caron_equiv/160928-Caron_suppl_table.xlsx"
    suppl_tbl = pd.read_excel(fpath, 'combined')

    filled_claw_dict = {i: suppl_tbl.query("total_claws == @i")[
        "claws_filled"].tolist() for i in np.unique(suppl_tbl.total_claws)}

    sim_conn = np.zeros((nkc, nglom))

    # assign all KCs (p = 1) or a fraction of KCs to different groups
    # KCs in the same group will connect to the same group of glomeruli
    kc_idx = group_division(int(p * nkc), ngroup) + \
        [ngroup] * (nkc - int(p * nkc))

    # assign glom to different groups
    assigned_glom = np.array(group_division(nglom, ngroup))

    for i in range(nkc):
        if kc_idx[i] == ngroup:
                # pick claw numbers from the claw table
                # ask KCs in the last group to just randomly pick any glomeruli
            col_idx = np.random.randint(
                nglom - 1, size=np.random.choice(suppl_tbl.total_claws, 1))
            for j in col_idx:
                sim_conn[i, j] += 1
        else:
                # Else, ask a fraction of the claws to exclusively connect to gloms of a group.
                # The remaining claws would randomly pick any glom.
                # Same as above, pick claw number from claw table
            claws = np.random.choice(suppl_tbl.total_claws, 1)[0]
            num_claws = int(round(claws * ratio))
            pick = np.random.choice(np.where(assigned_glom == kc_idx[i])[
                                    0], num_claws)
            rd = np.random.randint(nglom - 1, size=claws - num_claws)
            for j in np.concatenate((pick, rd)):
                sim_conn[i, j] += 1
        # sim_conn is the nglom(default 53) total KCs(default 2000) partially structured connectivity matrix

    p_4b = []
    p_4c = []
    for k in range(num_exp):
        # randomly select some KCs (sampled_kc) from the partially structured connectivity matrix
        select_kc = np.random.randint(nkc, size=sampled_kc)
        select_gk = sim_conn[select_kc, :]

        if unfilled_claw:
            new_select_gk = np.zeros(select_gk.shape)
            for row in range(select_gk.shape[0]):
                t6 = select_gk[row, :]
                conn_list = reduce(
                    lambda x, y: x + y, [int(t6[i]) * [i] for i in np.nonzero(t6)[0]])
                claws = len(conn_list)
                cols = np.random.choice(conn_list, np.random.choice(
                    filled_claw_dict[claws]), replace=False)
                for col in cols:
                    new_select_gk[row, col] += 1
            select_gk = new_select_gk
        print(k)

        t1 = simulated_vs_shuffle_simple(select_gk)
        p_4b.append(chi2_contingency(t1[0].transpose())[1])
        p_4c.append(chi2_contingency(t1[1].transpose())[1])
    r = pd.DataFrame({'p_4b': p_4b, 'p_4c': p_4c})
    # Shuffle the sampled matrix and compare 'observed' with shuffle data
    # for equivalent of F4b, F4c results
    return r

def get_random_conn(num_row, claw_set, prob):
    # given num_row, claw_set,
    # for example, claw_set=claw_nums, prob=mz_prob
    num_col = 2
    num_row = int(num_row)
    conn = np.zeros((num_row,num_col))
    for i in range(int(num_row)):
        for j in np.random.choice([0,1],size=int(np.random.choice(claw_set)),p=prob).astype(int):
            conn[i,j] += 1
    return  conn
