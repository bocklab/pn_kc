
# 180328 from build_connectivity_v3 in mushroom_py/v9

import os
import datetime
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from scipy.spatial import ConvexHull
from sqlalchemy import create_engine

from . import connect as cc
from . import neurogenesis


class DataMatrix(object):
    '''a 2D array to store pair-wise relationship
    cols - boutons, PNs, glomeruli
    rows - claws, KCs
    relationship - # syn, centroid dist, ect.'''
    def __init__(self, name, conn, col_ids, row_ids):
        self.name = name
        self.conn = conn
        self.col_ids = self._tolist(col_ids)
        self.row_ids = self._tolist(row_ids)

    def _tolist(self, ids):
        if not isinstance(ids, list):
            ids = [i for i in ids]
        return ids

    def retrieve(self, col_row_pairs):
        # col_row_pairs is a list in which each element is a two-element tuples
        # [(pn_skid_1, kc_skid_1), (pn_skid_2, kc_skid_2), etc.]
        r = [float(self.conn[self.row_ids.index(i[1]), self.col_ids.index(i[0])])
             for i in col_row_pairs]
        return r

    def to_sql(self, db_path, syn='1s'):
        conn = self.conn
        if isinstance(self.conn, dict):
            conn = self.conn[syn]
        tbl = pd.DataFrame(conn, columns=self.col_ids, index=self.row_ids)
        tbl.to_sql(self.name, create_engine(db_path))

    def subset(self, col_subsets=None, row_subsets=None, name=''):
        conn = self.conn.copy()
        if col_subsets is not None:
            col_ids = [self.col_ids[i] for i in col_subsets]
            conn = conn[:,col_subsets]
        else:
            col_ids = list(self.col_ids)
        if row_subsets is not None:
            row_ids = [self.row_ids[i] for i in row_subsets]
            conn = conn[row_subsets,:]
        else:
            row_ids = list(self.row_ids)
        return DataMatrix(name, conn, col_ids, row_ids)

class PairMatrix(DataMatrix):
    '''for now, col_ids need to be the same as row_ids'''
    def __init__(self, name, conn, ids):
        DataMatrix.__init__(self, name, conn, ids, ids)

    def reorder(self, reorder_idx, return_new=False):
        ids = [self.col_ids[i] for i in reorder_idx]
        # self.row_ids = self.col_ids
        conn = self.conn[np.ix_(reorder_idx, reorder_idx)]
        if return_new is not False:
            return PairMatrix(str(return_new), conn, ids)
        else:
            self.col_ids = ids
            self.row_ids = ids
            self.conn = conn

class GeoMatrix(DataMatrix):
    '''pair-wise geometric relationship'''
    def __init__(self, name, constructor_func,
                 col_ids, row_ids, col_neurons, row_neurons, op='neuron'):
        conn = constructor_func(col_ids, row_ids, col_neurons, row_neurons, op)
        DataMatrix.__init__(self, name, conn, col_ids, row_ids)

class PwMatrix(DataMatrix):
    ''''''
    def __init__(self, df, name):
        self._process_df(df)
        self.name = name

    def _process_df(self, df):
        ids = np.unique(df.iloc[:,:2]).tolist()
        self.conn = np.zeros((len(ids),len(ids)))
        for row in df.itertuples(index=False):
            self.conn[ids.index(row[0]),ids.index(row[1])]=row[2]
            self.conn[ids.index(row[1]),ids.index(row[0])]=row[2]
        ids = [int(i) for i in ids]
        self.col_ids = ids
        self.row_ids = ids

    def get_data(self, skids):
        data = []
        for i in pair_generator(skids):
            data.append(self.conn[self.col_ids.index(i[0]), self.col_ids.index(i[1])])
        return data

class ConnectivityMatrix(DataMatrix):
    '''pair-wise connectivity data'''
    def __init__(self, name, conn, col_ids, row_ids):
        DataMatrix.__init__(self, name, {'1s': conn}, col_ids, row_ids)
        self._process_conn_ci()
        self.ci_path_counts = None
        self._process_ci_source_matrix()
        self._process_conn_co()

    def _process_conn_ci(self, syn_threshold=3):
        # threshold and store connectivity matrix (e.g. PN - KC) with threshold of syn_threshold synapses
        # note that although the code say "5s", synaptic threshold has been changed to 3 (default) to be consistent with T1+ protocol. Todo: fix the "5s" keys.
        conn_5s = np.copy(self.conn['1s'])
        conn_5s[conn_5s < syn_threshold] = 0
        self.conn['5s'] = conn_5s
        self.ci_matrix = {'1s':get_ci_matrix(self.conn['1s']),
                          '5s':get_ci_matrix(self.conn['5s'])}

    def _process_conn_co(self):
        self.co_matrix = {'1s':get_ci_matrix(self.conn['1s'].transpose()),
                          '5s':get_ci_matrix(self.conn['5s'].transpose())}

    def _process_ci_source_matrix(self, syn='5s'):
        _, self.ci_source_matrix = get_ci_matrix(self.conn[syn], True)

    def _get_ci_similarity(self, syn='1s'):
        self.similarity = get_ci_similarity(self.ci_matrix[syn], self.conn['1s'])

    def skid_to_indices(self, kc_pair_skids):
    # given a list of pair tuple KC skeleton ids, return indices in rows (columns) in ci_matrix
        idx = [(self.row_ids.index(x), self.row_ids.index(y)) for x,y in kc_pair_skids]
        return idx

    def get_ci_list(self, kc_pair_skids, syn='5s'):
        # syn should be '5s' or '1s'
        kc_pair_indices = self.skid_to_indices(kc_pair_skids)
        return [self.ci_matrix[syn][i,j] for i,j in kc_pair_indices]

    def get_ci_source_from_pairs(self, kc_pair_skids, syn='5s'):
        kc_pair_indices = self.skid_to_indices(kc_pair_skids)
        r = []
        for i,j in kc_pair_indices:
            if self.ci_source_matrix[i,j] is not None:
                r.append([self.col_ids[x] for x in self.ci_source_matrix[i,j]])
            else:
                r.append([])
        return r

    def get_all_connected_pairs(self, syn='5s'):
    # the returned dictionary - col_ids and row_ids are two lists of the same length in which
    # col_ids[n] is connected to row_ids[n]
        nz = np.nonzero(self.conn[syn])
        return {'col_ids' : self.col_ids[nz[1]], 'row_ids' : self.row_ids[nz[0]]}

    def get_col_inputs_to_row(self, claw_id, syn='5s'):
        # claw_id is a single string id
        nz = np.nonzero(self.conn[syn][(self.row_ids.index(claw_id),)])
        return [self.col_ids[i] for i in nz[0]]

    def get_connections_of_pairs(self, col_row_pairs, syn='5s'):
        # col_row_pairs is a list in which each element is a two-element tuples
        # [(pn_skid_1, kc_skid_1), (pn_skid_2, kc_skid_2), etc.]
        results = []
        for i in [(self.row_ids.index(r), self.col_ids.index(c)) for c,r in col_row_pairs]:
            results.append(self.conn[syn][i])
        return results

    def count_ci_all_to_all_paths(self, max_len, syn='5s'):
        r = {}
        r['syn_threshold'] = syn
        r['max_len'] = max_len
        r['path_count_matrices'] = {}
        for i in range(1,max_len + 1):
            r['path_count_matrices'].update({i:np.linalg.matrix_power(self.ci_matrix[syn],i)})
        self.ci_path_counts = r

    def get_ci_paths_from_pairs(self, kc_pair_skids, max_len=5):
    # given a list of pair tuple KC skeleton ids (kc_pair_skids)
    # e.g. [(6486, 6618), (6486, 7118)]
    # the sequence in a tuple shouldn't matter because the path_count_matrix is symmetric
        if (self.ci_path_counts is None) or (len(self.ci_path_counts['path_count_matrices']) < 5):
            self.count_ci_all_to_all_paths(max_len)
        pc_matrices = self.ci_path_counts['path_count_matrices']
        skid_idx = [(self.row_ids.index(x), self.row_ids.index(y)) for x,y in kc_pair_skids]
        r = []
        for i in range(len(kc_pair_skids)):
            r.append([pc_matrices[j][skid_idx[i]] for j in range(1,max_len+1)])
        # n - number of KC pairs, m - max_len
        # [pair_1, pair_2, pair_3, ... pair_n] in which
        # pair_n = [paths_of_length_1,paths_of_length_2, ..., paths_of_length_m]
        return r

    def combine(self, id_map, output_name='combined', axis=1, conn='1s'):
        '''
        row - axis 0, column - axis 1'''
        cols = id_map.iloc[:,1]
        output_ids = pd.unique(cols)

        if isinstance(conn, str):
            conn = self.conn[conn]

        if axis == 1:
            axis_ids = self.col_ids
        elif axis == 0:
            conn = np.transpose(conn)
            axis_ids = self.row_ids

        combined = np.zeros((conn.shape[0],len(output_ids)))
        for i,e in enumerate(output_ids):
            idx = id_map[cols == e].iloc[:,0].tolist()
            combined[:,i] = np.sum(conn[:,np.array(find_elements(axis_ids, idx))],1)
        if axis == 1:
            data = ConnectivityMatrix(output_name, combined, output_ids, self.row_ids)
        elif axis == 0:
            combined = np.transpose(combined)
            data = ConnectivityMatrix(output_name, combined, self.col_ids, output_ids)
        return data

    def combine_new(self, id_map, output_name='combined', axis=1, conn='1s'):
        '''
        row - axis 0, column - axis 1'''
        cols = id_map.iloc[:,1]
        output_ids = pd.unique(cols)

        if isinstance(conn, str):
            conn = self.conn[conn]

        if axis == 1:
            axis_ids = self.col_ids
        elif axis == 0:
            conn = np.transpose(conn)
            axis_ids = self.row_ids

        combined = np.zeros((conn.shape[0],len(output_ids)))
        for i,e in enumerate(output_ids):
            idx = id_map[cols == e].iloc[:,0].tolist()
            axis_idx = find_elements(axis_ids, idx)
            if len(axis_idx) > 0:
                combined[:,i] = np.sum(conn[:,np.array(axis_idx)],1)
        if axis == 1:
            data = ConnectivityMatrix(output_name, combined, output_ids, self.row_ids)
        elif axis == 0:
            combined = np.transpose(combined)
            data = ConnectivityMatrix(output_name, combined, self.col_ids, output_ids)
        return data

    def fill_co_matrix_diagonal(self, id_map, co_matrix, syn='5s'):
    # related to Caron et al.
    #    e.g. co_matrix = ana.conn_data['glom_kc_contracted'].co_matrix['5s']

        t1 = self.conn[syn].copy()
        t1[t1>0] = 1

        t2 = self.combine(id_map=id_map, axis=0, conn=t1)

        for i in np.where(t2.conn['1s']>1)[1]:
            co_matrix[i,i] += 1
        return co_matrix

    def permutate(self, perm_func, geom_matrix=None):
        conn = self.conn['1s']
        if geom_matrix is not None:
            self.perm = perm_func(conn, geom_matrix)
        else:
            self.perm = perm_func(conn)

class ProbMatrix(ConnectivityMatrix):
    '''pair-wise common input probability'''
    def __init__(self, name, conn, col_ids, row_ids):
        DataMatrix.__init__(self, name, conn, col_ids, row_ids)

    def retrieve_from_pairs(self, kc_pairs):
        # kc_pairs is a list in which each element is a two-element tuples
        # [(pn_skid_1, kc_skid_1), (pn_skid_2, kc_skid_2), etc.]
        if not isinstance(self.row_ids,list):
            row_ids = self.row_ids.tolist()
        else:
            row_ids = self.row_ids
        r = [self.prob_matrix[row_ids.index(i[0]), row_ids.index(i[1])]
             for i in kc_pairs]
        return r

    def combine(self, id_map, output_name='combined', axis=1):
        return ConnectivityMatrix.combine(self, id_map, output_name, axis, self.conn)

class Mapping(object):
    """a class for the mapping in between Boutons, PNs, Glomeruli
    also between claws and KCs"""
    def __init__(self, connection, type_meta_anno, skids, segment_ids=None):
        self._mapping_skids_to_types(connection, type_meta_anno, skids)
        self._mapping_boutons_to_skids(segment_ids)
        self._mapping_boutons_to_types()

    def _mapping_skids_to_types(self, connection, type_meta_anno, skids):
        if type_meta_anno and skids is not None:
            type_ids = cc.retrieve_annotated_annotations(connection, type_meta_anno, 'id')
            self.long_skids = []
            skid_list = []
            type_list = []
            for i in skids:
                annos = cc.retrieve_annotations_for_skeleton(connection, i)
                anno_ids = [int(j) for j in list(annos.keys())]
                sk_anno_ids = [j for j in anno_ids if j in type_ids]
                self.long_skids.append('.'.join([str(i)] + [str(k) for k in sk_anno_ids]))
                skid_list.extend([i] * len(sk_anno_ids))
                type_list.extend(sk_anno_ids)
            self.skid_type = pd.DataFrame({'skids': skid_list, 'type_ids': type_list},
                                           columns=['skids', 'type_ids'])
        else:
            self.skid_type = None

    def _mapping_boutons_to_skids(self, segment_ids):
        if segment_ids is not None:
            skids = [int(i.split('.')[0]) for i in segment_ids]
            self.segment_skid = pd.DataFrame({'segments': segment_ids, 'skids': skids},
                                              columns=['segments', 'skids'])
        else:
            self.segment_skid = None

    def _mapping_boutons_to_types(self):
        segsk = self.segment_skid
        types = self.skid_type.type_ids
        seg_list = []
        type_list = []
        for i in pd.unique(types):
            per_type = pd.concat([segsk.loc[segsk.skids==j,'segments']
                for j in self.skid_type.loc[types==i,'skids'].tolist()])
            seg_list.extend(per_type.tolist())
            type_list.extend([i] * per_type.shape[0])
        self.segment_type = pd.DataFrame({'segments':seg_list,'types':type_list},
                                        columns=['segments', 'types'])

    def types_to_segments(self, types):
        '''
        Given a list of glom, return boutons from that glom.
        Note that gloms should be a list
        '''
        ids = self.skid_type.query('type_ids == @types').skids.tolist()
        return self.segment_skid.query('skids == @ids').segments.tolist()

    def types_to_skids(self, types):
        '''
        given a list of glom, return skids (e.g. PNs) of those glom categories
        '''
        return self.skid_type.query('type_ids == @types').skids.tolist()

    def skids_to_types(self, skids):
        '''
        same as above, given skids return types of those skids (e.g. glom)
        '''
        return self.skid_type.query('skids == @skids').type_ids.tolist()

    def to_sql(self, db_path, prefix=''):
        self.segment_skid.to_sql(prefix + 'segment_skid', create_engine(db_path))
        self.skid_type.to_sql(prefix + 'skid_type', create_engine(db_path))

def centroid_matrix_constructor(col_ids, row_ids, col_neurons, row_neurons, op='neuron'):
    conn=np.zeros((len(row_ids),len(col_ids)))
    it = np.nditer(conn, flags=['multi_index'], op_flags=['readwrite'])
    for x in it:
        row_id = row_ids[it.multi_index[0]]
        col_id = col_ids[it.multi_index[1]]
        if op == 'neuron':
            row_ctd = row_neurons[row_id].segments.centroid_of_all_segments()
            col_ctd = col_neurons[col_id].segments.centroid_of_all_segments()
            x[...] = np.linalg.norm(row_ctd - col_ctd)
        else:
            row_ctd = row_neurons[int(row_id.split('.')[0])].segments.id_to_ctd(row_id)
            col_ctd = col_neurons[int(col_id.split('.')[0])].segments.id_to_ctd(col_id)
            x[...] = np.linalg.norm(row_ctd - col_ctd)
    return conn


def ConvexHullVol_matrix_constructor(col_ids, row_ids, col_neurons, row_neurons, op='neuron'):
    conn=np.zeros((len(row_ids),len(col_ids)))
    it = np.nditer(conn, flags=['multi_index'], op_flags=['readwrite'])
    for x in it:
        row_id = row_ids[it.multi_index[0]]
        col_id = col_ids[it.multi_index[1]]
        if op == 'neuron':
            row_segs = row_neurons[row_id].segments
            col_segs = col_neurons[col_id].segments
            row_chv = row_segs.ConvexHullVol_of_all_segments()
            col_chv = col_segs.ConvexHullVol_of_all_segments()
            combined = np.concatenate((row_segs.all_segments_nodes(),
                                       col_segs.all_segments_nodes()),axis=0)
        else:
            row_segs = row_neurons[int(row_id.split('.')[0])].segments
            col_segs = col_neurons[int(col_id.split('.')[0])].segments
            row_chv = row_segs.id_to_ConvexHullVol(row_id)
            col_chv = col_segs.id_to_ConvexHullVol(col_id)
            combined = np.concatenate((row_segs.id_to_xyz(row_id),
                                       col_segs.id_to_xyz(col_id)),axis=0)
        intersect = row_chv + col_chv - neurogenesis.convex_hull_vol(combined)
        x[...] = intersect
    return conn

def find_elements(a_list, elements):
    '''given a list and some elements (a list), find the indices of the matched elements in the list'''
    return [i for i,e in enumerate(a_list) if e in elements]

def merge_dicts(args):
    ''' Given any number of dicts in a list, shallow copy and merge into a new dict '''
    result = {}
    for i in args:
        result.update(i)
    return result

def neurons_to_segment_ids(skids, neuron_dict):
    '''given skids and a neuron dict (keyed by skids), return a list of all segment ids'''
    sids = []
    for i in skids:
        sids.extend(neuron_dict[i].segments.ids)
    sids = [_f for _f in sids if _f]
    return sids

def get_connectivity_data(connection, col_neurons, row_neurons, c_skids, r_skids, btn_ids, claw_ids, info_file):
    conn_data = {}
    r = cc.get_pre_post_info(connection, c_skids, r_skids, info_file)
    # 0 "connector_id",
    # 1 "connector_xyz",
    # 2 "pre_node",
    # 3 "pre_skid",
    # 4 "pre_conf",
    # 5 "pre_creator",
    # 6 "pre_node_xyz",
    # 7 "post_node",
    # 8 "post_skid",
    # 9 "post_conf",
    # 10 "post_creator",
    # 11 "post_node_xyz"

    # bouton -> claw connectivity
    conn=np.zeros((len(claw_ids),len(btn_ids)))
    for i in r:
        sub_b = col_neurons[i[3]].node_to_segment(i[2])
        if sub_b:
            sub_c = row_neurons[i[8]].node_to_segment(i[7])
            if sub_c:
                conn[claw_ids.index(sub_c), btn_ids.index(sub_b)] += 1
    conn_data['bouton_claw'] = ConnectivityMatrix('bouton_claw', conn, btn_ids, claw_ids)

    # bouton -> KC connectivity
    bk_conn=np.zeros((len(r_skids),len(btn_ids)))
    for i in r:
        subs = col_neurons[i[3]].node_to_segment(i[2])
        if subs in btn_ids:
            bk_conn[r_skids.index(i[8]), btn_ids.index(subs)] += 1
    conn_data['bouton_kc'] = ConnectivityMatrix('bouton_kc', bk_conn, btn_ids, r_skids)

    # PN -> KC connectivity
    pn_kc_conn=np.zeros((len(r_skids),len(c_skids)))
    for i in r:
        pn_kc_conn[r_skids.index(i[8]), c_skids.index(i[3])] += 1
    conn_data['pn_kc'] = ConnectivityMatrix('pn_kc', pn_kc_conn, c_skids, r_skids)

"""
    # 020323 suppress here as KC -> PN connectivity is not relevant to the paper

    # KC -> PN connectivity
    r = cc.get_pre_post_info(connection, r_skids, c_skids)
    kc_pn = np.zeros((len(r_skids),len(c_skids)))
    for i in r:
        kc_pn[r_skids.index(i[3]), c_skids.index(i[8])] += 1
    conn_data['kc_pn'] = ConnectivityMatrix('kc_pn', kc_pn, c_skids, r_skids)

    # KC -> bouton connectivity
    kc_btn=np.zeros((len(r_skids),len(btn_ids)))
    for i in r:
        subs = col_neurons[i[8]].node_to_segment(i[7])
        if subs in btn_ids:
            kc_btn[r_skids.index(i[3]), btn_ids.index(subs)] += 1
    conn_data['kc_bouton'] = ConnectivityMatrix('kc_bouton', kc_btn, btn_ids, r_skids)
"""

    return conn_data


def get_ci_matrix(conn, ci_source=False):
    num_row,num_col=conn.shape
    ci_matrix = np.zeros((num_row,num_row))
    ci_source_matrix = np.empty((num_row,num_row),dtype=object)
    for i in range(num_col):
        nnz_idx = np.nonzero(conn[:,i])[0]
        if len(nnz_idx) >= 2:
            for j in [(x,y) for x in nnz_idx for y in nnz_idx if x != y]:
                ci_matrix[j] += 1
                if ci_source:
                    if ci_source_matrix[j] is None:
                        ci_source_matrix[j] = []
                    ci_source_matrix[j].append(i)
    if ci_source:
        return ci_matrix, ci_source_matrix
    else:
        return ci_matrix

def pair_generator(skids):
    return ((i,j) for i in skids for j in skids if (i!=j and i<j))

def contract_columns(conn, id_map):
    col_ids = pd.unique(id_map)
    contracted = np.zeros((conn.shape[0],len(col_ids)))
    for i in range(len(col_ids)):
        contracted[:,i] = np.sum(conn[:,np.array([col_ids[i] == x for x in id_map])],1)
    return (contracted, col_ids)

def process_pn_glom_table(fname):
    # process downloaded PN glom table
    raw_table = pd.read_excel(fname)
    pn_glom_map=raw_table.loc[:125,('skeleton id','cleaned_glom_type')]
    pn_glom_map.iloc[:,0] = [int(x) for x in pn_glom_map.iloc[:,0]]
    return pn_glom_map

def process_skid_pair_df(pair_df):
    '''process the skeleton pair dataframe into a list
    in which each element is a tuple of KC skeleton pairs
    note that each skid need to be an int'''
    kc1=[int(i) for i in pair_df.iloc[:,0].values.tolist()]
    kc2=[int(i) for i in pair_df.iloc[:,1].values.tolist()]
    return list(zip(kc1,kc2))

def process_skid_pair_csv(csv_fname):
    '''process the skeleton pair csv into a list
    in which each element is a tuple of KC skeleton pairs'''
    skid_pairs=pd.read_csv(csv_fname, header=None, delim_whitespace=True)
    kc1=skid_pairs.iloc[:,0].values.tolist()
    kc2=skid_pairs.iloc[:,1].values.tolist()
    return list(zip(kc1,kc2))

def get_ci_similarity(ci_matrix, conn):
    '''currently only work if conn is a glom_kc connnectivity matrix in which
    1 means 1 claw connects with a bouton. Two in the cell mean 2 claws have
    connections with their boutons.
    '''
    similarity = np.zeros(ci_matrix.shape)
    all_idx = list(range(ci_matrix.shape[0]))
    for col,row in [(x, y) for x in all_idx for y in all_idx if x>y]:
        n = ci_matrix[col, row]*2 / np.sum(conn[(col,row),:])
        similarity[col, row] = n
        similarity[row, col] = n
    return similarity

def segid_to_skid(segid):
    # segid should be a list
    return [int(i.split(".")[0]) for i in segid]

def get_gk_conn_in_claw_unit(gc_obj, mapping_obj):
    '''A glom kc connectivity matrix in which 1 means 1 claw connects with a glom.
    For example, a cell with 2 mean the KC (row) uses 2 claws to connect with the glom (col)
    gc_obj = ana.conn_data['glom_claw_contracted']
    mapping_obj = ana.kc_mapping.segment_skid

    '''

    conn = gc_obj.conn['5s'].copy()
    conn[conn>0]=1

    gc_name = 'glom_claw_in_claw_unit'
    gk_name = 'glom_kc_in_claw_unit'
    conn_data = {}
    conn_data[gc_name] = ConnectivityMatrix(gc_name, conn, gc_obj.col_ids, gc_obj.row_ids)
    conn_data[gk_name] = conn_data[gc_name].combine(mapping_obj, gk_name, axis=0)
    return conn_data
