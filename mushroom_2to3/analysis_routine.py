
# 180328 from analysis_routine_v1 in mushroom_py/v9

from . import neurogenesis
from .build_connectivity import *


class Analysis(object):
    '''An analysis class
    conn - connectivity matrix
    col_ids - id of each column of the connectivity matrix
    row_id - id of each row of the connectivity matrix

    col_ids and row_ids contains skeleton id that can be mapped to the neurons in
    col_neuron and row_neuron for further information of the neuron (name, synapses, tag, etc.)

    currently doesn't do claws of KCs in rows (can be added later)'''

    def __init__(self):
        self.conn_data = {}
        self.col_neurons = {}
        self.row_neurons = {}
        self.geo_data = {}

    @classmethod
    def init_connectivity(cls, connection, c_skids, r_skids):
        c=cls()

        c.col_neurons = neurogenesis.init_from_skid_list(connection, c_skids)
        c.row_neurons = neurogenesis.init_from_skid_list(connection, r_skids)

        btn_ids = neurons_to_segment_ids(c_skids, c.col_neurons)
        claw_ids = neurons_to_segment_ids(r_skids, c.row_neurons)

        c.conn_data = get_connectivity_data(connection, c.col_neurons, c.row_neurons,
                                            c_skids, r_skids, btn_ids, claw_ids)

        c.geo_data['centroids'] = GeoMatrix('segment_centroids', centroid_matrix_constructor,
                                            btn_ids, claw_ids, c.col_neurons, c.row_neurons, 'segments')

#       c.geo_data['intersect_volume'] = GeoMatrix('segment_ConvexHullVols', ConvexHullVol_matrix_constructor,
#                                            btn_ids, claw_ids, c.col_neurons, c.row_neurons, 'segments')

        c.combine_columns(connection, c_skids, r_skids, btn_ids, claw_ids)
        c.get_glom_kc_diagonal()
        c.get_conn_in_claw_unit()
        return c

    def SegmentId_to_data(self, segment_id):
        skid = int(segment_id.split('.')[0])
        if skid in self.col_neurons:
            r = self.col_neurons[skid].segments
        else:
            r = self.row_neurons[skid].segments
        return r

    def combine_columns(self, connection, c_skids, r_skids, btn_ids, claw_ids):
        self.pn_mapping = Mapping(connection, 'glom_class', c_skids, btn_ids)
        self.kc_mapping = Mapping(connection, 'kc_class', r_skids, claw_ids)
        self.conn_data['pn_kc_contracted'] = self.conn_data['bouton_kc'].combine(self.pn_mapping.segment_skid, 'pn_kc_contracted')
        self.conn_data['glom_kc_contracted'] = self.conn_data['pn_kc'].combine(self.pn_mapping.skid_type, 'glom_kc_contracted')

        self.conn_data['pn_claw_contracted'] = self.conn_data['bouton_claw'].combine(self.pn_mapping.segment_skid, 'pn_claw_contracted')

        self.conn_data['glom_claw_contracted'] = self.conn_data['pn_claw_contracted'].combine(self.pn_mapping.skid_type, 'glom_claw_contracted')

    def get_glom_kc_diagonal(self):
        # given an analysis object, process glom_claw and get diagonal in glom_kc matrix as in Caron et al.
        co_matrix = self.conn_data['glom_kc_contracted'].co_matrix['5s']
        co_matrix = self.conn_data['glom_claw_contracted'].fill_co_matrix_diagonal(self.kc_mapping.segment_skid, co_matrix, syn='5s')

    def get_conn_in_claw_unit(self):
        '''A glom kc connectivity matrix in which 1 means 1 claw connects with a glom.
        For example, a cell with 2 mean the KC (row) uses 2 claws to connect with the glom (col)'''

        gc = self.conn_data['glom_claw_contracted']
        conn = gc.conn['5s'].copy()
        conn[conn>0]=1

        gc_name = 'glom_claw_in_claw_unit'
        gk_name = 'glom_kc_in_claw_unit'
        self.conn_data[gc_name] = ConnectivityMatrix(gc_name, conn, gc.col_ids, gc.row_ids)
        self.conn_data[gk_name] = self.conn_data[gc_name].combine(self.kc_mapping.segment_skid, gk_name, axis=0)
        self.conn_data[gk_name]._get_ci_similarity()


    def skid_to_name(self, skids):
        '''given skids return neuron name,
        skids - a list of skeleton ids'''
        name = []
        for i in skids:
            if i in self.col_neurons:
                name.append(self.col_neurons[i].name[0])
            elif i in self.row_neurons:
                name.append(self.row_neurons[i].name[0])
            else:
                print(("couldn't find name for %s" % i))
        return name
