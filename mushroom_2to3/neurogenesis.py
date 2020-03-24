
import time
import numpy as np

from scipy.spatial import ConvexHull
from scipy.spatial import distance

from igraph import *

from . import connect_path as cc


url = cc.fafbServer[1] + '/?pid=1&xp=%s&yp=%s&zp=%s&tool=tracingtool&sid0=5&s0=1'

class Neuron(object):
    '''
    represent a Neuron
    graph - node-by-node skeleton structure
    group - Separated bouton (or claw) and non-bouton regions.
    reduced_graph - bouton graph in which each bouton(claw) or non-bouton region is represented by one node.
    '''
    def __init__(self, skid):
            self.graph = None
            self.skid = skid
            self.name = None
            self.tags = None
            self.synapse = None
            self.root = None
            self.compact_sk = None
            self.segments = NeuronSegments()

    @classmethod
    def init_from_json(cls, connection, skid):
        neuron=cls(skid)
        neuron.compact_sk = cc.get_compact_sk(connection, skid)
        # compact_skeleton in JSON, [[nodes], [connectors], {nodeID: [tags]}]
        # in [nodes], id, parent_id, user_id, location_x, location_y, location_z, radius, confidence
        # in [connectors],
        # from-treenode (Not neccessarily pre-synaptic!! The treenode in the current skeleton)
        # to-treenode
        # 1 for post-synaptic and 0 for pre-synaptic connections
        # x, y, z
        neuron._read_compact_sk(neuron.compact_sk)
        neuron.name = cc.skid_to_name(connection, [skid])
        neuron.root=cc.get_root_node(connection, skid)['root_id']

        if 'bouton border' in neuron.tags:
            group_border = 'bouton border'
        elif 'claw border' in neuron.tags:
            group_border = 'claw border'
        else:
            print(("couldn't find bouton or claw tags for %s" % skid))
        neuron._build_group(group_border)
        return neuron

    def incident_to_tagged_nodes(self, group_border):
        '''given a tag name,
        find all nodes associated with the tag,
        return incident edges of each of the node
        '''
        # all nodes associated with the tag
        nodes = self.graph.vs.select(name_in=self.tags[group_border]).indices

        # all incident edges to the tagged nodes
        incident_edges = [self.graph.incident(x, mode='IN') for x in nodes]
        incident_edges = [x[0] for x in incident_edges]
        return nodes, incident_edges

    def _read_compact_sk(self, compact_sk):
        '''parse the CATMAID compact skeleton json'''

        # construct a igraph graph from edge list in the compact skeleton json
        el = [(x[1],x[0]) for x in compact_sk[0]]
        self.graph = Graph.TupleList(el, directed=True)

        # store tag and synapse info
        self.tags = compact_sk[2]
        self.synapse = compact_sk[1]

    def _build_group(self, group_border):
        ''' main function to process skeleton and build reduced_graph (bouton or claw graph),
        group_border should be a string, e.g. 'bouton border',
        'claw border' works as well
        '''

        nodes, incident_edges = self.incident_to_tagged_nodes(group_border)

        # delete incident edges in a copy graph,
        # use the 'components' function to return each connected component as a cluster
        # http://igraph.org/python/doc/igraph.Graph-class.html#components
        temp_copy = self.graph.copy()
        temp_copy.delete_edges(incident_edges)

        # vertex_clusters - a VertexClustering object in igraph
        # http://igraph.org/python/doc/igraph.clustering.VertexClustering-class.html
        vertex_clusters = temp_copy.components(mode='weak')

        # use membership in vertex_clusters to contract the graph into bouton(or claw) graph
        reduced_graph = self.graph.copy()
        reduced_graph.contract_vertices(vertex_clusters.membership, combine_attrs=None)
        reduced_graph.simplify(multiple=False)

        # get root node index and figure out root_node containing components
        root_idx = self.graph.vs.select(name_eq=self.root).indices[0]
        root_component = vertex_clusters.membership[root_idx]

        # shortest paths of all other components to the root_component in reduced_graph
        len_to_root = reduced_graph.get_all_shortest_paths(root_component)

        # length of paths of all other components to root_comnponent in reduced_graph
        # root_component distance = 0
        dist_to_root = [len(x)-1 for x in len_to_root]

        # boolean list of whether a component is bouton (claw) or not
        subs_b = [1-len(x)%2 for x in len_to_root]

        # indices of all nodes in the reduced graph
        mg_idx = reduced_graph.vs.indices

        # generate an id for every component (bouton or not bouton)
        all_component_ids = ['%s.%s.%s.%s' % (x) for x in
                             zip([self.skid for x in mg_idx],subs_b,mg_idx,dist_to_root)]

        # set attributes ['group'] of each nodes to their group indices
        self.graph.vs['segment'] = vertex_clusters.membership

        # id for bouton (or claw) component, '' for non-bouton
        sub_id = [x*y for x,y in zip(all_component_ids,subs_b)]

        # fix the outgoing border node attribute
        for i in nodes:
            v = vertex_clusters.membership[i]
            if not sub_id[v]:
                if reduced_graph.neighbors(v, mode=IN):
                    self.graph.vs[i]['segment'] = reduced_graph.neighbors(v, mode=IN)[0]

        # along with the for loop below, generate a segment_type mapping
        # non-bouton (or claw) areas - 'skeleton', rout containing area - 'root', then 'bouton' or 'claw'
        seg_type = ['skeleton'] * len(all_component_ids)
        seg_type[root_component] = 'root'

        all_ids = [i[0] for i in self.compact_sk[0]]
        for idx,e in enumerate(sub_id):
            if e:
                # an example of compact_sk[0][0]
                # [3425104, 3425103, 26, 388236.0, 160616.0, 174405.0, -1.0, 5]
                nodes_id = self.graph.vs(segment_eq=idx)['name']
                nodes_xyz = np.array([self.compact_sk[0][all_ids.index(j)][3:6] for j in nodes_id])
                compact_connector = [i for i in self.compact_sk[1] if i[0] in nodes_id]
                self.segments.add(idx, e, nodes_id, nodes_xyz, compact_connector)
                seg_type[idx] = group_border.split(' ')[0]

        # sub_id - skid. group_boolean . ids in reduced_graph . distance to root
        # all component ids = dict(zip(reduced_graph.vs.indices,all_component_ids))
        self.segments.add_attr(reduced_graph=reduced_graph,
                               ids=sub_id, # only bouton (or claw) component has an id, others are ''
                               all_types=dict(list(zip(all_component_ids, seg_type))),
                               vc=vertex_clusters)

    def node_to_segment(self, treenode_id):
        '''given a treenode id, return which bouton (or None) it belongs to'''
        return self.segments.idx_to_id(self.graph.vs.select(name_eq=treenode_id)['segment'][0])


class NeuronSegments(object):
    class synapses:
        def __init__(self, compact_connectors):
            self.nodes = [i[0] for i in compact_connectors]
            self.connectors = [i[1] for i in compact_connectors]
            self.relation = [i[2] for i in compact_connectors]
            self.xyz = np.array([i[3:] for i in compact_connectors])
            # from-treenode (Not neccessarily pre-synaptic!! The treenode in the current skeleton)
            # to-treenode
            # relation: 1 for post-synaptic and 0 for pre-synaptic connections
            # x, y, z

    """docstring for """
    def __init__(self):
        self.segment_indices = []
        self.ids = []
        self.nodes_id = []
        self.nodes_xyz = []
        self.centroids = []
        self.ConvexHullVol = []
        self.attr = {}
        self.inputs = []
        self.outputs = []
        self.segment_urls = []

    def add(self, index, SegmentId, nodes_id, nodes_xyz, compact_connectors):
        ctd = nodes_xyz.mean(0)
        self.segment_indices.append(index)
        self.ids.append(SegmentId)
        self.nodes_id.append(nodes_id)
        self.nodes_xyz.append(nodes_xyz)
        self.centroids.append(ctd)
        self.ConvexHullVol.append(convex_hull_vol(nodes_xyz))

        # syns is the segment subset [connectors] content from compact_skeleton
        self.inputs.append(NeuronSegments.synapses(
                           [i for i in compact_connectors if i[2]==1]))
        self.outputs.append(NeuronSegments.synapses(
                           [i for i in compact_connectors if i[2]==0]))
        self.segment_urls.append(url % tuple(nodes_xyz[0,:]))

    def add_attr(self,**kwargs):
        self.attr.update(kwargs)

    def idx_to_id(self,index):
        if index not in self.segment_indices:
            r = None
        else:
            r = self.ids[self.segment_indices.index(index)]
        return r

    def id_to_ctd(self, SegmentId):
        return self.centroids[self.ids.index(SegmentId)]

    def id_to_ConvexHullVol(self,SegmentId):
        return self.ConvexHullVol[self.ids.index(SegmentId)]

    def ids_to_ctds(self):
        '''return a dict of ctds keyed by ids'''
        return dict(list(zip(self.ids,self.centroids)))

    def id_to_xyz(self, sid):
        return self.nodes_xyz[self.ids.index(sid)]

    def id_to_inputs(self, sid):
        return self.inputs[self.ids.index(sid)]

    def all_segments_nodes(self):
        return np.concatenate(self.nodes_xyz, axis=0)

    def ConvexHullVol_of_all_segments(self):
        return convex_hull_vol(self.all_segments_nodes())

    def centroid_of_all_segments(self):
        return self.all_segments_nodes().mean(0)

    def id_to_url(self, sid):
        return self.segment_urls[self.ids.index(sid)]


def convex_hull_vol(nodes_xyz):
    if nodes_xyz.shape[0] >= 4:
        r = ConvexHull(nodes_xyz).volume
    elif nodes_xyz.shape[0] < 2:
        r = 0
    elif nodes_xyz.shape[0] == 2:
        r = np.linalg.norm(nodes_xyz[(0,)] - nodes_xyz[(1,)])
    elif nodes_xyz.shape[0] == 3:
        #Heron's formula to calculate the area of triangle
        # calculate the semi-perimeter
        # s = (a + b + c) / 2
        # calculate the area
        # area = (s*(s-a)*(s-b)*(s-c)) ** 0.5
        sides = distance.pdist(nodes_xyz)
        s = sum(sides)/2
        r = (s*(s-sides[0])*(s-sides[1])*(s-sides[2])) ** 0.5
    return r


def init_from_skid_list(connection, skids):
    neurons = {i : Neuron.init_from_json(connection, i) for i in skids}
    return neurons

if __name__ == '__main__':
    main()
