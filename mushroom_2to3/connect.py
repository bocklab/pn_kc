# modified from script by Noah Nelson in 2015 summer

import requests
import json
import sys
from operator import and_
from functools import reduce

# use your CATMAID fafb Token
# %run /Users/zhengz11/myscripts/mushroom_v9/credential/fafb_tokens.py

fafbServer = (1, 'https://neuropil.janelia.org/tracing/fafb/v14')


class NoNeuron(Exception):
    pass


class Auth(requests.auth.AuthBase):
    def __init__(self, token):
        self.token = token

    def __call__(self, r):
        r.headers['X-Authorization'] = 'Token {}'.format(self.token)
        return r


class Connection:
    def __init__(self, server, projectID, token):
        self.server = server
        self.projectID = projectID
        self.token = token
        self.session = requests.Session()

    def addpath(self, path):
        self.url = self.server + path

    def get(self):
        return self.session.get(self.url, auth=Auth(self.token))

    def post(self, opts):
        return self.session.post(self.url, auth=Auth(self.token), data=opts)


def get_annotation_ids(connection, query):
    if not isinstance(query, list):
        print("input is not a list!")
    connection.addpath('/%s/annotations' % connection.projectID)
    r = connection.get()
    annotations = json.loads(r.text)
    id_list = []
    for anno_query in query:
        for annotation in annotations['annotations']:
            if anno_query == annotation['name']:
                id_list.append(annotation['id'])
    return id_list


def retrieve_annotations_query_targets(connection, annotation_id):
    connection.addpath('/%s/annotations/query-targets' % connection.projectID)
    opts = {'project_id': 1, 'annotated_with[0]': ','.join(
        map(str, annotation_id))}
    r = connection.post(opts)
    data = json.loads(r.text)
    results = [idx['skeleton_ids'][0] for idx in data['entities']]
    return results


def retrieve_annotated_annotations(connection, meta_anno, returns='name'):
    # returns can be 'name' (return annotation names) or 'id' (id of the annotation)
    response = get_annotation_ids(connection, [meta_anno])
    connection.addpath('/%s/annotations/query-targets' % connection.projectID)
    opts = {'project_id': 1, 'annotated_with[0]': response[0]}
    r = connection.post(opts)
    data = json.loads(r.text)
    results = [i[returns] for i in data['entities']]
    return results


def retrieve_annotations_for_skeleton(connection, skid):
    connection.addpath('/%s/annotations/forskeletons' % connection.projectID)
    opts = {'skeleton_ids[0]': skid}
    r = connection.post(opts)
    results = json.loads(r.text)
    return results['annotations']


def get_neuron_name(connection, skeleton_id):
    connection.addpath('/%s/skeleton/%s/neuronname' % (connection.projectID,
                                                       skeleton_id))
    r = connection.get()
    namedict = json.loads(r.text)
    if 'error' in namedict:
        raise NoNeuron("no neuron with id %s" % skeleton_id)
        return ""
    return namedict['neuronname']


def get_neurons_names(connection, sk_ids):
    # same as above except take a list of skeleton_ids as input and return a list of names as outputs
    return [get_neuron_name(connection, i) for i in sk_ids]


def get_root_node(connection, skid):
    connection.addpath('/%s/skeletons/%s/root' % (connection.projectID, skid))
    r = connection.get()
    data = json.loads(r.text)
    return data


def ids_to_annotations(connection, anno_ids):
    # anno_ids should be a list
    connection.addpath('/%s/annotations' % connection.projectID)
    r = connection.get()
    data = json.loads(r.text)
    all_annos = data['annotations']
    anno_dict = {i['id']: i['name'] for i in all_annos}
    annos = [anno_dict[i] for i in anno_ids]
    return annos


def skid_to_name(connection, skids):
    name_list = [get_neuron_name(connection, x) for x in skids]
    return name_list


def download_annotation_map(connection, anno_list):
    f_anno = open('annotation_ids.txt', 'w')
    for anno in anno_list:
        id_list = []
        f_anno.write(anno)
        anno_id = get_annotation_ids(connection, [anno])
        sk_ids = retrieve_annotations_query_targets(connection, anno_id)
        id_list = list(set(sk_ids))
        for skid in id_list:
            f_anno.write(',' + str(skid))
        f_anno.write('\n')
    f_anno.close()


def build_anno_sk_dict(connection, anno_list):
    ann_sk = {}
    for anno in anno_list:
        id_list = []
        anno_id = get_annotation_ids(connection, [anno])
        sk_ids = retrieve_annotations_query_targets(connection, anno_id)
        id_list = list(set(sk_ids))
        ann_sk[anno] = id_list
    return ann_sk


def retrieve_skids_from_anno(connection, anno_list):
    # can take a list of annotations,
    # output a list of skeletons that have any of the provided annotations
    r = []
    for i in anno_list:
        r.extend(retrieve_annotations_query_targets(connection,
                                                    get_annotation_ids(connection, [i])))
    return list(set(r))


def get_compact_sk(connection, skid, save_pickle=False):
    # "nodes", "connectors", "tags"
    # for each node,
    # 0: id
    # 1: parent_id
    # 2: user_id
    # 3: x
    # 4: y
    # 5: z
    # 6: radius
    # 7: confidence
    print(("downloading skeleton %s" % skid))
    connection.addpath('/%s/%s/1/1/compact-skeleton' % (connection.projectID,
                                                        skid))
    r = connection.get()
    if save_pickle is not False:
        with open("%s" % skid, 'wb') as f:
            pickle.dump(r, f)
    sk = json.loads(r.text)
    if 'error' in sk:
        raise NoNeuron("no neuron with skid %s" % skeleton_id)
        return ""
    return sk


def save_compact_sk(connection, skid):
    # "nodes", "connectors", "tags"
    # for each node,
    # 0: id
    # 1: parent_id
    # 2: user_id
    # 3: x
    # 4: y
    # 5: z
    # 6: radius
    # 7: confidence
    print(("downloading skeleton %s" % skid))
    connection.addpath('/%s/%s/1/1/compact-skeleton' % (connection.projectID,
                                                        skid))
    r = connection.get()
    with open("%s" % skid, 'w') as outfile:
        json.dump(json.loads(r.text), outfile)


def get_pre_post_info(connection, pre_skids, post_skids):
    # same as above, but this one accept different number in pre_skids and post_skids
    connection.addpath('/%s/connector/info' % connection.projectID)

    opts = {'pre[%s]' % i: pre_skids[i] for i in range(len(pre_skids))}
    opts.update({'post[%s]' % i: post_skids[i]
                 for i in range(len(post_skids))})

    r = connection.post(opts)
    data = json.loads(r.text)
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
    return data


def get_connector_detail(connection, connector_ids):
    # connector_ids should be a list
    data = []
    for connector_id in connector_ids:
        connection.addpath('/%s/connectors/%s' % (connection.projectID,
                                                  connector_id))
        r = connection.get()
        data.append(json.loads(r.text))
    return data


def remove_annotations(connection, entity_ids, annos):
    anno_ids = get_annotation_ids(connection, annos)
    connection.addpath('/%s/annotations/remove' % connection.projectID)
    opts = {'entity_ids[%s]' % i: entity_ids[i]
            for i in range(len(entity_ids))}
    opts.update({'annotation_ids[%s]' % i: anno_ids[i]
                 for i in range(len(anno_ids))})
    r = connection.post(opts)
    return r


def get_skids_from_annos(connection, and_annos, not_in_annos=None):
    # AND between each element of and_annos, each element is a list, which uses OR among its elements
    # e.g. [['Random Draw 1 KC', 'Random Draw 2 KC'], ['Complete']]
    # KCs from random draw 1 or 2, and also annotated 'Complete'
    skids = []
    for i in and_annos:
        skids.append(set(retrieve_skids_from_anno(connection, i)))
    skids = list(reduce(and_, skids))
    if not_in_annos is not None:
        not_in_skids = retrieve_skids_from_anno(connection, not_in_annos)
        skids = [i for i in skids if i not in not_in_skids]
    return skids


def get_labels_stats(connection):
    # [labelID, labelName, skeletonID, treenodeID]
    connection.addpath('/%s/labels/stats' % connection.projectID)
    r = connection.get()
    return json.loads(r.text)


def add_annotation(connection, skids, annos):
    # a list of skids, each one will be annotated with all anotations in annos(a list)
    connection.addpath('/%s/annotations/add' % connection.projectID)
    for anno in annos:
        opts = {'skeleton_ids[%s]' %
                cter: val for cter, val in enumerate(skids)}
        opts.update({'annotations[0]': anno})
        r = connection.post(opts)
    return r

# Basic fafb server initiation
def fafb_connection(token):
    projectID, server = fafbServer
    fafb_c = Connection(server, projectID, token)
    return fafb_c

def json_save(js, fname, save_path):
    f = fname + save_path
    with open(f, 'w+') as json_file:
        json.dump(js, json_file)
    json_file.close()
    print(f + ' saved')


if __name__ == '__main__':
    main()
