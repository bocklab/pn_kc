# modified from script by Noah Nelson in 2015 summer

import requests
import json
import sys
from operator import and_
from functools import reduce
import os

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


#------------------------------------------------------------------------------
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
    if isinstance(connection, str):
        connection = os.path.join(connection, "compact_sk")
        with open(os.path.join(connection, "%s") % skid) as outfile:
            sk = json.load(outfile)
        return sk
    else:
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

def save_compact_sk(connection, skid, path):
    # path = pn_kc/data/

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
    path = os.path.join(path, "compact_sk")
    if not os.path.exists(path):
        os.makedirs(path)
    with open(os.path.join(path, "%s" % skid), 'w+') as outfile:
        json.dump(json.loads(r.text), outfile)

def load_compact_sk(skid, path):
    path = os.path.join(path, "compact_sk")
    with open(os.path.join(path, "%s.json" % skid)) as outfile:
        r = json.load(outfile)
    return r



#------------------------------------------------------------------------------
def skid_to_name(connection, skids):
    if isinstance(connection, str):
        connection = os.path.join(connection, "neurons_names")
        with open(connection) as outfile:
                r = json.load(outfile)
        t = [str(j) for j in skids]
        return [r[i] for i in t]
    else:
        name_list = [get_neuron_name(connection, x) for x in skids]
        return name_list

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

def save_neurons_names(connection, sk_ids, path):
    path = os.path.join(path, "neurons_names")
    r = dict(zip(sk_ids, [get_neuron_name(connection, i) for i in sk_ids]))
    with open(path, 'w+') as outfile:
            json.dump(r, outfile)



#------------------------------------------------------------------------------
def get_root_node(connection, skid):
    if isinstance(connection, str):
        # right now it will load the root node file every time...
        # not ideal but...
        connection = os.path.join(connection, "root_node")
        with open(connection) as outfile:
            r = json.load(outfile)
        return r[str(skid)]
    else:
        connection.addpath('/%s/skeletons/%s/root' % (connection.projectID, skid))
        r = connection.get()
        data = json.loads(r.text)
        return data

def save_root_node(connection, skids, path):
    path = os.path.join(path, "root_node")
    f = {}
    for skid in skids:
        connection.addpath('/%s/skeletons/%s/root' % (connection.projectID, skid))
        r = connection.get()
        f[skid] = json.loads(r.text)
    with open(path, 'w+') as outfile:
        json.dump(f, outfile)
    print("root nodes saved")



#------------------------------------------------------------------------------
def retrieve_annotated_annotations(connection, meta_anno, returns='name'):
#    type_meta_anno: 'glom_class', 'kc_class'
    if isinstance(connection, str):
        # if connection is a pre-saved path,
        # returns can only be 'id'
        connection = os.path.join(connection, meta_anno)
        with open(connection) as outfile:
            r = json.load(outfile)
        return r
    else:
        # returns can be 'name' (return annotation names) or 'id' (id of the annotation)
        response = get_annotation_ids(connection, [meta_anno])
        connection.addpath('/%s/annotations/query-targets' % connection.projectID)
        opts = {'project_id': 1, 'annotated_with[0]': response[0]}
        r = connection.post(opts)
        data = json.loads(r.text)
        results = [i[returns] for i in data['entities']]
        return results

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

def save_annotated_annotations(connection, meta_anno, returns, path):
    response = get_annotation_ids(connection, [meta_anno])
    connection.addpath('/%s/annotations/query-targets' % connection.projectID)
    opts = {'project_id': 1, 'annotated_with[0]': response[0]}
    r = connection.post(opts)
    data = json.loads(r.text)
    results = [i[returns] for i in data['entities']]
    with open(os.path.join(path, meta_anno), 'w+') as outfile:
        json.dump(results, outfile)
    print(meta_anno + "saved")








#------------------------------------------------------------------------------
def retrieve_annotations_for_skeleton(connection, skid):
    if isinstance(connection, str):
        with open(os.path.join(connection, "annotations_for_skeleton/%s" % skid)) as file:
            r = json.load(file)
        return r
    else:
        connection.addpath('/%s/annotations/forskeletons' % connection.projectID)
        opts = {'skeleton_ids[0]': skid}
        r = connection.post(opts)
        results = json.loads(r.text)
        return results['annotations']

def save_annotations_for_skeleton(connection, skids, path):
    path = os.path.join(path, "annotations_for_skeleton/")
    if not os.path.exists(path):
        os.makedirs(path)
    for skid in skids:
        connection.addpath('/%s/annotations/forskeletons' % connection.projectID)
        opts = {'skeleton_ids[0]': skid}
        r = connection.post(opts)
        results = json.loads(r.text)
        with open(os.path.join(path, "%s" % skid), "w+") as file:
            json.dump(results['annotations'], file)
            print("%s saved" % skid)



#---------------------------------------------------------------------------
def get_pre_post_info(connection, pre_skids, post_skids, info_file="_"):
    if isinstance(connection, str):
        with open(os.path.join(connection, "pre_post_info", info_file)) as file:
            r = json.load(file)
        return r
    else:
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

def save_pre_post_info(connection, pre_skids, post_skids, path, info_file):
    # info_file is one of "RandomDraw", "t1p", or "all"
    connection.addpath('/%s/connector/info' % connection.projectID)

    opts = {'pre[%s]' % i: pre_skids[i] for i in range(len(pre_skids))}
    opts.update({'post[%s]' % i: post_skids[i]
                         for i in range(len(post_skids))})
    path = path + "pre_post_info/"
    if not os.path.exists(path):
        os.makedirs(path)
    r = connection.post(opts)
    data = json.loads(r.text)
    with open(os.path.join(path, info_file), "w+") as file:
        json.dump(data, file)
    print(info_file + " saved")






# Basic fafb server initiation
def fafb_connection(token):
    projectID, server = fafbServer
    fafb_c = Connection(server, projectID, token)
    return fafb_c


if __name__ == '__main__':
    main()
