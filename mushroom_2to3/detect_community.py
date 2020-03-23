
import pandas as pd
import numpy as np
import seaborn as sns
from .shuffle import *
from sklearn.cluster import KMeans

from .build_connectivity import *

#####
# load functions
def get_mean_std_binom(num_claws, prob):
    # usually prob=glom_prob
    # sd (sqrt of variance of binomial distribution)
    # mean of probability for each glom
    r = {"mean": np.array([num_claws * i for i in prob]),
         "std": np.array([np.sqrt(num_claws * i * (1 - i)) for i in prob])}
    return r

def get_mean_std_model_free(num_claws, avg, sd):
    return {"mean": num_claws * avg, "std": sd * num_claws}

def get_input_preference(conn, glom_prob, order=False, model=False):
    # for example, order could be order=glom_seq
    # could set order=False to suppress it
    results = []
    for i in range(conn.shape[1]):
        kc_idx = np.where(conn[:, i])[0]

        t1 = conn[kc_idx, :].sum(0)

        t1[i] = t1[i] - len(kc_idx)
        if model:
            r = get_glom_mean_std(sum(t1), glom_prob)
        else:
            r = get_mean_std_model_free(sum(t1))

        r["observed"] = t1
        r["std_delta"] = np.array(
            [(t1[i] - r['mean'][i]) / r['std'][i] for i in range(len(t1))])
        results.append(r)
    t2 = np.concatenate([i['std_delta'].reshape(
        1, len(i['std_delta'])) for i in results])
    if order is not False:
        t2 = t2[np.ix_(order, order)]
    return t2

def get_input_matrix(conn, model_fun):
    # for example, model_fun = lambda x: get_mean_std_model_free(x, r_mean, r_std)
    results = []
    for i in range(conn.shape[1]):
        kc_idx = np.where(conn[:, i])[0]

        t1 = conn[kc_idx, :].sum(0)

        t1[i] = t1[i] - len(kc_idx)
        r = model_fun(sum(t1))
        r["observed"] = t1
        r["std_delta"] = np.array(
            [(t1[i] - r['mean'][i]) / r['std'][i] for i in range(len(t1))])
        results.append(r)
    t2 = np.concatenate([i['std_delta'].reshape(
        1, len(i['std_delta'])) for i in results])
    return t2

def get_input_matrix(conn, model_fun, avg, sd):
    results = []
    for i in range(conn.shape[1]):
        kc_idx = np.where(conn[:, i])[0]

        t1 = conn[kc_idx, :].sum(0)

        t1[i] = t1[i] - len(kc_idx)
        r = model_fun(sum(t1), avg[i,:], sd[i,:])
        r["observed"] = t1
        r["std_delta"] = np.array(
            [(t1[i] - r['mean'][i]) / r['std'][i] for i in range(len(t1))])
        results.append(r)
    t2 = np.concatenate([i['std_delta'].reshape(
        1, len(i['std_delta'])) for i in results])
    return t2

def get_input_preference_model_free(conn, prob, ids, num_exp=1000):
    # conn: glom-kc connectivity in claw units
    # prob: probability for shuffle control of randomly picking gloms
    # ids: the ids for the gloms of picking probability prob
    # num_exp: number of repetitions
    # output a conditioned input data matrix object with its column and row ids
    avg, sd = shuffle_and_compute_mean_std(conn, prob, num_exp)
    data = get_input_matrix(
        conn, lambda x: get_mean_std_model_free(x, avg, sd))
    return PairMatrix('', data, ids)

def shuffle_and_compute_mean_std(conn, prob, num_exp=1000):
    stat = [get_input_probability(shuffle_glom_kc_w_prob(
        conn, prob)) for i in range(num_exp)]
    stat = np.array(stat)
#    stat[np.isnan(stat)] = 0
    sd = np.nanstd(stat, axis=0)
    avg = np.nanmean(stat, axis=0)
    return avg, sd

def get_input_probability(conn):
    results = []
    for i in range(conn.shape[1]):
        kc_idx = np.where(conn[:, i])[0]
        t1 = conn[kc_idx, :].sum(0)

        t1[i] = t1[i] - len(kc_idx)
        results.append(t1/np.sum(t1))
    return np.array(results)

def order_heatmap(m, op="cluster", save_name=False, new_fig=True):
    # check glom_seq match back to the original glomerular sequence
    # op will either be "cluster" or
    # a sequence of length 53, which would be used for ordering the matrix
    m[np.isnan(m)] = 0
    if op == "cluster":
        fi = KMeans(n_clusters=8).fit(m)
        sorted_idx = np.argsort(fi.labels_)
        data = m[np.ix_(sorted_idx, sorted_idx)]
    else:
        sorted_idx = op
        data = m[np.ix_(sorted_idx, sorted_idx)]
    if new_fig:
        plt.figure(figsize=(21, 11.970));
    f = sns.heatmap(data, xticklabels=sorted_idx, yticklabels=sorted_idx, vmin=-8.53, vmax=8.53)
    if save_name is not False:
        plt.savefig(save_name)
    return {'fig': f, 'data': data, 'order': sorted_idx}

def make_heatmap(data, idx, **kwargs):
    # reorder the matrix based on idx for X, idx for Y, and plot heatmap
    data = data[np.ix_(idx, idx)]
    return sns.heatmap(data, xticklabels=idx, yticklabels=idx, vmin=-8.53, vmax=8.53, **kwargs)

def km_cluster(m):
    m[np.isnan(m)] = 0
    fi = KMeans(n_clusters=4).fit(m)
    return np.argsort(fi.labels_)

def processing(ana_ip, op="cluster", shuffle=False):
    # check glom_seq match back to the original glomerular sequence
    # op will either be "cluster" or
    # a sequence of length 53, which would be used for ordering the matrix
    conn = ana_ip.conn_data['glom_kc_in_claw_unit'].conn['1s']
    if shuffle:
        conn = shuffle_glom_kc_w_prob(conn, glom_prob)
    r = get_results(conn)[np.ix_(glom_seq, glom_seq)]
    r[np.isnan(r)] = 0
    if op == "cluster":
        fi = KMeans(n_clusters=8).fit(r)
        sorted_idx = np.argsort(fi.labels_)
        data = r[np.ix_(sorted_idx, sorted_idx)]
    else:
        sorted_idx = op
        data = r[np.ix_(sorted_idx, sorted_idx)]
    plt.figure(figsize=(21, 11.970));
    f = sns.heatmap(data, xticklabels=sorted_idx, yticklabels=sorted_idx, vmin=-8.53, vmax=8.53)
    return {'fig': f, 'data': data, 'order': sorted_idx}

def permute(m):
    return np.random.permutation(m.flatten()).reshape(m.shape)

def sum_conditioned_inputs(conn, w, order, threshold=5):
    '''
    w=btns_per_glom or w=1
    order=glom_seq
    could set order=False to suppress it
    note that for now the inputs are conditioned on KCs with >5 synapses in the first place
    '''
    results = []
    if w is False:
        w = 1
    for i in range(conn.shape[1]):
        kc_idx = np.where(conn[:, i]>=threshold)[0]
        t1 = conn[kc_idx, :].sum(0) / w
        results.append(t1.reshape(1,len(t1)))

    t2 = np.concatenate(results)
    if order is not False:
        t2 = t2[np.ix_(order, order)]
    return t2

def reorder(ref_ids, target_ids):
    # given the inputs, output an ordering index list
    # that would reorder target_ids to match ref_ids,
    # with elements in target_ids not in ref_ids being padded at the back
    # e.g. [target_ids[i] for i in seq] would get the re-arranged target_ids that match ref_ids
    seq = [target_ids.index(i) for i in ref_ids if i in target_ids]
    seq.extend([i for i in range(len(target_ids)) if i not in seq])
    return seq

def df_lookup(ref_col, ref_vals, target_col, df):
    # ref_val should handle a list and return a list
    # for example,
    # ref_col = 'glom_class'
    # ref_val = 'VM7v'
    # target_col = 'glom_number'
    # df = glom_tbl
    # df[df['glom_class'].values == 'VM7v']['glom_number']
    # return the target_vals in a list
    return [df[df[ref_col].values == i][target_col].tolist()[0] for i in ref_vals]

def get_conditional_matrix(conn, prob, idx):
    '''
    Parameters
    ----------
    conn:                a connectivity matrix where columns are pre-synaptic
                         and rows are post-synaptic neurons
                         e.g. a glom-KC connectivity matrix
                         conn =  ana.conn_data['glom_kc_in_claw_unit'].conn['1s']
    prob :               A list of null-model probability for each column. All
                         claws randomly pick glomeruli (columns) based on this probability
    idx :                a list of index for both rows and columns (they are
                         the same) of the resulting conditional matrix

    Returns
    -------
    a PairMatrix object with a glom_vs_glom probability matrix (rows in the same order as columns). Also have row and column indeces.

    Comments
    -------
    Originally the function is called compare_to_random_conn(conn, prob, idx)
    '''
    stat = [get_input_probability(shuffle_glom_kc_w_prob(
        conn, prob)) for i in range(1000)]
    stat = np.array(stat)
    sd = np.nanstd(stat, axis=0)
    avg = np.nanmean(stat, axis=0)

    data = get_input_matrix(
        conn, lambda x, y, z: get_mean_std_model_free(x, y, z), avg, sd)

    return PairMatrix('', data, idx)

def get_raw_inputs(conn):
    results = []
    for i in range(conn.shape[1]):
        kc_idx = np.where(conn[:, i])[0]
        t1 = conn[kc_idx, :].sum(0)
        t1[i] = t1[i] - len(kc_idx)
        results.append(t1)
    return np.array(results)
