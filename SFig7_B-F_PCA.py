
# try run PCA on FAFB data
#-----------------------------------------------------
from sklearn.decomposition import PCA
num_exp = 10000
ana = ana_all_rd
conn_data = ana.conn_data['glom_kc_in_claw_unit']
ob_conn, glom_prob, glom_idx_ids = get_conn_prob_idx(conn_data)

results = []
for i in range(num_exp):
    sfl_conn = shuffle_glom_kc_w_prob(ob_conn, glom_prob)
    pca = PCA()
    pca.fit_transform(sfl_conn)
    results.append(pca.explained_variance_ratio_)

rd_results = np.array(results)
rd_r = {}
rd_r['mean'] = np.mean(rd_results, axis=0)
rd_r['std'] = np.std(rd_results, axis=0)

pca = PCA()
t1 = pca.fit_transform(rd_conn)
rd_exp_ratio = pca.explained_variance_ratio_

xa = np.arange(1, len(rd_exp_ratio)+1)

import matplotlib
matplotlib.rcParams.update({'font.size': 20})
fig, ax = plt.subplots()
ax.errorbar(xa, rd_r['mean'], yerr=rd_r['std'], fmt='ko')
ax.plot(xa, rd_exp_ratio, 'ro')
ax.set_xlabel('components')
ax.set_ylabel('Fraction variance explained')
ax.set_ylim(top=0.07)
ax.set_title("All KCs observed vs. random bouton model", pad = 20)


import matplotlib.ticker as plticker
# this locator puts ticks at regular intervals
ax.xaxis.set_major_locator(plticker.MultipleLocator(base=5))
fig.set_size_inches(12, 8)

plt.tick_params(axis='both', which='both', bottom=True, left=True)
# fig.savefig(save_path + '210201-allKCs_ObvsRandomBouton_pca.png')


# ob vs. random glom
n = len(glom_prob)
prob = [1/n]*n
results = []
for i in range(num_exp):
    sfl_conn = shuffle_glom_kc_w_prob(ob_conn, prob)
    pca = PCA()
    pca.fit_transform(sfl_conn)
    results.append(pca.explained_variance_ratio_)

rd_results = np.array(results)
rd_r = {}
rd_r['mean'] = np.mean(rd_results, axis=0)
rd_r['std'] = np.std(rd_results, axis=0)

pca = PCA()
t1 = pca.fit_transform(rd_conn)
rd_exp_ratio = pca.explained_variance_ratio_

xa = np.arange(1, len(rd_exp_ratio)+1)

import matplotlib
matplotlib.rcParams.update({'font.size': 20})
fig, ax = plt.subplots()
ax.errorbar(xa, rd_r['mean'], yerr=rd_r['std'], fmt='ko')
ax.plot(xa, rd_exp_ratio, 'ro')
ax.set_xlabel('components')
ax.set_ylabel('Fraction variance explained')
ax.set_ylim(top=0.07)
ax.set_title("All KCs observed vs. random glom model", pad = 20)


import matplotlib.ticker as plticker
# this locator puts ticks at regular intervals
ax.xaxis.set_major_locator(plticker.MultipleLocator(base=5))
fig.set_size_inches(12, 8)

plt.tick_params(axis='both', which='both', bottom=True, left=True)
# fig.savefig(save_path + '210201-allKCs_ObvsRandomGlom_pca.png')



# Observed vs. random claw model
results = []
sfl_list = shuffle_glom_kc_iterate(ob_conn, num_exp)
for i in range(num_exp):
    pca = PCA()
    pca.fit_transform(sfl_list[i])
    results.append(pca.explained_variance_ratio_)

rd_results = np.array(results)
rd_r = {}
rd_r['mean'] = np.mean(rd_results, axis=0)
rd_r['std'] = np.std(rd_results, axis=0)

pca = PCA()
t1 = pca.fit_transform(rd_conn)
rd_exp_ratio = pca.explained_variance_ratio_

xa = np.arange(1, len(rd_exp_ratio)+1)

import matplotlib
matplotlib.rcParams.update({'font.size': 20})
fig, ax = plt.subplots()
ax.errorbar(xa, rd_r['mean'], yerr=rd_r['std'], fmt='ko')
ax.plot(xa, rd_exp_ratio, 'ro')
ax.set_xlabel('components')
ax.set_ylabel('Fraction variance explained')
ax.set_ylim(top=0.07)
ax.set_title("All KCs observed vs. random claw model", pad = 20)


import matplotlib.ticker as plticker
# this locator puts ticks at regular intervals
ax.xaxis.set_major_locator(plticker.MultipleLocator(base=5))
fig.set_size_inches(12, 8)

plt.tick_params(axis='both', which='both', bottom=True, left=True)
# fig.savefig(save_path + '210201-allKCs_ObvsRandomClaw_pca.png')



# observed vs local random model
# first run some script in 200212-Ob_vs_LocalBtnRdModel.py
results = []
for i in range(num_exp):
    sfl_conn = shuffle_glom_kc_w_nearest5(ob_conn, conn_data.row_ids, rows, glom_to_pick)
    pca = PCA()
    pca.fit_transform(sfl_conn)
    results.append(pca.explained_variance_ratio_)

rd_results = np.array(results)
rd_r = {}
rd_r['mean'] = np.mean(rd_results, axis=0)
rd_r['std'] = np.std(rd_results, axis=0)

pca = PCA()
t1 = pca.fit_transform(ob_conn)
rd_exp_ratio = pca.explained_variance_ratio_

xa = np.arange(1, len(rd_exp_ratio)+1)

import matplotlib
matplotlib.rcParams.update({'font.size': 20})
fig, ax = plt.subplots()
ax.errorbar(xa, rd_r['mean'], yerr=rd_r['std'], fmt='ko')
ax.plot(xa, rd_exp_ratio, 'ro')
ax.set_xlabel('components')
ax.set_ylabel('Fraction variance explained')
ax.set_ylim(top=0.07)
ax.set_title("All KCs observed vs. local random model", pad = 20)


import matplotlib.ticker as plticker
# this locator puts ticks at regular intervals
ax.xaxis.set_major_locator(plticker.MultipleLocator(base=5))
fig.set_size_inches(12, 8)

plt.tick_params(axis='both', which='both', bottom=True, left=True)
# fig.savefig(save_path + '210201-allKCs_ObvsLocalRandom_pca.png')
