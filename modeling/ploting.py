
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

tbl = pd.read_csv('data/190422-comm_fract_vs_performance.csv')

t1 = pd.DataFrame({"mean": tbl.mean(0), "ci95": tbl.sem(0)*1.96, "commfract": [int(i)/16 for i in tbl.mean(0).index.values]})


fig, ax1 = plt.subplots()
ax1.errorbar(x=t1['commfract'], y=t1['mean'], yerr=t1['ci95'], fmt='o', markersize=9)
plt.xlabel('fraction of community inputs', fontsize=15)
plt.ylabel('errror rates', fontsize=15)

ax1.xaxis.set_tick_params(labelsize=15)

ax1.yaxis.set_tick_params(labelsize=15)

fig.set_size_inches(12, 8)
# fig.savefig("data/190422-community_fraction_big.png", bbox_inches='tight')
