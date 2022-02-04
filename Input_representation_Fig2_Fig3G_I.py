

import pandas as pd

path = wd + "/tables/"

tbl = pd.read_excel(path + "191029_tbl.xlsx", index_col=0)

#  extract the x axis order
fafb_tbl = pd.read_excel(path + "210202-fafb_pn_counts.xlsx", index_col=0)
hemi_tbl = pd.read_excel(path + "210201-hemibrain_pn_counts.xlsx", index_col=0).rename(columns={"hemi_pn":"gloms"})
hemi_tbl.loc[52, "gloms"] = "VP1"
m_tbl = fafb_tbl.merge(hemi_tbl, how='outer', on="gloms")
m_tbl.drop(index=54, inplace=True)
m_tbl = m_tbl.sort_values("fafb_counts", ascending=False).reset_index()

# Fig2B-C
# number of boutons per PN, per glom
p_tbl = tbl
p_tbl["cat_glom"] = pd.Categorical(p_tbl.short_glom_name, categories=m_tbl.gloms, ordered=True)
p_tbl.sort_values("cat_glom", inplace=True)
# p_tbl.reset_index(inplace=True)

fig, ax = plt.subplots()
ind = np.array(p_tbl.index)
width = 0.8

r1 = ax.bar(ind + 0.4, p_tbl.num_boutons, width, color=sns.xkcd_palette(['windows blue']), align='center')

plt.xticks(ind + 0.4, p_tbl.short_glom_name)

ticks = ax.get_xticklabels()
for i,tick in enumerate(ticks):
        tick.set_color(tbl.query("short_glom_name==@tick.get_text()").color.iloc[0])
#        if tick.get_text() in comm_gloms:
#            tick.set_weight("extra bold")

ax.set_ylabel("number of boutons")
ax.tick_params(axis='both', which='major', labelsize=14)
plt.xlim(-0.4,113.5)
fig.set_size_inches(70,10)
# fig.savefig(save_path + "210209-FAFB_boutonsPerPNs.png", bbox_inches="tight")




# boutons per glom
p_tbl = glom_btn_table
p_tbl["cat_glom"] = pd.Categorical(p_tbl.short_glom_name, categories=m_tbl.gloms, ordered=True)
p_tbl = p_tbl.sort_values("cat_glom")
p_tbl = p_tbl.reset_index(inplace=True)


#----------------------------------------------------------------------
# 210210 colorful bars

p_tbl = tbl
p_tbl["cat_glom"] = pd.Categorical(p_tbl.short_glom_name, categories=t_tbl.gloms, ordered=True)
p_tbl = p_tbl.sort_values("cat_glom")
p_tbl.reset_index()

fig, ax = plt.subplots()
ind = np.array(p_tbl.index)
width = 0.8

r1 = ax.bar(ind + 0.4, p_tbl.num_boutons, width, align='center')

plt.xticks(ind + 0.4, p_tbl.cat_glom, rotation=90)

col_list = []
ticks = ax.get_xticklabels()
for i,tick in enumerate(ticks):
    t1 = tbl.query("short_glom_name==@tick.get_text()").color.iloc[0]
    tick.set_color(t1)
    col_list.append(t1)
#        if tick.get_text() in comm_gloms:
#            tick.set_weight("extra bold")

r1 = ax.bar(ind + 0.4, p_tbl.num_boutons, width, color=col_list, align='center')
ax.set_ylabel("number of boutons", size=30)
ax.tick_params(axis='both', which='major', labelsize=30, left=True)
plt.xlim(-0.4,113.5)
ax.set_title("Number of boutons per PN", pad = 15)
fig.set_size_inches(70,10)
# fig.savefig(save_path + "210210-FAFB_boutonsPerPNs_ordered.png", bbox_inches="tight")





# boutons per glom
p_tbl = glom_btn_table
p_tbl["cat_glom"] = pd.Categorical(p_tbl.short_glom_name, categories=m_tbl.gloms, ordered=True)
p_tbl = p_tbl.sort_values("cat_glom")
p_tbl = p_tbl.reset_index()

fig, ax = plt.subplots()
ind = np.array(p_tbl.index)
width = 0.8

r1 = ax.bar(ind + 0.4, p_tbl.bouton_count, width, align='center')
plt.xticks(ind + 0.5, p_tbl.gloms)

ticks = ax.get_xticklabels()
col_list = []
for i,tick in enumerate(ticks):
    t1 = tbl.query("short_glom_name==@tick.get_text()").color.iloc[0]
    tick.set_color(t1)
    col_list.append(t1)
#        if tick.get_text() in comm_gloms:
#            tick.set_weight("extra bold")


r1 = ax.bar(ind + 0.4, p_tbl.bouton_count, width, colors=col_list, align='center')
plt.xlim(-0.4,54)
ax.set_ylabel("number of boutons")
ax.tick_params(axis='both', which='major', labelsize=16)
ax.set_title("Number of boutons per glom", pad = 15)
fig.set_size_inches(40,10)
# fig.savefig(save_path + "210210-FAFB_boutonsPerGlom.png", bbox_inches="tight")


# Figure 2A
# PNs per glom
fafb_tbl = pd.read_excel(path + "210202-fafb_pn_counts.xlsx", index_col=0)
hemi_tbl = pd.read_excel(path + "210201-hemibrain_pn_counts.xlsx", index_col=0).rename(columns={"hemi_pn":"gloms"})
hemi_tbl.loc[52, "gloms"] = "VP1"
m_tbl = fafb_tbl.merge(hemi_tbl, how='outer', on="gloms")
m_tbl.drop(index=54, inplace=True)

p_tbl = m_tbl

p_tbl = p_tbl.sort_values("fafb_counts", ascending=False).reset_index()
fig, ax = plt.subplots()
ind = np.array(p_tbl.index)
width = 0.8

r1 = ax.bar(ind + 0.8, p_tbl.fafb_counts, width, align='center')

plt.xticks(ind + 0.8, p_tbl.gloms, rotation=90)

col_list = []
ticks = ax.get_xticklabels()
for i,tick in enumerate(ticks):
    t1 = tbl.query("short_glom_name==@tick.get_text()").color.iloc[0]
    tick.set_color(t1)
    col_list.append(t1)
#        if tick.get_text() in comm_gloms:
#            tick.set_weight("extra bold")
r1 = ax.bar(ind + 0.8, p_tbl.fafb_counts, width, color=col_list, align='center')

plt.xlim(-0.4,54.4)
ax.set_ylabel("number of PNs", size=30)
ax.set_title("Number of PNs per PN subtype", pad = 15)
ax.tick_params(axis='both', which='major', labelsize=30, left=True)
fig.set_size_inches(40,10)
# fig.savefig(save_path + "210211-FAFB_PNsPerGlom.png", bbox_inches="tight")


# Figure 2C
# average number of boutons per glom
fafb_tbl = pd.read_excel(path + "210202-fafb_pn_counts.xlsx", index_col=0)

# boutons per glom
m_tbl = glom_btn_table
m_tbl = m_tbl.rename(columns = {"short_glom_name":"gloms"})
m_tbl = fafb_tbl.merge(m_tbl, how='outer', on="gloms")
m_tbl = m_tbl.assign(avg_btn = lambda x: x.bouton_count / x.fafb_counts).sort_values("avg_btn", ascending=False).reset_index()
p_tbl = m_tbl

fig, ax = plt.subplots()
ind = np.array(p_tbl.index)
width = 0.8

r1 = ax.bar(ind + 0.8, p_tbl.avg_btn, width, align='center')

plt.xticks(ind + 0.8, p_tbl.gloms, rotation=90)

col_list = []
ticks = ax.get_xticklabels()
for i,tick in enumerate(ticks):
    t1 = tbl.query("short_glom_name==@tick.get_text()").color.iloc[0]
    tick.set_color(t1)
    col_list.append(t1)
#        if tick.get_text() in comm_gloms:
#            tick.set_weight("extra bold")
r1 = ax.bar(ind + 0.8, p_tbl.avg_btn, width, color=col_list, align='center')

plt.xlim(-0.4,54.4)
ax.set_ylabel("average number of boutons", size=30)
ax.set_title("Number of boutons per PN", pad = 15)
ax.tick_params(axis='both', which='major', labelsize=30, left=True)
fig.set_size_inches(40,10)
# fig.savefig(save_path + "210211-AvgBoutonPerPN.png", bbox_inches="tight")




# Fig1E
stbl = fafb_tbl
stbl["significance"] = df_lookup("short_glom_name",fafb_tbl.gloms, "significance", tbl)
p_tbl = stbl.significance.value_counts()

fig, ax = plt.subplots()
r1 = ax.bar(p_tbl.index, p_tbl, color=df_lookup("significance", p_tbl.index, "color", tbl), align='center')

ax.tick_params(axis='both', which='major', labelsize=12, left=True)
ax.set_ylabel("number of glomeruli", size=12)
fig.set_size_inches(7,7)
# fig.savefig(save_path + "210328-NumGlomPerCat.png", bbox_inches="tight")

# fig.savefig(save_path + "210210-NumGlomPerCat.png", bbox_inches="tight")








# Figure 3I
# 210215
# number of claws per glomerulus
gloms = []
claws = []
for glom in pd.unique(tbl.short_glom_name):
    gloms.append(glom)
    claws.append(tbl.query("short_glom_name==@glom")["num_claws"].sum())


p_tbl = pd.DataFrame({"gloms":gloms, "num_claws": claws})

p_tbl = p_tbl.sort_values("num_claws", ascending=False)
p_tbl = p_tbl.reset_index()

fig, ax = plt.subplots()
ind = np.array(p_tbl.index)
width = 0.8

r1 = ax.bar(ind + 0.8, p_tbl.num_claws, width, align='center')

plt.xticks(ind + 0.8, p_tbl.gloms, rotation=90)

col_list = []
ticks = ax.get_xticklabels()
for i,tick in enumerate(ticks):
    t1 = tbl.query("short_glom_name==@tick.get_text()").color.iloc[0]
    tick.set_color(t1)
    col_list.append(t1)
#        if tick.get_text() in comm_gloms:
#            tick.set_weight("extra bold")
r1 = ax.bar(ind + 0.8, p_tbl.num_claws, width, color=col_list, align='center')

plt.xlim(-0.4,54.4)
ax.set_ylabel("number of claws", size=30)
ax.set_title("Number of claws per glom", pad = 15)
ax.tick_params(axis='both', which='major', labelsize=30, left=True)
fig.set_size_inches(40,10)
# fig.savefig(save_path + "210214-ClawsPerGlom.png", bbox_inches="tight")


# Figure 3G
# 210217
# claws per glom divided by boutons per glom

# boutons per glom
b_tbl = glom_btn_table

 # number of claws per glomerulus
gloms = []
claws = []
for glom in pd.unique(tbl.short_glom_name):
    gloms.append(glom)
    claws.append(tbl.query("short_glom_name==@glom")["num_claws"].sum())

c_tbl = pd.DataFrame({"gloms":gloms, "num_claws": claws})

b_tbl = glom_btn_table
b_tbl = b_tbl.rename(columns={"short_glom_name":"gloms"})
cb_tbl = b_tbl.merge(c_tbl, how="inner", on="gloms")
cb_tbl = cb_tbl.assign(claws_p_boutons = lambda x:x.num_claws/x.bouton_count)

p_tbl = cb_tbl
p_tbl = p_tbl.sort_values("claws_p_boutons", ascending=False)
p_tbl = p_tbl.reset_index()

fig, ax = plt.subplots()
ind = np.array(p_tbl.index)
width = 0.8

r1 = ax.bar(ind + 0.8, p_tbl.claws_p_boutons, width, align='center')

plt.xticks(ind + 0.8, p_tbl.gloms, rotation=90)

col_list = []
ticks = ax.get_xticklabels()
for i,tick in enumerate(ticks):
    t1 = tbl.query("short_glom_name==@tick.get_text()").color.iloc[0]
    tick.set_color(t1)
    col_list.append(t1)
#        if tick.get_text() in comm_gloms:
#            tick.set_weight("extra bold")
r1 = ax.bar(ind + 0.8, p_tbl.claws_p_boutons, width, color=col_list, align='center')

plt.xlim(-0.4,54.4)
ax.set_ylabel("number of claws", size=30)
ax.set_title("Number of claws per bouton", pad = 15)
ax.tick_params(axis='both', which='major', labelsize=30, left=True)
fig.set_size_inches(40,10)
# fig.savefig(save_path + "210218-ClawsPerBoutonForGlom.png", bbox_inches="tight")
