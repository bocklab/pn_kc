
library(ggplot2)
library(neuprintr)
library(hemibrainr)
library(tidyverse)
library(dendextend)
library(nat.nblast)
library(doParallel)
library(doMC)

library(dplyr)
numCores <- detectCores()
doMC::registerDoMC(numCores)

wd = "pn_kc/"

comm_gloms = c('DL2v', 'DM1', 'DM2', 'DM3', 'DM4', 'DP1m', 'VA2', 'VA4', 'VM2', 'VM3')

fpath = paste0(wd,"/200715-tbl.csv")
tbl = read.csv(fpath)

t2 = distinct(tbl, short_glom_name, .keep_all=T)
col_set = as.character(t2$rs_group_color)

col_set[col_set=="red"]="orange"
col_set[col_set=="purple"]="yellowgreen"
col_set[col_set=="blue"]="khaki4"
#"red", "orange"; "purple" "deeppink3"; "blue", "khaki4"

names(col_set) = t2$short_glom_name
col_set["VP4"] = "black"
col_set["VP1m"] =  "black"
col_set[col_set == "green"] = "green3"

# left side PNs
##-----------------------------------------------
library(xlsx)
library(catmaid)
library(elmr)
source(su)
# save_path = "/Users/zhengz11/myscripts/data_results/200710-PNKCrevision_optimal_cluster/"

t1 = read.xlsx(paste0(wd,"/210529-left_PN_tbl.xlsx"), sheetIndex=1)

pns = read.neurons.catmaid(t1$skeleton_id, conn=fafb_conn)

pns.ca = prune_in_volume(x = pns, surf = FAFB14NP.surf, neuropil = "MB_CA_L")

pns.ca.dps = dotprops(pns.ca/1000, k = 5, resample = 1)
# pns.ca.dps = dotprops(pns.ca, k = 5, resample = 1000)

names(pns.ca.dps) = t1$type

aba = nblast_allbyall(pns.ca.dps, .parallel = TRUE)

ls.pns.hc <- nhclust(scoremat=aba, method = "ward.D2")

# plot(pns.hc)
  
dend <- as.dendrogram(ls.pns.hc) %>% set("labels_col", col_set[labels(.)])
plot(dend)

t5 = get_leaves_attr(dend, "label")
n_pch = rep(NA, length(t5))
n_col = rep(NA, length(t5))
n_pch[t5 %in% comm_gloms] = 19
n_col[t5 %in% comm_gloms] = 3

ls_dend = set(dend, "leaves_pch", n_pch) %>% set("leaves_col", n_col) 

plot(ls_dend, main = "Left-side PNs")

# 210529_NBLAST_left-side_FAFB_pns.png


# right side PNs
##-----------------------------------------------

library(xlsx)
library(catmaid)
library(elmr)
source(su)

t1 = read.xlsx(paste0(wd,"/200721_right_PNs.xlsx"), sheetIndex=1)

pns = read.neurons.catmaid(t1$skeleton_id, conn=fafb_conn)

pns.ca = prune_in_volume(x = pns, surf = FAFB14NP.surf, neuropil = "MB_CA_R", OmitFailures = TRUE)

pns.ca.dps = dotprops(pns.ca/1000, k = 5, resample = 1)
# pns.ca.dps = dotprops(pns.ca, k = 5, resample = 1000)


names(pns.ca.dps) = t1$type

aba = nblast_allbyall(pns.ca.dps, .parallel = TRUE)

rs.pns.hc <- nhclust(scoremat=aba, method = "ward.D2")
# plot(pns.hc)

dend <- as.dendrogram(rs.pns.hc) %>% set("labels_col", col_set[labels(.)]) 

t5 = get_leaves_attr(dend, "label")
n_pch = rep(NA,length(t5))
n_col = rep(NA, length(t5))
n_pch[t5 %in% comm_gloms] = 19
n_col[t5 %in% comm_gloms] = 3

rs_dend = set(dend, "leaves_pch", n_pch) %>% set("leaves_col", n_col) 

plot(rs_dend, main = "Right-side PNs")

# 210529_NBLAST_right-side_FAFB_pns.png


# hemibrain
#--------------------------------
ca.mesh = neuprintr::neuprint_ROI_mesh(roi = "CA(R)")

pn_ids = intersect(as.numeric(upn.ids), neuprint_bodies_in_ROI("CA(R)")$bodyid)

# remove VP3, VP4, VP1m
pn_ids = setdiff(pn_ids,c(5813056072,663432544,543010474, 634759240))

pns.t = neuprint_read_skeletons(pn_ids, heal = F)

cant_stitch = c(5813055048, 542634818, 1944507292)

# using prune_in_volume
t2 = nlapply(pns.t, resample, 125, OmitFailures=FALSE)

# t2 = prune_twigs(t2, 375, OmitFailures = T)

t3 = nlapply(t2, prune_in_volume, as.hxsurf(ca.mesh), OmitFailures=FALSE)

pns.pr.dps = dotprops(t3/125, k = 5, resample = 1)
# pns.pr.dps = dotprops(t3, k = 5, resample = 125)
pn_types = sapply(neuprint_get_meta(names(pns.t))$type, function(x) sub("\\_.*", "", x))

names(pns.pr.dps) = pn_types

aba = nblast_allbyall(pns.pr.dps, .parallel = TRUE)
# saveRDS(aba, "/Users/zhengz11/myscripts/data_results/200710-PNKCrevision_optimal_cluster/200715-hemibrain_upn_calyx_nblast_aba.rds")

colnames(aba) = pn_types
rownames(aba) = pn_types

# write.csv(dm, paste0(save_path,"200723-hemibrain_pn_CalyxSkeleton_DistMatrix.csv"))
# write.csv(pn_types[1:110], paste0(save_path,"200723-hemibrain_pn_CalyxSkeleton_DistMatrix_names.csv"), row.names = F)

hm.pns.hc <- nhclust(scoremat=aba, method = "ward.D2")
# splot(pns.hc)

dend <- as.dendrogram(hm.pns.hc) %>% set("labels_col", col_set[labels(.)])
plot(dend)

t5 = get_leaves_attr(dend, "label")
n_pch = rep(NA, length(t5))
n_col = rep(NA, length(t5))
n_pch[t5 %in% comm_gloms] = 19
n_col[t5 %in% comm_gloms] = 3

set(dend, "leaves_pch", n_pch) %>% set("leaves_col", n_col) %>% plot(main = "Hemibrain PNs")

# 210529_NBLAST_hemibrain_pns.png
# old: 1650x450 -> new: 1580x350
