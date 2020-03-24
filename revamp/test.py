
from revamp import connect_path as ct

path = "/Users/zhengz11/myscripts/git_clone/pn_kc/test/"

t_pn_skids = pn_skids
t_kc_skids = rd[:30]

t_skids = t_pn_skids + t_kc_skids


for i in t_skids:
    save_compact_sk(fafb_c, i, path)
    save_annotations_for_skeleton(fafb_c, skids, path)


save_neurons_names(fafb_c, t_skids, path)

save_root_node(fafb_c, t_skids, path)

save_annotated_annotations(fafb_c, 'glom_class', 'id', path)

save_annotated_annotations(fafb_c, 'kc_class', 'id', path)

save_pre_post_info(fafb_c, t_pn_skids, t_kc_skids, 'testing')



# run it for
'glom_class'
'kc_class'
save_annotated_annotations



# need to call for specific pairing of pre_skids, and post_skids
save_pre_post_info(connection, pre_skids, post_skids):
