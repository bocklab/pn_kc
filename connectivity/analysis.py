
# local_path = "/Users/zhengz11/myscripts/git_clone/pn_kc/pn_kc/"
# path = local_path + "test/"

exec(open(local_path + "/connectivity/load_pn_metadata_v2.py").read())


pn_skids = load_json(local_path + "data/skids/pn")
rd = load_json(local_path + "data/skids/RandomDraw")
t1p = load_json(local_path + "data/skids/t1p")
bundle = load_json(local_path + "data/skids/bundle")

ana_all_rd = ar.Analysis.init_connectivity(local_path + "data/", pn_skids, rd + t1p, 'pn_all_kc')

ana_rd = ar.Analysis.init_connectivity(local_path + "data/", pn_skids, rd, 'pn_rd_kc')

exec(open(local_path + "/connectivity/load_pn_tbl.py").read())
