
local_path = "/Users/zhengz11/myscripts/git_clone/pn_kc/"
path = local_path + "test/"

exec(open(local_path + "/connectivity/load_pn_metadata_v2.py").read())
exec(open(local_path + "/connectivity/load_pn_tbl.py").read())

pn_skids = load_json(data_path + "pn")
rd = load_json(data_path + "RandomDraw")
t1p = load_json(t1p, data_path + "t1p")
bundle = load_json(bundle, data_path + "bundle")

ana_all_rd = ar.Analysis.init_connectivity(local_path + "data/", pn_skids, rd + t1p, 'pn_all_kc')
