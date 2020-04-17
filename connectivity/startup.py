import sys
import os
sys.path.append('/Users/zhengz11/myscripts/git_clone/pn_kc/')

# local path to the code
local_path = "/Users/zhengz11/myscripts/git_clone/pn_kc/"
# data_path = local_path + "test/skids/"

# exec(open(local_path + "/connectivity/startup.py").read())




# import things in the mushroom module
# I use it so ubiquitously that I just import everything
import mushroom_2to3.connect as cc
import mushroom_2to3.analysis_routine as ar

from mushroom_2to3.shuffle import *
from mushroom_2to3.build_connectivity import *
from mushroom_2to3.detect_community import *



# Only if you have token and want to connect to CATMAID
# import your personal token
# token = "xxxxxxxxxx"
# sys.path.append('/Users/zhengz11/myscripts/mushroom_v9/credential/')
# from fafb_tokens import token
# fafb_c = cc.fafb_connection(token)
