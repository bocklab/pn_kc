
# need to replace with your path to the mushroom package (flyconnectome/mushroom/)
mushroom_path = '/Users/zhengz11/myscripts/mushroom_v9/mushroom/'

import sys
sys.path.append(mushroom_path)

# import things in the mushroom module
# I use it so ubiquitously that I just import everything
import mushroom_2to3.connect as cc
import mushroom_2to3.analysis_routine as ar

from mushroom_2to3.shuffle import *
from mushroom_2to3.build_connectivity import *
from mushroom_2to3.detect_community import *

# import your personal token
# token = "xxxxxxxxxx"
sys.path.append('/Users/zhengz11/myscripts/mushroom_v9/credential/')
from fafb_tokens import token
fafb_c = cc.fafb_connection(token)
