# -*- coding: utf-8 -*-
"""
Created on Fri May  5 10:48:59 2023

@author: sunney


PQSM_alpha/
│
├── algorithms/
│   ├── __init__.py
│   ├── QSM.py
│   └── ...
│
├── datasets/
│   ├── Instances/
│   │   ├── 1.csv
│   │   ├── 2.csv
│   │   └── ...
│   ├── __init__.py
│   └── dataset.py
│   └── Constant_data.csv
├── config.yaml
│
└── main.py

"""
import yaml
import copy
import time
from datasets.dataset import Dataset
from common_utils.heuristic import starting_heuristic as heur
from algorithms.TRAP import solve_TRAP,remove_points
from algorithms.QSM import qsm


# Load configuration file
with open('config.yaml', 'r') as f:
    config = yaml.safe_load(f)
f.close()

for instance in range(6,7):
    d_o = Dataset(instance,Delta=[0.1,2,50],parms=config['dataset_parameters'])
    d_o.load_data()
    
    
    start_time = time.time()
    S_int = heur(d_o)
    dir_path= "C:\\Users\\sunney\\Desktop\\PaperV\\PQSM_alpha"
    
    # First compute the minimum value for f_3 and then solve TRAP
    try:
        m = solve_TRAP(d_o, S_int, bi=False)        # \hat{D}
        m.optimize()
        TRAP_time = copy.deepcopy(time.time()) - start_time
        print("TRAP_time:", TRAP_time)
    except Exception as e:
        print("Error occurred during TRAP optimization:", e)
    
    
    NDP_FirstStage = remove_points(m.ndp) 
    stock_limit=1.5
    QSM_instance = qsm(d_o,stock_limit,NDP_FirstStage,S_int, dir_path,TRAP_time )
    
    NDP_unique, FinalTime =QSM_instance.optimize()

