# -*- coding: utf-8 -*-
"""
Created on Fri May  5 11:39:38 2023

@author: sunney
"""
import numpy as np
from csv import reader
from pathlib import Path
from typing import Dict

def find_new_product(QualifiedMachines,Demand, PartType, T):
    '''
    finds a new product
    finds the first time period it is ordered
    '''
    new_part = next((i for (i,j) in QualifiedMachines if not QualifiedMachines[i,j]), None)
    time = next((t for t in range(1, T+1) if Demand[new_part, t] > 0), None)
    return new_part, time

class Dataset:
    """Loads data for a production planning problem instance."""
    
    def __init__(self,instance:int,Delta:list,parms: Dict[str, float] = None):   
        """Initializes the Dataset instance.
    
        Args:
            instance: The name of the instance to load.
            parms: Dictionary of optional parameters to use when loading the data. Extract from config files
                T: The number of time periods.The actual data can have more time period but you choose 1,..,T
                Tau: Number of set-ups in each time period
                Threshold: This is the threshold used for loading levels
                Gamma: Maximum number of qualifications
        """
        self.newPart = expando()  # store the part that is newly introduced
        self.time    = expando()  # store the time a new part is introduced
        self.Work_Center_ID = []       # Machines
        self.JobTypeData = {}     # Operations for each part
        self.PartType = []
        self.Capacity = {}
        self.FeasibleMachines = {}
        self.QualifiedMachines = {}
        self.Pijk = {}
        self.Demand = {}
        self.QualCost = {}
        self.Delta = Delta 
        self.instance = instance
        self.T = int(parms.get('T',8))
        self.Tau = int(parms.get('Tau',3))
        self.Threshold = parms.get('Threshold',0.7)
        self.Gamma = int(parms.get('Gamma',4))
        
    def load_data(self):    
        """Loads data from files.""" 
        data_path = Path('datasets') / ('Instances//{}.csv'.format(str(self.instance)))       
        
        with open(Path('datasets') / 'constant_data.csv', 'r') as read_obj:
            csv_reader = reader(read_obj)
            switch=0
            for row in csv_reader:
                if len(row)!=0:
                    if row[0][0]=='#':
                        if row[0][1:5]== 'Capa':
                            switch=1
                            
                            continue
                        elif row[0][1:5]=='Feas':
                            switch=2
                            
                            continue
                        elif row[0][1:5]=='Part':
                            switch=3
                            
                        else:
                            switch =4              # Qualified Machines
                            
                            continue
                    else:
                        if switch==3:
                            self.PartType = [int(r) for r in row]
                            
                        elif switch==2:
                            self.FeasibleMachines[int(row[0]),int(row[1]) ] =  [int(r) for r in row[2:] ] 
                        elif  switch==4:
                            self.QualifiedMachines[int(row[0]),int(row[1])]=   [int(r) for r in row[2:] ]  
                        elif  switch==1:
                                   self.Capacity[int(row[0]),int(row[1])]=     float([r for r in row[2:]  ][0])  
                                      
        with open(data_path , 'r') as read_obj:      
            csv_reader = reader(read_obj)
            switch=0
            for row in csv_reader:
                if len(row)!=0:
                    if row[0][0]=='#':
                        if row[0][1:3]== 'Pr':
                            switch=1
                        elif row[0][1:3]=='De':
                            switch=2
                        else:
                            switch=3
                    else:
                        if switch==1:
                            self.Pijk[int(row[0]), int(row[1]), int(row[2])] = float(row[3])
                        elif switch==2:
                            self.Demand[int(row[0]), int(row[1])] = float(row[2]) 
                        else:
                            if (int(row[0]), int(row[1]),int(row[2])) in self.QualCost.keys():
                                self.QualCost[int(row[0]), int(row[1]),int(row[2]) ]= min(float(row[3]), self.QualCost[int(row[0]), int(row[1]),int(row[2]) ])
                            else:
                                self.QualCost[int(row[0]), int(row[1]),int(row[2]) ] = float(row[3])
                 
        self.Work_Center_ID =  np.sort(np.unique([k for (i,j,k) in self.Pijk.keys()])) 
        
        self.JobTypeData = {i:[] for i in self.PartType}
        
        for (i,j,k) in self.Pijk.keys():
            self.JobTypeData[i]=self.JobTypeData[i] + [j]
        self.JobTypeData = {i: (list(dict.fromkeys(v))) for i, v in self.JobTypeData.items()}
        
        for i in self.JobTypeData:
            self.JobTypeData[i].sort()
        new,tt = find_new_product(self.QualifiedMachines,self.Demand, self.PartType, self.T)
        self.newPart = new 
        self.time = tt
        
class expando(object):
    pass       
        

