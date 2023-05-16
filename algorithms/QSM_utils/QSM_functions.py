# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 13:56:34 2021

@author: sunney
"""
from pytictoc import TicToc
import math
import gurobipy as grb
import numpy as np
import copy
from algorithms.load_model.model_file import full_model
# Helper functions for Quadrant Shrinking Method (QSM)




def map_keys(u_r,u_t,upper):
    d=upper
    key_r = None
    key_t = None
    min_val_r = np.inf
    min_val_t = np.inf

    # Loop through dictionary
    for key, value in d.items():
        # Check if value is less than minimum non-negative value seen so far for u_r
        if value[1] >= u_r[1] and (value[1] - u_r[1]) < min_val_r :
            min_val_r = value[1] - u_r[1]
            key_r = key
        # Check if value is less than minimum non-negative value seen so far for u_t
        if value[1] >= u_t[1] and (value[1] - u_t[1]) < min_val_t:
            min_val_t = value[1] - u_t[1]
            key_t = key

    # Print the keys with the minimum non-negative values
    # print(f"Key for u_r: {key_r}")
    # print(f"Key for u_t: {key_t}")
    return key_r,key_t

#






def checkClosest(upper,z,i,Delta):
    # for a dictionary upper identify to which quadrant z belongs to and is closest on the second axis i.e. g2 for proj=0 and g3 for proj=1
    # return u_r : upper bound on the right boundary
    # return u_t : upper bound on the top boundary
    proj=0
    u_r = np.array([upper[i][0],z[1]-Delta[1] if proj==0 else z[2]-Delta[2]])
    u_t = np.array([z[0]-Delta[0],upper[i][1]] )
    key_r,key_t= map_keys(u_r,u_t,upper)
    u_r_return = True if key_r==i else False
    u_t_return = True if key_t==i else False  

    
    return u_r_return, u_t_return


def termination_callback(model,where):
    """
    if time 1000 seconds and there has been no improvement in objective value for 50 nodes then exit
    """
    
    if where == grb.GRB.Callback.MIP:
        nodecnt = model.cbGet(grb.GRB.Callback.MIP_NODCNT)
        if nodecnt - model._lastNode >= 100:
            
            model._lastNode = nodecnt
            runtime = model.cbGet(grb.GRB.Callback.RUNTIME)
            if ((abs(model.cbGet(grb.GRB.Callback.MIP_OBJBST) - model._curr_obj)<=0.01) & (runtime >=1000)):
                model.terminate()
                
            else:
                model._curr_obj = model.cbGet(grb.GRB.Callback.MIP_OBJBST)
class Node:
    def __init__(self, data):
        self.item = data
        self.nref = None
        self.pref = None

class DoublyLinkedList:
    
    def __init__(self):
        self.start_node = None
        
    def insert_in_emptylist(self, data):
        if self.start_node is None:
            new_node = Node(data)
            self.start_node = new_node
        else:
            print("list is not empty")
            
    def insert_at_start(self, data):
        if self.start_node is None:
            new_node = Node(data)
            self.start_node = new_node
            print("node inserted")
            return
        new_node = Node(data)
        new_node.nref = self.start_node
        self.start_node.pref = new_node
        self.start_node = new_node
        
    def insert_at_end(self, data):
        if self.start_node is None:
            new_node = Node(data)
            self.start_node = new_node
            return
        n = self.start_node
        while n.nref is not None:
            n = n.nref
        new_node = Node(data)
        n.nref = new_node
        new_node.pref = n
        
    def delete_at_start(self):
        if self.start_node is None:
            print("The list has no element to delete")
            return 
        if self.start_node.nref is None:
            self.start_node = None
            return
        self.start_node = self.start_node.nref
        self.start_prev = None
        
    def delete_at_end(self):
        if self.start_node is None:
            print("The list has no element to delete")
            return 
        if self.start_node.nref is None:
            self.start_node = None
            return
        n = self.start_node
        while n.nref is not None:
            n = n.nref
        n.pref.nref = None
    

    def reverse_linked_list(self):
        if self.start_node is None:
            print("The list has no element to delete")
            return 
        p = self.start_node
        q = p.nref
        p.nref = None
        p.pref = q
        while q is not None:
            q.pref = q.nref
            q.nref = p
            p = q
            q = q.pref
        self.start_node = p    
    
    def traverse_list(self):
        if self.start_node is None:
            print("List has no element")
            return
        else:
            n = self.start_node
            while n is not None:
                print(n.item , " ")
                n = n.nref





def three_dimensional(d_o,stock_limit,u):
    
    ###########################################################################
    # Demand
    ###########################################################################
    
    tss= TicToc()
    tss.tic()
    Model_solved_count=0
    
    
    BM = {} # big M values for each product and time
    
    for (i,j,k,t) in [(i,j,k,t) for i in d_o.PartType for j in d_o.JobTypeData[i] for k in d_o.FeasibleMachines[i,j] for t in range(1,d_o.T+1)]:
        BM[i,j,k,t] = min(sum(d_o.Demand[i,t] for t in range(t,d_o.T+1) ), math.floor(d_o.Capacity[k,t]/d_o.Pijk[i,j,k])) 
    
    r_bar={};e_bar={}
    for (i,j) in [(i,j) for i in d_o.PartType for j in [0] + d_o.JobTypeData[i] ]:
        r_bar[i,j] =0
    for (i,j) in r_bar.keys():
        ind= d_o.JobTypeData[i].index(j) if j!=0 else 0
        if ind!=0:
            e_bar[i,j] = sum(r_bar[i,k] for k in d_o.JobTypeData[i][ind:]) 
        else:
            
            e_bar[i,j] = sum(r_bar[i,k] for k in d_o.JobTypeData[i]) +r_bar[i,0]
    
    ###########################################################################
    relax='false'
    model_full =full_model(d_o,BM,e_bar,relax,stock_limit)
    
    func3 = grb.quicksum(model_full._e[i,0,t] for i in d_o.PartType for t in range(1,d_o.T+1)) 

    model_full.setObjective(func3,grb.GRB.MINIMIZE)
    model_full.addConstr(grb.quicksum(model_full._o[t] for t in range(1,d_o.T+1)) <= u[0])

    model_full.addConstr(grb.quicksum(d_o.QualCost[i,j,k]*model_full._vv[i,j,k,t] for (i,j,k,t) in model_full._vv) <= u[1])
    model_full.update()
    model_full._lastNode = -grb.GRB.INFINITY
    model_full._curr_obj = grb.GRB.INFINITY
    model_full.Params.Timelimit = 1500
    model_full.optimize(termination_callback); 
    status = model_full.status
    Model_solved_count = Model_solved_count + (1 if status==2 else 0)
    sol_count = model_full.SolCount
    time=model_full.Runtime
    
    if (status==4 or status == 3 or sol_count < 1):
        z=None 
        
        return z,model_full.Runtime,Model_solved_count
    else:
        model_full.setObjective(grb.quicksum(model_full._o[t] for t in range(1,d_o.T+1)) +  grb.quicksum(d_o.QualCost[i,j,k]*model_full._vv[i,j,k,t] for (i,j,k,t) in model_full._vv), grb.GRB.MINIMIZE)
        model_full.addConstr(grb.quicksum(model_full._e[i,0,t] for i in d_o.PartType for t in range(1,d_o.T+1)) <= copy.deepcopy(model_full.objval))

        model_full._lastNode = -grb.GRB.INFINITY
        model_full._curr_obj = grb.GRB.INFINITY
        model_full.Params.Timelimit = 1500
        model_full.optimize(termination_callback); 
        Model_solved_count = Model_solved_count + (1 if status==2 else 0)
        time=model_full.Runtime+time
        func2 = round(sum(d_o.QualCost[i,j,k]*model_full.getVarByName('vv[{},{},{},{}]'.format(i,j,k,t)).x for (i,j,k,t) in model_full._vv),3)
        func3 = sum(model_full.getVarByName('e[{},0,{}]'.format(i,t)).x for i in d_o.PartType for t in range(1,d_o.T+1))
        func1 = sum(model_full.getVarByName('o[{}]'.format(t)).X  for t in range(1,d_o.T+1))
        z=[]
        z.append(func1); z.append(func2) ; z.append(func3);
        time= tss.tocvalue()
        return z, time,Model_solved_count







