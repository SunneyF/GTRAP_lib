# -*- coding: utf-8 -*-
"""
Created on Tue May  9 14:41:36 2023

@author: sunney
"""
import copy
import gurobipy as grb


# -----------------------------------------------------------------------------

def findone_F(S):
    ones = {}
    for key in S.keys():
        ones[key] = S[key].X
    return ones

def starting_heuristic(d_o):

    '''
    Starting heuristic as mentioned in https://doi.org/10.1111/itor.13180 Alg.1
    '''    
    S_0 = {}
    obj_0=0
    
    NQualCost= copy.deepcopy(d_o.QualCost)
    
    for t in range(1,d_o.T+1):
        if t== d_o.time:
            longe = len(d_o.JobTypeData[d_o.newPart])
            gamma=longe
            S_t,V_t,obj_t =solve_t(t,NQualCost,gamma,d_o)
        
        S_t,V_t,obj_t =solve_t(t,NQualCost,d_o.Gamma,d_o)
        
        if len(S_t)==0:
            return S_0
        else:
            for (i,j,k) in S_t:
                S_0[i,j,k,t] = S_t[i,j,k]
            for (i,j,k) in V_t:
                if V_t[i,j,k]==1:
                    NQualCost[i,j,k] =0
            obj_0 = obj_0 + obj_t
            S_new = {key: value for key, value in S_0.items() if value > 0.001}

    return S_new


def solve_t(t,NQualCost,gamma,d_o):
    
    model_t = grb.Model("TRAP_time_t")
    
    x_t= {} # integer variable for the number of orders
    s_t={}  # binary variable 
    vv_t={} # qualification variable
    
    # variables
    x_t = model_t.addVars(grb.tuplelist([(i,j,k) for i in d_o.PartType for j in d_o.JobTypeData[i] for k in d_o.FeasibleMachines[i,j] if d_o.Demand[i,t] !=0]), vtype=grb.GRB.INTEGER,name='x_t')
    s_t = model_t.addVars(grb.tuplelist([(i,j,k) for (i,j,k) in x_t]), vtype=grb.GRB.BINARY, name='s_t')
    vv_t = model_t.addVars(grb.tuplelist([(i,j,k) for i in d_o.PartType for j in d_o.JobTypeData[i] for k in d_o.FeasibleMachines[i,j] if k not in d_o.QualifiedMachines[i,j]]), vtype=grb.GRB.BINARY, name='vv_t')
    o = model_t.addVar(vtype=grb.GRB.CONTINUOUS,lb=0.0, ub=1-d_o.Threshold, name="o")
    model_t.setObjective(o + sum(NQualCost[i,j,k]*vv_t[i,j,k] for (i,j,k) in vv_t),grb.GRB.MINIMIZE)    

    # constraints #############################################################
    
    # 1: Demand of each job type
    model_t.addConstrs((grb.quicksum(x_t[i,j,k] for k in d_o.FeasibleMachines[i,j])  == d_o.Demand[i,t]  for i in d_o.PartType for j in d_o.JobTypeData[i] if d_o.Demand[i,t] !=0), name='flow')
    
    # 2: Setting the auxillary variables 's'
    model_t.addConstrs(d_o.Demand[i,t]*s_t[i,j,k] >= x_t[i,j,k] for (i,j,k) in s_t)
    
    # 3:  Tau constraints
    model_t.addConstrs(grb.quicksum(s_t[i,j,k] for k in d_o.FeasibleMachines[i,j] if (i,j,k) in s_t) <= d_o.Tau for i in d_o.PartType for j in d_o.JobTypeData[i]) 
       
    # 4: Min-max excess resource loading
    model_t.addConstrs((1/d_o.Capacity[k,t])*grb.quicksum(d_o.Pijk[i,j,k]*x_t[i,j,k] for i in d_o.PartType for j in d_o.JobTypeData[i] if (i,j,k) in x_t ) - d_o.Threshold <= o for k in d_o.Work_Center_ID)
    
    # 5: # of qualifications
    model_t.addConstr(sum(vv_t[i,j,k] for i in d_o.PartType for j in d_o.JobTypeData[i] for k in d_o.FeasibleMachines[i,j] if (i,j,k) in vv_t.keys())  <= gamma, name = "Budget")
    
    # -------------------------------------------------------------------------
    
    model_t.Params.Method=2
    model_t.Params.MIPGap = 0.005
    model_t.Params.OutputFlag=1
    model_t.setParam('TimeLimit', 500)   
    
    model_t.optimize()
    status = model_t.status
    
    if model_t.solCount==0:
        return {},{},0
    
    if status!=4: # neither infeasible or unbounded 
        S_0 = findone_F(s_t)
        V_0= findone_F(vv_t)
    else:
        print('error')
    return S_0,V_0,model_t.objval
