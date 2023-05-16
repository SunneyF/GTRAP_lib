import gurobipy as grb

def full_model(d_o,BM, e_bar,relax='false',stock_limit=1.5):
    
    # define the full model to find feasible solution and compare lower bounds
   
    
    model = grb.Model()
    
    x = model.addVars(grb.tuplelist([(i,j,k,t) for i in d_o.PartType for j in d_o.JobTypeData[i] for k in d_o.FeasibleMachines[i,j] for t in range(1,d_o.T+1) ]  ), vtype=grb.GRB.INTEGER if relax!='true' else grb.GRB.CONTINUOUS ,name='x')                                                                   # x variables
    s = model.addVars(grb.tuplelist([(i,j,k,t) for i in d_o.PartType for j in d_o.JobTypeData[i] for k in d_o.FeasibleMachines[i,j] for t in range(1,d_o.T+1) ]  ), vtype=grb.GRB.INTEGER if relax!='true' else grb.GRB.CONTINUOUS ,name='s')                                                                   # s variables
    vv = model.addVars(grb.tuplelist([(i,j,k,t) for i in d_o.PartType for j in d_o.JobTypeData[i] for k in d_o.FeasibleMachines[i,j] for t in range(1,d_o.T+1) if k not in d_o.QualifiedMachines[i,j]]), vtype=grb.GRB.BINARY,name='vv')                                                                        # z variables
    
    e =  model.addVars(grb.tuplelist([(i,j,t) for i in d_o.PartType for j in [0]+d_o.JobTypeData[i]  for t in range(0,d_o.T+1) ]  ), vtype=grb.GRB.INTEGER if relax!='true' else grb.GRB.CONTINUOUS ,name='e')                                                                                                  # x variables

    o = model.addVars(grb.tuplelist(t  for t in range(1,d_o.T+1)), vtype=grb.GRB.CONTINUOUS, ub=1-d_o.Threshold,name='o')
    
    m = model.addVars(grb.tuplelist([(i,t) for i in d_o.PartType for t in range(1,d_o.T+1) ] ), vtype=grb.GRB.INTEGER,name='m')
    
    # constrant upper limit on number of raw materials that can be ordered
    #model.addConstrs(m[i,t] <= 1.5*sum(d_o.Demand[i,t] for t in range(t,d_o.T+1) if (i,t) in d_o.Demand.keys()) for i in d_o.PartType for t in range(1,d_o.T+1))
    model.addConstrs(m[i,t] <= stock_limit*d_o.Demand[i,t]  for i in d_o.PartType for t in range(1,d_o.T+1))
    
    
    # constraint specific to QSM
    #model.addConstr(sum(o[t] for t in range(1,T+1)) >= min_x)

    # constraint 2(a): flow balance constraint    
    model.addConstrs((grb.quicksum(x[i,j,k,t] for k in d_o.FeasibleMachines[i,j]) + e[i,j,t-1] == d_o.Demand[i,t] + e[i,j,t]  for i in d_o.PartType for j in d_o.JobTypeData[i] for t in range(1,d_o.T+1) ), name='flow')
    
    # constraint 2(b): flow balance for raw material
    
    model.addConstrs((m[i,t] + e[i,0,t-1] == d_o.Demand[i,t] + e[i,0,t]  for i in d_o.PartType for t in range(1,d_o.T+1) ), name='flow_raw')

    
    # Constraint 2(c): constraint on r^{l}_{jt} \geq 0
    
    for i in d_o.PartType:
     cs=0
     for j in d_o.JobTypeData[i]:
        if j!=d_o.JobTypeData[i][-1]:
            model.addConstrs(e[i,j,t] >= e[i,d_o.JobTypeData[i][cs+1],t] for t in range(0,d_o.T+1) )
        else:
            break
        cs =cs +1
    
    for i in d_o.PartType:
        model.addConstrs(e[i,0,t] >= e[i,d_o.JobTypeData[i][0],t] for t in range(0,d_o.T+1) )

        
    
    # Constraint 2(e) : fixing initial inventory
    model.addConstrs(e[i,j,0]==e_bar[i,j]  for i in d_o.PartType for j in [0]+d_o.JobTypeData[i] )
    
    
    # Constraint 1(e): big M
    model.addConstrs((x[i,j,k,t] <= BM[i,j,k,t]*s[i,j,k,t] for i in d_o.PartType for j in d_o.JobTypeData[i] for k in d_o.FeasibleMachines[i,j] for t in range(1,d_o.T+1) ))
   
   
    # Constraint 1(f): tau constraint
    model.addConstrs((grb.quicksum(s[i,j,k,t] for k in d_o.FeasibleMachines[i,j] if (i,j,k,t) in s) <= d_o.Tau for i in d_o.PartType for j in d_o.JobTypeData[i] for t in range(1,d_o.T+1) ), name='route')
    

       
    #  constraint 1(k): qualification cost variable
    model.addConstrs((grb.quicksum(vv[i,j,k,l] for l in range(1,t+1))>= s[i,j,k,t] for (i,j,k,t) in s if k not in d_o.QualifiedMachines[i,j]), name='z_var')
    
    
    # constraint 1(g)
    
    model.addConstrs( (1/d_o.Capacity[k,t])*sum(d_o.Pijk[i,j,k]*x[i,j,k,t] for i in d_o.PartType for j in d_o.JobTypeData[i] if (i,j,k,t) in x) - d_o.Threshold <= o[t] for k in d_o.Work_Center_ID for t in range(1,d_o.T+1) )
    

    # constraint 1(i)
    model.addConstrs(sum(vv[i,j,k,t] for i in d_o.PartType for j in d_o.JobTypeData[i] for k in d_o.FeasibleMachines[i,j] if (i,j,k,t) in vv) <= d_o.Gamma for t in range(1,d_o.T+1))

    
    model.setParam('TimeLimit', 1500)
    model.setParam('MIPGap',1e-4)
    
    model.setParam('MIPFocus',1)

    
    model._x=x
    model._s=s
    model._vv = vv
    model._o = o
    model._e=e            
            
    return model