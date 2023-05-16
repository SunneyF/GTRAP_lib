import copy
import gurobipy as grb
import numpy as np



def remove_points(ndp_preProcessed):
    
    """ removes dominated points and any duplicates """
    ndp_preProcessed = np.round(ndp_preProcessed, 2)
    NDP = np.unique(ndp_preProcessed, axis=0)
    if ndp_preProcessed[0][0] == -1 and ndp_preProcessed[0][1] == -1:
        NDP = np.delete(NDP, 0, axis=0)
    NDP = keep_efficient(-NDP)
    NDP = -NDP
    return NDP

def keep_efficient(pts):
    'returns Pareto efficient row subset of pts'
    # sort points by decreasing sum of coordinates
    pts = pts[pts.sum(1).argsort()[::-1]]
    # initialize a boolean mask for undominated points
    # to avoid creating copies each iteration
    undominated = np.ones(pts.shape[0], dtype=bool)
    for i in range(pts.shape[0]):
        # process each point in turn
        n = pts.shape[0]
        if i >= n:
            break
        # find all points not dominated by i
        # since points are sorted by coordinate sum
        # i cannot dominate any points in 1,...,i-1
        undominated[i+1:n] = (pts[i+1:] > pts[i]).any(1)  ##
        # keep points undominated so far
        pts = pts[undominated[:n]]
    return pts


# -----------------------------------------------------------------------------

def findone_F(S):
    ones = {}
    for key in S.keys():
        ones[key] = S[key].X
    return ones

def termination_callback(model,where):
    """
    if time 1000 seconds and there has been no improvement in objective value for 50 nodes then exit
    """
    
    if where == grb.GRB.Callback.MIP:
        nodecnt = model.cbGet(grb.GRB.Callback.MIP_NODCNT)
        if nodecnt - model._lastNode >= 10:
            
            model._lastNode = nodecnt
            runtime = model.cbGet(grb.GRB.Callback.RUNTIME)
            if ((abs(model.cbGet(grb.GRB.Callback.MIP_OBJBST) - model._curr_obj)<=0.05) & (runtime >=1000)):
                model.terminate()
                
            else:
                model._curr_obj = model.cbGet(grb.GRB.Callback.MIP_OBJBST)
                
class expando(object):
    pass

class solve_TRAP:
    
    def __init__(self, d_o, S_int={},bi = True):  # bi-> bidirectional  or False
    
        self.d_o=d_o
        self.S_int = S_int
        self.variables = expando()
        self._build_model()
        self.sol_times=0  # time to solve the TRAP i.e B^{g_3}_eff(0)
        self.ndp = np.array([[-1,-1]])
        self.count_models = 0 # records # of scalarizations solved including infeasible ones
        self.count_solved_models =0
        self.bi = bi          # if True use bi-directional epsilon constraint otherwise just apply constraints on g_2
        
    def _build_model(self):
        
        self.model = grb.Model()
        self._build_variables()
        self._build_objective()
        self._build_constraints()
        self.model.Params.PrePasses=1
        self.model.Params.Method=3
        self.model.Params.MIPGap = .0001
        self.model.Params.MIPGapAbs=.005
        self.model.Params.Timelimit = 1500
        self.model.Params.OutputFlag=1
        
        
        
    def _build_variables(self):
        
        self.variables.x = self.model.addVars(grb.tuplelist([(i,j,k,t) for i in self.d_o.PartType for j in self.d_o.JobTypeData[i] for k in self.d_o.FeasibleMachines[i,j] for t in range(1,self.d_o.T+1) if (i,t) in self.d_o.Demand]), vtype=grb.GRB.INTEGER,name='x')
        self.variables.s = self.model.addVars(grb.tuplelist([(i,j,k,t) for (i,j,k,t) in self.variables.x]),vtype= grb.GRB.BINARY,name='s') 
        self.variables.vv = self.model.addVars(grb.tuplelist([(i,j,k,t) for i in self.d_o.PartType for j in self.d_o.JobTypeData[i] for k in self.d_o.FeasibleMachines[i,j] for t in range(1,self.d_o.T+1) if k not in self.d_o.QualifiedMachines[i,j]]), vtype=grb.GRB.BINARY,name='vv')
        self.variables.o = self.model.addVars(grb.tuplelist([t for t in range(1,self.d_o.T+1)]), vtype=grb.GRB.CONTINUOUS, ub=1-self.d_o.Threshold, lb=0, name='o') # make ub = 1-threshold
        self.model.update()
        
    def _build_objective(self,switch=0): 
        
        """
        Parameters
        ----------
        switch : INT, optional
            DESCRIPTION. The default is 0. S
            switch=['g1','g2'] --> priority g1->g2->g3
            switch = ['g2','g1'] -->priority g2->g1->g3

        Returns
        -------
        None.

        """
        
        max_g1 = self.d_o.T*(1-self.d_o.Threshold) ; max_g2 = max(self.d_o.QualCost.values())*self.d_o.Gamma*self.d_o.T
        func1 =  grb.quicksum(self.variables.o[t] for t in range(1,self.d_o.T+1)) + (.01/max_g2)*grb.quicksum(self.d_o.QualCost[i,j,k]*self.variables.vv[i,j,k,t] for (i,j,k,t) in self.variables.vv) if switch==0 else grb.quicksum(self.d_o.QualCost[i,j,k]*self.variables.vv[i,j,k,t] for (i,j,k,t) in self.variables.vv) + (1/max_g1)*grb.quicksum(self.variables.o[t] for t in range(1,self.d_o.T+1))

        self.model.setObjective(func1 , grb.GRB.MINIMIZE)
        
        
    def _build_constraints(self):
        
        # Constraint 2(a): Demand of each job type must be met
        self.model.addConstrs((grb.quicksum(self.variables.x[i,j,k,t] for k in self.d_o.FeasibleMachines[i,j]) == self.d_o.Demand[i,t] for i in self.d_o.PartType for j in self.d_o.JobTypeData[i] for t in range(1,self.d_o.T+1) if (i,t) in self.d_o.Demand))
        
        # Constraint 2(b): setting s variable
        self.model.addConstrs((self.d_o.Demand[i,t]*self.variables.s[i,j,k,t] >= self.variables.x[i,j,k,t] for (i,j,k,t) in self.variables.s))        
        
        # Constraint 2(c): tau constraints
        self.model.addConstrs((grb.quicksum(self.variables.s[i,j,k,t] for k in self.d_o.FeasibleMachines[i,j] if (i,j,k,t) in self.variables.s) <= self.d_o.Tau for i in self.d_o.PartType for j in self.d_o.JobTypeData[i] for t in range(1,self.d_o.T+1) if (i,t) in self.d_o.Demand), name='route')
        
        # Constraint 2(d) : min-max 
        self.model.addConstrs(((1/self.d_o.Capacity[k,t])*grb.quicksum(self.d_o.Pijk[i,j,k]*self.variables.x[i,j,k,t] for i in self.d_o.PartType for j in self.d_o.JobTypeData[i] if (i,j,k,t) in self.variables.x ) - self.d_o.Threshold <= self.variables.o[t] for k in self.d_o.Work_Center_ID for t in range(1,self.d_o.T+1)) )
        
       
        # constraint 2(e): gamma limitation 
        self.model.addConstrs((grb.quicksum(self.variables.vv[i,j,k,t] for i in self.d_o.PartType for j in self.d_o.JobTypeData[i] for k in self.d_o.FeasibleMachines[i,j] if (i,j,k,t) in self.variables.vv.keys())  <= self.d_o.Gamma for t in range(1,self.d_o.T+1)))

       
        #  constraint 2(f): qualification cost variable
        self.model.addConstrs((grb.quicksum(self.variables.vv[i,j,k,l] for l in range(1,t+1))>= self.variables.s[i,j,k,t] for (i,j,k,t) in self.variables.s if k not in self.d_o.QualifiedMachines[i,j]))
 
    def optimize(self):
        
        # g^Top--------------------------------------------------------
        
        self._build_objective(switch=0)
        self.model.update()
        self.model._curr_obj = grb.GRB.INFINITY
        self.model._lastNode = -grb.GRB.INFINITY
        self.model.Params.Timelimit = 1500
        
        self._starting_solution(self.S_int)
        
        self.model.optimize(termination_callback)
        
        status = self.model.status
        
        self.count_models=copy.deepcopy(self.count_models)+1 ; self.count_solved_models=copy.deepcopy(self.count_solved_models)+(1 if  status==2 else 0)
        
        if ((status!=2) & (abs(self.model.objval - self.model.objbound) > .009 )):
            
            print('inside here----------')
            print('Top')
            self.model.optimize(termination_callback); 
            self.sol_times = self.model.Runtime + self.sol_times
            status = self.model.status
            self.count_models=copy.deepcopy(self.count_models)+1 ; self.count_solved_models=copy.deepcopy(self.count_solved_models)+(1 if  status==2 else 0)

            assert self.model.solcount ==0, 'No z_top found'
            
        else:
            
            self.sol_times = self.model.Runtime + self.sol_times
        
        g_2 = sum(self.d_o.QualCost[i,j,k]*self.variables.vv[i,j,k,t].X for (i,j,k,t) in self.variables.vv)
        
        g_1 = round(sum(self.variables.o[t].X for t in range(1,self.d_o.T+1)),4) 
        
        self.ndp = np.concatenate((self.ndp, np.array([[g_1,g_2]])), axis=0)
        
        S_Top = findone_F(self.variables.s);
        
        z={}; z['T'] = [g_1, g_2]
        
        # g^Bot --------------------------------------------------
        
        self._build_objective(switch=1) # now the second objective function is to be minimized
        self.model.update()
        self.model.optimize(termination_callback);           
        status = self.model.status
        self.count_models=copy.deepcopy(self.count_models)+1 ; self.count_solved_models=copy.deepcopy(self.count_solved_models)+(1 if  status==2 else 0)

        self.sol_times = self.model.Runtime + self.sol_times
        if status!=2:
            #print('Bottom')
            self.model.optimize(termination_callback); 
            status = self.model.status

            self.count_models=copy.deepcopy(self.count_models)+1 ; self.count_solved_models=copy.deepcopy(self.count_solved_models)+(1 if  status==2 else 0)

            self.sol_times = self.model.Runtime + self.sol_times
            status = self.model.status
            assert self.model.solcount ==0, 'No feasible z_bottom found'
        else:
            self.sol_times = self.model.Runtime + self.sol_times
        g_2 = sum(self.d_o.QualCost[i,j,k]*self.variables.vv[i,j,k,t].X for (i,j,k,t) in self.variables.vv)
        g_1 = round(sum(self.variables.o[t].X for t in range(1,self.d_o.T+1)),4) 
        self.ndp = np.concatenate((self.ndp, np.array([[g_1,g_2]])), axis=0)
        S_Bot = findone_F(self.variables.s);  z['B'] = [g_1, g_2]
        self.epsilon_constraint(z,S_Top,S_Bot)
        
    def _add_new_constraints(self,z,func):
        if func=='Top':
            self.model.addConstr(grb.quicksum(self.d_o.QualCost[i,j,k]*self.variables.vv[i,j,k,t] for (i,j,k,t) in self.variables.vv) <= z['T'][1] -self.d_o.Delta[1] )
        else:
            
            self.model.addConstr(grb.quicksum(self.variables.o[t] for t in range(1,self.d_o.T+1)) <= z['B'][0]-self.d_o.Delta[0])
    
    def _starting_solution(self,S_int):
        for (i,j,k,t) in self.variables.s:
            if (i,j,k,t) in self.S_int:
                self.variables.s[i,j,k,t].start = self.S_int[i,j,k,t]
            else:
                self.variables.s[i,j,k,t].start = 0 
    
    def epsilon_constraint(self,z,S_Top={},S_Bot={}):
        
        while True:
            if (abs(z['T'][1]-z['B'][1]) <=self.d_o.Delta[1] or abs(z['B'][0]-z['T'][0]) <=self.d_o.Delta[0]):
                break
            
            func='Top' # add constraint in g_2
            
            self._add_new_constraints(z,func) # add epsilon bounds from top for g_2
            self._build_objective(switch=0)   # minimizing g_1
            self.model.update()
            self._starting_solution(S_Bot)    # starting feasible solution
            self.model._curr_obj = grb.GRB.INFINITY
            self.model._lastNode = -grb.GRB.INFINITY
            self.model.optimize(termination_callback); 
            status = self.model.status
            self.count_models=copy.deepcopy(self.count_models)+1; self.count_solved_models=copy.deepcopy(self.count_solved_models)+(1 if status==2 else 0) 
            
            self.sol_times =self.model.Runtime + self.sol_times
            if status!=2:
                self.model.optimize(termination_callback); 
                status = self.model.status

                self.count_models=copy.deepcopy(self.count_models)+1 ; self.count_solved_models=copy.deepcopy(self.count_solved_models)+(1 if  status==2 else 0)
                self.sol_times = self.model.Runtime + self.sol_times
            g_2 = sum(self.d_o.QualCost[i,j,k]*self.variables.vv[i,j,k,t].X for (i,j,k,t) in self.variables.vv)
            g_1 = round(sum(self.variables.o[t].X for t in range(1,self.d_o.T+1)),4) 
            
            S_Top = findone_F(self.variables.s);
            if  (z['B'][1]==g_2):  # repeat solution
                break
            z['T'] = [g_1, g_2]
            self.ndp = np.concatenate((self.ndp, np.array([[g_1,g_2]])), axis=0)
            
            if self.bi== True: 
                func='Bot'
                self._add_new_constraints(z,func) # add epsilon bounds from bottom
                self._build_objective(switch=1) # minimize functional_sequence[1] i.e. next bottom point
                self._starting_solution(S_Top)
                self.model.update()
                self.model._curr_obj = grb.GRB.INFINITY
                self.model._lastNode = -grb.GRB.INFINITY
                self.model.optimize(termination_callback);  
    
                status = self.model.status
    
                self.count_models=copy.deepcopy(self.count_models)+1 ; self.count_solved_models=copy.deepcopy(self.count_solved_models)+(1 if  status==2 else 0)
    
                self.sol_times = self.model.Runtime + self.sol_times 
                if status!=2:
                    self.model.optimize(termination_callback); 
                    status = self.model.status
    
                    self.count_models=copy.deepcopy(self.count_models)+1 ; self.count_solved_models=copy.deepcopy(self.count_solved_models)+(1 if  status==2 else 0)
    
                    
                    self.sol_times = self.model.Runtime + self.sol_times
                g_2 = sum(self.d_o.QualCost[i,j,k]*self.variables.vv[i,j,k,t].X for (i,j,k,t) in self.variables.vv)
                g_1 = round(sum(self.variables.o[t].X for t in range(1,self.d_o.T+1)),4) 
                
                S_Bot = findone_F(self.variables.s);
                if (z['T'][1]==g_2):  # repeat solution
                    break
                z['B'] = [g_1, g_2]
    
                self.ndp = np.concatenate((self.ndp, np.array([[g_1,g_2]])), axis=0)






#------------------------------------------------------------------------------
