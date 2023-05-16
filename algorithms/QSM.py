# -*- coding: utf-8 -*-
"""
Created on Fri May  5 10:59:51 2023

@author: sunney
"""

import numpy as np 
import copy
import os
from pytictoc import TicToc
from algorithms.QSM_utils.QSM_functions import DoublyLinkedList,three_dimensional,checkClosest
from algorithms.TRAP import remove_points


def write_a_txt(data,instance,types,Delta,dir_path):
    data_unique = remove_points(data) if type(data) != tuple else data
    file_path = "Results\\new_output_{}_{}_{}_{}.txt".format(instance,Delta[0], Delta[1],Delta[2])
    output_file = os.path.join(dir_path, file_path)
    with open(output_file,"{}".format(types)) as txt_file:
        if type(data)!=tuple:
            for line in data_unique:
                txt_file.write(str(line[0]) + ',' + str(line[1]) + ',' + str(line[2]) + "\n") # works with any number of elements in a line
        else:
            txt_file.write("u= " + str(data[0]) + ',' + str(data[1]) + "\n")


def find_localUpper(NDP_GTRAP_Delta,E, g_1_max, g_2_max):
    """ Given a set of NDP \delta optimal with g_3=0 find local upper bounds"""
    u = {}
    for i in E:
        if i==1:
           u[1] = [round(NDP_GTRAP_Delta[i-1][0],2),g_2_max ] 
        elif i<=np.shape(NDP_GTRAP_Delta)[0]:
            u[i] = [round(NDP_GTRAP_Delta[i-1][0],2),NDP_GTRAP_Delta[i-2][1] ] 
        else:
            u[i] = [g_1_max,NDP_GTRAP_Delta[i-2][1] ] 
    return u





class qsm:

    def __init__(self, params,stock_limit,NDP_FirstStage,S_int={}, dir_path ="C:\\Users\\sunney\\Desktop\\PaperV\\PQSM_alpha",TRAP_time=0 ):
        self.newPart = params.newPart  # store the part that is newly introduced
        self.time    = params.time  # store the time a new part is introduced
        self.Work_Center_ID = params.Work_Center_ID       # Machines
        self.JobTypeData = params.JobTypeData     # Operations for each part
        self.PartType = params.PartType
        self.Capacity = params.Capacity
        self.FeasibleMachines = params.FeasibleMachines
        self.QualifiedMachines = params.QualifiedMachines
        self.Pijk = params.Pijk
        self.Demand = params.Demand
        self.QualCost = params.QualCost
        self.Delta = params.Delta 
        self.instance = params.instance
        self.T = params.T
        self.Tau = params.Tau
        self.Threshold = params.Threshold
        self.Gamma = params.Gamma
        self.stock_limit = stock_limit
        self.S_int = S_int # intial feasible solution guess
        self.d_o = params
        self.NDP_FirstStage = NDP_FirstStage
        self.dir_path= dir_path
        self.TRAP_time = TRAP_time
    
    def optimize(self):
        file_path = "Results\\new_output_{}_{}_{}_{}.txt".format(self.instance,self.Delta[0], self.Delta[1],self.Delta[2])
        output_file = os.path.join(self.dir_path, file_path)

        # Optimization code goes here
                
        g_1_max = self.T*(1-self.Threshold) ; g_2_max = max(self.QualCost.values())*self.T*self.Gamma ; 
    
        bb = np.zeros((np.shape(self.NDP_FirstStage)[0],np.shape(self.NDP_FirstStage)[1]+1))
        bb[:,:-1]=copy.deepcopy(self.NDP_FirstStage)
        NDP_GTRAP = np.zeros((np.shape(self.NDP_FirstStage)[0],np.shape(self.NDP_FirstStage)[1]+1)) 
        NDP_GTRAP[:,:-1]=self.NDP_FirstStage 
        NDP_GTRAP= NDP_GTRAP[np.lexsort((NDP_GTRAP[:,0], -NDP_GTRAP[:,1] ))] # B^{g_3}_{eff}(0)
        
        # -- Find Delta optimal efficient solutions------------ #
        
        NDP_GTRAP = NDP_GTRAP - np.array([[self.Delta[0], self.Delta[1], 0]]) 
        f = lambda x: round(max(0,x),2)
        NDP_GTRAP_Delta = np.array([list(map(f, NDP_GTRAP[i,:])) for i in range(0,np.shape(NDP_GTRAP)[0]) ])
        NDP_GTRAP_Delta = remove_points(NDP_GTRAP_Delta) # NDPs for the GTRAP (2-dimensional)
        NDP_GTRAP_Delta= NDP_GTRAP_Delta[np.lexsort((NDP_GTRAP_Delta[:,0], -NDP_GTRAP_Delta[:,1] ))] 
    
        K = np.shape(NDP_GTRAP_Delta)[0]
        E = range(2,K+2) if NDP_GTRAP_Delta[0][0]==0 else range(1,K+2)
        upper=find_localUpper(NDP_GTRAP_Delta,E, g_1_max, g_2_max) 
        u_infeas = [] # set of Q(u) which are empty or previous attempt has not resulted in infeasiblity
        Time_QSM ={}  ; tss= TicToc()
        NDP_QSM ={}
        tss.tic()
        
        for i in upper.keys():
            
            D = DoublyLinkedList()
            D.insert_at_start((upper[i][0],upper[i][1]))
            L = []                              # initialize an empty list of non dominated points
            exits=0
            L_copy =[]
            while D.start_node is not None and exits ==0:
                right_boundary_n_treated = True
                while right_boundary_n_treated == True and exits==0:
                    if D.start_node is not None:
                        u = copy.deepcopy(D.start_node.item)
                        print('--------u={}--------'.format(u))
                        if u in u_infeas:
                            D.delete_at_start()
                            break
                    else:
                        break
                    D.delete_at_start()    # pop the first element of D and denote it by u
                    z,time,num = three_dimensional(self.d_o,self.stock_limit,u)  # find the ndp 
                    write_a_txt(u,self.d_o.instance,'a+',self.Delta,self.dir_path) if os.path.isfile(output_file)==True else write_a_txt(u,self.instance,'w',self.Delta,self.dir_path)
                    
                    if z is None:
                        right_boundary_n_treated= False
                        u_infeas.append(u)
                    else:
                        sd = [z]+ L; write_a_txt([z],self.instance,'a+',self.Delta,self.dir_path)
                        if ((z not in L) and (z in sd)):
                            L.append(z)
                            if z not in L_copy:
                             z_new = [round(z[0],2), z[1], z[2]]
                             L_copy.append(z_new)
                             
                            # CheckClosest
                            
                            # if proj=0 i.e. axis g3
                            # check closes in g2
                        u_r_return,u_t_return   = checkClosest(upper,z,i,self.Delta)
                        
                        if (u_t_return ==True and z[0] -self.Delta[0] >=0):
                            if ((D.start_node is None) or (D.start_node.item[0] < z[0] -self.Delta[0])):
                                D.insert_at_start((z[0]-self.Delta[0], u[1]))
                        if (u_r_return ==True and z[1] -self.Delta[1] >=0):
                            D.insert_at_start((u[0], z[1]-self.Delta[1])) 
                        else:
                            right_boundary_n_treated= False
    
                top_boundary_n_treated = True
    
                while top_boundary_n_treated == True and exits==0:
                    D.reverse_linked_list()
                    if D.start_node is not None:
                        u = D.start_node.item
                        if u in u_infeas:
                            D.reverse_linked_list()
                            D.delete_at_end()
                            break
                    else:
                        break
                    D.reverse_linked_list()
                    D.delete_at_end()
                    
    
                    z,time,num = three_dimensional(self.d_o,self.stock_limit,u) # find the ndp
                    write_a_txt(u,self.instance,'a+',self.Delta,self.dir_path) if os.path.isfile(output_file)==True else write_a_txt(u,self.instance,'w',self.Delta,self.dir_path)
    
                    if z is  None:
                        top_boundary_n_treated = False
                        u_infeas.append(u)
                    else:
                            # update_u is already False now
                        sd = [z]+ L; write_a_txt([z],self.instance,'a+',self.Delta,self.dir_path)
                        if ((z not in L) and (z in sd)):
                        
                            L.append(z)
                            if z not in L_copy:
                                z_new = [round(z[0],2), z[1], z[2]]
                                L_copy.append(z_new)
                        D.reverse_linked_list()
                        if D.start_node is not None:
                            back_node = copy.deepcopy(D.start_node.item) 
                        D.reverse_linked_list()
                        mm = z[1]-self.Delta[1] 
                        
                        u_r_return,u_t_return   = checkClosest(upper,z,i,self.Delta)
                        if (u_r_return ==True and z[1] -self.Delta[1] >=0):
                            
                            if ((D.start_node is None) or (back_node[1] < mm)):
                                D.insert_at_end((u[0], z[1]-self.Delta[1])) 
                                
                        if (u_t_return ==True and z[0] -self.Delta[0] >=0):
                            if z[0]-self.Delta[0] > 0:
                                D.insert_at_end((z[0]-self.Delta[0], u[1])) 
                        else:
                            top_boundary_n_treated = False
                                        
            #Time_QSM = {}
            Time_QSM[i] =tss.tocvalue(restart='True')
            NDP_QSM[i]= L_copy
    
        NDP_unique=[]
        for j in NDP_QSM.keys():
          if len(NDP_QSM[j])!=0:  
            NDP_unique = NDP_unique + NDP_QSM[j] 
        NDP_unique = np.append(NDP_unique,bb,axis=0)
        NDP_unique = remove_points(NDP_unique) 
        FinalTime= max(Time_QSM.values()) +self.TRAP_time
        
        with open(output_file,'a+') as txt_file:
            txt_file.write("Time= " + str(FinalTime)  + "\n")
        write_a_txt(NDP_unique,self.instance,'a+',self.Delta,self.dir_path) if os.path.isfile(output_file)==True else write_a_txt(NDP_unique,self.instance,'w',self.Delta,self.dir_path)
    
        return NDP_unique, FinalTime
        

    
    
        