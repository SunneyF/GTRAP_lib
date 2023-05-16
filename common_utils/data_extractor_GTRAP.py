import numpy as np
import os


def data_convert(dir_path, types):

    text_select = {'PQSM':'new_output_','TDA': 'TDA_', 'GPBA':'GPBAA_'}
    start_sep = {'PQSM':2,'TDA': 1, 'GPBA':1}
    end_sep = {'PQSM':6,'TDA': 5, 'GPBA':5}
    # Create an empty dictionary to store the data
    data_dict= {}
    # Loop through all files in the directory that start with 'new_output_'
    for filename in os.listdir(dir_path):
        if filename.startswith(text_select[types]):
            # Extract the values of a, b, c, and d from the filename
            field1 = filename.split('.txt') 
            a, b, c, d = field1[0].split('_')[start_sep[types]: end_sep[types]]
            
            # Read the contents of the file
            with open(os.path.join(dir_path, filename), 'r') as f:
                lines = f.readlines()
            
            # Find the index of the line that starts with 'No. of models solved'
            index = 0
            for i, line in enumerate(lines):
                if line.startswith('Time') or line.startswith('Time='):                       
                    index = i + 1
                    break
  
            # Create a numpy array from the lines after the 'No. of models solved' line
            arrays = np.empty((0, 3))
            #print(lines)
            arrays = np.genfromtxt(lines[index:], delimiter=',')
            if np.ndim(arrays) ==1:
                arrays = [[arrays]]
            
            # Add the numpy array to the dictionary with keys a, b, c, and d
            data_dict[(a, b, c, d)] = -1*arrays if types=='GPBA' else arrays
    
    # Print the resulting dictionary
    #print(data_dict)
    return data_dict

# get solutions from all three algorithms 




def identify_pareto_front(array):
    # sort the array based on the first column
    sorted_array = array[np.argsort(array[:, 0])]
    pareto_front = []
    # iterate through the sorted array
    for i in range(len(sorted_array)):
        is_pareto_optimal = True
        # compare the current point with all previous points
        for j in range(i):
            if (sorted_array[i][1] >= sorted_array[j][1]) and (sorted_array[i][2] >= sorted_array[j][2]):
                is_pareto_optimal = False
                break
        # if the point is Pareto optimal, add it to the Pareto front list
        if is_pareto_optimal:
            pareto_front.append(sorted_array[i])
    return np.array(pareto_front)


def closest(A):
    # sort lexicographically first, second, and third
    # find the rows for which the first two satisfy desired coverage and chose the one with minimum third
    # if none satisfy the coverage restriction for first two. The output the first row
    # merge sort shorturl.at/cfIOW 
    A_new = np.append(A, np.array([range(np.shape(A)[0])]).T, axis=1)
    max_r = 10000; ind_r =-1
    for i in A_new:
        test = max(i[0:3])
        if test < max_r:
            max_r= test
            ind_r = i[3]
    return max_r,ind_r

def closest_new(max_values,min_values,alpha_temp_i_three, Delta,R_norm,check_val=True):
    # first add row index as a new column to alpha_temp_i_three

    row_indices = np.arange(alpha_temp_i_three.shape[0]).reshape(-1, 1)
    arr_with_indices = np.hstack((row_indices, alpha_temp_i_three))

    #print(arr_with_indices)
    
    # identify rows which satisfy delta condition
    
    cond = np.all(arr_with_indices <= [np.shape(arr_with_indices)[0],Delta[0], Delta[1], Delta[2]], axis=1)

    # Subset the array based on the condition
    subset = arr_with_indices[cond] if check_val==True else arr_with_indices 
    if len(subset)==0:
        return -1,-1
    # now we have the subset array with first column the row number in alpha_temp_i
    
    # now normalize the subset
    
    row_indices = np.arange(R_norm.shape[0]).reshape(-1, 1)
    R_norm_with_indices = np.hstack((row_indices, R_norm))

    R_subset_norm = R_norm_with_indices[list(map(int,subset[:,0])),:]
    
    
    max_valuess = np.max(R_subset_norm[:,1:], axis=1)
    
    # Find the row index with the minimum of the maximum values
    min_index = np.argmin(max_valuess)
    idx = int(R_subset_norm[min_index,0])
    return max_valuess[0],idx
    


def find_CoverageGap_new(Store_representation, NDP_max, Delta,normalize=True,check_value=True):
    '''
    Input: a dictionary with indices delta and algorithm
    Output: coverage gap for each . Store it in a dictionary
    '''
    Alpha=0; 
    Z_set = set(map(tuple,NDP_max)) # convert a numpy array into a set of tuples
    R_set = set(map(tuple,Store_representation)) # convert a numpy array into a set of tuples
    
    Z_int_R = Z_set.intersection(R_set)  
    Z_new_set = Z_set.difference(Z_int_R)
    Z_new_array = np.array(list(map(list,Z_new_set))) # remove elements which are in set R
    
    if len(Z_new_array) ==0:
        Alpha=0; return Alpha
    R= Store_representation
    
    max_values = np.max(NDP_max,axis=0)
    min_values = np.min(NDP_max,axis=0)
    Z_new_array_norm = (Z_new_array-min_values)/(max_values-min_values) if normalize== True else Z_new_array
    R_norm = (R-min_values)/(max_values-min_values) if normalize == True else R
    mapping = {i:0 for i in range(np.shape(Z_new_array_norm)[0])} # maps (index) NDP from Z_new_array to  (index) NDP from R
    coverage = {i:0 for i in range(np.shape(Z_new_array_norm)[0])}
    for i1 in range(np.shape(Z_new_array_norm)[0]):
        alpha_temp_i = np.zeros(np.shape(R))
        alpha_temp_i_three = np.zeros(np.shape(R))
        for j1 in range(3):
            alpha_temp_i[:,j1] =   R_norm[:,j1]-Z_new_array_norm[i1,j1]
            alpha_temp_i_three[:,j1] = R[:,j1]-Z_new_array[i1,j1]

        # find the index of R where the closest representation for i \in Z exists
        #value,arg_min =closest(alpha_temp_i)
        value,arg_min = closest_new(max_values,min_values,alpha_temp_i_three, Delta,R_norm,check_value)
        if value==-1:
            return -1
        coverage[i1] = value
        mapping[i1] = int(arg_min) # index of R
    Alpha = max(coverage.values())
    return Alpha

folder_path = "C:\\Users\\sunney\\Desktop\\PaperV\\PQSM_alpha\\Results"
NDPs_PQSM = data_convert(dir_path=folder_path, types='PQSM')
NDPs_TDA = data_convert(dir_path=folder_path, types='TDA')
NDPs_GPBA = data_convert(dir_path=folder_path, types='GPBA')

alpha_pqsm={}
alpha_tda={}
alpha_gpba={}

for instance in range(8,9):
    Delta = [0.1,2,50]; 
    rounded_test_pqsm = np.round(NDPs_PQSM['{}'.format(instance),'{}'.format(Delta[0]),'{}'.format(Delta[1]),'{}'.format(Delta[2])],2)
    rounded_test_gpba = np.round(NDPs_GPBA['{}'.format(instance),'{}'.format(Delta[0]),'{}'.format(Delta[1]),'{}'.format(Delta[2])],2)
    rounded_test_tda = np.round(NDPs_TDA['{}'.format(instance),'{}'.format(Delta[0]),'{}'.format(Delta[1]),'{}'.format(Delta[2])],2) if ('{}'.format(instance),'{}'.format(Delta[0]),'{}'.format(Delta[1]),'{}'.format(Delta[2])) in NDPs_TDA.keys() else []
    
    joined = np.append(rounded_test_pqsm,rounded_test_gpba,axis=0) if len(rounded_test_gpba)!=0 else rounded_test_pqsm
    joined = np.append(joined,rounded_test_tda,axis=0) if len(rounded_test_tda)!=0 else joined
    
    joined_unique = np.unique(joined,axis=0)
    pareto_front = identify_pareto_front(joined_unique)
    
    alpha_pqsm[instance] = find_CoverageGap_new(rounded_test_pqsm, pareto_front,Delta,check_value=True)
    alpha_tda[instance] = find_CoverageGap_new(rounded_test_tda, pareto_front,Delta,check_value=True)   if len(rounded_test_tda)!=0 else -1
    alpha_gpba[instance] = find_CoverageGap_new(rounded_test_gpba, pareto_front,Delta,check_value=True) if len(rounded_test_gpba)!=0 else -1

print('PQSM={}'.format(alpha_pqsm))
print('alpha_tda={}'.format(alpha_tda))
print('alpha_gpba={}'.format(alpha_gpba))
#%%

import numpy as np

# Define the array
arr = np.array([[ 3.20e-01, 0.00e+00, -4.20e+01],
                [-2.80e-01, 1.00e+00, -4.20e+01],
                [-6.90e-01, 2.00e+00, -4.20e+01],
                [-8.90e-01, 3.00e+00, -4.20e+01],
                [ 2.10e-01, 0.00e+00, -3.20e+01],
                [-7.90e-01, 2.00e+00, 3.63e+02]])

# Define the condition
cond = np.all(arr <= [0.1, 1, 100], axis=1)

# Subset the array based on the condition
subset = arr[cond]

# Print the subset
print(subset)

#%%
import numpy as np

# Define the array
arr = np.array([[ 0.26446281, 0.0, -0.1037037 ],
                [-0.23140496, 0.33333333, -0.1037037 ],
                [-0.57024793, 0.66666667, -0.1037037 ],
                [-0.73553719, 1.0, -0.1037037 ],
                [ 0.17355372, 0.0, -0.07901235],
                [-0.65289256, 0.66666667, 0.8962963 ]])

# Find the maximum value in each row
max_values = np.max(arr, axis=1)

# Find the row index with the minimum of the maximum values
min_index = np.argmin(max_values)

# Print the result
print("Row index with minimum of the maximum value:", min_index)

