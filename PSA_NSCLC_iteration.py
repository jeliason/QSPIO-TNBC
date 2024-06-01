print('Importing libraries...')
import matlab.engine
import pickle as pkl
import numpy as np
import argparse
import pandas as pd
import os


print('Parsing arguments...')
parser = argparse.ArgumentParser()
parser.add_argument('--task_index', type=int, default=0)
parser.add_argument('--n_param_sets', type=int, default=1)
parser.add_argument('--write_directory', type=str, default='.')
args = parser.parse_args()
print('Number of parameter sets: ',args.n_param_sets)
print('Write directory: ',args.write_directory)

# try:
# 	start_index = int(os.environ['SLURM_PROCID'])
# except KeyError:
# 	start_index = args.start_index
start_index = args.task_index * args.n_param_sets
print('Params row start index: ' + str(start_index))

params_df = pd.read_csv('params_PSA.csv')
param_names = params_df.columns.values

print('Starting MATLAB engine...')
eng = matlab.engine.start_matlab()
eng.addpath('scripts/')
eng.addpath('mini models/')
eng.addpath('model/')
eng.addpath('parameters/')
eng.addpath('PK_fitting/')
eng.addpath('postprocessing/')
eng.addpath('visualization/')

for i in range(args.n_param_sets):
	idx = i + start_index
	params = params_df.iloc[idx].values
	param_struct = dict(zip(param_names,params))

	print('Running simulation for set: '+ str(idx))
	success,index,M1,M2,Teff,T1exh,T0,Th,Thexh,C_total,VDT = eng.PSA_NSCLC_iteration(param_struct,nargout=11)

	print('Saving results for set: ' + str(idx))

	result = [M1,M2,Teff,T1exh,T0,Th,Thexh,C_total]
	result = [np.array(r) for r in result]
	result = [success,index] + result + [VDT]

	result_dict = dict(zip(['success','idx','M1','M2','Teff','T1exh','T0','Th','Thexh','C_total','VDT'],result))
	pkl.dump(result_dict,open(args.write_directory + '/PSA_NSCLC_iteration_' + str(idx) + '.pkl','wb'))

print('Finished.')
