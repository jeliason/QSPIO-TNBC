import pickle as pkl
import numpy as np
import os

variables = ['M1','M2','Teff','T1exh','T0','Th','Thexh','C_total','VDT']

# Load the results from the PSA simulations
# get all pkl files
files = os.listdir('/nfs/turbo/umms-ukarvind/joelne/PSA_simulations/')
files = [f for f in files if f.endswith('.pkl')]

# get length of M1 from first file
with open('/nfs/turbo/umms-ukarvind/joelne/PSA_simulations/' + files[0],'rb') as f:
	data = pkl.load(f)
	M1 = data['M1']
	f.close()

length_var = M1.shape[0]

# initialize numpy array to store results
results = np.zeros((len(files),len(variables),length_var))

# load results into numpy array
for i,f in enumerate(files):
	with open('/nfs/turbo/umms-ukarvind/joelne/PSA_simulations/' + f,'rb') as f:
		data = pkl.load(f)
		for j,var in enumerate(variables):
			results[i,j,:] = data[var]
		f.close()

# save results
np.save('/nfs/turbo/umms-ukarvind/joelne/PSA_simulations.npy',results)

# get simulation number from each
simulations = np.array([f.split('_')[-1].split('.')[0] for f in files])

# save simulation numbers
np.save('/nfs/turbo/umms-ukarvind/joelne/PSA_simulation_numbers.npy',simulations)
