import pickle as pkl
import numpy as np
import os

variables = ['M1','M2','Teff','T1exh','T0','Th','Thexh','C_total','VDT']

total_simulations = 600

# Load the results from the PSA simulations
# get all pkl files
files = os.listdir('/nfs/turbo/umms-ukarvind/joelne/PSA_simulations/')
files = [f for f in files if f.endswith('.pkl')]

# get length of M1 from first file
with open('/nfs/turbo/umms-ukarvind/joelne/PSA_simulations/' + files[0],'rb') as f:
	data = pkl.load(f)
	M1 = data['M1']
	# variables = data.keys()
	f.close()

length_var = M1.shape[0]

# initialize numpy array to store results
results = np.zeros((total_simulations,len(variables),length_var))

# load results into numpy array
for i,f in enumerate(files):
	sim_number = f.split('_')[-1].split('.')[0]
	with open('/nfs/turbo/umms-ukarvind/joelne/PSA_simulations/' + f,'rb') as f:
		data = pkl.load(f)
		for j,var in enumerate(variables):
			results[sim_number,j,:] = data[var]
		f.close()

# save results
np.save('/nfs/turbo/umms-ukarvind/joelne/PSA_simulations.npy',results)
