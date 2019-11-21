import numpy as np
def gen_sample(n):
	fvalues = [0.9,1.0,1.0,1.0,1.0,1.0]
	vvalues = [0,0,0,0,0,0]
	pvalues = [0.01,0.001,np.random.uniform(0, 1),0.1,np.random.normal(0, 1e-12),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1),np.random.normal(0, 1)]
	return n, fvalues, vvalues, pvalues
