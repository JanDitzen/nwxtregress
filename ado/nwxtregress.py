from scipy.sparse import csc_matrix,lil_matrix, identity, linalg as sla
from scipy.linalg import lu, logm

from sfi import Mata
import numpy as np


def nwpython_lud(data_n,n_n,rgrid_n,output_n):
	
	dta = np.array(Mata.get(data_n))	

	rows = dta[:,0]
	cols = dta[:,1]
	dta = dta[:,2]
	rgrid = np.array(Mata.get(rgrid_n))	

	n = np.array(Mata.get(n_n),dtype=np.intc)
		
	rows = rows.flatten()
	cols = cols.flatten()
	dta = dta.flatten()
	n = n.flatten()	
	rgrid = rgrid.flatten()

	results = np.zeros((len(rgrid),1))
	
	for i in range(len(rgrid)):
		ri = rgrid[i]

		dta_x = - np.multiply(dta,ri)

		data_sp = csc_matrix((dta_x,(rows,cols)),shape = (n))	

		data_sp.setdiag(1,k=0)		
		#lud = sla.spilu(data_sp)
		lud = sla.splu(data_sp)
		Ai = lud.L+lud.U-identity(n[0])	
		detmi = np.array(Ai.diagonal())
		detmi[detmi<=0] = 1
		detmil = np.sum(np.log(detmi.data))
		results[i,:] = detmil	

	Mata.store(output_n,results)	


def nwpython_bp(data_n,vals_n,max_n,sdm_n,output_n):

	dta = np.array(Mata.get(data_n))
	vals = np.array(Mata.get(vals_n))
	max = np.array(Mata.get(max_n),dtype=np.intc)	
	sdm = np.array(Mata.get(sdm_n),dtype=np.intc)

	rows = dta[:,0]
	cols = dta[:,1]
	dta = dta[:,2]
	
	rows = rows.flatten()
	cols = cols.flatten()
	dta = dta.flatten()
	max = max.flatten()
	
	len = vals.shape
	gridl = len[1]
	nobs = len[0]

	sparse = csc_matrix((dta,(rows,cols)),shape = (nobs,nobs))

	wjjju = vals

	tracew = np.zeros((max[0],1))

	for i in range(0,max[0]):
		wjjju = sparse.dot(wjjju)
		tracew[i] = np.mean(np.multiply(np.array(vals),np.array(wjjju)))

	
	if (sdm == 0):
		tracew[1] = np.sum(((sparse.T).multiply(sparse)))/nobs

	Mata.store(output_n,tracew)