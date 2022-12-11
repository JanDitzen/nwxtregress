capture program drop nwxtregress_estat
program define nwxtregress_estat, rclass

	gettoken subcmd rest: 0 
	if regexm("`subcmd'","impact*")==1 {
		`0'
	}
	else {
		noi disp "no cmd"
	}
	return add
end

program define impact, rclass
	syntax [anything] , [trace seed(string) array(string) python]
		if "`anything'" == "" {
			local varlist "`e(indepvar)'"
		}
		else {
			local varlist "`anything'"
		}

		if "`trace'" != "" {
			local trace noi
		}
		else {
			local trace "qui"
		}

		if "`array'" == "" local array `e(output_mat)'	

		cap mata asarray(`array',"W")

		if _rc != 0 {
			noi disp in smcl as error "mata array `array' with saved contents from {cmd:nwxtregress} not found."
			exit
		}

		if "`seed'" != "" {
			set seed `seed'
		}
		
		
		mata asarray(`array',"usepython",("`python'"!=""))
		`trace' mata nwxtreg_calc_DIEffects(`array')
		

		/// Output
		local maxline = c(linesize)	

		**allow max linesize of 100
		if `maxline' > 100 {
			local maxline = 100
		}
		local abname = 14
		local col_i = `abname' + 1
		
		local level = c(level)
		local maxline = `abname' + 66
		
		/// display
		di as text "Average Impacts" _col(`=`maxline'-30') "Number of obs =" %12.0f e(N)

		/// Header
		local maxline = `maxline' - 2
		di as text "{hline `col_i'}{c TT}{hline `=`maxline'-`col_i''}"
		di as text %`col_i's  abbrev(e(depvar),`abname') "{c |}" _c
		local col = `col_i' + 1 + 6
		di as text _col(`col') "dy/dx" _c
		local col = `col' + 5 + 3
		di as text _col(`col') "Std. Err."  _c
		local col = `col' + 9 + 6
		di as text _col(`col') "z"  _c
		local col = `col' + 1 + 4
		di as text _col(`col') "P>|z|"  _c
		local col = `col' + 5 + 5
		di as text _col(`col') "[`level'% Conf. Interval]"    
		di as text "{hline `col_i'}{c +}{hline `=`maxline'-`col_i''}"

		mata nwxtreg_output_matrix((asarray(`array',"direct")[.,1]),(asarray(`array',"direct")[.,2]),(asarray(`array',"direct")[.,3]),(asarray(`array',"direct")[.,4]),"`e(indepvar)'","direct",`col_i',0,"`varlist'",1::rows((asarray(`array',"direct"))))
		di as text "{hline `col_i'}{c +}{hline `=`maxline'-`col_i''}"

		mata nwxtreg_output_matrix((asarray(`array',"indirect")[.,1]),(asarray(`array',"indirect")[.,2]),(asarray(`array',"indirect")[.,3]),(asarray(`array',"indirect")[.,4]),"`e(indepvar)'","indirect",`col_i',0,"`varlist'",1::rows((asarray(`array',"indirect"))))
		di as text "{hline `col_i'}{c +}{hline `=`maxline'-`col_i''}"

		mata nwxtreg_output_matrix((asarray(`array',"total")[.,1]),(asarray(`array',"total")[.,2]),(asarray(`array',"total")[.,3]),(asarray(`array',"total")[.,4]),"`e(indepvar)'","total",`col_i',0,"`varlist'",1::rows((asarray(`array',"total"))))
		di as text "{hline `col_i'}{c BT}{hline `=`maxline'-`col_i''}"

		tempname b_direct b_indirect b_total V_direct V_indirect V_total mVnames

		mata `mVnames' = tokens("`varlist'")'
		mata `mVnames' = J(rows(`mVnames' ),1,""),`mVnames'

		mata nwxtreg_post_matrix(asarray(`array',"direct")[.,1]',"`b_direct'",`mVnames',0)
		mata nwxtreg_post_matrix(asarray(`array',"indirect")[.,1]',"`b_indirect'",`mVnames',0)
		mata nwxtreg_post_matrix(asarray(`array',"total")[.,1]',"`b_total'",`mVnames',0)

		mata nwxtreg_post_matrix(asarray(`array',"directVcov"),"`V_direct'",`mVnames',1)
		mata nwxtreg_post_matrix(asarray(`array',"indirectVcov"),"`V_indirect'",`mVnames',1)
		mata nwxtreg_post_matrix(asarray(`array',"totalVcov"),"`V_total'",`mVnames',1)

		mata asarray_remove(`array',"direct")
		mata asarray_remove(`array',"indirect")
		mata asarray_remove(`array',"total")

		mata asarray_remove(`array',"directVcov")
		mata asarray_remove(`array',"indirectVcov")
		mata asarray_remove(`array',"totalVcov")

		return clear
		

		return matrix b_direct = `b_direct'
		return matrix b_indirect = `b_indirect'
		return matrix b_total = `b_total'

		return matrix V_direct = `V_direct'
		return matrix V_indirect = `V_indirect'
		return matrix V_total = `V_total'

		*return add
end


// -------------------------------------------------------------------------------------------------
// _dotspct
// -------------------------------------------------------------------------------------------------
cap program drop _dotspct
program _dotspct
        version 8.2
        syntax [anything] [, title(string) reps(integer -1) NODOTS ///
                DOTS(numlist >0 integer max=1)]
        if "`dots'"=="" local dots 1
        tokenize `anything'
        args i rc
        if "`i'" == "0" {
                if `"`title'"' != "" {
                        if `reps' > 0 {
                                di as txt "`title' ("   ///
                                   as res `reps' as txt ") " 
                        }
                        else {
                                di as txt "`title'" 
                        }
                }
                if "`nodots'" == "" {
                        di as txt "{hline 4}{c +}{hline 3} 10 "  ///
                                  "{hline 2}{c +}{hline 3} 20 "  ///
                                  "{hline 2}{c +}{hline 3} 30 "  ///
                                  "{hline 2}{c +}{hline 3} 40 "  ///
                                  "{hline 2}{c +}{hline 3} 50   %"
                }
                exit
        }	
	else if `reps' < 100 {
		local start = floor((`i'-1)*100/`reps')	+1
		local ende = min(`start' + floor(100/`reps')-1,100) 
		local ende = ((floor((`i')*100/`reps')+1 - `ende') >1) + `ende'
		forvalues ii = `start'(1)`ende' {
			_dotspct `ii' `rc' , reps(100) 
		}			
	}
	else if `reps' > 100 {
		 local check = floor(`i'*100/`reps') - floor((`i'-1)*100/`reps')
		 if `check' > 0 _dotspct `=floor(`i'*100/`reps')' `rc' , reps(100) 
		  else {
                        if `i' == `reps' _dotspct 100 `rc' , reps(100) 
                 }
	}
	else {
		_dots 1 `rc'
		if `i' == 50 di as txt "" _col(51) %5.0f 50
		if `i' == 100 di as txt "" _col(51) %5.0f 100
	}		
end


cap program drop output_table
program define output_table
	syntax anything , [printnoz ]

	tokenize `anything'
	local varname `1'
	local coeff `2'
	local se `3'
	local tval `4' 
	local lower  `5'
	local upper `6'
	local col `7'

	if "`printnoz'" != "" {
		local tval ""
		local pval ""
	}
	else {
		local pval= 2*(1 - normal(abs(`tval' )))
	}
	di as text %`col's abbrev("`varname' ",`=`col'-1') "{c |}"  _continue
	local col = `col' + 3
	di as result _column(`col') %9.8g `coeff' _continue
	local col = `col' + 8 + 3
	di as result _column(`col') %9.8g `se' _continue
	local col = `col' + 8 + 3
	di as result _column(`col') %6.2f `tval' _continue
	
	local col = `col' + 10
	di as result _column(`col') %5.3f `pval' _continue
	local col = `col' + 10
	di as result _column(`col') %9.7g `lower' _continue
	local col = `col' + 11
	di as result _column(`col') %9.7g `upper'
end


// -------------------------------------------------------------------------------------------------
// Python programs
// -------------------------------------------------------------------------------------------------

if c(version) >= 16 {

python:
from scipy.sparse import csc_matrix, identity, linalg as sla
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
		lud = sla.spilu(data_sp)
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

end



}
