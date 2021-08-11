capture program drop nwxtregress_estat
program define nwxtregress_estat, rclass

	gettoken subcmd rest: 0 
	if regexm("`subcmd'","impact*")==1 {
		impact `rest'
	}
	else {
		noi disp "no cmd"
	}
	return add
end

program define impact, rclass
	syntax [anything] , [trace seed(string) array(string)]

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
		
		

		`trace' mata calc_DIEffects(`array')
		

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

		mata output_matrix((asarray(`array',"direct")[.,1]),(asarray(`array',"direct")[.,2]),(asarray(`array',"direct")[.,3]),(asarray(`array',"direct")[.,4]),"`e(indepvar)'","direct",`col_i',0,"`varlist'",1::rows((asarray(`array',"direct"))))
		di as text "{hline `col_i'}{c +}{hline `=`maxline'-`col_i''}"

		mata output_matrix((asarray(`array',"indirect")[.,1]),(asarray(`array',"indirect")[.,2]),(asarray(`array',"indirect")[.,3]),(asarray(`array',"indirect")[.,4]),"`e(indepvar)'","indirect",`col_i',0,"`varlist'",1::rows((asarray(`array',"indirect"))))
		di as text "{hline `col_i'}{c +}{hline `=`maxline'-`col_i''}"

		mata output_matrix((asarray(`array',"total")[.,1]),(asarray(`array',"total")[.,2]),(asarray(`array',"total")[.,3]),(asarray(`array',"total")[.,4]),"`e(indepvar)'","total",`col_i',0,"`varlist'",1::rows((asarray(`array',"total"))))
		di as text "{hline `col_i'}{c BT}{hline `=`maxline'-`col_i''}"

		tempname b_direct b_indirect b_total V_direct V_indirect V_total mVnames

		mata `mVnames' = tokens("`varlist'")'
		mata `mVnames' = J(rows(`mVnames' ),1,""),`mVnames'

		mata post_matrix(asarray(`array',"direct")[.,1]',"`b_direct'",`mVnames',0)
		mata post_matrix(asarray(`array',"indirect")[.,1]',"`b_indirect'",`mVnames',0)
		mata post_matrix(asarray(`array',"total")[.,1]',"`b_total'",`mVnames',0)

		mata post_matrix(asarray(`array',"directVcov"),"`V_direct'",`mVnames',1)
		mata post_matrix(asarray(`array',"indirectVcov"),"`V_indirect'",`mVnames',1)
		mata post_matrix(asarray(`array',"totalVcov"),"`V_total'",`mVnames',1)

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
// Aux files
// -------------------------------------------------------------------------------------------------
** auxiliary file with auxiliary programs
findfile "_spatial_models_aux.ado"
include "`r(fn)'"
findfile "_nwxtreg_output.ado"
include "`r(fn)'"
findfile "_sparse_matrix_op.ado"
include "`r(fn)'"
findfile "_direct_indirect_BP.ado"
include "`r(fn)'"
