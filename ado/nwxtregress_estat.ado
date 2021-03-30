capture program drop nwxtregress_estat
program define nwxtregress_estat, rclass

	gettoken subcmd rest: 0, parse(",")
	if "`subcmd'" == "impact" {
		noi disp "`rest'"
		Impact `rest'
	}
	else {
		noi disp "no cmd"
	}

end

program define Impact
	syntax [anything]

		if "`anything'" == "" {
			local varlist "`e(indepvar)'"
		}
		else {
			local varlist "`anything'"
		}

		local output_mat = e(output_mat)
		qui mata calc_DIEffects(`output_mat')

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

		mata output_matrix((asarray(`output_mat',"direct")[.,1]),(asarray(`output_mat',"direct")[.,2]),(asarray(`output_mat',"direct")[.,3]),(asarray(`output_mat',"direct")[.,4]),"`e(indepvar)'","direct",`col_i',0,"`varlist'")
		di as text "{hline `col_i'}{c +}{hline `=`maxline'-`col_i''}"

		mata output_matrix((asarray(`output_mat',"indirect")[.,1]),(asarray(`output_mat',"indirect")[.,2]),(asarray(`output_mat',"indirect")[.,3]),(asarray(`output_mat',"indirect")[.,4]),"`e(indepvar)'","indirect",`col_i',0,"`varlist'")
		di as text "{hline `col_i'}{c +}{hline `=`maxline'-`col_i''}"

		mata output_matrix((asarray(`output_mat',"total")[.,1]),(asarray(`output_mat',"total")[.,2]),(asarray(`output_mat',"total")[.,3]),(asarray(`output_mat',"total")[.,4]),"`e(indepvar)'","total",`col_i',0,"`varlist'")
		di as text "{hline `col_i'}{c BT}{hline `=`maxline'-`col_i''}"

end


// -------------------------------------------------------------------------------------------------
// Aux files
// -------------------------------------------------------------------------------------------------
** auxiliary file with auxiliary programs
findfile "_spatial_models_aux.do"
include "`r(fn)'"
findfile "_output.do"
include "`r(fn)'"
findfile "_sparse_matrix_op.do"
include "`r(fn)'"
findfile "_direct_indirect_BP.do"
include "`r(fn)'"
