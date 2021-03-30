// -------------------------------------------------------------------------------------------------
// Output Table
// -------------------------------------------------------------------------------------------------
capture mata mata drop post_matrix
mata:
	function post_matrix(real matrix mat, string scalar as, string matrix varnames, real scalar rows)
	{
		varname = tokens(varnames)
		rcNames = J(cols(varname),1,""),varname'

		st_matrix(as,mat)
		st_matrixcolstripe(as,rcNames)

		 if (rows==1) {
		 	st_matrixrowstripe(as,rcNames)
		 }

	}
end

cap program drop output_table_main
program define output_table_main, eclass
	syntax anything(name=matname) , varnames(string) [col(real 0)] spatialnames(string) touse(string) idvar(string) tvar(string) [cmdname(string)] weightname(string) 

	/// varnames
	local y = word("`varnames'",1)
	local varnames : list varnames - y

	/// prepare data for output
	tempname b bp varcov pvarcov
	mata post_matrix(asarray(`matname',"bhat"),"`b'","`varnames'",0)
	mata post_matrix(asarray(`matname',"phat"),"`bp'","`spatialnames'",0)
	mata post_matrix(asarray(`matname',"bvarcov"),"`varcov'","`varnames'",1)
	mata post_matrix(asarray(`matname',"pvarcov"),"`pvarcov'","`spatialnames'",1)
	

	mata st_local("K",strofreal(asarray(`matname',"dimensions")[3]))
	mata st_local("K",strofreal(asarray(`matname',"dimensions")[3]))


	/// Add spatial coeff 
	matrix `b' = `b', `bp'
	matrix `varcov' = (`varcov' , J(rowsof(`varcov'),colsof(`pvarcov'),0)) \ (J(rowsof(`pvarcov'),colsof(`varcov'),0), `pvarcov' )

	matrix colnames `varcov' = `varnames' `spatialnames'
	matrix rownames `varcov' = `varnames' `spatialnames'

	/// Post
	ereturn clear 
	ereturn post `b' `varcov' , depname(`y') esample(`touse')

	mata st_numscalar("e(N)",(asarray(`matname',"dimensions")[1]))
	mata st_numscalar("e(N_g)",(asarray(`matname',"dimensions")[2]))
	mata st_numscalar("e(K)",(asarray(`matname',"dimensions")[3]))
	mata st_numscalar("e(T)",(asarray(`matname',"dimensions")[4]))
	mata st_numscalar("e(Tmin)",(asarray(`matname',"dimensions")[5]))
	mata st_numscalar("e(Tmax)",(asarray(`matname',"dimensions")[6]))
	mata st_numscalar("e(Tavg)",(asarray(`matname',"dimensions")[7]))

	mata st_numscalar("e(r2)",(asarray(`matname',"r2")[1]))
	mata st_numscalar("e(r2_a)",(asarray(`matname',"r2")[2]))
	mata st_numscalar("e(MCdraws)",(asarray(`matname',"MCMC")[1]))

	ereturn local idvar "`idvar'"
	ereturn local tvar "`tvar'"
	ereturn local indepvar "`varnames'"

	local maxline = c(linesize)	
	**allow max linesize of 100
	if `maxline' > 100 {
		local maxline = 100
	}
	local abname = 14
	local col_i = `abname' + 1
	
	local level = c(level)
	local maxline = `abname' + 66
	di ""

	di in gr "`cmdname'" ///
		_col(`=`maxline'-80+50') in gr "Number of obs" _col(`=`maxline'-80+68') "=" ///
		_col(`=`maxline'-80+71') in ye %9.0f e(N)

	di in gr "Panel Variable (i): " in ye abbrev(e(idvar),`abname') in gr ///
		_col(`=`maxline'-80+50') "Number of groups" _col(`=`maxline'-80+68') "=" ///
		_col(`=`maxline'-80+71') in ye %9.0g e(N_g) 


	if e(T) == e(Tavg) {
		di in gr "Time Variable (t): " in ye abbrev(e(tvar),`abname') in gr
				_col(`=`maxline'-80+50') "Obs. of group" _col(`=`maxline'-80+68') "=" ///
				_col(`=`maxline'-80+71') in ye %9.0g e(T) 
	}
	else {
		di in gr "Time Variable (t): " in ye abbrev(e(tvar),`abname') in gr ///
				_col(`=`maxline'-80+50') "Obs. of group:" _col(`=`maxline'-80+68') ///
				_col(`=`maxline'-80+71') in ye %9.0g e(T) 
		di in gr _col(`=`maxline'-80+68-4') in gr "min = " ///
				_col(`=`maxline'-80+71') in ye %9.0f e(Tmin) 	
		di in gr _col(`=`maxline'-80+68-4') in gr "avg = " ///
				_col(`=`maxline'-80+71') in ye %9.0f e(Tavg) 	
		di in gr _col(`=`maxline'-80+68-4') in gr "max = " ///
				_col(`=`maxline'-80+71') in ye %9.0f e(Tmax) 	
	}
				  
	di in gr _col(`=`maxline'-80+50') in gr "R-squared" _col(`=`maxline'-80+68') "=" ///
				_col(`=`maxline'-80+71') in ye %9.2f e(r2) 	
	di in gr _col(`=`maxline'-80+50') in gr "Adj. R-squared" _col(`=`maxline'-80+68') "=" ///
				_col(`=`maxline'-80+71') in ye %9.2f e(r2_a)	

	/// Header
	local maxline = `maxline' - 2
	di as text "{hline `col_i'}{c TT}{hline `=`maxline'-`col_i''}"
	di as text %`col_i's  abbrev(e(depvar),`abname') "{c |}" _c
	local col = `col_i' + 1 + 6
	di as text _col(`col') "Coef." _c
	local col = `col' + 5 + 3
	di as text _col(`col') "Std. Err."  _c
	local col = `col' + 9 + 6
	di as text _col(`col') "z"  _c
	local col = `col' + 1 + 4
	di as text _col(`col') "P>|z|"  _c
	local col = `col' + 5 + 5
	di as text _col(`col') "[`level'% Conf. Interval]"    
	di as text "{hline `col_i'}{c +}{hline `=`maxline'-`col_i''}"

	/// Loop over variables
	mata output_matrix(asarray(`matname',"bhat"),asarray(`matname',"bvarcov"),asarray(`matname',"blower"),asarray(`matname',"bupper"),"`varnames'","",`col_i',0,"")
	di as text "{hline `col_i'}{c +}{hline `=`maxline'-`col_i''}"

	/// display W
	mata output_matrix(asarray(`matname',"phat"),asarray(`matname',"pvarcov"),asarray(`matname',"plower"),asarray(`matname',"pupper"),"`varnames'","`weightname'",`col_i',0,"")
	di as text "{hline `col_i'}{c +}{hline `=`maxline'-`col_i''}"

	/// display sigma
	mata output_matrix(asarray(`matname',"shat"),asarray(`matname',"svarcov"),asarray(`matname',"slower"),asarray(`matname',"supper"),"\sigma_u","",`col_i',1,"")
	di as text "{hline `col_i'}{c BT}{hline `=`maxline'-`col_i''}"

end

cap mata mata drop output_matrix
mata:
	function output_matrix(real matrix b, real matrix var, real matrix lower, real matrix upper, string matrix vnames, string scalar header, real scalar col, real scalar printNoZ,string scalar restrict)
	{
		
		if (header!="") {
			printf("{txt}%-s", header)
			printf("{txt}{space %s}{c |} \n",strofreal(col-udstrlen(header)))  

		}
		if (printNoZ==1) {
			printZi = "printnoz"
		}
		else {
			printZi = ""
		}
		if (rows(var)==cols(var)) {
			sd = sqrt(diagonal(var))
		}
		else {
			sd = sqrt(var)
		}
		vnamesi = tokens(vnames)
		
		i = 1
		Ks = max((cols(b),rows(b)))
		while (i<=Ks) {
			
			if (restrict != "") {			
				isin = sum(vnamesi[i]:==tokens(restrict))
			}
			else {
				isin = 1
			}
			if (isin >= 1) {
				t = b[i]:/sd[i]
				cmd = sprintf("output_table %s %s %s %s %s %s %s , %s",vnamesi[i],strofreal(b[i]),strofreal(sd[i]),strofreal(t),strofreal(lower[i]),strofreal(upper[i]),strofreal(col),printZi)
			
				stata(cmd)
			}
			i++
		}
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
