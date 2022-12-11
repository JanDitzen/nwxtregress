*! nwxtregress
*! v. 0.12
*! 11.12.2022 - see https://janditzen.github.io/nwxtregress/
/*
Change Log:
- 11.12.2022 - bug in BarryPace trick solved
- 01.12.2022 - fix if spatial weight matrix was empty and python used
*/

/*
This file contains all Stata programs.
mata programs are contained in lnwxtregress
*/

capture program drop nwxtregress
program define nwxtregress, eclass
syntax varlist(ts min=2) [if] 	, 	[	///
		/// spatial weight matrix settings
		dvarlag(string)					/// name of spatial weights matrix for dependent var
		seed(string)					/// set seed	
		* 								/// ivlag weights
		update							/// update from github
		version							/// shows version
		*								/// remaining options
		] 
		
		version 14.2

		if "`update'" != "" {
			qui nwxtregress, version
			local v_installed "`r(version)'"	
			cap net uninstall nwxtregress
			if _rc == 0 {
				noi disp "Version `v_installed' removed."
				local rm = 0
			}
			else {
				cap ssc uninstall nwxtregress
				if _rc == 0 {
					noi disp "Version `v_installed' (from SSC) removed."
				}
				else {
					noi disp "No version from net install or ssc install found!"
				}
			}
			noi disp "Updating from Github ... ", _c
			cap qui net install nwxtregress , from("https://janditzen.github.io/nwxtregress/") force replace
			if _rc == 0 {
				noi disp " Update successfull."
				qui nwxtregress, version
				noi disp "New version is `r(version)'"
			}
			else {
				noi disp "Update not successfull!"
			}
			exit
		}

		if "`version'" != "" {
			local version 0.1
			noi disp "This is version `version' - xx.xx.2022"
			nwxtregress , version
			return local version "`version'"
			exit
		}

		**** Main Program

		*noi disp as text ""
		*noi disp as error "This is an alpha version! Results may change.", _c
		noi disp as text ""

		*** necessary checks
		* moremata
		*** check that moremata is installed
		qui{
			cap findfile moremata.hlp
			if _rc != 0 {
				noi disp "Moremata is required. Please install:"
				noi disp in smcl "{stata ssc install moremata}"
				error 198
			}
		}


		* set seed
		if "`seed'" != "" {
			set seed `seed'
		}

		_nwxtreg `varlist' `if' , `options' dvarlag(`dvarlag')		
		

end


capture program drop _nwxtreg
program define _nwxtreg, eclass

	syntax varlist(ts min=2) [if] 	, 	///
		dvarlag(string)					/// weights, can be Stata or mata
		///ivarlag(string)					/// weights, can be Stata or mata
		[								/// optional options
		trace							/// trace for seeing all output
		/// spatial option
		NOSParse 						/// do not use sparse matrix internally
		/// sampling options			
		draws(integer 2000)				/// number of griddy gibs draws
		gridlength(integer 1000)		/// grind length
		NOMIT(integer 500)				/// number of omitted draws
		BArrypace(numlist)				/// settings for BarryPace Trick, iterations, maxorder default: 50 100
		usebp							/// use BP instead of LUD
		/// FE, constant etc.
		fe 								/// add fixed effects
		STANDardize						/// standardize data
		NOCONStant						/// suppress constant
		/// computational options
		python 							/// use pyhton for LUD
		/// output
		asarray(string)					/// name of array
		* 								/// ivlag weights
		]		
		
		*** Get Array Name
		GetArrayName `asarray'	

		*** Check if option python is used that python is installed
		if `c(stata_version)' >= 16 {
			if "`python'" != "" nwxtreg_check_python
		}
		else {
			if "`python'" != "" {
				noi disp as smcl "Option {cmd:python} requires Stata 16 or higher."
				local python = ""
			}
		}
		*** check that mata library is installed
		qui mata mata mlib index
		tempname mataversion
		cap mata st_numscalar("`mataversion'",nwxtregress_m_version())
		if _rc != 0 {
			noi disp "mata library for nwxtregress not found."
			exit
 		}
 		else {
			if `mataversion' < 1.01 {
				noi disp "mata library outdated. Pleae update nwxtregress:"
				noi disp as smcl "{stata:nwxtregress, update}"
				exit
			}
		}			

		if "`trace'" == "" {
			local trace "qui"
		}
		else {
			local trace noi
		}
		noi di ""
		`trace' {
			

			*** mark sample
			tempname touse
			marksample touse

			*** get id and t var
			_xt
			local tvar_o `r(tvar)'
			local idvar `r(ivar)'	

			*** make sure data is sorted
			issorted `idvar' `tvar_o'
			
			*** generate new tvar from 1...T
			tempvar tvar
			local tvar `tvar_o'
			
			tsunab indepdepvars :  `varlist'			
			local spatial =word("`indepdepvars'",1)		
			
			**** BP options 
			if "`usebp'" == "" {
				local uselud "uselud"
			}
			else {
				local uselud ""
			}

			*** adjust touse
			markout `touse' `indepdepvars' `spatial'

			*** process spatial weigth matrices
			tempname nwxtreg_Warray 
			mata `nwxtreg_Warray' = asarray_create()

			*** Depvar
			gettoken 1 2: dvarlag , parse(",") 
			local w_dep `1'
			spmat_set `1' , `=subinstr("`2'",",","",.)' varlist(`spatial') arrayname("`nwxtreg_Warray'") isdep
			
			local sar = 1
			local sdmsar = "SAR"
			
			*** Indepvar
			local rest "`options'"
			while  "`rest'" != "" {
				tempname wname

				if "`wname'" != "`w_dep'" {

					gettoken cur rest: rest, bind

					local 0 ", `cur'"
					syntax [anything], ivarlag(string) 

					gettoken 1 2: ivarlag , parse(",") 
					local 2 = subinstr("`2'",",","",.)
					gettoken tmp1 tmp2 : 1, parse(":")
					local tmp2 = subinstr("`tmp2'",":","",.)
					spmat_set `tmp1' , `2' varlist(`tmp2') arrayname("`nwxtreg_Warray'") wname(`wname')
				}
				else {
					mata asarray(`nwxtreg_Warray',"Wdep_indep",1)	
				}
				local sar = 0
				local sdmindic sdm
			}
			
			*** Use sparse matrices internally
			if "`nosparse'" == "" {
				local nosparse = 1
			}
			else {
				local nosparse = 0
			}


			*** barry pace options
			if "`uselud'" != "" {
				local bp1 = 0
				local bp2 = 0
			}
			else {
				if "`barrypace'" != "" {
					tokenize `barrypace'
					local bp1 = `1'
					local bp2 = `2'
				}
				else {
					local bp1 = 50
					local bp2 = 100
				}
			}
			
			*** clear ereturn
			tempname DataTrans
			mata `DataTrans' = ("`fe'"!=""),("`noconstant'"!=""),("`standardize'"!="")
			mata `asarray' = nwxtreg_estimate("`indepdepvars'","`touse'","`idvar' `tvar'",`sar',`nwxtreg_Warray',`nosparse',`gridlength',`draws',`nomit',(`bp1',`bp2'),`DataTrans',("`python'"!=""),("`trace'"!=""))
			
		}
		
		output_table_main `asarray', varnames("`indepdepvars'") touse(`touse') idvar(`idvar') tvar(`tvar_o') cmdname("Spatial `sdmsar'") warray("`nwxtreg_Warray'") `sdmindic'

		ereturn local estat_cmd "nwxtregress_estat"
		ereturn local predict "nwxtregress_predict"

		ereturn hidden local output_mat "`asarray'"
		ereturn hidden scalar nomit = `nomit'
		tempname sparse
		mata st_numscalar("`sparse'",asarray(asarray(`nwxtreg_Warray',"Wy"),"sparse"))
		ereturn hidden scalar issparse = `sparse'
		
		ereturn local cmd "nwxtregress"
		ereturn hidden local type "`sdmsar'"

		mata mata drop `nwxtreg_Warray'
		cap mata mata drop  nwxtreg_* __*
		
end

// issorted checks if data is sorted, if not then sorts it
capture program drop issorted
program define issorted
	syntax	varlist 
	
	local sorted : dis "`:sortedby'"
	if "`sorted'" != "`varlist'" {
	    noi disp "sort data"
	    sort `varlist'
	}

end

// -------------------------------------------------------------------------------------------------
// spatial weight settings
// -------------------------------------------------------------------------------------------------
capture program drop spmat_set
program define spmat_set
	syntax anything(name=wname_init) , [mata frame(string) id(string) sparse timesparse isdep NORMalize(string)] varlist(string) arrayname(string) [wname(string)] [zero(string)] [addzero]

	if `c(stata_version)' < 16 & "`frame'" != "" {
		noi disp as error in smcl "Option {it:frame} requires Stata 16 or higher."
		error 199
	}


	if "`mata'" != "" & "`frame'" != "" {
		noi disp as error in smcl "Options {it:mata} and {it:frame} cannot be combined."
		error 184
	}
	if "`wname'" == "" {
		tempname wname
	}
	tempname w_ary 
	mata: `w_ary' = asarray_create()

	if "`sparse'" != "" & "`timesparse'" != "" {
	*** time sparse implies sparse
			local sparse ""
	}

	/// default: non sparse remove zeros, sparse treat them as 0.0001
	if "`zero'" == "" & "`sparse'`timesparse'" != "" local zero 0.0001 
	else if "`zero'" == "" & "`sparse'`timesparse'" == "" local zero .

	if "`addzero'" != "" {
		local addzero = 1
	}
	else {
		local addzero = 0
	}

	
	if "`mata'" == "" & "`frame'" == "" {
		/// spmat
		tempname w_id 
		spmatrix matafromsp `wname' `w_id' = `wname_init'
		mata asarray(`w_ary',"W_id",`w_id')
		mata asarray(`w_ary',"W",`wname')
	}
	else if "`mata'" == "" & "`frame'" != "" {
		/// frame
		if "`sparse'`timesparse'" == "" {
			frame `frame': mata asarray(`w_ary',"W",st_data(.,"`wname_init'"))
			frame `frame': mata asarray(`w_ary',"W_id",1::rows(st_data(.,"`wname_init'")))
		}
		else {
			if "`id'" == "" {
				noi disp as error in smcl "Option {it:id()} required."
				error 198
			}
		}
		tempname data idt
		frame `frame': mata `data' = st_data(.,"`wname_init'")
		frame `frame': mata `idt' = st_data(.,"`id'")

		mata asarray(`w_ary',"W",(`idt',`data'))

		mata mata drop `data' `idt'
		
		tempname w_id
		frame `frame': mata `w_id' = st_data(.,"`id'")
		
		local wname_init = word("`wname_init'",1)
		local wname  "`wname_init'"

		if "`timesparse'" != "" mata asarray(`w_ary',"W_id",uniqrows(`w_id'[.,2]))
		if "`sparse'" != "" 	mata asarray(`w_ary',"W_id",uniqrows(`w_id'[.,1]))		
 		
	}
	else if "`mata'" != "" & "`frame'" == "" {
		mata asarray(`w_ary',"W",`wname_init')
		mata `wname' = `wname_init'	

		if "`id'" == "" {
			if "`sparse'`timesparse'" == "" {
				mata asarray(`w_ary',"W_id",(1::rows(`wname_init')))
			}
			else if "`sparse'" != "" & "`timesparse'" == "" {
				mata asarray(`w_ary',"W_id",uniqrows(`wname_init'[.,1]))
			}
			else if "`sparse'" == "" & "`timesparse'" != "" {
				mata asarray(`w_ary',"W_id",uniqrows(`wname_init'[.,2]))
			} 
		}
		else {
			mata asarray(`w_ary',"W_id",`id')
		}
	}

	*** Sparse indicator
	if "`sparse'`timesparse'" == "" mata asarray(`w_ary',"sparse",0)
 	else if "`sparse'" != "" 	mata asarray(`w_ary',"sparse",1)
	else if "`timesparse'" != ""  	mata asarray(`w_ary',"sparse",2)
	

	*** which normalisation
	if "`normalize'" == ""  		mata asarray(`w_ary',"norm",1)
	else if "`normalize'" == "none"		mata asarray(`w_ary',"norm",0)
	else if "`normalize'" == "row" 		mata asarray(`w_ary',"norm",1)
 	else if regexm("`normalize'","col*") 	mata asarray(`w_ary',"norm",2)
	else if regexm("`normalize'","spec*")  	mata asarray(`w_ary',"norm",3)
	else if regexm("`normalize'","minmax")  mata asarray(`w_ary',"norm",4)	
	
	mata asarray(`w_ary',"zero",`zero')
	mata asarray(`w_ary',"addzero",`addzero')

	mata asarray(`w_ary',"vars","`varlist'")
	mata asarray(`w_ary',"isdepvar",("`isdep'"!=""))

	if "`isdep'" != "" {
		mata asarray(`w_ary',"Wyname","`wname_init'")
		mata asarray(`arrayname',"Wy",`w_ary')
	}
	else {
		mata asarray(`w_ary',"Wyname","`wname_init'")
		mata asarray(`arrayname',"`wname'",`w_ary')
	}
	cap mata mata drop `w_ary' 
	disp "Spatial Weight Set"
	*cap mata mata drop `wname'
end

// -------------------------------------------------------------------------------------------------
// GetArrayName
// -------------------------------------------------------------------------------------------------
cap program drop GetArrayName
program define GetArrayName
syntax [anything]

	if "`anything'" == "" {
		local anything "_NWXTREG_OBJECT"
	}

 	cap mata asarray_contains(`anything', "W")

	if _rc != 0 {
		c_local asarray `anything'
	}
	else {
		local i = 0
		while _rc == 0 {
			local i = `i' + 1
			cap mata asarray_contains(`anything'`i', "W")			
		}
		c_local asarray `anything'`i'
	}
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

// -------------------------------------------------------------------------------------------------
// parse spatial lags
// -------------------------------------------------------------------------------------------------

*** load weight matrix
*** sparse indicates that spatial weight matrix is already sparse
*** time sparse indicates that spatial weitght matrix is sparse and different in time
*** if sparse matrix is used definition is: |t|,i, j,, value
*** id is a mata matrix indicating the id


cap program drop parse_spatial_lag
program define parse_spatial_lag
	syntax anything(name=wname) , [mata id(string) sparse timesparse] name(string)

	tempname w_id 

	if "`sparse'" != "" & "`timesparse'" != "" {
		*** time sparse implies sparse
		local sparse ""
	}

	if "`mata'" == "" {
		spmatrix matafromsp `wname' `w_id' = `wname'
	}

	if "`id'" == "" {
		if "`sparse'`timesparse'" == "" {
			/// here correct!
			mata `w_id' = (1::rows(`wname'))
		}
		else if "`sparse'" != "" & "`timesparse'" == "" {
			mata `w_id' = uniqrows(`wname'[.,1])
		}
		else if "`sparse'" == "" & "`timesparse'" != "" {
			mata `w_id' = uniqrows(`wname'[.,2])
		} 
	}
	else {	
		mata `w_id' = `id'
	}

	*** 1 if W is sparse
	if "`sparse'`timesparse'" == "" {
		local sparse = 0
	}
	else if "`sparse'" != "" {
		local sparse = 1
	}
	else if "`timesparse'" != ""  {
		local sparse = 2
	}

	*** Use sparse matrices internally
	if "`nosparse'" == "" {
		local nosparse = 1
	}
	else {
		local nosparse = 0
	}

	*** Return
	
	c_local nosparse `nosparse'

end

// -------------------------------------------------------------------------------------------------
// Output Table
// -------------------------------------------------------------------------------------------------
cap program drop output_table_main
program define output_table_main, eclass
	syntax anything(name=matname) , varnames(string) [col(real 0)] touse(string) idvar(string) tvar(string) [cmdname(string)] warray(string)  [sdm]

	/// varnames
	local y = word("`varnames'",1)
	local varnames : list varnames - y

	/// prepare data for output
	tempname b bp varcov pvarcov

	mata nwxtreg_post_matrix(asarray(`matname',"bhat"),"`b'", asarray(`matname',"Wx_order"),0)
	mata nwxtreg_post_matrix(asarray(`matname',"phat"),"`bp'",asarray(`matname',"Wy_order"),0)
	mata nwxtreg_post_matrix(asarray(`matname',"bvarcov"),"`varcov'",asarray(`matname',"Wx_order"),1)
	mata nwxtreg_post_matrix(asarray(`matname',"pvarcov"),"`pvarcov'",asarray(`matname',"Wy_order"),1)	

	mata st_local("K",strofreal(asarray(`matname',"dimensions")[3]))
	
	/// Add spatial coeff 
	matrix `b' = `b', `bp'
	matrix `varcov' = (`varcov' , J(rowsof(`varcov'),colsof(`pvarcov'),0)) \ (J(rowsof(`pvarcov'),colsof(`varcov'),0), `pvarcov' )

	local tmp: colnames `b'
	matrix colnames `varcov' = `tmp'
	matrix rownames `varcov' = `tmp'
	local tmp: coleq `b'
	matrix coleq `varcov' = `tmp'
	matrix roweq `varcov' = `tmp'

	/// Post
	ereturn clear 
	ereturn post `b' `varcov' , depname(`y') esample(`touse')

	mata st_numscalar("e(N)",(asarray(`matname',"dimensions")[1]))
	mata st_numscalar("e(N_g)",(asarray(`matname',"dimensions")[2]))
	mata st_numscalar("e(K)",(asarray(`matname',"dimensions")[3]))
	mata st_numscalar("e(Kfull)",(asarray(`matname',"dimensions")[4]))
	mata st_numscalar("e(T)",(asarray(`matname',"dimensions")[5]))
	mata st_numscalar("e(Tmin)",(asarray(`matname',"dimensions")[6]))
	mata st_numscalar("e(Tmax)",(asarray(`matname',"dimensions")[7]))
	mata st_numscalar("e(Tavg)",(asarray(`matname',"dimensions")[8]))

	mata st_numscalar("e(r2)",(asarray(`matname',"r2")[1]))
	mata st_numscalar("e(r2_a)",(asarray(`matname',"r2")[2]))
	mata st_numscalar("e(MCdraws)",(asarray(`matname',"MCMC")[1]))

	mata st_numscalar("e(HasCons)",(asarray(`matname',"DataTrans")[2]:==0),"hidden")

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

	di in gr _col(`=`maxline'-80+50') in gr "Number of obs" _col(`=`maxline'-80+68') "=" ///
		_col(`=`maxline'-80+71') in ye %9.0f e(N)

	di in gr "Panel Variable (i): " in ye abbrev(e(idvar),`abname') in gr ///
		_col(`=`maxline'-80+50') "Number of groups" _col(`=`maxline'-80+68') "=" ///
		_col(`=`maxline'-80+71') in ye %9.0g e(N_g) 


	if e(T) == e(Tavg) {
		di in gr "Time Variable (t): " in ye abbrev(e(tvar),`abname') in gr ///
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
	tempname index
	mata `index' = selectindex(asarray(`matname',"Wx_order")[.,1]:=="`y'")
	mata nwxtreg_output_matrix(asarray(`matname',"bhat"),asarray(`matname',"bvarcov"),asarray(`matname',"blower"),asarray(`matname',"bupper"),asarray(`matname',"Wx_order")[.,2],"",`col_i',0,"",`index')
	di as text "{hline `col_i'}{c +}{hline `=`maxline'-`col_i''}"

	/// display W, start with dependent variable W
	mata nwxtreg_output_matrix(asarray(`matname',"phat"),asarray(`matname',"pvarcov"),asarray(`matname',"plower"),asarray(`matname',"pupper"),asarray(`matname',"Wy_order")[.,2],asarray(`matname',"Wy_order")[.,1],`col_i',0,"",1)

	/// now display other variables with same spatial weight matrix (that is rho W)
	mata `index' = selectindex(asarray(`matname',"Wx_order")[.,1]:==asarray(`matname',"Wy_order")[1,1])
	mata nwxtreg_output_matrix(asarray(`matname',"bhat"),asarray(`matname',"bvarcov"),asarray(`matname',"blower"),asarray(`matname',"bupper"),asarray(`matname',"Wx_order")[.,2],"",`col_i',0,"",`index')
	di as text "{hline `col_i'}{c +}{hline `=`maxline'-`col_i''}"

	if "`sdm'"!= "" {
	
		/// now loop over remaining
		mata st_local("ListW",invtokens(uniqrows(asarray(`matname',"Wx_order")[.,1])'))
		mata st_local("Wy",asarray(`matname',"Wy_order")[1,1])
		local ListW : list ListW - Wy
		local ListW : list ListW - y

		foreach wi in `ListW' {
			mata `index' = selectindex(asarray(`matname',"Wx_order")[.,1]:=="`wi'")
			mata nwxtreg_output_matrix(asarray(`matname',"bhat"),asarray(`matname',"bvarcov"),asarray(`matname',"blower"),asarray(`matname',"bupper"),asarray(`matname',"Wx_order")[.,2],"`wi'",`col_i',0,"",`index')
			di as text "{hline `col_i'}{c +}{hline `=`maxline'-`col_i''}"
		}
	}
	
	/// display sigma
	mata nwxtreg_output_matrix(asarray(`matname',"shat"),asarray(`matname',"svarcov"),asarray(`matname',"slower"),asarray(`matname',"supper"),"/sigma_u","",`col_i',1,"",1)
	di as text "{hline `col_i'}{c BT}{hline `=`maxline'-`col_i''}"

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
// python check
// -------------------------------------------------------------------------------------------------
cap program drop nwxtreg_check_python
program define nwxtreg_check_python
	qui {
		python query
		if r(initialized) == 1 {
			cap python which numpy
			local HasNumpy = _rc	

			cap python which scipy
			local HasScipy = _rc

			cap python which sfi
			local HasSfi = _rc

			if `=`HasNumpy'+`HasSfi'+`HasScipy'' > 0 {
				noi disp as smcl "{cmd:nwxtregress} option {it:python} requires the following Python packages:"
				if `HasNumpy' != 0 noi disp "  numpy"
				if `HasScipy' != 0 noi disp "  scipy"
				if `HasSfi' != 0 noi disp "  sfi"
				noi disp "Please install them before using the option {it:python}."
			}
			exit
		}
		else {
			noi disp as error "Python not initialized."
			exit
		}
	}	


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
