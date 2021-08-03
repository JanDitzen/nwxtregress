capture program drop _nwxtregsar
program define _nwxtregsar, eclass

	syntax varlist(ts min=2) [if] 	, 	///
		dvarlag(string)					/// weights, can be Stata or mata
		[								/// optional options
		trace							/// trace for seeing all output
		/// spatial option
		NOSParse 						/// do not use sparse matrix internally
		/// sampling options			
		draws(integer 2000)				/// number of griddy gibs draws
		gridlength(integer 1000)		/// grind length
		NOMIT(integer 500)				/// number of omitted draws
		BArrypace(numlist)				/// settings for BarryPace Trick, iterations, maxorder default: 50 100
		usebp							/// use LUD instead of Barry-Pace
		/// output
		direct 							/// show total, indirect and direct effects
		]

		if "`usebp'" == "" {
			local uselud "uselud"
		}
		else {
			local uselud ""
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
			*egen `tvar' = group(`tvar_o') if `touse'
			local tvar `tvar_o'
			*** tsset varlist 
			tsunab indepdepvars :  `varlist'
			
			local spatial =word("`indepdepvars'",1)		
			
			*** adjust touse
			markout `touse' `indepdepvars' `spatial'

			*** load weight matrix
			*** sparse indicates that spatial weight matrix is already sparse
			*** time sparse indicates that spatial weitght matrix is sparse and different in time
			*** if sparse matrix is used definition is: |t|,i, j,, value
			*** id is a mata matrix indicating the id
			/*
			local 0 `dvarlag'
			syntax anything(name=wname) , [mata id(string) sparse timesparse]

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
			*/
			*** process spatial weigth matrices
			tempname nwxtreg_Warray 
			mata `nwxtreg_Warray' = asarray_create()

			gettoken 1 2: dvarlag , parse(",") 
			spmat_set `1' , `=subinstr("`2'",",","",.)' varlist(`spatial') arrayname("`nwxtreg_Warray'") isdep

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

timer on 1
			mata nwxtreg_output = estimatesar("`indepdepvars'","`touse'","`idvar' `tvar'",`nwxtreg_Warray',`nosparse',`gridlength',`draws',`nomit',(`bp1',`bp2'))
timer off 1
			
			
			

			*mata mata drop `output'


		}		
		
		output_table_main nwxtreg_output, varnames("`indepdepvars'") touse(`touse') idvar(`idvar') tvar(`tvar_o') cmdname("Spatial SAR") warray("`nwxtreg_Warray'")

		ereturn local estat_cmd "nwxtregress_estat"

		ereturn hidden local output_mat "nwxtreg_output"
		ereturn hidden scalar nomit = `nomit'
		*ereturn hidden scalar issparse = `sparse'
		ereturn local cmd "nwxtregress"
		ereturn hidden local type "sar"
		
		mata mata drop `nwxtreg_Warray'

		*ereturn local indepvar 
end

// -------------------------------------------------------------------------------------------------
// Estimate SAR Mata program
// -------------------------------------------------------------------------------------------------
capture mata mata drop estimatesar()
mata:
	function estimatesar(	///
							string scalar Varnames,			/// list of variable names, order y x
							string scalar tousename,		///
							string scalar idtname,			///
							/// W options
							transmorphic matrix nwxtreg_Warray,	///
							real scalar W_make_sparse,		/// use internally sparse matrices
							/// griddy gibbs options
							real scalar gridlength,			///
							real scalar ndraws,				///
							real scalar nomit,				///
							real matrix BarryPace 			///
							)

	{
timer_on(10)
Varnames
		/// load data
		LoadData(Varnames,idtname,tousename,Y=.,X=.,idt=.,N=.,T=.,K=.,PanelSetup=.,Unbalanced=.,uniqueid=.,uniquet=.)
		
		NT = rows(X)

		_studentizeT(Y,T,idt,uniquet)
		_studentizeT(X,T,idt,uniquet)

		Wii = asarray(nwxtreg_Warray,"Wy")
		W_id = asarray(Wii,"W_id")
		W_sparse = asarray(Wii,"sparse")

		Wy_order = (asarray(Wii,"Wyname"),asarray(Wii,"vars"))
		Wx_order = (J(cols(tokens(Varnames))-1,1,tokens(Varnames)[1]) , (tokens(Varnames)[2..cols(tokens(Varnames))])')

		if (W_make_sparse == 1 & asarray(Wii,"sparse") == 0 ) {
			"def. sparse"
			W = SparseDefine(asarray(Wii,"W"),W_id,W_sparse)
		}
		else {
			W = asarray(Wii,"W")
		}
		/// order spatial weight matrix ( only row normalise now!)
		/// here matrix is restricted to correct time periods
		"start order"
		W = OrderW(W,W_id,idt,Unbalanced,N,T,uniqueid,1,W_sparse)
		asarray_keys(W)
		///asarray(W,1)[1..100,.]
		if (W_sparse == 2) {
			W_sparse = 1
		}

		"order done"
timer_off(10)
timer_on(11)
		/// Y and X are in long format, i.e. X is N*T(n) x K matrix, Y is N*T(n) x 1 vector
		/// construct X'X, X'Y and Wy
		XX = quadcross(X,X)
		XY = quadcross(X,Y)

		Wy = spatialLag(Y,W,idt,T,W_id,uniqueid,uniquet,W_sparse)

		XWy = quadcross(X,Wy)

		/// inital estimations (needed?)
		olsp(Y,X,b0=.,e0=.,XY,XX)

		olsp(Wy,X,b0=.,ed=.,XWy,XX)
timer_off(11)
timer_on(12)
timer_on(13)
		/// inital grid
		initgrid(idt,W,detval=.,gridlength/10,uniquet,N,T,W_sparse,BarryPace,trace=.,0)
timer_off(13)
"grid done"
		///fgrid = (-(1-1/gridlength)*gridlength..(1-1/gridlength)*gridlength):/gridlength
		fgrid = range(-0.99,0.99,1/gridlength)
		detinit = detval

timer_on(14)
		lndet = spline3eval(spline3(detval[.,1],detval[.,2]),fgrid)
		///_editmissing(lndet,0)
timer_off(14)
		detval = fgrid, lndet

		/// griddy gibbs sampler
		iter = 1
		rho0 = 0.1

		sige = 1
		TI = I(K) * 1e+12
		TI = invsym(TI) * I(K)  
		TIc = TI * J(K,1,0)
		nu = 0
		d0 = 0
timer_off(12)
		/// do GriddyGibbs
gridr = asarray_create("real",1)
asarray(gridr,1,(J(ndraws+nomit,3,.)))

		GriddyGibbs(ndraws,&Y,&X,&XX,&XY,&XWy,&Wy,&detval,rho0,K,NT,sige,TI,TIc,nu,d0,nomit,sdraws=.,bdraws=.,pdraws=.,e0=.,ed=.,ng=.,gridr)



		/// Coefficients, Var/Cov, T-Stats
		smean = quadmeanvariance(sdraws)
		svarcov = smean[2..rows(smean),.]
		shat = smean[1,.] 
"here"
		pmean = quadmeanvariance(pdraws)		
		pvarcov = pmean[2..rows(pmean),.]
		phat = pmean[1,.] 
"here2"	
		bhat = quadmeanvariance(bdraws)
		varcov = bhat[2..K+1,.]
		bhat = bhat[1,.]
		sds = sqrt(diagonal(varcov))
		tstat = bhat :/ sds'
"her"

		/// Statistics
		residuals = e0 - phat * ed
"res"
		sigma = quadcross(residuals,residuals) / (N-K)
		R2 = R2_calc(residuals,Y,N,K)
"dinbe"
		/// Upper and lower vales for CI Interval
		bupper = mm_quantile(bdraws,1,st_numscalar("c(level)"))
		blower = mm_quantile(bdraws,1,1-st_numscalar("c(level)"))

		supper = mm_quantile(sdraws,1,st_numscalar("c(level)"))
		slower = mm_quantile(sdraws,1,1-st_numscalar("c(level)"))

		pupper = mm_quantile(pdraws,1,st_numscalar("c(level)"))
		plower = mm_quantile(pdraws,1,1-st_numscalar("c(level)"))

		/// Direct and Indirect Effects
		
		
"DI done"
		/// Output
		output = asarray_create()
		asarray(output,"bhat",bhat)
		asarray(output,"phat",phat)
		asarray(output,"shat",shat)

		asarray(output,"sds",sds)

		asarray(output,"blower",blower)
		asarray(output,"bupper",bupper)
		asarray(output,"slower",slower)
		asarray(output,"supper",supper)
		asarray(output,"plower",plower)
		asarray(output,"pupper",pupper)

		asarray(output,"bvarcov",varcov)
		asarray(output,"pvarcov",pvarcov)
		asarray(output,"svarcov",svarcov)
		asarray(output,"varnames",Varnames)
		asarray(output,"BarryPace",BarryPace)
		asarray(output,"Wx_order",Wx_order)
		asarray(output,"Wy_order",Wy_order)
		/// draws for (in)direct effects
		asarray(output,"pdraws",pdraws)
		asarray(output,"sdraws",sdraws)
		asarray(output,"bdraws",bdraws)
		if (BarryPace==0) {
			trace = 0			
			asarray(output,"idt",idt)
		}
		else {
			trace = 1
		}
		asarray(output,"trace",trace)
		asarray(output,"W",W)
		asarray(output,"sdm",0)

		/// to remove
		asarray(output,"initgrid",detinit)
		asarray(output,"detval",detval)
		/// Stats
		Tminmax = PanelSetup[.,2] - PanelSetup[.,1] :+ 1
		Tminmax = min(Tminmax),max(Tminmax),mean(Tminmax)

		asarray(output,"r2",R2)
		asarray(output,"dimensions",(rows(Y),N,K,K,T,Tminmax))
		asarray(output,"MCMC",(ndraws,nomit))

		/// for testing
		asarray(output,"Y",(Y,idt))
		asarray(output,"X",(X,idt))
		asarray(output,"Wy",(Wy,idt))
		asarray(output,"gridr",gridr)
		return(output)
		

	}

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
