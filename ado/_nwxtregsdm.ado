capture program drop _nwxtregsdm
program define _nwxtregsdm, eclass

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
		/// output
		asarray(string)					/// name of array
		* 								/// ivlag weights
		]

		if "`usebp'" == "" {
			local uselud "uselud"
		}
		else {
			local uselud ""
		}

		GetArrayName `asarray'

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

			*** process spatial weigth matrices
			tempname nwxtreg_Warray 
			mata `nwxtreg_Warray' = asarray_create()

			*** Depvar
			gettoken 1 2: dvarlag , parse(",") 
			spmat_set `1' , `=subinstr("`2'",",","",.)' varlist(`spatial') arrayname("`nwxtreg_Warray'") isdep

			
			*** Indepvar
			local rest "`options'"
			while  "`rest'" != "" {
				tempname wname

				gettoken cur rest: rest, bind

				local 0 ", `cur'"
				syntax [anything], ivarlag(string) 

				gettoken 1 2: ivarlag , parse(",") 
				local 2 = subinstr("`2'",",","",.)
				gettoken tmp1 tmp2 : 1, parse(":")
				local tmp2 = subinstr("`tmp2'",":","",.)
				spmat_set `tmp1' , `2' varlist(`tmp2') arrayname("`nwxtreg_Warray'") wname(`wname')
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
timer on 1
			mata `asarray' = estimatesdm("`indepdepvars'","`touse'","`idvar' `tvar'",`nwxtreg_Warray',`nosparse',`gridlength',`draws',`nomit',(`bp1',`bp2'),`DataTrans')
timer off 1
			
		}		
		
		output_table_main `asarray', varnames("`indepdepvars'") touse(`touse') idvar(`idvar') tvar(`tvar_o') cmdname("Spatial SDM") warray("`nwxtreg_Warray'") sdm

		ereturn local estat_cmd "nwxtregress_estat"
		ereturn local predict "nwxtregress_predict"

		ereturn hidden local output_mat "`asarray'"
		ereturn hidden scalar nomit = `nomit'
		tempname sparse
		mata st_numscalar("`sparse'",asarray(asarray(`nwxtreg_Warray',"Wy"),"sparse"))
		ereturn hidden scalar issparse = `sparse'
		
		ereturn local cmd "nwxtregress"
		ereturn hidden local type "sdm"

		mata mata drop `nwxtreg_Warray'

end

// -------------------------------------------------------------------------------------------------
// Estimate SAR Mata program
// -------------------------------------------------------------------------------------------------
capture mata mata drop estimatesdm()
mata:
	function estimatesdm(	///
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
							real matrix BarryPace, 			///
							real matrix DataTrans 			///
							)

	{
timer_on(10)
"start SDM"
		/// load data
		LoadData(Varnames,idtname,tousename,Y=.,X=.,idt=.,N=.,T=.,K=.,PanelSetup=.,Unbalanced=.,uniqueid=.,uniquet=.)
		
		NT = rows(X)
		
		_DataTransform(Y,X,idt,uniquet,DataTrans,T,Varnames,K)

		///X = X1
		Wx_order = (J(cols(tokens(Varnames))-1,1,tokens(Varnames)[1]) , (tokens(Varnames)[2..cols(tokens(Varnames))])')

		/// loop over all weights
		nW = rows(asarray_keys(nwxtreg_Warray)) 
		
		while (nW > 0) {
			Wii = asarray(nwxtreg_Warray,asarray_keys(nwxtreg_Warray)[nW])
			"Weight"
			asarray_keys(nwxtreg_Warray)[nW]
			W_id = asarray(Wii,"W_id")

			if (W_make_sparse == 1 & asarray(Wii,"sparse") == 0 ) {
				"def. sparse"
				Wi = SparseDefine(asarray(Wii,"W"),W_id,0)
				
			}
			else {
				Wi = asarray(Wii,"W")
			}
			"start order"
			if (asarray_keys(nwxtreg_Warray)[nW]:=="Wy") {
				/// default is for Wy
				"Dep var"
				sum(Wi)
				W_sparse = asarray(Wii,"sparse")

				W = OrderW(Wi,W_id,idt,Unbalanced,N,T,uniqueid,asarray(Wii,"norm"),W_sparse)
				"sum is"
				sum(asarray(W,asarray_keys(W)[1]))
				Wy = spatialLag(Y,W,idt,T,W_id,uniqueid,uniquet,W_sparse)

				Wy_order = (asarray(Wii,"Wyname"),asarray(Wii,"vars"))
			}
			else {
				"Indep var"
				sum(Wi)
				Wi = OrderW(Wi,W_id,idt,Unbalanced,N,T,uniqueid,asarray(Wii,"norm"),asarray(Wii,"sparse"))
				"sum is"
				asarray_keys(Wi)
				sum(asarray(Wi,asarray_keys(Wi)[1]))
				"order done, spatial lag"
				XW = spatialLag(st_data(.,asarray(Wii,"vars"),tousename),Wi,idt,T,W_id,uniqueid,uniquet,asarray(Wii,"sparse"))
				rows(X),rows(XW)
				X = X, XW
				"X done"
				Wx_order
				asarray(Wii,"vars")
				J(cols(tokens(asarray(Wii,"vars"))),1,asarray_keys(nwxtreg_Warray)[nW])
				///Wx_order = Wx_order \ (J(cols(tokens(asarray(Wii,"vars"))),1,asarray_keys(nwxtreg_Warray)[nW]),tokens(asarray(Wii,"vars"))')
				Wx_order = Wx_order \ (J(cols(tokens(asarray(Wii,"vars"))),1,asarray(Wii,"Wyname")),tokens(asarray(Wii,"vars"))')
			}
			
			"order done"
			nW--

		}
		/// order spatial weight matrix ( only row normalise now!)
		if (W_sparse == 2) {
			W_sparse = 1
		}
		/// update K
		Ksp = cols(X)

		XX = quadcross(X,X)
		XY = quadcross(X,Y)
		XWy = quadcross(X,Wy)
		
		/// inital estimations 
		olsp(Y,X,b0=.,e0=.,XY,XX)
		olsp(Wy,X,b0=.,ed=.,XWy,XX)

		/// inital grid
		initgrid(idt,W,detval=.,gridlength/10,uniquet,N,T,W_sparse,BarryPace,trace=.,1)


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
		TI = I(Ksp) * 1e+12
		TI = invsym(TI) * I(Ksp)  
		TIc = TI * J(Ksp,1,0)
		nu = 0
		d0 = 0
timer_off(12)
		/// do GriddyGibbs
gridr = asarray_create("real",1)
asarray(gridr,1,(J(ndraws+nomit,3,.)))

		GriddyGibbs(ndraws,&Y,&X,&XX,&XY,&XWy,&Wy,&detval,rho0,Ksp,NT,sige,TI,TIc,nu,d0,nomit,sdraws=.,bdraws=.,pdraws=.,e0=.,ed=.,ng=.,gridr)



		/// Coefficients, Var/Cov, T-Stats
		smean = quadmeanvariance(sdraws)
		svarcov = smean[2..rows(smean),.]
		shat = smean[1,.] 

		pmean = quadmeanvariance(pdraws)		
		pvarcov = pmean[2..rows(pmean),.]
		phat = pmean[1,.] 

		bhat = quadmeanvariance(bdraws)
		varcov = bhat[2..Ksp+1,.]
		bhat = bhat[1,.]
		sds = sqrt(diagonal(varcov))
		tstat = bhat :/ sds'

		/// Statistics
		residuals = e0 - phat * ed

		///if (DataTrans[2]:==0) {
		///	K = K + 1
		///	Ksp = Ksp + 1
		///}


		sigma = quadcross(residuals,residuals) / (N-K)
		R2 = R2_calc(residuals,Y,N,Ksp)

		/// Upper and lower vales for CI Interval
		bupper = mm_quantile(bdraws,1,st_numscalar("c(level)"))
		blower = mm_quantile(bdraws,1,1-st_numscalar("c(level)"))

		supper = mm_quantile(sdraws,1,st_numscalar("c(level)"))
		slower = mm_quantile(sdraws,1,1-st_numscalar("c(level)"))

		pupper = mm_quantile(pdraws,1,st_numscalar("c(level)"))
		plower = mm_quantile(pdraws,1,1-st_numscalar("c(level)"))

		/// constant
		///if (DataTrans[2]:==0) bhat = _EstConst(Y,(Wy,X),(phat,bhat),cols(Wy),sigma,Wx_order,Varnames,varcov,tstat)		
		
		
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
		///if (sum(BarryPace):==0) {
		///	trace = 0			
			asarray(output,"idt",idt)
		///}
		///else {
			///asarray(output,"trace",trace)
		///}
		asarray(output,"trace",trace)
		asarray(output,"W",W)
		asarray(output,"sdm",1)
		asarray(output,"DataTrans",DataTrans)

		/// to remove
		asarray(output,"initgrid",detinit)
		asarray(output,"detval",detval)

		/// Stats
		Tminmax = PanelSetup[.,2] - PanelSetup[.,1] :+ 1
		Tminmax = min(Tminmax),max(Tminmax),mean(Tminmax)

		asarray(output,"r2",R2)
		asarray(output,"dimensions",(rows(Y),N,K,Ksp,T,Tminmax))
		asarray(output,"MCMC",(ndraws,nomit))

		asarray(output,"residuals",(residuals,idt))

		asarray(output,"Wstuff",nwxtreg_Warray)

		/// for testing
		///asarray(output,"Y",(Y,idt))
		///asarray(output,"X",(X,idt))
		///asarray(output,"Wy",(Wy,idt))
		///asarray(output,"gridr",gridr)
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
