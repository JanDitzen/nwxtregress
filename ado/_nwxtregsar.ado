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
		gridlength(integer 10000)		/// grind length
		NOMIT(integer 500)				/// number of omitted draws
		BArrypace(numlist)				/// settings for BarryPace Trick, iterations, maxorder default: 50 100
		uselud							/// use LUD instead of Barry-Pace
		/// output
		direct 							/// show total, indirect and direct effects
		]


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
			egen `tvar' = group(`tvar_o') if `touse'
			
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
			
			mata nwxtreg_output = estimatesar("`indepdepvars'","`touse'","`idvar' `tvar'",`wname',`w_id',`sparse',`nosparse',`gridlength',`draws',`nomit',(`bp1',`bp2'))
timer off 1
			
			
			

			*mata mata drop `output'


		}		
		
		output_table_main nwxtreg_output, varnames("`indepdepvars'") spatialnames("`spatial'") touse(`touse') idvar(`idvar') tvar(`tvar_o') cmdname("Spatial SAR") weightname("`wname'") 

		ereturn local estat_cmd "nwxtregress_estat"

		ereturn hidden local output_mat "nwxtreg_output"
		ereturn hidden scalar nomit = `nomit'
		ereturn hidden scalar issparse = `sparse'
		ereturn local cmd "nwxtregress"
		
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
							real matrix Weight,				///
							real matrix W_id,				/// uniqe ordering of W
							real scalar W_sparse,			///
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

		if (W_make_sparse == 1 & W_sparse == 0 ) {
			W = SparseDefine(Weight,W_id)
			W_sparse = 1
		}
		else {
			W = Weight
		}
		/// order spatial weight matrix ( only row normalise now!)
		W = OrderW(W,W_id,idt,Unbalanced,N,T,uniqueid,1,W_sparse)

		if (W_sparse == 2) {
			W_sparse = 1
		}

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
		initgrid(idt,W,detval=.,gridlength/10,uniquet,N,T,W_sparse,BarryPace,trace=.)
timer_off(13)
"grid done"
"detval grid"
"trace"

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
		asarray(output,"varnames",varnames)
		asarray(output,"BarryPace",BarryPace)

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


		/// to remove
		asarray(output,"initgrid",detinit)
		asarray(output,"detval",detval)
		/// Stats
		Tminmax = uniqrows(idt[.,2],1)
		Tminmax = min(Tminmax),max(Tminmax),mean(Tminmax)

		asarray(output,"r2",R2)
		asarray(output,"dimensions",(rows(Y),N,K,T,Tminmax))
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
// griddy gibbs
// -------------------------------------------------------------------------------------------------
capture mata mata drop GriddyGibbs()
mata:
	function GriddyGibbs(			///
							real scalar ndraws,					///
							/// pointers to data
							pointer(real matrix) Y,				/// yvec in matlab
							pointer(real matrix) X,				/// is xpx in matlab 
							pointer(real matrix) XX,			/// x in matlab
							pointer(real matrix) XY,			///
							pointer(real matrix) XWy,			///					
							pointer(real matrix) Wy,			///	
							pointer(real matrix) detval,		///
							/// scalars
							real scalar rho,					///
							real scalar K,						///
							real scalar NT,						///
							real scalar sige,					///
							real matrix TI,						///
							real matrix TIc,					///
							real scalar nu,						///
							real scalar d0,						///
							/// options
							real scalar nomit, 					///
							/// Return
							real matrix sdraws,					///
							real matrix bdraws,					///
							real matrix pdraws,					///
							real matrix e0,						///
							real matrix ed,						///
							real scalar ng, 						///
							transmorphic  matrix gridr 					///
							)

	{
		"start griddy"
		saves = J(ndraws+nomit,1,.)
		bsave = J(ndraws+nomit,K,.)
		psave = J(ndraws+nomit,1,.)

		cmd = sprintf("noi _dotspct 0 , title(Griddy Gibbs) reps(%s )",strofreal(ndraws))	
		stata(cmd)



		i = 1
		while (i<=ndraws+nomit) {
		
			timer_on(80)			
			GriddyGibbsi(i,Y,X,XX,XY,XWy,Wy,detval,&bsave,&saves,&psave,rho,K,NT,sige,TI,TIc,nu,d0,e0=.,ed=.,gridr)		

			cmd = sprintf("noi _dotspct %s 0 , reps(%s)",strofreal(i),strofreal(ndraws+nomit))
			stata(cmd)
			i++
			timer_off(80)
		}
		/// remove nomit obs
		sdraws = saves[nomit+1..ndraws+nomit,.]
		bdraws = bsave[nomit+1..ndraws+nomit,.]
		pdraws = psave[nomit+1..ndraws+nomit,.]
		
		ng = max((rows(sdraws),cols(sdraws)))

	}

end

capture mata mata drop GriddyGibbsi()
mata:
	function GriddyGibbsi(										///
							real scalar iter,					///
							pointer(real matrix) Y,				/// yvec in matlab
							pointer(real matrix) X,				/// x in matlab
							pointer(real matrix) XX,			/// is xpx in matlab
							pointer(real matrix) XY,			/// is xpy in matlab
							pointer(real matrix) XWy,			///	normal matrix, not sparse!		
							pointer(real matrix) Wy,			///	normal matrix, not sparse!
							pointer(real matrix) detval,		///
							pointer(real matrix) bsave,			///
							pointer(real matrix) ssave,			///	
							pointer(real matrix) psave,			///					
							real scalar rho,					///
							real scalar K,						///
							real scalar NT,						///
							real scalar sige,					///
							real matrix TI,						///
							real matrix TIc,					///
							real scalar nu,						///
							real scalar d0,						///	
							/// return
							real matrix e0o,					///					
							real matrix edo, 					///
							transmorphic  matrix gridr 	///
						)
	{
		AI = m_xtdcce_inverter((*XX) + sige * TI)*I(K)	
		
		ys = (*Y) - rho * (*Wy)
		
		b = quadcross((*X),ys) + sige*TIc

		b0 = m_xtdcce_inverter(( (*XX) + sige * TI) ) * b	

		sigeAI = sige * AI
		bhat =  drawnorm(1,0,sigeAI)' + b0
		xb = (*X) * bhat

		(*bsave)[iter,.] = bhat'

		/// update sige
		nu1 = NT + 2*nu
		e = (ys - xb)
		d1 = 2*d0 + quadcross(e,e)

		k = asarray(gridr,1)
		k[iter,.] = d1,quadcross(e,e),nu1
		asarray(gridr,1,k)
		
		chi = rchi2(1,1,nu1)
		sige = d1/chi

		(*ssave)[iter,1] = sige

		/// Griddy Gibbs
		b0 = m_xtdcce_inverter((*XX) + sige*TI) * ((*XY) + sige * TIc)		
		bd = m_xtdcce_inverter((*XX) + sige*TI) * ((*XWy) + sige * TIc)
		e0 = (*Y) - (*X) * b0
		ed = (*Wy) - (*X)*bd
		
		epe0 = quadcross(e0,e0)
		eped = quadcross(ed,ed)
		epe0d = quadcross(ed,e0)

		nmk = NT/2
		nrho = max((rows((*detval)[.,1]),cols((*detval)[.,1])))
		
		iota = J(nrho,1,1)
		/// For each potential value of rho, calculate the log likelihood.
		z = epe0*iota - 2*(*detval)[.,1]*epe0d + (*detval)[.,1]:*(*detval)[.,1]*eped
		z = -nmk*log(z)
		den = (*detval)[.,2] + z
		nn = max((rows(den),cols(den)))
		yy = (*detval)[.,1]
		/// % Calculation of log likelihood missing a scaling value, and need an adjustment.

		adj = max(den)
		den = den :- adj
		xx = exp(den)

		/// Use the trapezoid rule to integrate the area under the pdf.
		isum = sum(((yy[2..nn,1] + yy[1..nn-1,1])) :* (xx[2..nn,1] - xx[1..nn-1,1])/2)
		z= abs(xx/isum)

		/// store cdf
		den = mm_colrunsum(z,0,1)

		/// Using a uniform distribution, randomly select a value of rho using the cdf (by inverting it).
		rnd = runiform(1,1,0,1)*sum(z)

		/// rnd scalar, so use ":"
		ind = selectindex(den :<= rnd)

		idraw = max(ind)

		if (idraw > 0 & idraw < nrho) {
			rho = (*detval)[idraw,1]
		}
		(*psave)[iter,1] = rho
		
		e0o = e0
		edo = ed
	}
end
// -------------------------------------------------------------------------------------------------
// init grid
// -------------------------------------------------------------------------------------------------
capture mata mata drop initgrid()
mata:
	function initgrid(	real matrix idt ,			///
						transmorphic matrix W,		///
						real matrix detval,			///
						real matrix gridlength,		///
						real matrix uniquet, 		///
						real scalar N,				///
						real scalar T,				///
						real scalar issparse,		///
						real matrix BarryPace, 		///
						real matrix trace 			///
						)
	{
		///rgrid = (-(1-1/gridlength)*gridlength..(1-1/gridlength)*gridlength):/gridlength
		rgrid = range(-0.99,0.99,1/gridlength)
		ngrid = rows(rgrid)
		detval = J(ngrid,2,.)	
		
		if (BarryPace[1,1]:== 0) {
			i=1
			while (i<=ngrid) {			
				detm = 0
				t=1				
				timer_on(90)
				while(t<=T) {
					detm = detm + LUD(W,rgrid,N,uniquet[t,1],i,issparse,idt)
					t++
				}
				timer_off(90)
				detval[i,.] = rgrid[i] , detm
				i++
			}
		}	
		else {
			i = 1
			while (i<=ngrid) {
				detm = 0
				t=1
				timer_on(90)
				///while(t<=T) {
				///	detm = detm + BarryPace(W,idt,T,N,BarryPace[1],BarryPace[2],issparse)
				"BP res"
				///BarryPace(W,idt,T,N,BarryPace[1],BarryPace[2],issparse)
				trace =  BarryPace(W,idt,T,N,BarryPace[1],BarryPace[2],issparse)
				detm = sum(log(trace))
				///}
				timer_off(90)
				detval[i,.] = rgrid[i] , detm
				i++
			}
		}
		///dasdas[1,1]+1
		"init grid done"
	}
end

// -------------------------------------------------------------------------------------------------
// rho W
// -------------------------------------------------------------------------------------------------
capture mata mata drop rhoW()
mata:
	function rhoW(	transmorphic matrix W,	///
					real scalar rho,		///
					real scalar t,			///	
					real scalar issparse	///				
					)

	{
		if (eltype(W) == "real") {
			ret = SparseMultiply(W,rho,issparse) 
		}
		else {
			ret = SparseMultiply(asarray(W,t),rho,issparse)
		}
		return(ret)
	}
end

capture mata mata drop rhoWIminus()
mata:
	function rhoWIminus(	transmorphic matrix W,	///
							real scalar rho,		///
							real scalar N,			///		
							real scalar t,			///
							real scalar issparse	///			
							)

	{
		if (eltype(W) == "real") {
			ret =  SparseMultiply(W,rho,issparse)
		}
		else {
			wi = asarray(W,t)
			ret = SparseMultiply(wi,rho,issparse)			
		}


		if (issparse == 1) {
			/// build I - rhoW; rhoW has no elements on diagonal, so multiply 3rd col of ret with -1 and add 
			/// diagonal elements to ret.
			ret[.,3] = -ret[.,3]
			/// add I(N) to bottom
			/// get uniqelements
			el1 = uniqrows(ret[.,1])
			el2 = uniqrows(ret[.,2])
			el = uniqrows((el1\el2))
			ret = ret \ (el,el,J(rows(el),1,1))
		}
		else {
			ret = I(N) - ret 
		}
		return(ret)
	}
end



// -------------------------------------------------------------------------------------------------
// OLS program
// -------------------------------------------------------------------------------------------------
capture mata mata drop olsp()
mata:
	function olsp( ///
					real matrix Y,		///
					real matrix X,		///
					real matrix betaP,	///
					real matrix eps,|	///
					real matrix XY,		///
					real matrix XX		///
					)
	{
		if (args() == 4) {
			betaP = m_xtdcce_inverter(quadcross(X,X)) * quadcross(X,Y)
		}
		else {
			betaP = m_xtdcce_inverter(XX) * (XY)
		}
		eps = Y - X * betaP
	}
end


// -------------------------------------------------------------------------------------------------
// spatial lag
// -------------------------------------------------------------------------------------------------
/// produces T(n)*N x K matrix with spatial lags. 
capture mata mata drop spatialLag()
mata:
	function spatialLag(real matrix mat,				///
						transmorphic matrix W,			///
						real matrix idt,				///
						real matrix T,					///
						real matrix Wuniq,				///
						real matrix uniqueid,			///
						real matrix uniquet, 			///
						real scalar issparse			///
						)
	{
		eltype(W)
		if (eltype(W) == "real") {		
			///ret = W * mat
			///Wuniqi = uniqrows(Wuniq,1)
			ret = SparseMultiply(W,mat,issparse,uniqueid)
		}
		else {
			"unbal case"
			/// unbalanced case, W is N(t)xN(t) array; T arrays in total
			ret = J(rows(mat),cols(mat),.)
			t = 1
			while (t<=T) {
				index = selectindex(idt[.,2]:==uniquet[t])
				wt = asarray(W,uniquet[t])
				
				idtt = uniqrows(idt[index,1])
				
				ret[index,.] = SparseMultiply(wt, mat[index,.],issparse,idtt)
				///ret[index,.] = asarray(W,t) * mat[index,.]
				t++
			}
		}
		return(ret)
	}
end



//-------------------------------------------------------------------------------------------------
// order spatial weights mat
// -------------------------------------------------------------------------------------------------
/// in unbalanced case OrderW returns W as an array with W for each period.
/// if W is sparse, then W is returned as sparse (if unablanced then as array;
/// if balanced then as sparse matrix
capture mata mata drop OrderW()
mata:
	function OrderW(	///
						real matrix W,			///
						real matrix W_id,		///
						real matrix idt, 		///
						real matrix unbalanced,	///
						real scalar N,			///
						real scalar T,			///
						real matrix uniqueid,	///
						real scalar rownorm,	///
						real scalar issparse	///
					)
	{
		if (issparse==2) unablanced = 1

		if (unbalanced == 0 ) {
			Wfin = OrderWt(W,W_id,N,T,uniqueid,issparse,0)
		}
		else {
			Wfin = OrderWtUn(W,W_id,N,T,idt,issparse)
		}
		
		/// correct W_sparse to 1 if W_sparse == 2, because now it is normal sparse format
		if (issparse == 2) issparse = 1
	
		Wfin = Wnorm(Wfin,rownorm,issparse)
		
		return(Wfin)
	}

end

/// unbalanced case
capture mata mata drop OrderWtUn()
mata:
		function OrderWtUn(									///
							real matrix W,					///
							real matrix W_id,				///
							real scalar N,					///
							real scalar T,					///
							real matrix idt, 				///
							real scalar issparse			///
							)
		{
			Wfin = asarray_create("real",1)
			asarray_notfound(Wfin, 0)
			Tuniq = uniqrows(idt[.,2])
			tt = 1
			while (tt<=T) {
				/// loop over all periods and use OrderWt for each period				
				indext = selectindex(Tuniq[tt]:==idt[.,2])				
				idtt = idt[indext,1]
				Ni = rows(idtt)
				OrderWtt = OrderWt(W,W_id,Ni,1,idtt,issparse,tt)
				///_makesymmetric(OrderWtt)		
				asarray(Wfin,Tuniq[tt], OrderWtt)
				
				tt++
			}

			return(Wfin)

		}
		
end
/// return N*T x N spatial weight matrices
capture mata mata drop OrderWt()
mata:
	function OrderWt(							///
						real matrix W,			///
						real matrix W_id,		///	ids in W
						real scalar N,			/// number of obs (N) from idt
						real scalar T,			/// number of obs (T) from idt
						real matrix uniqueid,	/// uniqueids in t from data
						real scalar issparse,	///
						real scalar timeperiod	/// use if issparse = 2, i.e. sparse matrix is already in time format.
					)
	{
		
		real matrix SortW
		pointer(matrix) sorter
		pointer(matrix)	Wi


		if (issparse > 0) {
			"sparse case, order of W not important"		
			issparse	
			if (issparse == 1) {
				/// sparse matrix contains on first two rows id. Here we only need first row and we alter W_id to the one from sparse matrix				
				///W_id = uniqrows(W[.,1])
				W_id = W[.,(1,2)]
				Wi = &W

				sortsi = (1,2,3)
							
				
			}
			else if (issparse == 2) {
				/// select only rows for given period in time
				Wt = uniqrows(W[.,1])
				indext = selectindex(Wt[timeperiod]:==W[.,1])
				Wadj = W[indext,.]
				W_id = Wadj[.,(2,3)]				
				Wi = &Wadj
				sortsi = (2,3,4)				
			}
			/// now select rows which are permitted, W_id1 and W_id2 need to be elements of uniqueid
			W_id1 = J(rows((*Wi)),1,0)
			W_id2 = W_id1

			W_index1 = W_id1
			W_index2 = W_id2

			i = 1
			while (i<=N) {				
				ii= uniqueid[i]
				W_id1 = W_id1 :+ (W_id[.,1]:==ii)
				W_id2 = W_id2 :+ (W_id[.,2]:==ii)

				tmp1 = selectindex(W_id[.,1]:==ii)
				W_index1[tmp1,1] = J(rows(tmp1),1,i)

				tmp1 = selectindex(W_id[.,2]:==ii)
				W_index2[tmp1,1] = J(rows(tmp1),1,i)

				i++
			}
			index = W_id1 :* W_id2
			index = selectindex(index)
			W_sorted = (*Wi)[index,sortsi]
			/// add index
			W_sorted = W_sorted, W_index1[index], W_index2[index]

		}
		else if (issparse == 0) {
			SortW = J(N,1,.)
			i =1 
			while (i<=rows(uniqueid)) {
				ii = uniqueid[i]
				sorti = selectindex((ii:==W_id))
				if (sum(sorti) > 0) {
					SortW[i,1] = sorti
				}
				i++
			}
			SortW = SortW[selectindex(SortW:!=.)]
			W_sorted = W[SortW,SortW]
			W_sorted = J(T,1,1)#W_sorted
		}

		return(W_sorted)
	}
end
// -------------------------------------------------------------------------------------------------
// LoadData
// -------------------------------------------------------------------------------------------------
capture mata mata drop LoadData()
mata:
	function LoadData(	/// Input
						string scalar varnames,				///
						string scalar idtname,				///
						string scalar tousename,			///
						/// Returns
						real matrix y,						///
						real matrix x,						///
						real matrix idt,					///
						real scalar N,						///
						real scalar T, 						///
						real scalar K,						///
						real matrix panelsetup,				///
						real scalar unbalanced,				///
						real matrix uniqueid,				///
						real matrix uniquet 				///
						)
	{
		
		varname = st_tsrevar(tokens(varnames))
		y = st_data(.,varname[1],tousename)
		x = st_data(.,varname[1,2..cols(varname)],tousename)

		///real matrix X
		///real matrix Y
		///real matrix idt

		///st_view(y,.,varname[1],tousename)
		///st_view(x,.,varname[1,2..cols(varname)],tousename)
		///st_view(idt,.,idtname,tousename)

		K = cols(x)
		
		idt = st_data(.,idtname,tousename)

		uniquet = uniqrows(idt[.,2])
		T = rows(uniquet)	

		panelsetup = panelsetup(idt,1)

		stats = panelstats(idt)
		
		uniqueid = uniqrows(idt[.,1])
		N = rows(uniqueid)
		unbalanced = (stats[3]!=stats[4])
	}
end

// -------------------------------------------------------------------------------------------------
// Normalisation
// -------------------------------------------------------------------------------------------------
// subtract mean and divide by SD, do for each colum from matrix
capture mata mata drop Wnorm()
mata:
	function Wnorm(	 transmorphic matrix mat, 	/// 	
					real scalar type,			/// 1: for row normalisation
					real scalar issparse		/// 1: sparse matrix 
					)


	{
		if (issparse == 1) {
			"normalise sparse matrix"
			if (type == 1)	{
				loc = 1
			}
			else {
				loc = 2
			}

			if (eltype(mat) == "real") {
				mat = Wnormi(mat,loc)
			} 
			else {
				keys = sort(asarray_keys(mat),1)				
				nkeys = rows(keys)
				i = 1
				while (i<= nkeys) {
					Wi = asarray(mat,keys[i])
					Wi = Wnormi(Wi,loc)
					asarray(mat,keys[i],Wi)
					i++
				}
			}

		}
		else {
			pointer(function) fn

			if (type == 1) {
				fn = &myrowsum()
			}
			else {
				fn = &mycolsum()
			}
				
			if (eltype(mat) == "real") {
				mat = mat:/(*fn)(mat)
			}
			else {
				keys = asarray_keys(mat)
				nkeys = rows(keys)
				i = 1
				while (i<= nkeys) {
					Wi = asarray(mat,keys[i])
					Wi = Wi :/ (*fn)(Wi)
					asarray(mat,keys[i],Wi)
					i++
				}
			}
		}
		return(mat)
	}

end

capture mata mata drop Wnormi()
mata:
	function Wnormi(real matrix mat, real scalar loc)
	{
		uniq = uniqrows(mat[.,loc])
		Nuniq = rows(uniq)
		i = 1
		while (i<=Nuniq) {
			uniqi = uniq[i]
			index = selectindex(mat[.,loc]:==uniqi)
			rcsum = quadsum(mat[index,3])
			mat[index,3] = mat[index,3]:/rcsum
			i++
		}
		_editmissing(mat,0)
		return(mat)
	}
end

// -------------------------------------------------------------------------------------------------
// studentize
/// take T*K variable
// -------------------------------------------------------------------------------------------------
capture mata mata drop studentize()
mata:
	function studentize(real matrix mat)
	{
		nvar = cols(mat)
		nobs = rows(mat)

		if (nvar == 1) {
			var = quadmeanvariance(mat)
			mean = var[1]
			sd = sqrt(var[2])
	
			if (sd == 0) {
				res = J(nobs,1,0)
			}
			else {
				res = (mat :- mean) :/ sd
			}
		}
		else {
			var = quadmeanvariance(mat)
			mean = var[1,.]
			sd = sqrt(diagonal(var[(2..nvar+1),.]))
			res = (mat :- mean) :/ sd'
			///_editmissing(res,0)
		}
		return(res)
	}
end

/// TS function, loops over all time periods and studentizes
capture mata mata drop _studentizeT()
mata:
	function _studentizeT(real matrix mat, real matrix T, real matrix idt,  real matrix uniqueT)
	{
		
		t = 1
		while (t<=T) {
			index = selectindex(idt[.,2]:==uniqueT[t])
			mat[index,.] = studentize((mat[index,.]))
			t++
		}
	}

end



// -------------------------------------------------------------------------------------------------
// pointer functions to build in functions
// -------------------------------------------------------------------------------------------------
capture mata mata drop myrowsum()
mata: function myrowsum(x) return(rowsum(x))

capture mata mata drop mycolsum()
mata: function mycolsum(x) return(colsum(x))


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
		    local mult = floor(`i'/`reps'*100)-floor((`i'-1)/`reps'*100)
			forvalues s=1(1)`mult' {				
			     _dots 1 0
				 if (`i'-1)*`mult'+`s' >= 50 & (`i'-1)*`mult'+`s'-1 < 50 {
					 di as txt "" _col(51) %5.0f 50
				}
			}
			
		}
		else {
		   local check = round(`i'/`reps',0.01) - round((`i'-1)/`reps',0.01)
		 
		   if `check' > 0 {
		       _dots 1 0
		   }	   
			if `i'/`reps' >= 0.5 & (`i'-1)/`reps' < 0.5 {
				di as txt "" _col(51) %5.0f 50
			}
		}
		 
		if `i'/`reps' == 1 {
			di as txt "" _col(51) %5.0f 100
        }
end



// -------------------------------------------------------------------------------------------------
// R2 program
// -------------------------------------------------------------------------------------------------
/// Returns R2 and R2_adjusted
capture mata mata drop R2_calc()
mata:
	function R2_calc(real matrix residuals,real matrix Y, real scalar N, real scalar K)
	{
		rsqr1 = quadcross(residuals,residuals)		
		yt = Y :- mean(Y)

		rsqr2 = quadcross(yt,yt)
		
		R2 = 1- rsqr1/rsqr2
		
		rsqr1 = rsqr1/(rows(Y)-K)
		rsqr2 = rsqr2/(rows(Y)-1)
		R2adj = 1 - (rsqr1/rsqr2)

		return(R2,R2adj)
	}
	
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
