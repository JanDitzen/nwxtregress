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
	syntax anything(name=wname_init) , [mata id(string) sparse timesparse isdep] varlist(string) arrayname(string)

	tempname w_ary wname
	mata: `w_ary' = asarray_create()


	if "`sparse'" != "" & "`timesparse'" != "" {
	*** time sparse implies sparse
			local sparse ""
	}

	if "`mata'" == "" {
		tempname w_id 
		spmatrix matafromsp `wname' `w_id' = `wname_init'
		mata asarray(`w_ary',"W_id",`w_id')
		mata asarray(`w_ary',"W",`wname')
	}
	else {
		mata asarray(`w_ary',"W",`wname_init')
		mata `wname' = `wname_init'
	}

	if "`id'" == "" {
		if "`sparse'`timesparse'" == "" {
			/// here correct!
			mata asarray(`w_ary',"W_id",(1::rows(`wname')))
		}
		else if "`sparse'" != "" & "`timesparse'" == "" {
			mata asarray(`w_ary',"W_id",uniqrows(`wname'[.,1]))
		}
		else if "`sparse'" == "" & "`timesparse'" != "" {
			mata asarray(`w_ary',"W_id",uniqrows(`wname'[.,2]))
		} 
	}
	else {
		mata asarray(`w_ary',"W_id",`id')
	}

	*** 1 if W is sparse
	if "`sparse'`timesparse'" == "" {
		mata asarray(`w_ary',"sparse",0)
	}
 	else if "`sparse'" != "" {
 		mata asarray(`w_ary',"sparse",1)
	}
	else if "`timesparse'" != ""  {
		mata asarray(`w_ary',"sparse",2)
	}

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
	mata mata drop `w_ary' `wname'
end


// -------------------------------------------------------------------------------------------------
// drawnorm
// -------------------------------------------------------------------------------------------------
capture mata mata drop drawnorm()
mata:
	function  drawnorm(real scalar n,real mean, real matrix variance) 
		{
			res = invnormal(uniform(n,cols(variance)))*cholesky(variance)' :+ mean
			return(res)	
			
		}

end


// -------------------------------------------------------------------------------------------------
// Inverter Programs
// -------------------------------------------------------------------------------------------------
// Mata utility for sequential use of solvers
// Default is cholesky;
// if that fails, use QR;
// if overridden, use QR.
// By Mark Schaffer 2015
capture mata mata drop cholqrsolve()
mata:
	function cholqrsolve (  numeric matrix A,
							numeric matrix B,
						  | real scalar useqr)
	{
			
			if (args()==2) useqr = 0
			
			real matrix C

			if (!useqr) {
					C = cholsolve(A, B)
					if (C[1,1]==.) {
							C = qrsolve(A, B)
					}
			}
			else {
					C = qrsolve(A, B)
			}
			return(C)

	};
end


capture mata mata drop cholqrinv()
mata:
	function cholqrinv (  numeric matrix A,
						  | real scalar useqr)
	{
			if (args()==2) useqr = 0

			real matrix C

			if (!useqr) {
					C = cholinv(A)
					if (C[1,1]==.) {
							C = qrinv(A)
					}
			}
			else {
					C = qrinv(A)
			}
			return(C)

	};
end

///Program for matrix inversion.
///Default is cholesky
///if not full rank use invsym (Stata standard) 
///and obtain columns to use
///options: 
///1. if columns are specified, force use invsym
///2. allow for old method (cholinv, if fails qrinv)
///output
///return: inverse
///indicator for rank (1x2, rank and rows), which method used and variables used

capture mata mata drop m_xtdcce_inverter()
mata:
	function m_xtdcce_inverter(	numeric matrix A,
								| real scalar useold,
								real matrix rank,
								real matrix coln,
								string scalar method)
								
	{
		real matrix C
		
		if (args() == 1) {
			useold = 0
			coln = 0
		}
		if (args() == 2){
			coln = 0
		}
		if (useold == 2) {
			coln = (1..cols(A))
		}
		
		if (useold == 1) {			
			C = cholqrinv(A)
			qrinv(A,rank)
			method = "cholqr"		
		}
		else {
			if (coln[1,1] == 0) {
				/// calculate rank seperate. if A is not full rank, cholinv still produces results
				/// 1..cols(A) makes sure variables from left are not dropped
				C = invsym(A,(1..cols(A)))
				rank = rows(C)-diag0cnt(C)
				
				if (rank < rows(A)) {	
					/// not full rank, use invsym
					method = "invsym"
					coln = selectindex(colsum(C:==0):==rows(C):==0)			
				}
				else {
					/// full rank use cholsolve
					C = cholinv(A)
					method = "chol"
				}				
			}
			else {
				C = invsym(A,coln)
				rank = rows(C)-diag0cnt(C)
				method = "invsym"
			}			
		}
		rank = (rank, rows(C))
		return(C)
	}

end
/// same as inverter, rank is for matrix A (which is inverted) 
capture mata mata drop m_xtdcce_solver()
mata:
	function m_xtdcce_solver(	numeric matrix A,
								numeric matrix B,
								| real scalar useold,
								real matrix rank,
								real matrix coln,
								string scalar method)
								
	{
		real matrix C
		real scalar A1
		
		if (args() == 2) {
			useold = 0
			coln = 0
		}
		if (args() < 5){
			coln = 0
		}		
		
		if (useold == 2) {
			coln = (1..cols(A))
		}
		
		if (useold == 1) {			
			C = cholqrsolve(A,B)
			qrinv(A,rank)
			method = "cholqr"
			rank = (rank, rows(C))
		}
		else {
			if (coln[1,1] == 0) {
				
				/// calculate rank seperate. if A is not full rank, cholsolve still produces results
				/// 1..cols(A) makes sure variables from left are not dropped
				A1 = invsym(A,(1..cols(A)))
				rank = rows(A1)-diag0cnt(A1)
				
				if (rank < rows(A)) {	
					/// not full rank, solve by hand
					C = A1 * B
					method = "invsym"
					coln = selectindex(colsum(A1:==0):==rows(A1):==0)			
				}
				else {
					/// full rank use cholsolve
					C = cholsolve(A,B)
					method = "chol"
					coln = 0
				}
			}
			else {
				/// coln is defined, use invsym on specified columns
				A1 = invsym(A,coln)
				C = A1 * B
				method = "invsym"
				rank = rows(A1)-diag0cnt(A1)
			}
			rank = (rank, rows(A1))
		}		
		return(C)		
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
						real matrix trace, 			///
						real scalar sdm				///
						)
	{
		///rgrid = (-(1-1/gridlength)*gridlength..(1-1/gridlength)*gridlength):/gridlength
		rgrid = range(-0.99,0.99,1/gridlength)
		ngrid = rows(rgrid)
		detval = J(ngrid,2,.)	
		
		cmd = sprintf("noi _dotspct 0 , title(Initialise Grid) reps(%s)",strofreal(ngrid))	
		stata(cmd)

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

				cmd = sprintf("noi _dotspct %s 0 , reps(%s)",strofreal(i),strofreal(ngrid))
				stata(cmd)

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
				trace =  BarryPace(W,idt,T,N,BarryPace[1],BarryPace[2],issparse,sdm)
				detm = sum(log(trace))
				///}
				timer_off(90)
				detval[i,.] = rgrid[i] , detm
				cmd = sprintf("noi _dotspct %s 0 , reps(%s)",strofreal(i),strofreal(ngrid))
				stata(cmd)
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
			add = (el,el,J(rows(el),1,1))
			if (cols(ret) > 3) {
				el1 = uniqrows(ret[.,4])
				el2 = uniqrows(ret[.,5])
				el22 = uniqrows((el1\el2))
				add = add , el22,el22
			}
			ret = ret \ add
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

		///if (unbalanced == 0 ) {
		///	Wfin = OrderWt(W,W_id,N,T,uniqueid,issparse,0)
		///}
		///else {
			"W id"
			W_id
			Wfin = OrderWtUn(W,W_id,N,T,idt,issparse)
		///}
		
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
			Tuniq
			Tuniq = uniqrows(idt[.,2])
			tt = 1
			
			while (tt<=T) {			
				/// loop over all periods and use OrderWt for each period
				indext = selectindex(Tuniq[tt]:==idt[.,2])				
				idtt = idt[indext,1]
				Ni = rows(idtt)				
				OrderWtt = OrderWt(W,W_id,Ni,1,idtt,issparse,Tuniq[tt])				
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
		"is sparse"
		issparse
		cols(W),rows(W)
		if (issparse > 0) {
			"sparse case, order of W not important"		
			if (issparse == 1) {
				/// sparse matrix contains on first two rows id. Here we only need first row and we alter W_id to the one from sparse matrix				
				///W_id = uniqrows(W[.,1])
				W_idi = W[.,(1,2)]
				Wi = &W

				sortsi = (1,2,3)
							
				
			}
			else if (issparse == 2) {
				/// time sparse!
				/// select only rows for given period in time
				///Wt = uniqrows(W[.,1])
				///indext = selectindex(Wt[timeperiod]:==W[.,1])
				W[1..20,1]
				indext = selectindex(timeperiod:==W[.,1])

				Wadj = W[indext,.]
				W_idi = Wadj[.,(2,3)]				
				Wi = &Wadj
				sortsi = (2,3,4)				
			}
			else {
				W_idi = W_id
			}
			"start ordering"
			rows((*Wi)),cols((*Wi))
			///W_id1
			///W_id2
			/// now select rows which are permitted, W_id1 and W_id2 need to be elements of uniqueid
			W_id1 = J(rows((*Wi)),1,0)
			W_id2 = W_id1

			W_index1 = W_id1
			W_index2 = W_id2
			
			i = 1
			while (i<=N) {				
				ii= uniqueid[i]
				W_id1 = W_id1 :+ (W_idi[.,1]:==ii)
				W_id2 = W_id2 :+ (W_idi[.,2]:==ii)

				tmp1 = selectindex(W_idi[.,1]:==ii)
				W_index1[tmp1,1] = J(rows(tmp1),1,i)

				tmp1 = selectindex(W_idi[.,2]:==ii)
				W_index2[tmp1,1] = J(rows(tmp1),1,i)

				i++
			}			
			index = W_id1 :* W_id2
			index = selectindex(index)
			W_sorted = (*Wi)[index,sortsi]
			/// add index
			W_sorted = W_sorted, W_index1[index], W_index2[index],W_sorted[.,3]

		}
		else if (issparse == 0) {
			"not sparse!"
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
			"sort"
			
			SortW = SortW[selectindex(SortW:!=.)]
cols(W),rows(W)
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
				"row norm"
				loc = (1,2)
				
			}
			else {
				"col norm"
				loc = (2,1)
			}

			if (eltype(mat) == "real") {
				"is real"
				mat = Wnormi(mat,loc)
			} 
			else {
				"non real"
				keys = sort(asarray_keys(mat),1)
				nkeys = rows(keys)
				i = 1
				while (i<= nkeys) {
					Wi = asarray(mat,keys[i])
					
					/// remove diagonals
					Wi = Wi[selectindex(Wi[.,1]:!=Wi[.,2]),.]
					Wi = Wnormi(Wi,loc,ki=.)
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
	function Wnormi(real matrix mat, real matrix loc,real matrix kk)
	{
		

		mat = sort(mat,(loc))
		index = panelsetup(mat,loc[1])
		i = panelstats(index)[1]
kk = .,.,.
		while (i>0) {			
			wi = panelsubmatrix(mat,i,index)
			ws = quadsum(wi[.,3])
			ws = ws + (ws==0)	
			kk = kk \ (i,ws,rows(wi))
			///i
			
			
			mat[|index[i,1],3 \ index[i,2],3|] = wi[.,3] :/ ws
			///
			i--
		}

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
