// -------------------------------------------------------------------------------------------------
// DirectIndirect Effects
// -------------------------------------------------------------------------------------------------
capture mata mata drop calc_DIEffects()
mata:
	function calc_DIEffects(string scalar arrayname)
	{
		
		bdraws = asarray(arrayname,"bdraws")
		pdraws = asarray(arrayname,"pdraws")
		trace = asarray(arrayname,"trace")

		N = asarray(arrayname,"dimensions")[1]
		T = asarray(arrayname,"dimensions")[4]
		K = asarray(arrayname,"dimensions")[3]

		ndraws = asarray(arrayname,"MCMC")[1]
		///nomit = asarray(arrayname,"MCMC")[2]
		nomit = 0

		issparse = st_numscalar("e(issparse)")
(N,T,K,ndraws,nomit,issparse)
		trace
		if (trace==0) {

			W = asarray(arrayname,"W")
			idt = asarray(arrayname,"idt")
			trace =  BarryPace(W,idt,T,N,50,100,issparse)
		}
		"trace"
trace
		trs = (1\ trace:/T)

		ntrs = length(trs)
		trbig = trs
		ree = (0..(ntrs-1))

		rmat = J(1,ntrs,.)
		
		/// arrays, each coeff (K) has one key in array
		totalm = asarray_create("real",1)
		direct = asarray_create("real",1)
		indirect = asarray_create("real",1)

		rmat = (pdraws#J(1,ntrs,1)):^ree

		j=1
		while (j<=K) {			
			total_j = bdraws[.,j]:*rmat
			direct_j = total_j * trbig
			indirect_j= total_j :- direct_j

			asarray(totalm,j,total_j)
			asarray(direct,j,direct_j)
			asarray(indirect,j,indirect_j)

			j++

		}
		
		totalo = calc_Effect(totalm,ndraws,nomit,K,ntrs)
		direct = calc_Effect(direct,ndraws,nomit,K,ntrs)
		indirect = calc_Effect(indirect,ndraws,nomit,K,ntrs)
		
		/// write to array
		asarray(arrayname,"total",totalo)
		asarray(arrayname,"indirect",indirect)
		asarray(arrayname,"direct",direct)


	}
end

capture mata mata drop  calc_Effect()
mata:
	function calc_Effect(transmorphic matrix mat, real scalar ndraws, real scalar nomit, real scalar K, real scalar ntrs)
	{
		outmat = J(K,4,0)
		save = J(ndraws-nomit,K,0)
		i=1
		while (i<=K) {

			tmp = asarray(mat,i)
			rows(tmp),cols(tmp)
			///total_mean = mean(tmp)
			"colsum"
			///total_std = sqrt(diagonal(quadvariance(tmp)))
			total_sum = (quadcolsum(tmp'))'
			rows(total_sum),cols(total_sum)
			///cum_mean = mm_colrunsum(total_mean)
			///cum_std = mm_colrunsum(total_std)

			///save[.,i] = total_sum

			cmean = mean(total_sum)

			smean = quadvariance(total_sum)
			

			upper = mm_quantile(total_sum,1,st_numscalar("c(level)"))
			lower = mm_quantile(total_sum,1,1-st_numscalar("c(level)")) 

			outmat[i,.] = (cmean, smean,lower, upper)
			i++
		}

		return(outmat)
	}

end

// -------------------------------------------------------------------------------------------------
// Barry Pace Trick
// -------------------------------------------------------------------------------------------------
capture mata mata drop BarryPace()
mata:
	function BarryPace(transmorphic matrix W, real matrix idt, real scalar T, real scalar N,real scalar uiter, real scalar maxorder,real scalar issparse)
	{
		"start BP"
		if (eltype(W) == "real") {
			Wi = W
		}
		

		tracew = J(maxorder,1,0)
		traces = J(maxorder,1,0)
		Tuniq = uniqrows(idt[.,2],1)
		
		t=1
		timer_on(91)
		while(t<=T) {			
			if (eltype(W) != "real") {
				Wi = asarray(W,Tuniq[t,1])
			}

			nobs = Tuniq[t,2]
			rv = rnormal(nobs,uiter,0,1)
			wjjju = rv
			indext = selectindex(Tuniq[t,1]:==idt[.,2])
			indext = idt[indext,1]
			j = 1
			timer_on(92)
			while (j<=maxorder) {				
				timer_on(93)
				wjjju = SparseMultiply(Wi,wjjju,issparse,indext) 
				timer_off(93)
				tracew[j] = mean(mean(rv:*wjjju)')
				j++
			}
			timer_off(92)
			timer_on(94)
			tracew
			tracew[2] = sum(sum(SparseXSparse(SparseTranspose(Wi),Wi,issparse,1))) :/ nobs
			timer_off(94)
			
			traces = traces + tracew

			traces[1] = 0
			t++
		}
		timer_off(91)

		"BP done"
		return(traces)
	}

end

// -------------------------------------------------------------------------------------------------
// LUD program
// -------------------------------------------------------------------------------------------------

capture mata mata drop LUD()
mata:
	function LUD(transmorphic matrix W, real matrix rgrid,  real scalar N, real scalar t, real scalar i, real scalar issparse, real matrix idt)
	{
			timer_on(91)
				B = rhoWIminus(W,rgrid[i],N,t,issparse)
			timer_off(91)
			timer_on(92)
				B = SparseUndefine(B,idt,issparse)
			timer_off(92)
			timer_on(93)
				lud(B,L=0,U=0,p=0)
			timer_off(93)
				/// LU decomposition, check if results are same
				lu1 = L+U-I(rows(B))
				lu1 = diagonal(lu1)
				lu1 = sum(log(lu1))
				return(lu1)
			
	}
end
