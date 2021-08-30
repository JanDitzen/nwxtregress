/*
Programs for sparse matrices
Change Log:
30.08.2021: - bug in sparseXsparse fixed
*/

// -------------------------------------------------------------------------------------------------
// SparseDefine
// -------------------------------------------------------------------------------------------------
/// create sparse matrix with rows, cols and value. only non-zero values are selected
/// input: mat: matrix; uniquer: ids of mat (rows); uniquec: ids of mat (cols)
/// output: N*T(i) x 3 matrix, with id rows, id cols and value.

capture mata mata drop SparseDefine()
mata:
	function SparseDefine(real matrix mat, real matrix uniquer, real scalar W_sparse)
	{
		timer_on(21)
		/// fastest way to convert matrix to missings and then use rownonmissing		
		"sparse define prog. start"
		mat2 = editvalue(mat,0,.)
		row_nonmiss = rownonmissing(mat2)		
		nonmiss = sum(row_nonmiss)

		if (nonmiss > 0) {			
			"matrix converted to sparse"
			
			index_row_nonmiss = selectindex(row_nonmiss:>0)			
			num_rows = rows(index_row_nonmiss)			

			/// initalise result sparse matrix
			
			mat_sp = J(nonmiss,3,.)
			
			i = 1
			counter = 1
			while (i <= num_rows) {	

				idi = uniquer[i]
				
				element_i = index_row_nonmiss[i,1]
				

				numi = row_nonmiss[element_i,1]

				index = selectindex(mat2[element_i,.]:!=.)

				counteri = counter
				counter = counter + numi - 1
				

				ti = uniquer[index]

				mat_sp[|counteri,1 \ counter,3|] = J(numi,1,idi) , ti , mat2[element_i,index]'
				
				counter++
				i++
			}
			W_sparse = 1
			return(mat_sp)			
		}
		else {
			return(mat)
		}

		timer_off(21)
	}

end

// -------------------------------------------------------------------------------------------------
// Sparse Matrices - Multiplication
// -------------------------------------------------------------------------------------------------
/// Sparse is N(s)x3, X is Nxk, Nuniq is Nx1 matrix with unique identifier of cross-sections, 
/// N(s) is non zero elements in W; Sparse has 
/// row, column, value
/// program determines if matrix is sparse automatically
capture mata mata drop SparseMultiply()
mata:
	function SparseMultiply(real matrix Sparse, real matrix X, real scalar issparse,| real matrix XNuniq)
	{
		timer_on(22)
		if (issparse == 0) {
			result = Sparse*X 
		}
		else {
			if (args() == 3) {
				/// case that X is a scalar and W sparse, return sparse
				if (cols(Sparse)==3) {
					result = Sparse[.,(1,2)], Sparse[.,3]:*X
				}
				else {
					result = Sparse[.,(1,2)], Sparse[.,3]:*X , Sparse[.,(4,5)]
				}
			}
			else {
				timer_on(26)
				Nuniq = uniqrows(Sparse[.,1],1)	
				/// X is not a scalar, return non sparse
				result = J(rows(XNuniq),cols(X),0)
				
				i = 1
				ii = 1
				
				max_xW = max( (max(XNuniq[.,1]), max(Nuniq[.,1]) ) ) 
				
				K = cols(X)
				res_matX = J(max_xW,K,.)
				res_matX[XNuniq[.,1],.] = X
				timer_off(26)

				while (i<= rows(Nuniq)) {	
				
					timer_on(27)
					iii = ii
					ii = ii  + Nuniq[i,2] - 1
					
					Wi = Sparse[|iii,1 \ ii,3|]	
							
					posi = selectindex(Nuniq[i,1] :==XNuniq[.,1])
					posX = J(rows(Wi[.,2]),1,.) 
					timer_off(27)
					ix = 1
					timer_on(29)
					while (ix <= rows(Wi[.,2])) {
						///posi = selectindex(Wi[ix,2]:==XNuniq[.,1])
						
						index = selectindex(Wi[ix,2]:==XNuniq[.,1])
						
						if (sum(index) > 0 ) {
							posX[ix] = index
						} 
						else {
							"nosum!"

						}
						ix++
					}
					timer_off(29)
					///if (sum(posX) > 0) {				
						posX = posX[selectindex(posX:!=.)]
						result[posi,.] = Wi[.,3]' * X[posX,.]
				
					


					ii++
					i++
				}
			}
		}
		timer_off(22)
		return(result)
		
	}
end

/// Transpose of sparse matrix
capture mata mata drop SparseTranspose()
mata:
	function SparseTranspose(transmorphic matrix spmat) 
	{
		timer_on(23)
		if (eltype(spmat) == "real") {
			spmatn = spmat[.,2],spmat[.,1],spmat[.,3]
		}
		else {
			spmatn = asarray_create("real",1)
			keys = asarray_keys(mat)
			nkeys = rows(keys)
			i = 1 
			while (i<=nkeys) {
				spmati = asarray(spmat,i)
				spmati = spmati[.,2],spmati[.,1],spmati[.,3]
				asarray(spmatn,keys[i],spmati)
			}
		}
		timer_off(23)
		return(spmatn)
		
	} 

end

/// multiply two sparse matrices. Always assumes that W * W2 or W :* W2
/// program messy, can be improved!
capture mata mata drop SparseXSparse()
mata:
	function SparseXSparse(transmorphic matrix Sparse, transmorphic matrix Sparse2, real scalar issparse,real scalar pointbypoint)
	{
		timer_on(24)
		if ( (eltype(Sparse) == "real") & ( eltype(Sparse2) == "real") )  {
			returni = SparseXSparsei(Sparse,Sparse2,issparse,pointbypoint)
		}
		else  {
			returni = asarray_create("real",1)
			if ((eltype(Sparse) == "real") & (eltype(Sparse2) != "real")) {
				keys = asarray_keys(Sparse2)
				nkeys = rows(keys)
				i = 1
				while (i<= nkeys) {
					sparsei = asarray(Sparse2,keys[i])
					sparsei = SparseMultiply(Sparse,sparsei,issparse,pointbypoint)
					asarray(returni,keys[i],sparsei)
					i++
				} 
			}
			if ((eltype(Sparse) != "real") & (eltype(Sparse2) == "real")) {
				keys = asarray_keys(Sparse)
				nkeys = rows(keys)
				i = 1
				while (i<= nkeys) {
					sparsei = asarray(Sparse,keys[i])
					sparsei = SparseMultiply(sparsei,Sparse2,issparse,pointbypoint)
					asarray(returni,keys[i],sparsei)
					i++
				} 
			}
			if ((eltype(Sparse) != "real") & (eltype(Sparse2) != "real")) {
				keys = asarray_keys(Sparse)
				keys2 = asarray_keys(Sparse2)
				nkeys = rows(keys)
				nkeys2 = rows(keys2)
				if (nkeys:==nkeys2) {
					i = 1
					while (i<=nkeys) {
						sparsei = SparseMultiply(asarray(Sparse,keys[i]),asarray(Sparse2,keys2[i]),issparse,pointbypoint)
						asarray(returni,keys[i],sparsei)
						i++
					}
				}
				else {
					"keys mismatch"
				}
			}
		}
timer_off(24)
		return(returni)
		
	}
end

capture mata mata drop SparseXSparsei()
mata:
	function SparseXSparsei(transmorphic matrix Sparse, transmorphic matrix Sparse2, real scalar issparse,real scalar pointbypoint)
	{
		timer_on(25)
		if (issparse == 0) {
			if (pointbypoint == 0) {
				result = Sparse * Sparse2
			}
			else {
				result = Sparse :* Sparse2
			}
		}
		else {
			if (pointbypoint == 0 ) {
				"not available"
			}
			else {
				/// point by point. loop over each element and find corresponding element
				/// elements in neiter will remain empty
				SparseN = J(rows(Sparse),3,.)
				i = 1
				while (i<=rows(Sparse)) {
					eli = Sparse[i,.]
					indexi = selectindex((eli[1]:==Sparse2[.,1] :*eli[2]:==Sparse2[.,2]))
					if (sum(indexi) > 0) {
						SparseN[i,.] = Sparse[i,(1,2)] , Sparse[i,3] * Sparse2[indexi,3]
					}
					i++
				}
				if (sum(SparseN)>0) {
					result = SparseN[selectindex(SparseN[.,1]:!=.),.]
				}
			}
		}
		timer_off(25)
		return(result)
		
	}
end

// -------------------------------------------------------------------------------------------------
// make non sparse
// -------------------------------------------------------------------------------------------------
capture mata mata drop SparseUndefine() 
mata:
	function SparseUndefine(real matrix sparse, real matrix idt, real scalar issparse)
	{
		if (issparse == 1) {
			
			uniid1 = panelsetup(idt[.,1],1)
			N = rows(uniid1)
		
			result = J(N,N,0)

			sparseSorted = sort(sparse,(4,5))
			indexBig = panelsetup(sparseSorted,4)

			i = rows(indexBig)
			while (i>0) {

				col_i = panelsubmatrix(sparseSorted,i,indexBig)

				if (rows(col_i) > 0) {
					result[i,col_i[.,5]] = col_i[.,3]'				
				}
				
				i--
			}
			
			return(result)
		}
		else {
			return(sparse)
		}
	}
end
