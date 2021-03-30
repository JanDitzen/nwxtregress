clear all
// -------------------------------------------------------------------------------------------------
// Sparse Matrices - convert to sparse
// -------------------------------------------------------------------------------------------------
capture mata mata drop SparseDefine()
mata:
	function SparseDefine(real matrix mat)
	{
	
		/// fastest way to convert matrix to missings and then use rownonmissing		
		timer_on(80)
		mat2 = editvalue(mat,0,.)
		timer_off(80)
		timer_on(81)
		row_nonmiss = rownonmissing(mat2)		
		timer_off(81)
		timer_on(82)
		nonmiss = sum(row_nonmiss)
		timer_off(82)
		if (nonmiss > 0) {			
			
			timer_on(83)
			index_row_nonmiss = selectindex(row_nonmiss:>0)			
			num_rows = rows(index_row_nonmiss)			

			/// initalise result sparse matrix
			mat_sp = J(nonmiss,3,.)
			timer_off(83)
			i = 1
			counter = 1
			timer_on(84)
			while (i <= num_rows) {	
						
				element_i = index_row_nonmiss[i,1]
				numi = row_nonmiss[element_i,1]

				index = selectindex(mat2[element_i,.]:!=.)
			
				counteri = counter
				counter = counter + numi - 1
				
				mat_sp[|counteri,1 \ counter,3|] = J(numi,1,element_i) , index' , mat2[element_i,index]'
				
				counter++
				i++
			}
			
			
			timer_off(84)
			return(mat_sp)
			
		}
		else {
			return(mat)
		}


	}

end

// -------------------------------------------------------------------------------------------------
// Sparse Matrices - Multiplication
// -------------------------------------------------------------------------------------------------
/// Sparse is N(s)x3, X is Nxk, Nuniq is Nx2 matrix with unique identifier of cross-sections and frequence, ReturnSparse indicator is return is sparse matrix or not
/// N(s) is non zero elements in W; Sparse has 
/// row, column, value
/// Returns 
capture mata mata drop SparseMultiply()
mata:
	function SparseMultiply(real matrix Sparse, real matrix X, real matrix Nuniq, real matrix XNuniq, real scalar ReturnSparse)
	{
		/// initalise result (result)
		if (ReturnSparse == 0) {
			result = J(rows(XNuniq),cols(X),0)
			i = 1
			ii = 1
			timer_on(70)
			while (i<= rows(Nuniq)) {
				iii = ii
				ii = ii  + Nuniq[i,2] - 1
				Wi = Sparse[|iii,1 \ ii,3|]
				posi = selectindex(Nuniq[i,1] :==XNuniq[.,1])
				result[posi,.] = Wi[.,3]' * X[Wi[.,2],.]
				ii++
				i++
			}
			timer_off(70)
			return(result)
		}
		

	}
end

* edit options here
mata N =  10000
* share of non zero items for % of cross-sections
mata share = 0.01

**** no changes after here needed
mata W = J(N,N,0)
timer on 99
mata W1 = editvalue(W,0,.)
mata xx = vec(W1)
timer off 99
timer list
mata:
	i = 1
	while (i<=N*share) {
		W[1,i] = 1
		W[i,1] = 1
		i++
	} 
end



	timer on 19
	timer on 11
	mata Sw = SparseDefine(W)
	timer off 11

	mata Wid = uniqrows(Sw[.,1],1)
	mata X = rnormal(N,3,0,1)
	mata Xn = (1::N) , J(N,1,1)

	timer on 12
	mata Wid
	mata s2 = SparseMultiply(Sw,X,Wid,Xn,0)
	timer off 12
	timer off 19

	timer on 20
	mata Sx = W*X
	timer off 20

	mata Sx == s2

}

timer list

