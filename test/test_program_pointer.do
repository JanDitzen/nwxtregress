capture mata mata drop poin()
capture mata mata drop testabcded()
mata:
	function poin(pointer(real matrix) poin) 
	{
		res = mean((*poin))
		return(res)
	}


	X = rnormal(10,1,0,1)
	
	mean(X)

	poin(&X)

	Xp = &X

	Xp
	(*Xp)

	function testabcded(pointer(real matrix) poin , real scalar ind,|real scalar ttt)
	{
		v = (*poin)[ind] * 2
		(*poin)[ind] = v
		if (args() == 2 ) {
			testabcded(poin,ind+1,1)
		}
	}

	s = (1,2,3)
	testabcded(&s,1)
	s
end
