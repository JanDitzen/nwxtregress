*** testing all the sub programs
clear all

cd "C:\Users\jan\Dropbox (Personal)\Spatial Models in Stata\Stata\ado\"

do "C:\Users\jan\Dropbox (Personal)\Spatial Models in Stata\Stata\ado\estimatesar.ado"


local N = 5
local T = 4

set obs `=`N'*`T''

egen id = seq(), block(`T')
by id, sort: gen t = _n

drop if id == 1 & t == 2
drop if id == 3 & t == 4

mata W = rnormal(`N',`N',0,1)
mata W = W'W
mata W_id = (2,1,4,5,3)

mata idt = st_data(.,"id t")

mata N = 5 
mata T = 4
mata uniq = uniqrows(idt[.,1]) 


mata OrderWt(W,W_id,N,T,uniq,0)
mata Ws= OrderWtUn(W,W_id,N,T,idt,0)

mata (asarray(Ws,2))
mata Wnorm(asarray(Ws,2),1,0)

mata Ws = OrderW(W,W_id,idt,1,N,T,uniq,1,0)
mata asarray(Ws,2)

drawnorm x1 x2

gen y = x1 + x2 + rnormal()
gen touse = 1
mata LoadData("y x1 x2","id t","touse",y=.,x=.,idt=.,N=.,T=.,K=.,ps=.,uinbal=.,uN=.,uT =.)
mata y
mata _studentizeT(y,T,idt,uT)
mata y
mata _studentizeT(x,T,idt,uT)
mata ccep(y,x,b=.,e=.)

** take only first array
///mata Ws = asarray(Ws,1)
mata W_id = uniqrows(W_id',1)
mata Wy = spatialLag(y,Ws,idt,T,W_id,uN,uT,0)
mata XX = quadcross(x,x)
mata XWy = quadcross(x,Wy)
mata ccep(y,x,b=.,e=.,XWy,XX)
mata:
	rho = (-99..99)/100
	eltype(Ws)
	asarray(Ws,1)
	rhoW(Ws,rho[1],1,0)

end

/// now use sparse mata matrix
mata W_id = (2,1,4,5,3)
mata Ws = SparseDefine(W,W_id,W_id)
mata Ws= OrderW(Ws,W_id,idt,1,N,T,W_id,0,1)
mata Wy = spatialLag(y,Ws,idt,T,W_id,uN,uT,0)
