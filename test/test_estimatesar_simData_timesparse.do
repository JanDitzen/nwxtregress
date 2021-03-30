clear all
cls
cap program drop estimatesar

adopath + "C:\Users\jan\Dropbox (Personal)\Spatial Models in Stata\Stata\ado"

include "C:\Users\jan\Dropbox (Personal)\Spatial Models in Stata\Stata\ado\spatial_models_aux.ado"

*** create spatial weight matrix
mata W = J(11,11,0)

forvalues i = 2(1)11 {
	disp `i'
	mata W[`i',`=`i'-1']=1
}
mata W = W:+W'
mata W_id = 1::11

mata Wss = SparseDefine(W,W_id)
mata Ws = J(0,4,.)
forvalues i = 1(1)10 {
	mata Ws = Ws \ (J(rows(Wss),1,`i') , Wss)
}
mata Ws

set obs 100
egen id = seq(), block(10)
by id, sort: gen t = _n
drop if inlist(t,2,4,7) & id == 1
drop if id == 3 & t > 10
drop if id == 7

drawnorm x x2

gen y = 2*x + 3*x2 + rnormal()
xtset id t
timer on 99
estimatesar y x , spatial(y) weights(Ws,mata timesparse) trace
timer off 99

estimatesar y x x2 , spatial(y) weights(W,mata) trace


 
