clear all
cls
cap program drop estimatesar

adopath + "C:\Users\jan\Dropbox (Personal)\Spatial Models in Stata\Stata\ado"

*** create spatial weight matrix
mata W = J(11,11,0)

forvalues i = 2(1)11 {
	disp `i'
	mata W[`i',`=`i'-1']=1
}
mata W = W:+W'
mata W_id = 1::11

set obs 100
egen id = seq(), block(10)
by id, sort: gen t = _n
drop if inlist(t,2,4,7) & id == 1
drop if id == 3 & t > 10
drop if id == 7

drawnorm x x2

gen y = 2*x + 3*x2 + rnormal()
xtset id t
estimatesar y x , spatial(y) weights(W,mata id(W_id)) trace
estimatesar y x x2 , spatial(y) weights(W,mata) trace


 
