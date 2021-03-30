clear all
cls
cap program drop estimatesar
cap program drop nwxtregress_estat

adopath + "C:\Users\\`c(username)'\\Dropbox (Personal)\Spatial Models in Stata\Stata\ado"

** get W
import delim "C:\Users\\`c(username)'\\Dropbox (Personal)\Spatial Models in Stata\Data\TNIC2017.csv" 
putmata w = (*), replace
 
use "C:\Users\\`c(username)'\\Dropbox (Personal)\Spatial Models in Stata\Data\small", clear
///import delim "C:\Users\\`c(username)'\\Dropbox (Personal)\Spatial Models in Stata\Data\compustat_density", clear 

xtset gvkey fyear

cap log using "D:\temp\com", replace smcl

timer on 2
nwxtreg inv  l_logsale l_age l_cash zscore roa q l_leverage ,  dvarlag(w, mata timesparse) draws(2000) gridlength(10000) nomit(500) uselud
timer off 2


timer on 3
estimatesar inv  l_logsale l_age l_cash zscore roa q l_leverage , spatial(inv) weights(w, mata timesparse) draws(2000) gridlength(500) nomit(500) 
timer off 3
estat impact 

/*
timer on 3
estimatesar l_rd  l_logsale l_age l_cash roa q l_leverage , spatial(l_rd) weights(w, mata timesparse) draws(110) gridlength(100)  uselud
timer off 3
log close
timer list
*/
