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


estimatesar inv  l_logsale l_age l_cash zscore roa q l_leverage , spatial(inv) weights(w, mata timesparse) draws(2000) gridlength(500) nomit(500) uselud

gen smpl = e(sample)

mata W = asarray(nwxtreg_output,"W")
mata asarray(W,7) 
** element 2 is zero!; in Matlab it is 1!
/*

val =

  (13,1)        1
  (51,19)       1
  (48,30)       1
   (5,46)       1
  (19,51)       1
*/

*** next, check Y and X
mata y = asarray(nwxtreg_output,"Y")
mata y[selectindex(y[.,3]:==1),.]

mata x = asarray(nwxtreg_output,"X")
mata x[selectindex(x[.,cols(x)]:==1),.]
