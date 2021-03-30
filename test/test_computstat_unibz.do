clear all
cls
cap program drop estimatesar

adopath + "C:\Users\\`c(username)'\\Dropbox\Spatial Models in Stata\Stata\ado"

** get W
import delim "C:\Users\\`c(username)'\\Dropbox\Spatial Models in Stata\Data\TNIC2017.csv" 
putmata w = (*), replace

import delim "C:\Users\\`c(username)'\\Dropbox\Spatial Models in Stata\Data\compustat_density", clear 

xtset gvkey fyear

estimatesar l_rd  l_logsale l_age l_cash roa q l_leverage , spatial(l_rd) weights(w, mata timesparse) draws(50) gridlength(400) trace
 
