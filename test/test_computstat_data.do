clear all
cls
cap program drop estimatesar

adopath + "C:\Users\\`c(username)'\\Dropbox (Personal)\Spatial Models in Stata\Stata\ado"

import delim "C:\Users\\`c(username)'\\Dropbox (Personal)\Spatial Models in Stata\Data\compustat_density", clear 

xtset gvkey fyear

egen id = group(gvkey)
keep if id < 100
save "C:\Users\\`c(username)'\\Dropbox (Personal)\Spatial Models in Stata\Data\small", replace
	
export excel  "C:\Users\\`c(username)'\\Dropbox (Personal)\Spatial Models in Stata\Data\small.xlsx", replace firstrow(var)

 export delimited  "C:\Users\\`c(username)'\\Dropbox (Personal)\Spatial Models in Stata\Data\small.csv", replace
