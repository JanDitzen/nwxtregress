clear

cd "C:\Users\jan\Dropbox (Personal)\Spatial Models in Stata\Stata\test\Stata USM"

use homicide_1960_1990

xtset _ID year

spset

spmatrix create contiguity W if year == 1990

spxtregress hrate ln_population ln_pdensity gini i.year, re dvarlag(W) errorlag(W)
