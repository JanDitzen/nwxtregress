clear all

use https://janditzen.github.io/nwxtregress/examples/IO.dta
keep if Year == 1998
replace sam = 0 if sam < 0
replace sam = 0 if ID1==ID2
keep ID1 ID2 sam
reshape wide sam, i(ID1) j(ID2)
spset ID1
spmatrix fromdata WSpmat = sam* , replace


use https://janditzen.github.io/nwxtregress/examples/VA.dta
nwxtregress cap_cons compensation net_surplus , dvarlag(WSpmat) seed(1234)
estat impact

frame create IO
frame IO: use https://janditzen.github.io/nwxtregress/examples/IO.dta

nwxtregress cap_cons compensation net_surplus ,dvarlag(sam, frame(IO) id(Year ID1 ID2)  timesparse) seed(1234) 

frame IO: putmata Wt = (Year ID1 ID2 sam), replace
nwxtregress cap_cons compensation net_surplus , dvarlag(Wt, mata timesparse) seed(1234)

nwxtregress cap_cons compensation net_surplus , dvarlag(Wt,mata timesparse) ivarlag(Wt: compensation,mata timesparse )  seed(1234)


mata: Wt2 = Wt[selectindex(Wt[.,4]:>2601.996),.]
nwxtregress cap_cons compensation net_surplus , dvarlag(Wt, mata timesparse) ivarlag(Wt: net_surplus, mata timesparse) ivarlag(Wt2: compensation, mata timesparse) seed(1234)

estat impact


predict xb
predict residuals, residual

