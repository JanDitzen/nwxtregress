{smcl}
{hline}
{hi:help nwxtregress}{right: v. 0.01 - 03. August 2021}

{hline}
{title:Title}

{p 4 4}{cmdab:nwxtreg:ress} - Network Regressions in Stata with unbalanced panel data and time varying network structures or spatial weight matrices..{p_end}

{title:Syntax}

{p 4}{ul:Spatial Autocorrelation Model (SAR)}{p_end}

{p 4 13}{cmdab:nwxtreg:ress} {depvar} [{indepvars}]  [if] {cmd:,}
 {cmd:dvarlag(W1[,}{it:options1}{cmd:])}
 [{it:mcmcoptions}
 {cmd:nosparse}]
 {p_end}

{p 4}{ul:Spatial Durbin Model (SAR)}{p_end}

{p 4 13}{cmdab:nwxtreg:ress} {depvar} [{indepvars}]  [if] {cmd:,}
    {cmd:dvarlag(W1[,}{it:options1}{cmd:])}
    {cmd:ivarlag(W2[,}{it:options1}{cmd:])}
    [{it:mcmcoptions}
    {cmd:nosparse}]
{p_end}

{p 4 13}where {cmd:W1} and {cmd:W2} are spatial weight matrices. Default is {help Sp} object.
{cmd:dvarlag()} and {cmd:ivarlag()} define the spatial lag of the dependent
and independent variables. 
{cmd:ivarlag()} is repeatable with different spatial weights.{p_end}

{p 4}{ul: Maintenance}{p_end}

{p 4 13}{cmd:nwxtregress , [update version]}{p_end}

{p 8}{cmd:version} displays version. 
{cmd:update} updates {cmd:nwxtregress} from {browse "https://janditzen.github.io/nwxtregress/":GitHub}
using {help net install}.{p_end}


{p 4}{ul:General Options}{p_end}
{synoptset 20}{...}
{synopt:option}Description{p_end}
{synoptline}
{synopt:{it:nosparse}}do not convert weight matrix internally to a sparse matrix{p_end}
{synoptline}
{p2colreset}{...}

{p 4}{ul:{it:Options1}}{p_end}
{synoptset 20}{...}
{synopt:Options1}Description{p_end}
{synoptline}
{synopt:{it:mata}}declares weight matrix is mata matrix.{p_end}
{synopt:{it:sparse}}if weight matrix is sparse.{p_end}
{synopt:{it:timesparse}}weight matrix is sparse and varying over time.{p_end}
{synopt:{it:id(string)}}vector of IDs if W is a non sparse mata matrix{p_end}
{synoptline}
{p2colreset}{...}

{p 4}{ul:{it:mcmcoptions}}{p_end}
{synoptset 20}{...}
{synopt:mcmcoptions}Description{p_end}
{synoptline}
{synopt:{it:draws()}}number of griddy gibs draws, default 2000{p_end}
{synopt:{it:gridlength()}}grid length, default 1000{p_end}
{synopt:{it:nomit()}}number of omitted draws, default 500{p_end}
{synopt:{it:barrypace(numlist)}}settings for BarryPace Trick. Order is iterations, maxorder. Default is 50 and 100{p_end}
{synopt:{it:usebp}}use BarryPace trick instead of LUD for inverse of (I−ρW).{p_end}
{synopt:{it:seed(#)}}sets the seed{p_end}
{synoptline}
{p2colreset}{...}

{p 4 4}Data has to be xtset before using {cmd:nwxtregress}; see {help xtset}. {depvars} and {indepvars}  may contain time-series operators, see {help tsvarlist}.{p_end}

{title:Contents}
{p 4}{help nwxtregress##description:Description}{p_end}
{p 4}{help nwxtregress##options:Options}{p_end}
{p 4}{help nwxtregress##DIE:Direct, Indirect and Total Effects}{p_end}
{p 4}{help nwxtregress##saved_vales:Saved Values}{p_end}
{p 4}{help nwxtregress##examples:Examples}{p_end}
{p 4}{help nwxtregress##references:References}{p_end}
{p 4}{help nwxtregress##about:About, Authors and Version History}{p_end}

{marker description}{title:Description}

{p 4 4}{cmd:nwxtregress} estimates Spatial Autoregressive (SAR) or Spatial Durbin (SDM) models. The spatial weight matrices are allowed to be time varying and the dataset can be unbalanced.{p_end}

{p 4 4}The SAR is:{p_end}

{col 8}Y = rho W1 Y + beta X + eps

{p 4 4}The SDM is:{p_end}

{col 8}Y = rho W1 Y + beta X + gamma W2 X + eps

{p 4 4}where W1 and W2 are spatial weight matrices, Y the dependent and X the independent variables.{p_end}


{marker options}{...}
{title:Options}

{phang}
{opt mata} declares weight matrix is mata matrix. Default is to use a spatial weight matrix from the Sp environment.

{phang}
{opt sparse} if weight matrix is in sparse format. Sparse format implies that the first two column define the origin and the destination of the flow, the third column the value of the flow.

{phang}
{opt timesparse} weight matrix is sparse and varying over time. As sparse but first column includes the time period.

{phang}
{opt id(string)}  vector of IDs if W is a non sparse mata matrix.

{phang}
{opt nosparse} not convert weight matrix internally to a sparse matrix.

{phang}
{opt draws()} number of griddy gibs draws, default 2000.

{phang}
{opt gridlength()} grid length, default 1000.

{phang}
{opt nomit()} number of omitted draws, default 500.

{phang}
{opt barrypace(numlist)} settings for BarryPace Trick. Order is iterations, maxorder. Default is 50 and 100.

{phang}
{opt usebp} use BarryPace trick instead of LUD for inverse of (I−ρW).

{phang}
{opt seed(#)} sets the {help seed}.

{phang}
{opt version} display version.

{phang}
{opt update} update from Github.


{marker DIE}{...}
{title:Direct, indirect and total Effects}

{p 4 4}Direct, indirect and total effects can be calculated using {cmd:estat impact}.{p_end}

{marker saved_vales}{title:Saved Values}

{p 4}{cmd:nwxtregress test} stores the following in {cmd:e()}:{p_end}

{col 4} Matrices
{col 8}{cmd: e(b)}{col 27} Coefficient Matrix. 
{col 8}{cmd: e(V)}{col 27} Variance-Covariance Matrix. 

{col 4} Scalars
{col 8}{cmd: e(N)}{col 27} Number of observations
{col 8}{cmd: e(N_g)}{col 27}  Number of groups
{col 8}{cmd: e(T)}{col 27}  Number of time periods
{col 8}{cmd: e(Tmin)}{col 27}  Minimum number of time periods
{col 8}{cmd: e(Tavg)}{col 27} Average number of time periods
{col 8}{cmd: e(Tmax)}{col 27} Maximum number of time periods
{col 8}{cmd: e(K)}{col 27} Number of regressors excluding spatial lags
{col 8}{cmd: e(Kfull)}{col 27} Number of regressors including spatial lags
{col 8}{cmd: e(r2)}{col 27} R-squared
{col 8}{cmd: e(r2_a)}{col 27} adjusted R-squared
{col 8}{cmd: e(MCdraws)}{col 27} Number of MCMC draws

{col 4} Macros
{col 8}{cmd: e(sample)}{col 27} Sample

{marker examples}{...}
{title:Examples}

{p 4 4}An example dataset with USE/MAKE table data from the BEA’s website and links between industries is available {browse "https://github.com/JanDitzen/nwxtregress/tree/main/examples":GitHub}. 
The dataset {it:W.dta} contains the linkages (spatial weights) and the dataset {it:SAM.dta} the firm data. 
We want to estimate capital consumption by using compensation and net surplus as explanatory variables.{p_end}

{p 4 4}First we load the data from the W dataset and convert into a SP object for the year 1998.{p_end}

{col 8}{stata use W}
{col 8}{stata keep if Year == 1998}
{col 8}{stata spmatrix fromdata W = sam* ,replace}

{p 4 4}Next, we load the dataset with the firm data and estimate a SAR with a time constant spatial weight matrix. 
We also obtain the total, direct and indirect effects using {cmd:estat impact}. 
For reproducibility we set a {it:seed}.{p_end}

{col 8}{stata use SAM, clear}
{col 8}{stata nwxtregress cap_cons compensation net_surplus , dvarlag(W) seed(1234)}
{col 8}{stata estat impact}

{p 4 4}The disadvantage is that the spatial weight are constant across time. 
To allow for time varying spatial weights, we load the W dataset again and save it into mata:{p_end}

{col 8}{stata use W}
{col 8}{stata putmata W = (ID1 ID2 sam)}

{p 4 4}Using the SAM dataset again, we can estimate the SAR model with time varying spatial weights. 
To do so we use the option timesparse and mata.{p_end}

{col 8}{stata nwxtregress cap_cons compensation net_surplus , dvarlag(W,mata timesparse) seed(1234)}

{p 4 4}Finally, we want to estimate an SDM by adding the option {cmd:ivarlag()}:{p_end}

{col 8}{stata "nwxtregress cap_cons compensation net_surplus , dvarlag(W,mata timesparse) ivarlag(W: compensation,mata timesparse )  seed(1234)"}

{marker references}{title:References}


{marker about}{title:How to install}

{p 4 4}The latest version of the nwxtregress package can be obtained by typing in Stata:{p_end}

{col 8}{stata "net from https://janditzen.github.io/nwxtregress/"}


{title:Acknowledgments}

{p 4 4}We are grateful to James LeSage for guidance and for sharing code that was instrumental to our ability to create the {cmd:nwxtreg} command. 
This article also benefited from comments from Zack Liu, Ioannis Spyridopoulos, and Adam Winegar. 
All remaining errors are our own.{p_end}

{p 4 4}Jan Ditzen acknowledges financial support from Italian Ministry MIUR under 
the PRIN project Hi-Di NET - Econometric Analysis of High Dimensional Models
 with Network Structures in Macroeconomics and Finance (grant 2017TA7TYC).{p_end}


{title:Authors}

{p 4}Jan Ditzen (Free University of Bozen-Bolzano){p_end}
{p 4}Email: jan.ditzen@unibz.it{p_end}
{p 4}Web: {browse www.jan.ditzen.net}{p_end}

{p 4}William Grieser (Texas Christian University){p_end}
{p 4}Email: w.grieser@tcu.edu{p_end}
{p 4}Web: {browse www.williamgrieser.com/}{p_end}

{p 4}Morad Zekhnini (Michigan State University){p_end}
{p 4}Email: zekhnini@msu.edu{p_end}
{p 4}Web: {browse "https://sites.google.com/view/moradzekhnini/home"}{p_end}

{title:Changelog}
{p 4 4}{ul:Version 0.01}{p_end}
