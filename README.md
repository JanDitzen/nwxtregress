# nwxtregress
Network Regressions in Stata with unbalanced panel data and time varying network structures or spatial weight matrices.

Current Version: 0.13 (17.01.2023)

[nwxtregress](https://github.com/janditzen/nwxtregress) | ![version](https://img.shields.io/github/v/release/janditzen/nwxtregress) | ![release](https://img.shields.io/github/release-date/janditzen/nwxtregress) | |

__Table of Contents__
1. [Syntax](#1-syntax)
2. [Description](#2-description)
3. [Options](#3-options)
4. [Postestimation (Predict, Direct, indirect and total effects)](#4-postestimation)
5. [Saved Values](#5-saved-values)
6. [Examples](#6-examples)
7. [References](#7-references)
8. [How to install](#8-how-to-install)
9. [Questions?](#9-questions?)
10. [About](#10-authors)


# 1. Syntax

### SAR

```
nwxtregress depvar  indepvars [if], 
        ivarlag(W1[, sparse timesparse mata id(string)])
        [mcmcoptions absorb(varlist, keepsingeltons) transform(transfrom_varlist, transform_options)]
```

### SDM

```
nwxtregress depvar  indepvars [if], 
        ivarlag(W1[, sparse timesparse mata id(string)])
        dvarlag(Ws:varlist[, sparse timesparse mata id(string)]
        [mcmcoptions absorb(varlist, keepsingeltons) transform(transfrom_varlist, transform_options)]
```

Data has to be ```xtset``` before use. W1 and Ws define the spatial weight matrix, default is ***Sp*** object.
```dvarlag()``` and ```ivarlag()``` define the spatial lag of the dependent and independent variables.
 ```ivarlag()``` is repeatable and multiple spatial weight matrices are supported.

```nwxtregress``` requires Stata 14.2 or higher. 
```python``` and ```frame``` can only be used with Stata 16 or higher.

#### Options for ```ivarlag()``` and ```dvarlag()```

option | Description
--- | ---
**mata** | declares weight matrix is ```mata``` matrix.
**sparse** | if weight matrix is sparse.
**timesparse** | weight matrix is sparse and varying over time.
**id(string)** | vector of IDs if W is a non sparse mata matrix
**normalize(string)** | which normalization to use.
**zero(real)** | how to treat zeros in spatial weight matrix.

#### General Options

options | Description
--- | ---
**nosparse** | not convert weight matrix internally to a sparse matrix
**asarray(name)** | change name of array estimation results and info
**standardize** | standardizes all variables, short for ***transform(_all, by(idvar))***

#### MCMC Options

mcmcoptions | Description
--- | ---
**draws()** | number of griddy gibs draws, default 2000
**gridlength()** | grid length, default 1000
**nomit()** | number of omitted draws, default 500
**barrypace(numlist)** | settings for BarryPace Trick. Order is iterations, maxorder. Default is 50 and 100
**usebp** | use BarryPace trick instead of LUD for inverse of (I−ρW).
**python** | use Python to calculate LUD or Barry Pace trick.
**seed(#)** | sets the seed

#### Transform Options

transformoptions | Description
--- | ---
**transform_varlist** | variables to be transformed. ***_all*** transformes all dependent and independent variables. If not specified, ***cmd:_all*** assumed.
**by(varname)** | variable defining level of transformation.
**after** | transform variables after spatial lags are calcuted
**wy** | transform spatial lag of dependent variable
**wx** | transform spatial lag of independent variables as defined by ***varlist***
**nom:ean** | do not demean data
**nosd** | do not standardize data (standard deviation of 1)

#### Maintenance:

```
nwxtregress , [update version]
```

***nwxtregress, version*** displays the current version.
***nwxtregress, update*** updates ***nwxtregress*** from GitHub.


# 2. Description

```nwxtregress``` estimates Spatial Autoregressive (SAR) or Spatial Durbin (SDM) models. The spatial weight matrices are allowed to be time varying and the dataset can be unbalanced.

The SAR is:

```
y = rho W1 y + beta X + eps

```

The SDM is:

```
Y = rho W1 Y + beta X + gamma W2 X + eps

```

where **W1** and **W2** are spatial weight matrices, Y the dependent and X the independent variables.

```nwxtregress``` can handle spatial weights in three formats: 1. square matrix, 2. sparse and 3. time sparse.
Sparse matrices have the advantage that they save space and thus computational time
and allow for time varying weights.
The [Sp environment](http://www.stata.com/manuals/sp.pdf) only supports the square matrix format. 
```nwxtregress``` can read **square**, **sparse** and **time sparse** formats if the 
data for the weights is in ``mata`` or saved in a ``frame``.{p_end}

#### 1. Square matrix format

The spatial weights are a matrix with dimension N_g x N_g. It is time constant. An Example with a 5 x 5 matrix is:

        0    0.1  0.2  0
        0    0    0.1  0.2
        0.3  0.1  0    0
        0.2  0    0.2  0

#### 2. Sparse matrix format

The sparse matrix format is a **v x 3** matrix, where **v** is the number of non-zero elements in the spatial weight matrix. The weight matrix is time constant. The first column indicates the destination, the second the origin of the flow.
A sparse matrix of the matrix from above is:

           Destination  Origin    Flow
           1            2         0.1
           1            3         0.2
           2            3         0.1
           2            4         0.2
           3            1         0.3
           3            2         0.1
           4            1         0.2
           4            3         0.2

#### 3. Time-Sparse format

The time sparse format can handle time varying spatial weights. The first column indicates the time period, the remaining are the same as for the sparse matrix. 
For example, if there are two time periods and we have the matrix from 
above for the first and the square for the second period:


           Time    Destination    Origin    Flow
           1       1              2         0.1
           1       1              3         0.2
           1       2              3         0.1
           1       2              4         0.2
           1       3              1         0.3
           1       3              2         0.1
           1       4              1         0.2
           1       4              3         0.2
              (next time period)
           2       1              2         0.1
           2       1              3         0.4
           2       2              3         0.1
           2       2              4         0.4
           2       3              1         0.9
           2       3              2         0.1
           2       4              1         0.4
           2       4              3         0.4

Internally, nextregress will always use the time sparse format. This ensures that unbalanced panels do not pose a problem.  nextregress comes with functions for creating sparse matrices, coplying a sparse matrix into a squared format, and functions for mathematical operations (transpose and multiplication).


# 3. Options

#### Options

Option | Description
 --- | --- 
**frame(name)** | declares weight matrix is saved in a ```frame```. Default is to use a spatial weight matrix from the **Sp** environment. If a frame is used, data can be in sparse, timesparse or square matrix format.
**mata** | declares weight matrix is ```mata``` matrix. Default is to use a spatial weight matrix from the **Sp** environment. If a mata matrix is used, data can be in sparse, time sparse or square matrix format.
**sparse** | if weight matrix is in sparse format. Sparse format implies that the first two column define the origin and the destination of the flow, the third column the value of the flow.
**timesparse** | weight matrix is sparse and varying over time. As **sparse** but first column includes the time period.
**id(string)** | vector of IDs if W is a non sparse mata matrix. If a frame is used, then **id()** contains the varible names of the time indicator (if applicable), the origin and destination of the flows.
**normalize(string)** |  which normalization to use for spatial weight matrix.  Default is row normalisation.  Can be none, row (default), column, spectral or minmax, see normalisation option of [spmat creat](http://www.stata.com/manuals/spspmatrixcreate.pdf). The normalisation is done for each time period individually.
**zero(real)** | defines how to treat zeros in spatial weight matrics.  Default is to remove zero entries for non-sparse matrices and to set zeros to 0.0001 if weight matrix is (time)sparse.
**nosparse** | not convert weight matrix internally to a sparse matrix. Option is not recommended to use.
**asarray(name)** | nwxtregress saves intermediate results such as the spatial weight matrix in an internal time sparse format, residuals and results from the MCMC in an array, see stored values.  It is not recommended to change contents of the array and the option to change the name should only be rarely used. The default name is NWXTREG_OBJECT#, where # is a counter if the array already existed.
**draws()** | number of griddy gibs draws, default 2000.
**gridlength()** | grid length, default 1000.
**nomit()** | number of omitted draws, default 500.
**barrypace(numlist)** | settings for BarryPace Trick. Order is iterations, maxorder. Default is 50 and 100.
**usebp** | use BarryPace trick instead of LUD for inverse of (I−ρW).
**python** | use Python to calculate the LU Decomposition or BarryPace trick.  Requires installation of Python, scipy, sfi and numpy. Using Python to calculate the LUD is faster by a factor 4-10.
**seed(#)** | sets the seed.
**version** | display version.
**update** | update from Github.

## 3.1 High Dimensional Fixed Effects

***nwxtregress*** can remove high dimensional fixed effects using [reghdfe](http://scorreia.com/software/reghdfe/). The fixed effects are partialled out before spatial lags are cacluated. Constant is automatically removed when ```cmd:absorb()``` is used. The syntax is:

``
asorb(varlist, keepsingeltons)
``

Option | Description
 --- | --- 
varlist | categorical variables that identify the fixed effects to be absorbed. 
veepsingelton | keep singelton units.

## 3.2 Transform Data

```nwxtregress``` can demean and standardize dependent and independent variables, before or after the calculation of the spatial lags. Spatial lags can be transformed as well. The syntax is:

``
transform([varlist] [, by(varname)) after nomean nosd wy wx])
``

Option | Description
 --- | --- 
varlist | variables to be transformed. ```_all``` implies all dependent and independent variabkes. If left empty, ```_all``` assumed.
by(varname) | variable defining transformation. Default is ```by(ID)```, where ID identifies the cross-sections. ```by(_all)``` transforms data across all cross-sections.
after | transform data after caculation of spatial lags. Default is to transform data first.
nomean | do not demean data.
nosd | do not standardize data.
wy | transform spatial lag of dependent variable. Implies ```after```.
wx | transform spatial lags of independent variables as defined in ```it:varlist```. Implies ```cmd:after```.
transform | short for transform(_all).
standardize | short for transform(_all).

# 4. Postestimation

## 4.1 Direct, indirect and total effects.

Direct, indirect and total effects. can be calculated using ```estat impact```. The syntax is

```
 estat impact [varlist] [, options]
```

Option | Description
 --- | --- 
seed(#) | set seed for Barry Pace matrix inversion.
array(name) | name of array with saved contents from nwxtregress, see stored results.

``varlist`` defines the variables for which the direct, indirect and total effects are displayed.  If not specified, then estat impact will calculate the effects for all explanatory variables (indepvars).

``estat impact`` saves the following in r():



Matrix | Description
 --- | --- 
**r(b_direct)** | Coefficient Matrix of direct effects
**r(V_direct)** | Variance covariance matrix of direct effects
**r(b_indirect)** | Coefficient Matrix of indirect effects
**r(V_indirect)** |  Variance covariance matrix of indirect effects 
**r(b_total)** | Coefficient Matrix of total effects
**r(V_total)** | Variance covariance matrix of total effects

## 4.2 Predict

``predict`` can be used after nwxtregress. The syntax for predict is:

```
predict [type] varname [, options]
```

Option | Description
 --- | --- 
xb | calculate linear prediction.
res | calculate residuals.
replace | replace if varname exists.
array(name) | name of array with saved contents from nwxtregress, see stored results.


# 5. Saved Values

***nwxtregress*** saves the following in ***e()***

#### Matrices

Matrices | Description
---|---
***b*** | Coefficient Matrix
***V*** | Variance-Covariance Matrix

#### Scalars

Scalars | Description
---|---
N | Number of observations
N_g | Number of groups
T | Number of time periods
Tmin | Minimum number of time periods
Tavg | Average number of time periods
Tmax | Maximum number of time periods
K | Number of regressors excluding spatial lags
Kfull | Number of regressors including spatial lags
r2 | R-squared
r2_a | adjusted R-squared
MCdraws | Number of MCMC draws

#### Macros

Macro | Description
---|---
sample | sample


#### mata arrays

 In addition to e() and r() nwxtregress saves informations about the estimation in a mata array.  The contents are the weight matrix in time sparse format, residuals and results from the MCMC.  Storing those saves time for ``estat impact`` and ``predict``.  The name default name of the array is ```_NWXTREG_OBJECT#```, but can be set with the option **asarray()**.  In general it is not recommended to change this setting.


# 6. Examples

An example dataset with USE/MAKE table data from the BEA’s website and links between industries is available [GitHub](https://github.com/JanDitzen/nwxtregress/tree/main/examples). The dataset IO.dta contains the linkages (spatial weights) and the dataset VA.dta the firm data.  We want to estimate capital consumption by using compensation and net surplus as explanatory variables.


First we load the data from the W dataset and convert into a **SP** object for the year 1998.

```
use https://janditzen.github.io/nwxtregress/examples/IO.dta
keep if Year == 1998
replace sam = 0 if sam < 0
replace sam = 0 if ID1==ID2
keep ID1 ID2 sam
reshape wide sam, i(ID1) j(ID2)
spset ID1
spmatrix fromdata WSpmat = sam* , replace
```

Next, we load the dataset with the firm data and estimate a SAR with a time constant spatial weight matrix.  We also obtain the total, direct and indirect effects using estat impact.  For reproducibility we set a seed.

```
use https://janditzen.github.io/nwxtregress/examples/VA.dta
nwxtregress cap_cons compensation net_surplus , dvarlag(WSpmat) seed(1234)
estat impact

```

The disadvantage is that the spatial weight are constant across time and we had to get rid of all negative numbers.  To allow for time varying spatial weights, we load the W dataset again and but load it into the frame IO:

```
frame create IO
frame IO: use https://janditzen.github.io/nwxtregress/examples/IO.dta
```

Using the VA dataset again, we can estimate the SAR model with time varying spatial weights.  To do so we use the options frame(name), where name indicates the frame and the weight matrix name corresponds to the variable names.  The data is in timesparse format so we need to use the option timesparse.  Finally it is nessary to define the year identifier and the origin and destination of the flows using the id() option:


```
nwxtregress cap_cons compensation net_surplus ,
dvarlag(sam, frame(IO) id(Year ID1 ID2)  timesparse) 
seed(1234) 
```

Alternatively we can load the spatial weight matrix into mata:

```
frame IO: putmata Wt = (Year ID1 ID2 sam), replace
nwxtregress cap_cons compensation net_surplus , 
dvarlag(Wt, mata timesparse) seed(1234)
```

If we want to estimate an SDM by adding the option ivarlag():

```
nwxtregress cap_cons compensation net_surplus , 
dvarlag(Wt,mata timesparse) ivarlag(Wt: compensation,mata timesparse )  
seed(1234)
```

Use Python (requires Stata 16 or later) to improve speed of calculating the LUD:

```
nwxtregress cap_cons compensation net_surplus , dvarlag(Wt,mata timesparse) ivarlag(Wt: compensation,mata timesparse )  seed(1234) python
```

Transform data by demeaning and standardising it:

```
nwxtregress cap_cons compensation net_surplus , dvarlag(Wt,mata timesparse) ivarlag(Wt: compensation,mata timesparse )  seed(1234) transform(_all, by(ID)
```

or

```
nwxtregress cap_cons compensation net_surplus , dvarlag(Wt,mata timesparse) ivarlag(Wt: compensation,mata timesparse )  seed(1234) standardize
```

Partial out firm and year fixed effects (requires reghdfe):

```
nwxtregress cap_cons compensation net_surplus , dvarlag(Wt,mata timesparse) ivarlag(Wt: compensation,mata timesparse )  seed(1234) absorb(ID Year)
```

We can also define two different spatial weight matrices:

```
mata: Wt2 = Wt[selectindex(Wt[.,4]:>2601.996),.]
nwxtregress cap_cons compensation net_surplus , 
dvarlag(Wt, mata timesparse) 
ivarlag(Wt: net_surplus, mata timesparse) 
ivarlag(Wt2: compensation, mata timesparse) seed(1234)
```

Total, direct and indirect effects can be calculated using estat impact:

```
estat impact
```

To predict fitted values and residuals predict can be used:

```
predict xb
predict residuals, residual
```

# 7. References

# 8. How to install


The latest version of the ***nwxtregress*** package can be obtained by typing in Stata:

```
net from https://janditzen.github.io/nwxtregress/
``` 

or 

```
net install nwxtregress , from(https://janditzen.github.io/nwxtregress/)

```

# 9. Questions?

Questions? Feel free to write us an email, open an [issue](https://github.com/JanDitzen/nwxtregress/issues) or [start a discussion](https://github.com/JanDitzen/nwxtregress/discussions).


# 10. Authors

#### Jan Ditzen (Free University of Bozen-Bolzano)

Email: jan.ditzen@unibz.it

Web: www.jan.ditzen.net

#### William Grieser (Texas Christian University)

Email: w.grieser@tcu.edu

Web: https://www.williamgrieser.com/

#### Morad Zekhnini (Michigan State University)

Email: zekhnini@msu.edu

Web: https://sites.google.com/view/moradzekhnini/home

## About:
This version 0.13 as of 17.01.2023

## Changelog:
Version 0.13
- added options absorb() and transform()
- bugfixes when using fixed effects
- Python support for BarryPace trick
Version 0.12
- Support for Python to calculated LUD
Version 0.03 (alpha)
- Bugs in sparse matrix multiplication and return if non sparse matrix is used fixed.
