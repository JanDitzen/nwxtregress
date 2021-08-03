# nwxtregress
Network Regressions in Stata with unbalanced panel data and time varying network structures or spatial weight matrices.

__Table of Contents__
1. [Syntax](#1-syntax)
2. [Description](#2-description)
3. [Options](#3-options)
4. [Direct and Indirect Effects](#4-direct-and-indirect-effects)
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
        [mcmcoptions nosparse]
```

### SDM

```
nwxtregress depvar  indepvars [if], 
        ivarlag(W1[, sparse timesparse mata id(string)])
        dvarlag(Ws:varlist[, sparse timesparse mata id(string)]
        [mcmcoptions nosparse]
```

Data has to be ```xtset``` before use. W1 and Ws define the spatial weight matrix, default is ***Sp*** object.
```dvarlag()``` and ```ivarlag()``` define the spatial lag of the dependent and independent variables.
 ```dvarlag()``` is repeatable and multiple spatial weight matrices are supported.

#### Options for ```ivarlag()``` and ```dvarlag()```

option | Description
--- | ---
**mata** | declares weight matrix is ```mata``` matrix.
**sparse** | if weight matrix is sparse.
**timesparse** | weight matrix is sparse and varying over time.
**id(string)** | vector of IDs if W is a non sparse mata matrix

#### General Options

options | Description
--- | ---
**nosparse** | not convert weight matrix internally to a sparse matrix

mcmcoptions | Description
--- | ---
**draws()** | number of griddy gibs draws, default 2000
**gridlength()** | grid length, default 1000
**nomit()** | number of omitted draws, default 500
**barrypace(numlist)** | settings for BarryPace Trick. Order is iterations, maxorder. Default is 50 and 100
**usebp** | use BarryPace trick instead of LUD for inverse of (I−ρW).
**seed(#)** | sets the seed

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

# 3. Options

#### Options

Option | Description
 --- | --- 
**mata** | declares weight matrix is ```mata``` matrix. Default is to use a spatial weight matrix from the **Sp** environment.
**sparse** | if weight matrix is in sparse format. Sparse format implies that the first two column define the origin and the destination of the flow, the third column the value of the flow.
**timesparse** | weight matrix is sparse and varying over time. As **sparse** but first column includes the time period.
**id(string)** | vector of IDs if W is a non sparse mata matrix.
**nosparse** | not convert weight matrix internally to a sparse matrix.
**draws()** | number of griddy gibs draws, default 2000.
**gridlength()** | grid length, default 1000.
**nomit()** | number of omitted draws, default 500.
**barrypace(numlist)** | settings for BarryPace Trick. Order is iterations, maxorder. Default is 50 and 100.
**usebp** | use BarryPace trick instead of LUD for inverse of (I−ρW).
**seed(#)** | sets the seed.
**version** | display version.
**update** | update from Github.

 # 4. Direct and Indirect Effects

Direct and indirect effects can be calculated using ```estat impact```.

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

# 6. Examples

An example dataset with USE/MAKE table data from the BEA’s website and links between industries is available [GitHub](https://github.com/JanDitzen/nwxtregress/tree/main/examples). The dataset W.dta contains the linkages (spatial weights) and the dataset SAM.dta the firm data. We want to estimate capital consumption by using compensation and net surplus as explanatory variables.

First we load the data from the W dataset and convert into a **SP** object for the year 1998.

```
use W
keep if Year == 1998
spmatrix fromdata W = sam* ,replace
```

Next, we load the dataset with the firm data and estimate a SAR with a time constant spatial weight matrix. We also obtain the total, direct and indirect effects using ```estat impact```. For reproducibility we set a seed.

```
use SAM, clear
nwxtregress cap_cons compensation net_surplus , dvarlag(W) seed(1234)
estat impact

```

The disadvantage is that the spatial weight are constant across time. To allow for time varying spatial weights, we load the W dataset again and save it into mata:

```
use W
putmata W = (ID1 ID2 sam)

```

Using the SAM dataset again, we can estimate the SAR model with time varying spatial weights. To do so we use the option ***timesparse*** and ***mata***.

```
nwxtregress cap_cons compensation net_surplus , 
dvarlag(W,mata timesparse) seed(1234)

```

Finally, we want to estimate an SDM by adding the option ***ivarlag()***

```
nwxtregress cap_cons compensation net_surplus , 
dvarlag(W,mata timesparse) 
ivarlag(W: compensation,mata timesparse )  seed(1234)

```

# 7. References

# 8. How to install


The latest version of the ***nwxtregress*** package can be obtained by typing in Stata:

```
net from https://janditzen.github.io/nwxtregress/
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
