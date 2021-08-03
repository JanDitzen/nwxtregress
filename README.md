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

Data has to be ```xtset``` before use. W1 and Ws define the spatial weight matrix, default is ***Sp*** object. ```dvarlag()``` is repeatable and can have multiple spatial weight matrices.

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
***nwxtregress, update*** updates ***xtbreak*** from GitHub.


# 2. Description

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

 # 4. Saved Values

 # 5. Examples

 # 6. References

# 7. How to install


# 8. Questions?

Questions? Feel free to write us an email, open an [issue](https://github.com/JanDitzen/xtbreak/issues) or [start a discussion](https://github.com/JanDitzen/xtbreak/discussions).


# 9. Authors

#### Jan Ditzen (Free University of Bozen-Bolzano)

Email: jan.ditzen@unibz.it

Web: www.jan.ditzen.net

#### William Grieser (Texas Christian University)

Email: w.grieser@tcu.edu

Web: https://www.williamgrieser.com/

#### Morad Zekhnini (Michigan State University)

Email: zekhnini@msu.edu

Web: https://sites.google.com/view/moradzekhnini/home
