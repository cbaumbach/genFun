Installation
============
(tested under Ubuntu 16.04)

Install "devtools" package and its dependencies
-----------------------------------------------

At the shell prompt:

```shell
sudo apt-get install libcurl4-openssl-dev libssl-dev
```

At the R prompt:

```R
install.packages("devtools")
```

Install "genFun" package and its dependencies "miscFun" and "treeFun"
---------------------------------------------------------------------

At the R prompt:

```R
devtools::install_github("cbaumbach/miscFun")
devtools::install_github("cbaumbach/treeFun")
devtools::install_github("cbaumbach/genFun")
```
