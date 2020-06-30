# Pre-SISMID computational instructions 

Some of the exercises in this course will be computational, all using R https://www.r-project.org/. Even if you already have R installed, it may be worth checking if an updated version is available.

SISMID recommends all participants who are not familiar with R to work through some R tutorials in advance of the course. For example, https://cran.r-project.org/doc/contrib/Paradis-rdebuts_en.pdf. We also recommend the use of R Studio https://rstudio.com/, a free IDE (integrated development environment) for R, though this is optional. 

## devtools
To save time during the course, we suggest you install devtools in advance, if you do not already have it. Instructions are available here:
https://www.r-project.org/nosvn/pandoc/devtools.html
under the heading “Updating to the latest version of devtools”. Note that this involves downloading RTools for Windows, Xcode for Mac or a compiler and dev libraries for Linux. 

## R packages
We will be using several R packages in this module, you can install them in advance by entering the following commands in the R console: 

 ```{r }
devtools::install_github('xavierdidelot/TransPhylo')
install.packages("outbreaker2")
install.packages("adegenet", dep = TRUE)
install.packages("phangorn", dep = TRUE)
```
