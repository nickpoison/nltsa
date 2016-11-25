# first, open R in this directory
rm(list=ls())
source('sourceDir.R')
sourceDir('R')
  # require(testthat)
  # sourceDir('inst/statistical_tests')  # don't think I need these tests??
  # sourceDir('inst/tests')
# get the data 
 files=list.files('data')
 setwd('data')
 for (i in 1:length(files)) load(files[i])
 setwd('..')
#
rm(i)
rm(files) 
rm(sourceDir)
#rm(exampleCPS)  # don't need this anymore
#rm(arch)        # don't need this ??? not sure
# version number
nltsa.version=Sys.Date()
#
save(list=ls(), file='nltsa.rda')
#clean
rm(list=ls())
#
cat("now do this: load('nltsa.rda')", "\n")


