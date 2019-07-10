#! /usr/bin/Rscript
library(tidyverse)

#args <- commandArgs(trailingOnly = TRUE)
args<-c("o-tshaw")
load(paste("./outputs/104/", args[1], ".rda", sep=""))

