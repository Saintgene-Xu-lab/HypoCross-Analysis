library(clustree)
library(ape)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(qs)


meta <- qread("meta.qs")
p <- clustree(meta, prefix = "L")
qsave(p,file = "/tree.plot.qs")
p1 <- clustree(meta, prefix = "L",return = "layout") # layout graph
qsave(p1,file = "tree.layout.qs")
p2 <- clustree(meta, prefix = "L",return = "graph") # layout graph
qsave(p2,file = "tree.graph.qs")

rm(list=ls())

