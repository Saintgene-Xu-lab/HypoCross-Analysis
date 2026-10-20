############### 
library(clustree)
library(tidygraph)
library(ape)
library(tidyverse)
library(reshape2)
library(stringr)
library(ggplot2)
path <- "/data/analysis/hy_ION2/data/neu_data/sub_neuron"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  readRDS(x)})

data2 <- lapply(data,function(x){
  x <- x[,c("anno",paste0("anno",c(2:7)),paste0("integrated_snn_res.",seq(0.01,0.09,0.01)),paste0("integrated_snn_res.",seq(0.1,1,0.1)))]
})

re <- c(paste0("integrated_snn_res.",seq(0.1,1,0.1)),paste0("integrated_snn_res.",seq(0.01,0.09,0.01)))
for (i in 1:length(data2)){
  for (j in 1:length(re)){
    data2[[i]][,re[j]] <- paste0(data2[[i]]$anno6,"_",as.character(data2[[i]][,re[j]]))
    
  }
  data2[[i]]$id <- rownames(data2[[i]])
}

data_all <- do.call(rbind,data2)
rownames(data_all) <- data_all$id

# 0.05 128
# 0.1 215
data_all2 <- data_all[,c("integrated_snn_res.0.05","integrated_snn_res.0.1")]
p <- clustree(data_all2, prefix = "integrated_snn_res.") # layout graph
p
### 
data_all_set <- data_all
data_all_set$integrated_snn_res.0.05_m <- as.character(data_all_set$integrated_snn_res.0.05)
length(table(data_all_set$integrated_snn_res.0.05_m))
if(T){
  data_all_set[data_all_set$integrated_snn_res.0.05 %in% c("IN9_1"),"integrated_snn_res.0.05_m"] <- "IN9_0"
  data_all_set[data_all_set$integrated_snn_res.0.05 %in% c("IN8_3"),"integrated_snn_res.0.05_m"] <- "IN8_0"
  data_all_set[data_all_set$integrated_snn_res.0.05 %in% c("IN31_1"),"integrated_snn_res.0.05_m"] <- "IN31_0"
  data_all_set[data_all_set$integrated_snn_res.0.05 %in% c("IN29_1"),"integrated_snn_res.0.05_m"] <- "IN29_0"
  data_all_set[data_all_set$integrated_snn_res.0.05 %in% c("IN2_1"),"integrated_snn_res.0.05_m"] <- "IN2_0"
  data_all_set[data_all_set$integrated_snn_res.0.05 %in% c("IN17_1"),"integrated_snn_res.0.05_m"] <- "IN17_0"
  data_all_set[data_all_set$integrated_snn_res.0.05 %in% c("EN9_1"),"integrated_snn_res.0.05_m"] <- "EN9_0"
  data_all_set[data_all_set$integrated_snn_res.0.05 %in% c("EN4_2"),"integrated_snn_res.0.05_m"] <- "EN4_0"
  data_all_set[data_all_set$integrated_snn_res.0.05 %in% c("EN3_3"),"integrated_snn_res.0.05_m"] <- "EN3_0"
  data_all_set[data_all_set$integrated_snn_res.0.05 %in% c("EN29_1"),"integrated_snn_res.0.05_m"] <- "EN29_0"
  data_all_set[data_all_set$integrated_snn_res.0.05 %in% c("EN24_1"),"integrated_snn_res.0.05_m"] <- "EN24_0"
  data_all_set[data_all_set$integrated_snn_res.0.05 %in% c("EN21_1"),"integrated_snn_res.0.05_m"] <- "EN21_0"
  data_all_set[data_all_set$integrated_snn_res.0.05 %in% c("EN18_1"),"integrated_snn_res.0.05_m"] <- "EN18_0"
  data_all_set[data_all_set$integrated_snn_res.0.05 %in% c("EN16_1"),"integrated_snn_res.0.05_m"] <- "EN16_0"
  data_all_set[data_all_set$integrated_snn_res.0.05 %in% c("EN11_1"),"integrated_snn_res.0.05_m"] <- "EN11_0"
}

length(table(data_all_set$integrated_snn_res.0.05_m))

data_all_set$integrated_snn_res.0.05 <- data_all_set$integrated_snn_res.0.05_m
data_all2 <- data_all_set[,c("integrated_snn_res.0.05","integrated_snn_res.0.1")]
p2 <- clustree(data_all2, prefix = "integrated_snn_res.") # layout graph

meta3 <- readRDS("")

other_cell <- meta3[setdiff(rownames(meta3),rownames(data_all_set)),]

other_cell$anno7 <- as.character(other_cell$anno7)
other_cell$anno6 <- as.character(other_cell$anno6)
other_cell$anno5 <- as.character(other_cell$anno5)
other_cell$anno4 <- as.character(other_cell$anno4)
other_cell$anno3 <- as.character(other_cell$anno3)
other_cell$anno2 <- as.character(other_cell$anno2)
other_cell$integrated_snn_res.0.05 <- paste0(other_cell$anno6,"_","0")
other_cell$integrated_snn_res.0.1 <- paste0(other_cell$anno6,"_","0")



data_all2 <- data_all_set[,c("anno",paste0("anno",c(2:7)),"integrated_snn_res.0.05","integrated_snn_res.0.1")]

neu_all <- data_all2
neu_all <- rbind(neu_all,other_cell[,c("anno",paste0("anno",c(2:7)),"integrated_snn_res.0.05","integrated_snn_res.0.1")])

#### transfer to levels 

neu_all$L1 <- "Neuron"
neu_all$L2 <- neu_all$anno3
neu_all$L3 <- NULL
neu_all$L4 <- neu_all$anno6
neu_all$L5 <- paste0("L5_",neu_all$integrated_snn_res.0.05)
neu_all$L6 <- paste0("L6_",neu_all$integrated_snn_res.0.1)

## L3
if(T){
  meta <- readRDS("")
  meta$anno6 <- meta3$anno6
  meta$anno9 <- as.character(meta$anno6)
  meta[meta$anno6 %in% c("CHAT"),]$anno9 <- c("9")
  meta[meta$anno6 %in% c("EN13","EN16"),]$anno9 <- c("7")
  meta[meta$anno6 %in% c("EN17"),]$anno9 <- c("7")
  meta[meta$anno6 %in% c("IN13"),]$anno9 <- c("7_1")
  meta[meta$anno6 %in% c("EN4"),]$anno9 <- c("6")
  meta[meta$anno6 %in% c("EN12","EN26"),]$anno9 <- c("5")
  meta[meta$anno6 %in% c("EN14","EN6","EN8"),]$anno9 <- c("4")
  
  meta[meta$anno6 %in% c("IN4","IN27","IN22","IN14"),]$anno9 <- c("3")
  meta[meta$anno6 %in% c("EN10","EN2","EN21","EN7"),]$anno9 <- c("2")
  
  meta[meta$anno6 %in% c("IN26"),]$anno9 <- c("14")
  meta[meta$anno6 %in% c("EN22"),]$anno9 <- c("13")
  meta[meta$anno6 %in% c("EN19"),]$anno9 <- c("12")
  meta[meta$anno6 %in% c("EN18"),]$anno9 <- c("11")
  meta[meta$anno6 %in% c("IN23"),]$anno9 <- c("10")
  meta[meta$anno6 %in% c("EN9","EN5","EN30","EN3","EN29","EN28","EN25","EN24","EN23","EN15","EN11"),]$anno9 <- c("1")
  
  meta[meta$anno6 %in% c("IN15","IN9","IN8","IN7","IN6","IN5","IN31","IN30","IN3","IN29","IN28","IN25",
                         "IN24","IN21","IN20","IN2","IN19","IN18","IN17","IN16","IN12","IN11","IN10","IN1"),]$anno9 <- c("0")
  
  meta[meta$anno6 %in% c("EN1","EN20","EN27"),]$anno9 <- c("0_1")
}
dim(meta)
dim(neu_all)
meta <- meta[rownames(neu_all),]
neu_all$L3 <- paste0("L3_",meta$anno9)

saveRDS(neu_all,file = "neu_L1_L6_2.rds")






