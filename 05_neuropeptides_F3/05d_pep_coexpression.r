
if(T){
  library(tidyverse)
  library(qs)
  library(reshape2)
  library(ggplot2)
  library(ggpubr)
}

path <- ""

########################################################################################################################
############################# Part 1 pep stat caculation 
########################################################################################################################

mouse_cell <- qread(paste0(path,"mou_pet_id.qs"))
macaque_cell <- qread(paste0(path,"mac_pet_id.qs"))
marmoset_cell <- qread(paste0(path,"mar_pet_id.qs"))

peti2 <- qread(file = paste0(path,"p_pep_stat_use_all_spe.qs"))
peti <- unique(as.character(peti2$features.plot))

mouse_cell2 <- mouse_cell[which(lapply(mouse_cell, length) > 100)]
macaque_cell2 <- macaque_cell[which(lapply(macaque_cell, length) > 100)]
marmoset_cell2 <- marmoset_cell[which(lapply(marmoset_cell, length) > 100)]

common_pep <- intersect(names(mouse_cell2),names(macaque_cell2))
common_pep <- intersect(common_pep,names(marmoset_cell2))

##### mouse
slot <- list()
for (i in 1:length(common_pep)) {
  use_pep <- common_pep[i]
  other_pep <- setdiff(common_pep,use_pep)
  use_pep_cell <- mouse_cell2[[use_pep]]
  use_pep_mat <- matrix(nrow = length(use_pep_cell),ncol = 37,data = 0)  
  colnames(use_pep_mat) <- other_pep
  rownames(use_pep_mat) <- 1:length(use_pep_cell)
  ####
  print(paste0("Doing",use_pep))
  for (j in 1:length(use_pep_cell)){
    for (k in 1:length(other_pep)){
      if (length(intersect(use_pep_cell[j],mouse_cell2[[other_pep[k]]])) > 0) {
        use_pep_mat[j,other_pep[k]] <- 1
      }else{
        use_pep_mat[j,other_pep[k]] <- 0
      }
      
    }
    
  }
  print(paste0("Done",use_pep))
  slot[[use_pep]] <- use_pep_mat
  
}
mou_slot <- slot

#### marmoset
slot <- list()
for (i in 1:length(common_pep)) {
  use_pep <- common_pep[i]
  other_pep <- setdiff(common_pep,use_pep)
  use_pep_cell <- marmoset_cell2[[use_pep]]
  use_pep_mat <- matrix(nrow = length(use_pep_cell),ncol = 37,data = 0)  
  colnames(use_pep_mat) <- other_pep
  rownames(use_pep_mat) <- 1:length(use_pep_cell)
  ####
  print(paste0("Doing",use_pep))
  for (j in 1:length(use_pep_cell)){
    for (k in 1:length(other_pep)){
      if (length(intersect(use_pep_cell[j],marmoset_cell2[[other_pep[k]]])) > 0) {
        use_pep_mat[j,other_pep[k]] <- 1
      }else{
        use_pep_mat[j,other_pep[k]] <- 0
      }
      
    }
    
  }
  print(paste0("Done",use_pep))
  slot[[use_pep]] <- use_pep_mat
  
}
mar_slot <- slot

#### macaque
slot <- list()
for (i in 1:length(common_pep)) {
  use_pep <- common_pep[i]
  other_pep <- setdiff(common_pep,use_pep)
  use_pep_cell <- macaque_cell2[[use_pep]]
  use_pep_mat <- matrix(nrow = length(use_pep_cell),ncol = 37,data = 0)  
  colnames(use_pep_mat) <- other_pep
  rownames(use_pep_mat) <- 1:length(use_pep_cell)
  ####
  print(paste0("Doing",use_pep))
  for (j in 1:length(use_pep_cell)){
    for (k in 1:length(other_pep)){
      if (length(intersect(use_pep_cell[j],macaque_cell2[[other_pep[k]]])) > 0) {
        use_pep_mat[j,other_pep[k]] <- 1
      }else{
        use_pep_mat[j,other_pep[k]] <- 0
      }
      
    }
    
  }
  print(paste0("Done",use_pep))
  slot[[use_pep]] <- use_pep_mat
  
}
mac_slot <- slot

save(mou_slot,mar_slot,mac_slot ,file = "stat_pep_co_expre.Rdata")
################################ -----------------------------------------------------------------------


########################################################################################################################
############################# Part 2 pep stat caculation 
########################################################################################################################

load("stat_pep_co_expre.Rdata")
color3  <- c("#e20612","#F57F17","#00b0eb")
names(color3) <-  c("mac", "mar", "mou")

mou_slot2 <- list()
mac_slot2 <- list()
mar_slot2 <- list()
pep_use <- names(mou_slot)
for (i in 1:length(pep_use)) {
  a <- data.frame(mou_slot[[pep_use[i]]])
  a$sum <- apply(a, 1, sum)
  a$id <- paste0(rownames(a),"_",pep_use[i])
  a <- a[,c("id","sum")]
  a$pep <- pep_use[i]
  a$spe <- "mou"
  mou_slot2[[pep_use[i]]] <- a
  
  a <- data.frame(mac_slot[[pep_use[i]]])
  a$sum <- apply(a, 1, sum)
  a$id <- paste0(rownames(a),"_",pep_use[i])
  a <- a[,c("id","sum")]
  a$pep <- pep_use[i]
  a$spe <- "mac"
  mac_slot2[[pep_use[i]]] <- a
  
  a <- data.frame(mar_slot[[pep_use[i]]])
  a$sum <- apply(a, 1, sum)
  a$id <- paste0(rownames(a),"_",pep_use[i])
  a <- a[,c("id","sum")]
  a$pep <- pep_use[i]
  a$spe <- "mar"
  mar_slot2[[pep_use[i]]] <- a
}

mou_slot3  <- do.call(rbind,mou_slot2)
mou_slot3$pep <- factor(mou_slot3$pep,levels = pep_use)
mac_slot3  <- do.call(rbind,mac_slot2)
mac_slot3$pep <- factor(mac_slot3$pep,levels = pep_use)
mar_slot3  <- do.call(rbind,mar_slot2)
mar_slot3$pep <- factor(mar_slot3$pep,levels = pep_use)

plot_slot <- rbind(mou_slot3,mar_slot3)
plot_slot <- rbind(plot_slot,mac_slot3)

p <- plot_slot %>% ggplot(aes(x = pep,y = sum,fill = spe))+
  geom_boxplot(outlier.shape = NA)+
  #scale_color_manual(values = color3)+
  scale_fill_manual(values = color3)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
ggsave(p,filename = paste0(path,"pep_co_expre_number.pdf"),width = 8,height = 4)

plot_slot %>% ggplot(aes(x = pep,y = sum,fill = spe))+
  geom_violin(outlier.shape = NA)+
  #scale_color_manual(values = color3)+
  scale_fill_manual(values = color3)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))

### plot　
mou_max_min <- data.frame(mou_slot3 %>% group_by(pep) %>% summarise(max_co = max(sum),
                                                                    min_co = min(sum)))
mou_max_min$spe <- "mou"

mar_max_min <- data.frame(mar_slot3 %>% group_by(pep) %>% summarise(max_co = max(sum),
                                                                    min_co = min(sum)))
mar_max_min$spe <- "mar"

mac_max_min <- data.frame(mac_slot3 %>% group_by(pep) %>% summarise(max_co = max(sum),
                                                                    min_co = min(sum)))
mac_max_min$spe <- "mac"


plot_min_max <- rbind(mou_max_min,mar_max_min)
plot_min_max <- rbind(plot_min_max,mac_max_min)

p <- plot_min_max %>% ggplot(aes(x = pep,y = max_co,color = spe))+
  geom_point(size = 2.5,aes(shape = spe))+
  theme_bw()+
  scale_color_manual(values = color3)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  scale_y_continuous(limits = c(0,13))
ggsave(p,filename = paste0(path,"pep_co_expre_max_number.pdf"),width = 8,height = 4)

#### plot 

get_co_expre_fre <- function(data){
  a <- data.frame(table(data$sum))
  return(a)
}

mou_slot2_fre <- lapply(mou_slot2, function(x){
  x <- get_co_expre_fre(data = x)
})
mac_slot2_fre <- lapply(mac_slot2, function(x){
  x <- get_co_expre_fre(data = x)
})
mar_slot2_fre <- lapply(mar_slot2, function(x){
  x <- get_co_expre_fre(data = x)
})

use_pep <- names(mou_slot2)
for (i in 1:length(use_pep)) {
  mou_slot2_fre[[use_pep[i]]]$pep <- use_pep[i]
  mou_slot2_fre[[use_pep[i]]]$spe <- "mou"
  mou_slot2_fre[[use_pep[i]]]$ratio <- mou_slot2_fre[[use_pep[i]]]$Freq / sum(mou_slot2_fre[[use_pep[i]]]$Freq)
  
  mar_slot2_fre[[use_pep[i]]]$pep <- use_pep[i]
  mar_slot2_fre[[use_pep[i]]]$spe <- "mar"
  mar_slot2_fre[[use_pep[i]]]$ratio <- mar_slot2_fre[[use_pep[i]]]$Freq / sum(mar_slot2_fre[[use_pep[i]]]$Freq)
  
  mac_slot2_fre[[use_pep[i]]]$pep <- use_pep[i]
  mac_slot2_fre[[use_pep[i]]]$spe <- "mac"
  mac_slot2_fre[[use_pep[i]]]$ratio <- mac_slot2_fre[[use_pep[i]]]$Freq / sum(mac_slot2_fre[[use_pep[i]]]$Freq)
}

mou_slot2_fre2 <- do.call(rbind,mou_slot2_fre)
mar_slot2_fre2 <- do.call(rbind,mar_slot2_fre)
mac_slot2_fre2 <- do.call(rbind,mac_slot2_fre)

plot_fre_slot <- rbind(mou_slot2_fre2,mar_slot2_fre2)
plot_fre_slot <- rbind(plot_fre_slot,mac_slot2_fre2)

plot_fre_slot$Var1 <- as.numeric(as.character(plot_fre_slot$Var1))
plot_fre_slot$ratio <- as.numeric(as.character(plot_fre_slot$ratio))


p <- plot_fre_slot %>% ggplot(aes(x = Var1,y =  ratio,color = spe))+
  geom_smooth(span = 0.5)+
  scale_color_manual(values = color3)+
  theme_bw()+
  scale_x_continuous(breaks = c(0,1:12),limits = c(0,12))+
  theme(legend.position = "none")+
  xlab("Co-expression number")+
  ylab("Frequence")
ggsave(p,filename = paste0(path,"pep_co_expre_Cumulative_distribution_frequency.pdf"),width = 8,height = 4)






