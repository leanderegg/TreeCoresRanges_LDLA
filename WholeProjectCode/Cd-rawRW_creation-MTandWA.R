######## Creating cleaned raw RW for BAI analysis ########

#=========== Montanta =========================
#++++++++++++++++++++++++++++++++++++++++++++++
# NOTE: code between the ++++ markes resuses
# nomenclature for different species. so must be rerun for each species
## low elevation
# +++++++++++++++++ TSHE ++++++++++++++++++++++++++++
csv1 <- "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/MT-TSHE-L-RWdata05_07_14.csv" 
data1 <- read.csv(csv1, sep = ",", header = T)
treedata1 <- core_average(data=data1)
## mid elevation
csv2 <- "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/MT-TSHE-M-RWdata05_12_14.csv"
data2 <- read.csv(csv2, sep = ",", header = T)
treedata2 <- core_average(data=data2)
## High elevation
csv3 <- "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/MT-TSHE-H-RWdata05_09_14.csv"
data3 <- read.csv(csv3, sep = ",", header = T)
colnames(data3)[1] <- "Transect" 
treedata3 <- core_average(data=data3)


ncols <- min(c(ncol(treedata1), ncol(treedata2), ncol(treedata3)))
rwdata_all <- rbind(treedata1[,1:ncols], treedata2[1:ncols], treedata3[1:ncols])
rwdata_all$Tree <- as.character(rwdata_all$Tree)
####### Need to change when switching nomenclature
rwdata_all$Tag <- paste(rwdata_all$Species, rwdata_all$Elev, rwdata_all$Tree, sep = "-")
#made unique tree tag
rwdata_alltag <- rwdata_all[,c(1:4,ncol(rwdata_all),5:(ncol(rwdata_all)-1))]
#reordered columns to get tag in right place
rwdata_alltag <- rwdata_alltag[order(rwdata_alltag$Tag),]

write.csv(rwdata_alltag, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/MT-TSHE-avgRingWidth_7-31-15.csv")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++






# +++++++++++++++++ PSME ++++++++++++++++++++++++++++
csv1 <- "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/MT-PSME-L-RWdata05_19_14.csv" 
data1 <- read.csv(csv1, sep = ",", header = T)
treedata1 <- core_average(data=data1)
## mid elevation
csv2 <- "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/MT-PSME-M-RWdata06_04_14.csv"
data2 <- read.csv(csv2, sep = ",", header = T)
treedata2 <- core_average(data=data2)
## High elevation
csv3 <- "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/MT-PSME-H-RWdata05_16_14.csv"
data3 <- read.csv(csv3, sep = ",", header = T)
colnames(data3)[1] <- "Transect" 
treedata3 <- core_average(data=data3)


ncols <- min(c(ncol(treedata1), ncol(treedata2), ncol(treedata3)))
rwdata_all <- rbind(treedata1[,1:ncols], treedata2[1:ncols], treedata3[1:ncols])
rwdata_all$Tree <- as.character(rwdata_all$Tree)
####### Need to change when switching nomenclature
rwdata_all$Tag <- paste(rwdata_all$Species, rwdata_all$Elev, rwdata_all$Tree, sep = "-")
#made unique tree tag
rwdata_alltag <- rwdata_all[,c(1:4,ncol(rwdata_all),5:(ncol(rwdata_all)-1))]
#reordered columns to get tag in right place
rwdata_alltag <- rwdata_alltag[order(rwdata_alltag$Tag),]

write.csv(rwdata_alltag, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/MT-PSME-avgRingWidth_7-31-15.csv")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# +++++++++++++++++ ABLA ++++++++++++++++++++++++++++
csv1 <- "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/MT-ABLA-L-RWdata06_05_14.csv" 
data1 <- read.csv(csv1, sep = ",", header = T)
treedata1 <- core_average(data=data1)
## mid elevation
csv2 <- "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/MT-ABLA-M-RWdata06_05_14.csv"
data2 <- read.csv(csv2, sep = ",", header = T)
treedata2 <- core_average(data=data2)
## High elevation
csv3 <- "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/MT-ABLA-H-RWdata06_05_14.csv"
data3 <- read.csv(csv3, sep = ",", header = T)
colnames(data3)[1] <- "Transect" 
treedata3 <- core_average(data=data3)


ncols <- min(c(ncol(treedata1), ncol(treedata2), ncol(treedata3)))
rwdata_all <- rbind(treedata1[,1:ncols], treedata2[1:ncols], treedata3[1:ncols])
rwdata_all$Tree <- as.character(rwdata_all$Tree)
####### Need to change when switching nomenclature
rwdata_all$Tag <- paste(rwdata_all$Species, rwdata_all$Elev, rwdata_all$Tree, sep = "-")
#made unique tree tag
rwdata_alltag <- rwdata_all[,c(1:4,ncol(rwdata_all),5:(ncol(rwdata_all)-1))]
#reordered columns to get tag in right place
rwdata_alltag <- rwdata_alltag[order(rwdata_alltag$Tag),]

write.csv(rwdata_alltag, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/MT-ABLA-avgRingWidth_7-31-15.csv")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





#=========== Washington =========================
#++++++++++++++++++++++++++++++++++++++++++++++
# NOTE: code between the ++++ markes resuses
# nomenclature for different species. so must be rerun for each species
## low elevation
# +++++++++++++++++ TSHE ++++++++++++++++++++++++++++
# csv1 <- "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/MT_RingWidths/MT_RW_csvs_05-14/MT-TSHE-L-RWdata05_07_14.csv" 
# data1 <- read.csv(csv1, sep = ",", header = T)
# treedata1 <- core_average(data=data1)
## mid elevation
csv2 <- "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/WA-TSHE-M-RWdata06_09_14.csv"
data2 <- read.csv(csv2, sep = ",", header = T)
treedata2 <- core_average(data=data2)
## High elevation
csv3 <- "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/WA-TSHE-H-RWdata06_09_14.csv"
data3 <- read.csv(csv3, sep = ",", header = T)
colnames(data3)[1] <- "Transect" 
treedata3 <- core_average(data=data3)


ncols <- min(c(ncol(treedata2), ncol(treedata3)))
rwdata_all <- rbind(treedata2[1:ncols], treedata3[1:ncols])
rwdata_all$Tree <- as.character(rwdata_all$Tree)
####### Need to change when switching nomenclature
rwdata_all$Tag <- paste(rwdata_all$Species, rwdata_all$Elev, rwdata_all$Tree, sep = "-")
#made unique tree tag
rwdata_alltag <- rwdata_all[,c(1:4,ncol(rwdata_all),5:(ncol(rwdata_all)-1))]
#reordered columns to get tag in right place
rwdata_alltag <- rwdata_alltag[order(rwdata_alltag$Tag),]

write.csv(rwdata_alltag, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/WA-TSHE-avgRingWidth_7-31-15.csv")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++






# +++++++++++++++++ PSME ++++++++++++++++++++++++++++
csv1 <- "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/WA-PSME-L-RWdata06_09_14.csv" 
data1 <- read.csv(csv1, sep = ",", header = T)
treedata1 <- core_average(data=data1)
## mid elevation
csv2 <- "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/WA-PSME-M-RWdata06_09_14.csv" 
data2 <- read.csv(csv2, sep = ",", header = T)
treedata2 <- core_average(data=data2)
## High elevation
csv3 <- "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/WA-PSME-H-RWdata06_10_14.csv" 
data3 <- read.csv(csv3, sep = ",", header = T)
colnames(data3)[1] <- "Transect" 
treedata3 <- core_average(data=data3)


ncols <- min(c(ncol(treedata1), ncol(treedata2), ncol(treedata3)))
rwdata_all <- rbind(treedata1[,1:ncols], treedata2[1:ncols], treedata3[1:ncols])
rwdata_all$Tree <- as.character(rwdata_all$Tree)
####### Need to change when switching nomenclature
rwdata_all$Tag <- paste(rwdata_all$Species, rwdata_all$Elev, rwdata_all$Tree, sep = "-")
#made unique tree tag
rwdata_alltag <- rwdata_all[,c(1:4,ncol(rwdata_all),5:(ncol(rwdata_all)-1))]
#reordered columns to get tag in right place
rwdata_alltag <- rwdata_alltag[order(rwdata_alltag$Tag),]

write.csv(rwdata_alltag, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/WA-PSME-avgRingWidth_7-31-15.csv")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# +++++++++++++++++ ABLA ++++++++++++++++++++++++++++
csv1 <- "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/WA-ABLA-L-RWdata06_10_14.csv" 
data1 <- read.csv(csv1, sep = ",", header = T)
treedata1 <- core_average(data=data1)
## mid elevation
csv2 <- "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/WA-ABLA-M-RWdata06_10_14.csv" 
data2 <- read.csv(csv2, sep = ",", header = T)
treedata2 <- core_average(data=data2)
## High elevation
csv3 <- "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/WA-ABLA-H-RWdata06_13_14.csv" 
data3 <- read.csv(csv3, sep = ",", header = T)
colnames(data3)[1] <- "Transect" 
treedata3 <- core_average(data=data3)


ncols <- min(c(ncol(treedata1), ncol(treedata2), ncol(treedata3)))
rwdata_all <- rbind(treedata1[,1:ncols], treedata2[1:ncols], treedata3[1:ncols])
rwdata_all$Tree <- as.character(rwdata_all$Tree)
####### Need to change when switching nomenclature
rwdata_all$Tag <- paste(rwdata_all$Species, rwdata_all$Elev, rwdata_all$Tree, sep = "-")
#made unique tree tag
rwdata_alltag <- rwdata_all[,c(1:4,ncol(rwdata_all),5:(ncol(rwdata_all)-1))]
#reordered columns to get tag in right place
rwdata_alltag <- rwdata_alltag[order(rwdata_alltag$Tag),]

write.csv(rwdata_alltag, "/Users/leeanderegg/Dropbox/TreeCoresRanges(Leander) (1)/WA_RingWidths/WA_RW_csvs/WA-ABLA-avgRingWidth_7-31-15.csv")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




