
#based on
#https://yulab-smu.github.io/treedata-book/chapter2.html#tidytree


## REQUIRED PACKAGES & LIBRARIES
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")
BiocManager::install("ggnewscale")
BiocManager::install("ggplot2")
BiocManager::install("treeio")
library(ape)
library(ggplot2)
library(ggtree)
library(ggimage)
library(treeio)
library(tibble)
library(tidyr)
library(tidytree)
library(dplyr)
library(qsub)
library(stringr)
library(ggnewscale)
library(readr)
library(tidyverse)
library(reshape2)
library(gtable)
library(gridExtra)
install.packages("cowplot")
library(cowplot)


## REQUIRED FUNCTIONS 
source("~/orkun/projects/ongoing/AD_monitoring_project/analysisBySeb/ork_analyses_results/finalScripts/func_getLegend.R")

### read gtdb taxa names and info
##Archea gtdb taxonomy
#gtdbFile <- readr::read_tsv("~/orkun/projects/ongoing/AD_monitoring_project/analysisBySeb/tree_of_MAGs_gtdb_June2020/gtdb/ar122_taxonomy.tsv", col_names=FALSE, 
               #             col_types = "cc");

##Bacteria gtdb taxonomy
gtdbFile <- readr::read_tsv("~/orkun/projects/ongoing/AD_monitoring_project/analysisBySeb/tree_of_MAGs_gtdb_June2020/gtdb/bac120_taxonomy.tsv", col_names=FALSE, 
                            col_types = "cc");

### read original tree file into a tree object and read tips' information file
##Archea
#tree_original<-read.tree("~/orkun/projects/ongoing/AD_monitoring_project/analysisBySeb/tree_of_MAGs_gtdb_June2020/fasttree_familly_tree_ar122.nwk")
#myTable_original <- readr::read_tsv("~/orkun/projects/ongoing/AD_monitoring_project/analysisBySeb/tree_of_MAGs_gtdb_June2020/diversitySet_gtdb_v95_classify/gtdbtk.ar122.summary.tsv", col_names=TRUE, 
                #                    col_types = "cccdcddccddccccdidc");

###Bacteria tree - GTDB
tree_original<-read.tree("~/orkun/projects/ongoing/AD_monitoring_project/analysisBySeb/tree_of_MAGs_gtdb_June2020/fasttree_familly_tree_bac120.nwk")

##Bacteria MAGs - GTDB classification summary
myTable_original <- readr::read_tsv("~/orkun/projects/ongoing/AD_monitoring_project/analysisBySeb/tree_of_MAGs_gtdb_June2020/diversitySet_gtdb_v95_classify/gtdbtk.bac120.summary.tsv", col_names=TRUE, 
                                    col_types = "cccdcddcdcddcccciidc",na = c("", "NA"));
myTable_original$red_value <- as.numeric(myTable_original$red_value)


#get only needed subset of gtdb taxonomy file and remove file which is very large
colnames(gtdbFile) <- c("gtdbID","taxa")
gtdbFile_tipsOnly <- gtdbFile[which(gtdbFile$gtdbID %in% tree_original$tip.label), ]
remove(gtdbFile)

#seperate taxonomy info on gtdb file
dum=data.frame(domain=rep(NA,length(gtdbFile_tipsOnly$taxa)),
               phylum=rep(NA,length(gtdbFile_tipsOnly$taxa)),
               class=rep(NA,length(gtdbFile_tipsOnly$taxa)),
               order=rep(NA,length(gtdbFile_tipsOnly$taxa)),
               family=rep(NA,length(gtdbFile_tipsOnly$taxa)),
               genus=rep(NA,length(gtdbFile_tipsOnly$taxa)),
               species=rep(NA,length(gtdbFile_tipsOnly$taxa)))
for (i in 1:length(gtdbFile_tipsOnly$taxa)) {
  for (j in 1:7) {
    x=str_split(gtdbFile_tipsOnly$taxa[i],";")[[1]][j]
    dum[i,j]=substring(x,4)
  }
}
gtdbFile_tipsOnly <- cbind(gtdbFile_tipsOnly,dum)
remove(dum)

#take classification info for each bin and seperate
dum=data.frame(domain=rep(NA,length(myTable_original$classification)),
               phylum=rep(NA,length(myTable_original$classification)),
               class=rep(NA,length(myTable_original$classification)),
               order=rep(NA,length(myTable_original$classification)),
               family=rep(NA,length(myTable_original$classification)),
               genus=rep(NA,length(myTable_original$classification)),
               species=rep(NA,length(myTable_original$classification)))
for (i in 1:length(myTable_original$classification)) {
  for (j in 1:7) {
    x=str_split(myTable_original$classification[i],";")[[1]][j]
    dum[i,j]=substring(x,4)
  }
}
myTable_original <- cbind(myTable_original,dum)
remove(dum)

#OPTIONAL ---- PICK SUB-CLADE OF THE ORIGINAL TREE
###WORKS WITH PLOT 1 BELOW
#going back 5 levels from a234_m20724
#tree_collapsed_clade<-tree_subset(tree_original,'a234_m20724',levels_back=5) 
##

#re-name tree file- historical remnant - so that rest of code works
tree_collapsed_clade<-tree_original

#get leaf info from tree, which is for user genome only, and add other leafs from gtdb file
Leafs_clade <- as.vector(tree_collapsed_clade$tip.label)
myTable_collapsed_clade <- myTable_original[which(myTable_original$user_genome %in% Leafs_clade), ]
myTable_collapsed_clade <- myTable_collapsed_clade[ , which(colnames(myTable_collapsed_clade) %in% c("user_genome","red_value","phylum","class","order","family","genus","species"))]
dum <- gtdbFile_tipsOnly[which(gtdbFile_tipsOnly$gtdbID %in% tree_collapsed_clade$tip.label),]
dum <- dum[ , which(colnames(dum) %in% c("gtdbID","domain","phylum","class","order","family","genus","species"))]
dum$domain <- "NA"
colnames(dum) <- c("user_genome","red_value","phylum","class","order","family","genus","species")
myTable_collapsed_clade <- rbind(dum,myTable_collapsed_clade)
remove(dum)

#re-assign tip labels using info from gtdb
for (i in 1:length(tree_collapsed_clade$tip.label)) {
  dum <- which(gtdbFile_tipsOnly$gtdbID %in% tree_collapsed_clade$tip.label[i])
  if (length(dum) > 0) tree_collapsed_clade$tip.label[i] <- gtdbFile_tipsOnly$species[dum]
  dumber <- which(gtdbFile_tipsOnly$gtdbID %in% myTable_collapsed_clade$user_genome[i])
  if (length(dumber) > 0) myTable_collapsed_clade$user_genome[i] <- gtdbFile_tipsOnly$species[dumber]
}

#OPTIONAL ---- PICK ONLY THE MAGS FROM THE ORIGINAL TREE AND DROP ANY GTDB SPECIES
tree_drop_clades<-drop.tip(tree_collapsed_clade,gtdbFile_tipsOnly$species)  #drop GTDB species
a<-myTable_collapsed_clade$user_genome[which(is.na(myTable_collapsed_clade$red_value))] #find MAGs with RED=NA
tree_drop_clades2<-drop.tip(tree_drop_clades,a) #further drop MAGs with RED=NA
myTable_drop_clade <- myTable_collapsed_clade[which(myTable_collapsed_clade$user_genome %in% tree_drop_clades2$tip.label), ]

#get files to pass to plot
#treeToPlot <- tree_collapsed_clade
#tableToPlot <- myTable_collapsed_clade
treeToPlot <- tree_drop_clades2
tableToPlot <- myTable_drop_clade

## PLOTTING PART ## 

#PLOT 1 -
#Tree PLOT - Whole bacteria tree
#group leafs for coloring purposes at level of choice
leafList_clade <- split(tableToPlot$user_genome,tableToPlot$phylum)  #use family for coloring clade
leafListToPlot <- leafList_clade

#first create layers of info to be added as heatmap to the tree plot
nSpeciesInfo <- data.frame(as.numeric(tableToPlot[,2])) #RED
rownames(nSpeciesInfo) <- tableToPlot[,1]
colnames(nSpeciesInfo) <- c("RED")

#create tree object
p <- ggtree(treeToPlot, layout='circular', branch.length="branch length") %<+% tableToPlot

#add median_RED info for leaves onto same tree
p2 <- p + new_scale_fill()
pTree<-gheatmap(p2, nSpeciesInfo, offset=0.05, width=0.03, colnames = FALSE) +
  scale_fill_viridis_c(option="D", direction=-1,
                       na.value="white", name="RED")

pTreeFull<-groupOTU(pTree, leafListToPlot, 'Phylum') +
  aes(color=phylum) +             #use with tree_collapsed_clade
  geom_tippoint(size=0.25) +
  #geom_tiplab(aes(angle=-10),size = 2,align = F,hjust=.5) +
  #geom_tiplab(aes(angle=0,label=label,subset = label %in% labToShow),size = 2.5,align = F,hjust=-0.01) +
  #geom_tiplab(aes(label=label, subset = label=="Bin_c1"),color="464545",size=2.75) +
  #geom_nodelab() +
  scale_color_discrete("Class",breaks = names(leafListToPlot)) +
  #  labs(tag="A") + 
  theme(plot.margin = margin(t=0, unit="cm"),
        plot.tag=element_text(family="Times", face="bold", size=20),
        legend.key.size = unit(0.2,"line"),
        legend.position = "right",legend.box = "horizontal",
        legend.margin = margin(t=0, unit="cm"),
        legend.title = element_text(size = 8,face="bold"))



#PLOT 2 -
#Tree PLOT - all of bacteria from above
#p <- ggtree(treeToPlot, layout="circular", branch.length="branch length") %<+% tableToPlot
p <- ggtree(treeToPlot, branch.length="branch length") %<+% tableToPlot

a<-which(treeToPlot$tip.label %in% leafList_clade$Myxococcota)
b<-which(treeToPlot$edge[,2] %in% a)
c<-min(treeToPlot$edge[b,1])
cladeID1 <- treeToPlot$edge[which(treeToPlot$edge[,2]==c),1]

a<-which(treeToPlot$tip.label %in% leafList_clade$Bacteroidota)
b<-which(treeToPlot$edge[,2] %in% a)
c<-min(treeToPlot$edge[b,1])
cladeID2 <- treeToPlot$edge[which(treeToPlot$edge[,2]==c),1]

a<-which(treeToPlot$tip.label %in% leafList_clade$Firmicutes)
b<-which(treeToPlot$edge[,2] %in% a)
c<-min(treeToPlot$edge[b,1])
cladeID3 <- treeToPlot$edge[which(treeToPlot$edge[,2]==c),1]

p2 <- scaleClade(p, cladeID1-4, .2) 
p3 <- scaleClade(p2, cladeID2-3, .3)
p4 <- scaleClade(p3, cladeID3-5, .3)
p5 <- collapse(p4, cladeID1, 'min', alpha=0.6, fill = "#AC88FF", clade_name="Myxococcota") %>%
  collapse(cladeID2-2, 'min', alpha=0.6, fill = "#E08B00" , clade_name="Bacteroidota") %>%
  collapse(cladeID3-2, 'min', alpha=0.6, fill = "#00BF7D" , clade_name="Firmicutes")
#note that the colors in the above lines are picked to match those used by 
#scale_color_discrete used below. The scale_color_discrete uses real colors that are accessible by
#scales::hue_pal()(n) where n is the number of breaks used, in this case 25

#add phylum coloring and plot tree 
pTreeFull<-groupOTU(p5, leafListToPlot, 'Phylum') +
  aes(color=phylum) +             #use with tree_collapsed_clade
  geom_tippoint(size=0.25) +
  geom_tiplab(aes(angle=-10),size = 2,align = F,hjust=.5) +
  scale_color_discrete("Phylum",breaks = names(leafListToPlot)) +
#  labs(tag="A") + 
  theme(plot.margin = margin(t=0, unit="cm"),
                        plot.tag=element_text(family="Times", face="bold", size=20),
                        legend.key.size = unit(0.5,"line"),
                        legend.position = "right",legend.box = "horizontal",
                        legend.margin = margin(t=0, unit="cm"),
                        legend.title = element_text(size = 8,face="bold"))
#  guides(color=guide_legend(order=1),fill=guide_legend(order=2))

#PLOT 3 -
#MAG analysis PLOT
pBar<-ggplot(df4,aes(phylumEdited)) + 
  geom_bar(aes(fill=phylumEdited, alpha=variable, y=value),colour="black",size=0.1,position="stack",stat="identity") + 
  scale_y_continuous() + scale_alpha_discrete(labels=c("<0.6","0.6 - 0.85","0.85 - 1.0","known species")) +
  scale_color_discrete("Phylum", breaks = names(leafList)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.margin = margin(t=0, unit="cm"),
        legend.key.size = unit(0.5,"line"),
        legend.margin = margin(t=0, unit="cm"),
        legend.position = c(0.85, 0.8),
        legend.title = element_text(size = 8,face="bold"),
        plot.tag=element_text(family="Times",face="bold", size=20),
        axis.text.x = element_blank()) +
  labs(x="",y="nDMAGs",fill="Phylum",alpha="RED Level") + guides(fill="none")
#  scale_fill_viridis_d() + 
#  coord_polar(theta="y")


## PROCESSS PLOT FOR FINAL DISPLAY ## 
#Remove the legend from the tree plot
legendTree <- func_getLegend(pTreeFull)
pTreeFull <- pTreeFull + theme(legend.position="none")


#arrange plots, using cowplot
ggdraw() + 
  draw_plot(pTreeFull, x = 0, y = 0.05, width =1, height = 0.95) +
  draw_plot(legendTree, x = 0.7, y = 0.09, width = .35, height = .15)




