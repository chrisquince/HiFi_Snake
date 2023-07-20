library(ggtree)
library(treeio)
library(ggnewscale)
library(tidyverse)
library(wesanderson)

library(ape)
library(tidyr)
library(stringr)
library(gtable)
library(gridExtra)
library(ggimage)
library(tidytree)
library(readr)
library(reshape2)
library(cowplot)
library(ggplot2)
library(tibble)
library(dplyr)

## REQUIRED FUNCTIONS 
func_getLegend <- function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }

propagate_from_child = function(node,tbl_tree,node_taxa){
  nb_taxa = length(node_taxa[toString(node),])
  parent = parent(tbl_tree,node)$node
  if (length(parent)!=0)
  {
    # if parent doesn't have taxonomy, just give it child taxonomy
    if (!(parent %in% rownames(node_taxa)))
    {
      node_taxa[toString(parent),] = node_taxa[toString(node),]
    } else { # take the consensus of childs
      match_nb = max(which(node_taxa[toString(node),]==node_taxa[toString(parent),]))+1
      if (match_nb<=nb_taxa)
      {
        node_taxa[toString(parent),match_nb:nb_taxa] = NA
      }
    }
    # propagate current annotation to parents
    node_taxa = propagate_from_child(parent,tbl_tree,node_taxa)
  }
  return(node_taxa)
}

propagate_from_parents = function(parent,tbl_tree,taxa,node_taxa,node_taxa_mags)
{
  for (node in child(tbl_tree,parent)$node)
  {
    if (!(node %in% rownames(node_taxa)))
    {
      node_taxa_mags[toString(node),] = taxa
      node_taxa_mags = propagate_from_parents(node,tbl_tree,taxa,node_taxa,node_taxa_mags)
    }
  }
  return(node_taxa_mags)
}

annotate_from_parents = function(node,tbl_tree,node_taxa,node_taxa_mags)
{
  nb_taxa = length(node_taxa[toString(node),])
  if (!(node %in% rownames(node_taxa_mags)))
  {
    parent = parent(tbl_tree,node)$node
    # while parent node is not annotated by refference, checking parent
    while (!(parent %in% rownames(node_taxa)))
    {
      parent = parent(tbl_tree,parent)$node
    }
    # use gd_parent and parent common annotation to decide on mag annotation
    gd_parent = parent(tbl_tree,parent)$node
    match_nb = max(which(node_taxa[toString(gd_parent),]==node_taxa[toString(parent),]))+1
    taxa = node_taxa[toString(gd_parent),]
    if (match_nb<=nb_taxa)
    {
      taxa[match_nb:nb_taxa] = NA
    }
    # then consider that all childs and parent have this taxonomic annotation
    node_taxa_mags = propagate_from_parents(parent,tbl_tree,taxa,node_taxa,node_taxa_mags)
  }
  return(node_taxa_mags)
}

node_annotation = function(tree_original,reff_taxa,mag_taxa,mag_metadata)
{
  ### --- deal with gtdb reffs --------
  # get a mapping of gtdb name to node name
  tbl_tree = as_tibble(tree_original)
  tbl_tree_gtdb = tbl_tree[tbl_tree$label %in% reff_taxa$genome,]

  # it may happens that red value is 1 even if no gtdb reff is around, for the simple fact that not all gtdb refs are present, if the red_value is 1, consider the mag as a gtdb reff.
  tbl_tree_gtdb = rbind(tbl_tree_gtdb,tbl_tree[tbl_tree$label %in% mag_taxa[mag_taxa$red_value==1,'genome'],])
  rownames(tbl_tree_gtdb) = tbl_tree_gtdb$label

  # build a node taxonomy dataframe
  node_taxa = data.frame(reff_taxa[,c("domain", "phylum", "class", "order", "family", "genus", "species")],row.names = tbl_tree_gtdb[reff_taxa$genome,]$node)
  # do the same with mags
  node_taxa_2 = data.frame(mag_taxa[mag_taxa$red_value==1,c("domain", "phylum", "class", "order", "family", "genus", "species")],row.names = tbl_tree_gtdb[mag_taxa[mag_taxa$red_value==1,'genome'],]$node)
  node_taxa = rbind(node_taxa, node_taxa_2)

  # iterate over all gtdb refference genomes and annotate intermediary nodes
  for (node in rownames(node_taxa))
  {
    node_taxa = propagate_from_child(as.integer(node),tbl_tree,node_taxa)
  }
  ### --- deal with mags --------
  # get a mapping of mag name to node name
  mags = tree_original$tip.label[!(tree_original$tip.label %in% reff_taxa$genome)]
  tbl_tree_mags = tbl_tree[tbl_tree$label %in% mags,]
  rownames(tbl_tree_mags) = tbl_tree_mags$label
  # build a node taxonomy dataframe
  node_taxa_mags = data.frame(reff_taxa[0,c("domain", "phylum", "class", "order", "family", "genus", "species")])

  # iterate over all mags and find annotation for mags and intermediary nodes
  for (mag in mags)
  {
    node = tbl_tree_mags[mag,]$node
    node_taxa_mags = annotate_from_parents(node,tbl_tree,node_taxa,node_taxa_mags)
  }

  # finaly! concatenate both node annotation
  node_annotation = rbind(node_taxa,node_taxa_mags)
  node_annotation$label = tbl_tree[rownames(node_annotation),]$label
  node_annotation$node = as.numeric(rownames(node_annotation))

  # use mag annotation to change internal node annotation
  for (mag in mags)
  {
    node = tbl_tree_mags[mag,]$node
    node_annotation = propagate_from_child(as.integer(node),tbl_tree,node_annotation)
  }

  ## append mag annotation to node_taxa_mags
  nb_additional_col = dim(mag_metadata[,-1])[2]
  # first mags
  tip_mags = node_annotation[node_annotation$label %in% mag_metadata$MAGs,]
  tip_mags = cbind(tip_mags,mag_metadata[match(tip_mags$label,mag_metadata$MAGs),-1])
  # the all other nodes
  other_nodes = node_annotation[!(node_annotation$label %in% mag_metadata$MAGs),]
  colnames_nodes = colnames(other_nodes)
  other_nodes = cbind(other_nodes,data.frame(matrix(NA, nrow = dim(other_nodes)[1], ncol = nb_additional_col)))
  colnames(other_nodes) = c(colnames_nodes,colnames(mag_metadata[,-1]))
  # merge again
  node_annotation = rbind(tip_mags,other_nodes)

  ## rename nodes so that we can pass taxonomic annotation of nodes when plotting
  is_tip = tbl_tree$label %in% tree_original$tip.label
  tbl_tree$label[!is_tip]=paste(tbl_tree$node[!is_tip],tbl_tree$label[!is_tip], sep = '_')
  is_tip = node_annotation$label %in% tree_original$tip.label
  node_annotation$label[!is_tip] = paste(rownames(node_annotation)[!is_tip],node_annotation$label[!is_tip], sep = '_')

  # package output
  result=list("tree"=as.phylo(tbl_tree),'annotation'=node_annotation)
  return(result)
}

gheatmap2 = function (p, data, offset = 0, width = 1)
{
    pal <- wes_palette("Zissou1", 100, type = "continuous")
    variable <- value <- lab <- y <- NULL
    width <- width * (p$data$x %>% range(na.rm = TRUE) %>% diff)/ncol(data)
    if (width==0){width = 1/ncol(data)}
    isTip <- x <- y <- variable <- value <- from <- to <- NULL
    df <- p$data
    df <- df[df$isTip, ]
    start <- max(df$x, na.rm = TRUE) + offset
    dd <- as.data.frame(data)
    i <- order(df$y)
    i <- i[!is.na(df$y[i])]
    lab <- df$label[i]
    dd <- dd[lab, , drop = FALSE]
    dd$y <- sort(df$y)
    dd$lab <- lab
    dd <- gather(dd, variable, value, -c(lab, y))
    i <- which(dd$value == "")
    if (length(i) > 0) {
        dd$value[i] <- NA
    }
    dd$variable <- factor(dd$variable, levels = colnames(data))
    V2 <- start + as.numeric(dd$variable) * width
    mapping <- data.frame(from = dd$variable, to = V2)
    mapping <- unique(mapping)
    dd$x <- V2
    dd$width <- width

    # draw the heatmap 
    p2 <- p + geom_tile(data = dd, aes(x, y, fill = value), 
            width = width, color = "white",inherit.aes = FALSE)
    # add annotation on tiles 
    p2 <- p2 + geom_text(data = dd, aes(x,y,label=value, fontface="bold"),size=4/sqrt(ncol(data)),inherit.aes = FALSE)
    # add color 
    p2 <- p2 + scale_fill_gradientn(colours = pal,na.value="white")
    p2 <- p2 + theme(legend.position = "right")
    attr(p2, "mapping") <- mapping
    return(p2)
}



## ------ Path to files ------
GTDB_TAXA = "/mnt/gpfs/chris/Installed/miniconda3/envs/STRONG/share/gtdbtk-1.2.0/db/taxonomy/gtdb_taxonomy.tsv"
BAC_TREE = "/mnt/gpfs/seb/Project/Anaerobic_Digester/WGS/Assemblies/MAGs/gtdb_v95/fasttree_familly_tree_bac120.nwk"
BAC_SUMMARY = "/mnt/gpfs/seb/Project/Anaerobic_Digester/WGS/Assemblies/MAGs/gtdb_v95/gtdbtk.bac120.summary.tsv"
METADATA = "/mnt/gpfs/seb/Project/Anaerobic_Digester/WGS/Assemblies/MAGs/MAGAnalysis/bin_analysis/dMAG_InfoSummary_19_10_20.tsv"

# --- interactive part, choose which TREE you're working with
TREE = BAC_TREE
SUMMARY = BAC_SUMMARY


## --------------- Load tree -------------------
tree_original<-read.tree(TREE)

## --------------- Load taxonomy metadata ---------------
##----- get gtdb reffs taxonomic annotation ----
gtdb_all_taxa <- read.table(GTDB_TAXA, sep='\t')
colnames(gtdb_all_taxa) <- c("genome","taxa")
reff_taxa = gtdb_all_taxa[which(gtdb_all_taxa$genome %in% tree_original$tip.label), ]

## deal with silly semicolon format
# split by semicolon
List_of_taxa_vector = strsplit(reff_taxa$taxa, ";")
# convert this strange object to dataframe : https://stackoverflow.com/questions/43662457/convert-list-of-vectors-to-data-frame
reff_taxa = cbind(reff_taxa$genome, as.data.frame(do.call(rbind, List_of_taxa_vector)))
colnames(reff_taxa) = c("genome","domain", "phylum", "class", "order", "family", "genus", "species")

##----- get taxonomic annotation from summary file-----
mag_summary = read.table(SUMMARY, header=TRUE, sep="\t")

# correct red value, first they are numeric, second N/A is 1
mag_summary$red_value[mag_summary$red_value=="N/A"]=1
mag_summary$red_value = as.numeric(mag_summary$red_value)
rownames(mag_summary) = mag_summary$user_genome
# my worst onliner in R so far, split, convert to dataframe, concatenate to the other data.frame
mag_taxa = cbind(mag_summary[,c("user_genome","red_value")],as.data.frame(do.call(rbind, strsplit(mag_summary$classification, ";"))))
colnames(mag_taxa) = c("genome","red_value","domain", "phylum", "class", "order", "family", "genus", "species")

##----- get 0.1% diversity, prevalence, contamination, completion -----
mag_metadata_all = read.table(METADATA, header=TRUE, sep="\t")
mag_metadata = mag_metadata_all[,c("MAGs","completion","contamination")]
mag_metadata$prevalence = rowSums(mag_metadata_all[,c("AD16.1","AD16.2", "AD18D", "AD18H",  "AD7", "AD20", "MINI","AD2.1", "AD2.2", "AD10", "AD11", "AD12", "AD3")]!=0)
mag_metadata$diversity = mag_metadata_all$nb_mags_0.001

## --------------- Tree! ---------------
# use de novo tree topology to reannotate mag and redefine novelty
result = node_annotation(tree_original,reff_taxa,mag_taxa,mag_metadata)

# tree annotated has inner node renamed in a unique way, allowing to make correspond node annotation
tree_annotated = result[["tree"]]
node_annotation = result[["annotation"]] 

# select novelty and plot subgraphs
novel_phylum = node_annotation[is.na(node_annotation[,"phylum"])&(node_annotation$label %in% mag_taxa$genome),"label"]
novel_class = node_annotation[is.na(node_annotation[,"class"])&(node_annotation$label %in% mag_taxa$genome),"label"]
novel_order = node_annotation[is.na(node_annotation[,"order"])&(node_annotation$label %in% mag_taxa$genome),"label"]

# for each mag to plot, find the 1rst parent node including a gtdb node and plot
Mag_to_plot = c(novel_phylum,novel_class)
mag_ploted = c()
pdf("temporary_trees.pdf")
for (mag in Mag_to_plot)
{
  if (!(mag %in% mag_ploted))
  {
    mags = c(mag)
    nb = 1
    while (sum(mags %in% reff_taxa$genome)==0)
    {
      clade_tree<-tree_subset(tree_annotated,mag,levels_back=nb)
      nb = nb + 1
      mags = clade_tree$tip.label
    }
    clade_tree2<-tree_subset(tree_annotated,mag,levels_back=nb)
    if (length(clade_tree2$tip.label)<=50)
    {
      clade_tree = clade_tree2
    }
    plot_tree(clade_tree,node_annotation)
    mag_ploted = c(mag_ploted,mags)
  }
} 
dev.off()

plot_tree = function(clade_tree,node_annotation)
{
  clade_tree_tbl = as_tibble(clade_tree)
  clade_data = node_annotation[match(clade_tree_tbl$label,node_annotation$label),]
  clade_data$node = clade_tree_tbl$node
  rownames(clade_data)=clade_data$label
  # assess level of novelty
  root_node = MRCA(clade_tree,clade_data$label[!(is.na(clade_data$contamination))])
  taxa = clade_data[root_node,][1:7]
  min_determined = max(which(!(is.na(taxa))))+1
  taxa_col = c("domain", "phylum", "class", "order", "family", "genus", "species")
  if (min_determined>7){min_determined=7}
  color = taxa_col[min_determined]
  clade_data$color = clade_data[,color]

  p = ggtree(clade_tree,aes(color=class)) %<+% clade_data
  p = p + geom_tiplab(align = TRUE,offset=0.1*max(p$data$x)) + theme_tree2(axis.text.x = element_text(angle = 45, hjust = 1))
  p = p + aes(color=class)
  offset = max(p$data$x)
  width = 0.5
  p2 = gheatmap2(p,clade_data[,c("completion","contamination","diversity","prevalence")],offset=offset,width=width)
  p2 = p2 + ggtitle(paste("new",color,"diversity : ",paste(taxa[1:min_determined],collapse=";")))
  print(p2)
}


p2 = gheatmap2(p,clade_data[,c("completion","contamination")],offset=offset,width=width)+new_scale_fill()
p3 = gheatmap2(p2,clade_data[,c("prevalence","diversity")],offset=offset+0.8*width,width=width)+new_scale_fill()

p3 = p3 + scale_x_ggtree(breaks,labels)

p4 = gheatmap2(p3,clade_data[,],offset=offset+2*width,width=width)+ggnewscale()


tmp = clade_data[clade_data$label %in% clade_tree$tip.label,c("contamination","completion")]


pTree<-gheatmap(p,tmp , offset=0.05, width=0.03, colnames = FALSE) 


+scale_fill_viridis_c(option="D", direction=-1, na.value="white")


p <- ggtree(tree_annotated, branch.length="branch length") %<+% node_annotation


p1 <- gheatmap(p, nSpeciesInfo, offset=0, width=0.05, colnames = FALSE) +
  scale_fill_viridis_c(breaks=my_breaks, labels=my_breaks,
                       trans = scales::pseudo_log_trans(sigma = 1),
                       option="B", direction=-1, na.value="white",
                       name="nSpecies_MAGs")



p2 <- p + new_scale_fill()
# nSpeciesInfo <- data.frame(tableToPlot$red_value) #RED
# rownames(nSpeciesInfo) <- tableToPlot[,1]
# colnames(nSpeciesInfo) <- c("RED")
pTreeFull<-groupOTU(p2, clade_data$label[clade_data$label %in% mag_taxa$genome], 'Phylum') +
aes(color=phylum)

"p__OLB16"

### read gtdb taxa names and info
##Archea gtdb taxonomy
#gtdbFile <- readr::read_tsv("~/orkun/projects/ongoing/AD_monitoring_project/analysisBySeb/tree_of_MAGs_gtdb_June2020/gtdb/ar122_taxonomy.tsv", col_names=FALSE, 
               #             col_types = "cc");

##Bacteria gtdb taxonomy

### read original tree file into a tree object and read tips' information file
##Archea
#tree_original<-read.tree("~/orkun/projects/ongoing/AD_monitoring_project/analysisBySeb/tree_of_MAGs_gtdb_June2020/fasttree_familly_tree_ar122.nwk")
#mag_taxa <- readr::read_tsv("~/orkun/projects/ongoing/AD_monitoring_project/analysisBySeb/tree_of_MAGs_gtdb_June2020/diversitySet_gtdb_v95_classify/gtdbtk.ar122.summary.tsv", col_names=TRUE, 
                #                    col_types = "cccdcddccddccccdidc");

###Bacteria tree - GTDB



#OPTIONAL ---- PICK SUB-CLADE OF THE ORIGINAL TREE
###WORKS WITH PLOT 1 BELOW
#going back 5 levels from a234_m20724
tree_collapsed_clade<-tree_subset(tree_original,'a234_m20724',levels_back=1) 
##

#re-name tree file- historical remnant - so that rest of code works
# tree_collapsed_clade<-tree_original

#get leaf info from tree, which is for user genome only, and add other leafs from gtdb file
Leafs_clade <- as.vector(tree_collapsed_clade$tip.label)
myTable_collapsed_clade <- mag_taxa[mag_taxa$gtdbID %in% Leafs_clade,]
selected_reff_taxa = reff_taxa[reff_taxa$gtdbID %in% tree_collapsed_clade$tip.label,]
selected_reff_taxa$red_value = NA
myTable_collapsed_clade <- rbind(selected_reff_taxa,myTable_collapsed_clade)


#OPTIONAL ---- PICK ONLY THE MAGS FROM THE ORIGINAL TREE AND DROP ANY GTDB SPECIES
tree_drop_clades<-drop.tip(tree_collapsed_clade,reff_taxa$species)  #drop GTDB species
a<-myTable_collapsed_clade$user_genome[which(is.na(myTable_collapsed_clade$red_value))] #find MAGs with RED=NA
tree_drop_clades2<-drop.tip(tree_drop_clades,a) #further drop MAGs with RED=NA
myTable_drop_clade <- myTable_collapsed_clade[which(myTable_collapsed_clade$user_genome %in% tree_drop_clades2$tip.label), ]

#get files to pass to plot
treeToPlot <- tree_collapsed_clade
tableToPlot <- myTable_collapsed_clade
# treeToPlot <- tree_drop_clades2
# tableToPlot <- myTable_drop_clade

## PLOTTING PART ## 

#PLOT 1 -
#Tree PLOT - Whole bacteria tree
#group leafs for coloring purposes at level of choice
leafList_clade <- split(tableToPlot$gtdbID,tableToPlot$class)  #use family for coloring clade
leafListToPlot <- leafList_clade

#first create layers of info to be added as heatmap to the tree plot
nSpeciesInfo <- data.frame(tableToPlot$red_value) #RED
rownames(nSpeciesInfo) <- tableToPlot[,1]
colnames(nSpeciesInfo) <- c("RED")

#create tree object
p <- ggtree(treeToPlot, branch.length="branch length") %<+% tableToPlot

#add median_RED info for leaves onto same tree
p2 <- p + new_scale_fill()
pTree<-gheatmap(p2, nSpeciesInfo, offset=0.05, width=0.03, colnames = FALSE) +scale_fill_viridis_c(option="D", direction=-1, na.value="white", name="RED")

test = p2+ geom_hilight(aes(node=treeToPlot$tip.label))



pTreeFull<-groupOTU(pTree, leafListToPlot, 'Phylum') +
  aes(color=order) +             #use with tree_collapsed_clade
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




