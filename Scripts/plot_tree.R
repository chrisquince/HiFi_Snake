library(ggtree)
library(ape)
library(treeio)
library(stringr)

#https://yulab-smu.github.io/treedata-book/chapter2.html#tidytree
# first create a taxonomic annotation of all nodes

### read original tree file into a tree object and read tips' information file
# load tree
tree_original<-read.tree("/mnt/gpfs/seb/Database/RefSeq_07_2020/KO_trees/K00925/K00925.tree")

# get metadata
taxa = read.table("/mnt/gpfs/seb/Database/RefSeq_07_2020/KO_trees/K00925/K00925_leaves_to_lineage.tsv",sep="\t")
taxa <- data.frame(lapply(taxa, function(x) {gsub(".__", "", x)}))
colnames(taxa) = c("tip_name",'domain',"phylum","class","order","family","genus","species")
rownames(taxa) = taxa$tip_name
archea = taxa[grep("Archaea", taxa[,2]),]

#get mag nodes of interest
tip_archea = archea[,1]
mag_archea = tip_archea[grep("a.+_m.+",tip_archea)]


########## first tree Methanosarcina #########################################
# there are 10 archea mags, let's find the minimum number of tree encompassing them
archea_tree1 = tree_subset(tree_original, mag_archea[1], levels_back = 7)
p = ggtree(archea_tree1,aes(color=class),size=2) %<+% taxa + geom_tiplab(align=TRUE)
is_archea = data.frame(taxa[archea_tree1$tip.label,]$domain)
rownames(is_archea)=archea_tree1$tip.label
p1 = gheatmap(p, is_archea,offset=0.8*max(p$data$x), width=0.1, colnames = FALSE)

########## second tree iainarchaeia/Nanoarchaeia #########################################
# still 6 archea mags 
left_to_plot = mag_archea[!(mag_archea %in% archea_tree1$tip.label)] 
archea_tree2 = tree_subset(tree_original, left_to_plot[1], levels_back = 3)
p = ggtree(archea_tree2,aes(color=class),size=2) %<+% taxa
p = ggtree(archea_tree2,aes(color=class),size=2) %<+% taxa + geom_tiplab(align=TRUE,offset=.05*max(p$data$x))
is_archea = data.frame(taxa[archea_tree2$tip.label,]$domain)
rownames(is_archea)=archea_tree2$tip.label
p2 = gheatmap(p, is_archea,offset=0.8*max(p$data$x), width=0.1, colnames = FALSE)
# still 6 archea mags 
archea_tree2 = tree_subset(tree_original, left_to_plot[1], levels_back = 2)
p = ggtree(archea_tree2,aes(color=class),size=2) %<+% taxa
p = ggtree(archea_tree2,aes(color=class),size=2) %<+% taxa + geom_tiplab(align=TRUE,offset=.05*max(p$data$x))
is_archea = data.frame(taxa[archea_tree2$tip.label,]$domain)
rownames(is_archea)=archea_tree2$tip.label
p2_bis = gheatmap(p, is_archea,offset=0.8*max(p$data$x), width=0.1, colnames = FALSE)
########## 3rd tree other Nanoarchaeota #########################################
# still 2 archea mags 
left_to_plot = left_to_plot[!(left_to_plot %in% archea_tree2$tip.label)] 
archea_tree3 = tree_subset(tree_original, left_to_plot[1], levels_back = 3)
p = ggtree(archea_tree3,aes(color=class),size=2) %<+% taxa
p = ggtree(archea_tree3,aes(color=class),size=2) %<+% taxa + geom_tiplab(align=TRUE,offset=.05*max(p$data$x))
is_archea = data.frame(taxa[archea_tree3$tip.label,]$domain)
rownames(is_archea)=archea_tree3$tip.label
p3 = gheatmap(p, is_archea,offset=0.8*max(p$data$x), width=0.1, colnames = FALSE)
########## last tree iainarchaeia/Nanoarchaeia #########################################
# still 1 archea mags 
left_to_plot = left_to_plot[!(left_to_plot %in% archea_tree3$tip.label)] 
archea_tree4 = tree_subset(tree_original, left_to_plot[1], levels_back = 3)
p = ggtree(archea_tree4,aes(color=class),size=2) %<+% taxa
p = ggtree(archea_tree4,aes(color=class),size=2) %<+% taxa + geom_tiplab(align=TRUE,offset=.05*max(p$data$x))
is_archea = data.frame(taxa[archea_tree4$tip.label,]$domain)
rownames(is_archea)=archea_tree4$tip.label
p4 = gheatmap(p, is_archea,offset=0.8*max(p$data$x), width=0.1, colnames = FALSE)

pdf("acetate_kinase_tree.pdf")
print(p1)
print(p2)
print(p2_bis)
print(p3)
print(p4)
dev.off()
