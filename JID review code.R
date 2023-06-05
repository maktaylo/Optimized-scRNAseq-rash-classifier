#JID Case review code

library(concaveman)
library(ggforce)
library(dplyr)
library(Seurat)
library(patchwork)
library(reshape2)
library(gplots)
require(gtools)
library(psych)
library(RColorBrewer)
library(mclust)
library(factoextra)
library(ggrepel)
library(GGally)
library(Hmisc)
library(corrplot)
library(network)
library(sna)
library(intergraph)
library(destiny)
library(ggbiplot)
library(monocle3)
library(viridis)
library(garnett)
library(org.Mm.eg.db)
library(Matrix)
library(igraph)
library(intergraph)
library(ggVennDiagram)
library(stringr)
library(topGO)
library(WGCNA)
library("openxlsx")
library(lsmeans)
library(emmeans)
library(multcomp)
library(ggcorrplot)
library(Category)
library("ggpubr")
library(scran)
library(CellChat)
library(ggplot2)                  
library(patchwork)
library(ggalluvial)
library(igraph)
library(dplyr)
options(stringsAsFactors = FALSE)
library(limma)
library(xgboost)
library(caTools)
library(raster)

options(future.globals.maxSize = 20000 * 1024^2) #CC uses this to bypass vector memory limits
#plan("multiprocess", workers = 4)

`%!in%` <- Negate(`%in%`) #Define this operator

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$Hours)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

format_summary = function(x) {
  output = paste(round(x[1],0)," (",round(x[2],0),", ",round(x[3],0),")",sep="")
  return(output)
}

cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

#summary function
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length(na.omit(xx[[col]])),
                     mean = mean   (na.omit(xx[[col]])),
                     sd   = sd     (na.omit(xx[[col]]))
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# Basic function to convert human to mouse gene names: https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

#------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------#

#Read in this object
humandat = readRDS("/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/data/Manuscript\ Object.rds")

#Seurat object --> CDS
exp_mat <- humandat@assays[["RNA"]]@data #pull out NORMALIZED counts from Seurat object
cell_metadat <- humandat@meta.data #pull out cell meta-data from Seurat object
gene_annot = data.frame(humandat@assays[["RNA"]]@counts@Dimnames[[1]])#pull out gene names from Seurat object
names(gene_annot) = "gene_short_name"
row.names(gene_annot) = gene_annot$gene_short_name #row.names of gene_metadata must be equal to row.names of expression_data

#Check various metadata factors
data.frame(names(cell_metadat))
data.frame(table(cell_metadat$ident))
data.frame(table(cell_metadat$donor))
data.frame(table(cell_metadat$dis))
data.frame(table(cell_metadat$treat))
data.frame(table(cell_metadat$chem))
data.frame(table(cell_metadat$version))
data.frame(table(cell_metadat$Site))
data.frame(table(cell_metadat$final_clustering))
data.frame(table(cell_metadat$clusters))
data.frame(table(cell_metadat$ID))
data.frame(table(cell_metadat$type))
data.frame(table(cell_metadat$Ident1)) #same as Ident2
data.frame(table(cell_metadat$Ident2))

#------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------#
#Make CDS object
human_human_big_cds <- new_cell_data_set(exp_mat,
                                         cell_metadata = cell_metadat,
                                         gene_metadata = gene_annot)

#------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------#

#create 1 giant ass cds object from all the cds's created above
rm(humandat,exp_mat)

#Subset to patients of interest
pts = c("170","198","230","231","232","233","236","165","173","194","199","211","222","234","235","202")
human_human_big_cds = human_human_big_cds[,colData(human_human_big_cds)$donor3 %in% pts] #subset to all normal patients and a diseased patient
colData(human_human_big_cds)$dis_updated = colData(human_human_big_cds)$donor3
colData(human_human_big_cds)$dis_updated = gsub("170","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("198","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("230","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("231","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("232","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("233","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("236","AD",colData(human_human_big_cds)$dis_updated)

colData(human_human_big_cds)$dis_updated = gsub("165","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("173","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("194","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("199","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("211","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("222","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("234","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("235","PV",colData(human_human_big_cds)$dis_updated)

colData(human_human_big_cds)$dis_updated = gsub("202","PV",colData(human_human_big_cds)$dis_updated) #previously unknown sample

#--------------------------------------------------------------------#
#--------------------------------------------------------------------#
#--------------------------------------------------------------------#
#--------------------------------------------------------------------#
#--------------------------------------------------------------------#
#Optimize manifold by looping through min_distances

human_human_big_cds <- preprocess_cds(human_human_big_cds, 
                             method="PCA",
                             num_dim = 100,
                             norm_method = "none") #these data have already been normalized by Seurat


min_distances = seq(0.05,1,by=0.05) #This is showing a lot of interesting bridgeness with min distanes bewteen 0.4 and 0.6

pdf(file=paste("/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/results/manifold_optimization/figs/min_distances_v1.pdf",sep=""),height=10,width=10)
for (i in 1:length(min_distances)) {
  #i=10
  print(paste("i=",i,sep=""))
  
  #4b. ?reduce_dimension 
  human_human_big_cds <- reduce_dimension(
    human_human_big_cds,
    max_components = 2,
    reduction_method = c("UMAP"),
    preprocess_method = "PCA", #must now used ALIGNEMNT PRE-PROCESSING METHOD bc it replaces the other pre-process shit
    umap.metric = "cosine",
    umap.min_dist = min_distances[i],
    umap.n_neighbors = 10L,
    umap.fast_sgd = FALSE,
    umap.nn_method = "annoy",
    cores = 1,
    verbose = FALSE
  )
  
  #5. Cluster and partition ?cluster_cells cells into supergroups community detection (Phenograph algorithm) Community detection of cells as part of Phenograph algorithm https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4508757/
  #5b.  ?cluster_cells
  human_human_big_cds <- cluster_cells(human_human_big_cds,
                                   reduction_method = c("UMAP"),
                                   k = 10, #a bigger k will result in lower resolution 
                                   cluster_method = c("leiden"),
                                   num_iter = 2,
                                   partition_qval = 0.05, #Numeric, the q-value cutoff to determine when to partition
                                   weight = FALSE,
                                   resolution = NULL,
                                   random_seed = NULL,
                                   verbose = F)
  print( 
    plot_cells(human_human_big_cds, color_cells_by="dis_updated", group_cells_by="cluster",
               group_label_size = 0,
               show_trajectory_graph=F,
               x=1,
               y=2,
               alpha=0.5,
               cell_size=0.4) +
      ggtitle(paste(" Minimal distance = ",min_distances[i],sep="")) + xlab("") + ylab("") +
      theme_classic() +
      theme(axis.text=element_text(size=15,color='black'),
            title=element_text(size=15),
            axis.ticks = element_line(size=0),
            legend.position="bottom") +
      scale_color_manual(values = c("AD" = "darkgoldenrod3",
                                    "PV" = "olivedrab4"))
     
  )
}
dev.off()
#--------------------------------------------------------------------#
#--------------------------------------------------------------------#
#--------------------------------------------------------------------#
#Optimal manifold, write down sample means, disease means

#optimal min distance is 0.45
human_human_big_cds <- reduce_dimension(
  human_human_big_cds,
  max_components = 2,
  reduction_method = c("UMAP"),
  preprocess_method = "PCA", #must now used ALIGNEMNT PRE-PROCESSING METHOD bc it replaces the other pre-process shit
  umap.metric = "cosine",
  umap.min_dist = 0.45, #update optimal min distance value here
  umap.n_neighbors = 10L,
  umap.fast_sgd = FALSE,
  umap.nn_method = "annoy",
  cores = 1,
  verbose = FALSE
)

#5. Cluster and partition ?cluster_cells cells into supergroups community detection (Phenograph algorithm) Community detection of cells as part of Phenograph algorithm https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4508757/
#5b.  ?cluster_cells
human_human_big_cds <- cluster_cells(human_human_big_cds,
                                     reduction_method = c("UMAP"),
                                     k = 10, #a bigger k will result in lower resolution 
                                     cluster_method = c("leiden"),
                                     num_iter = 2,
                                     partition_qval = 0.05, #Numeric, the q-value cutoff to determine when to partition
                                     weight = FALSE,
                                     resolution = NULL,
                                     random_seed = NULL,
                                     verbose = F)

#Now make new dataframe of 
vars_of_interest = c("ident","donor","dis","treat","chem","final_clustering","clusters","ID","type","Ident1","Ident2","dis_updated")
allcell_allgene = data.frame(human_human_big_cds@int_colData@listData$reducedDims@listData$UMAP)
allcell_allgene = data.frame(allcell_allgene, colData(human_human_big_cds)[,match(vars_of_interest,names(colData(human_human_big_cds)))]) #combine UMAP coords with metadata factors of interest
names(allcell_allgene)[c(1,2)] = c("pc1","pc2")

#Get average across pc1 and pc2
#UNLOAD DPLYR PACKAGE#
pc1_dfc= summarySE(allcell_allgene, measurevar="pc1", groupvars=c("dis_updated","ident"))
names(pc1_dfc)[c(5)] = c("pc1_se")
pc2_dfc= summarySE(allcell_allgene, measurevar="pc2", groupvars=c("dis_updated","ident"))
names(pc2_dfc)[c(5)] = c("pc2_se")

#--------------------------------------------------------------------#
#Now make hull plot: https://www.r-bloggers.com/2021/07/ggforce-make-a-hull-plot-to-visualize-clusters-in-ggplot2/
re_int = cbind(pc1_dfc[,c(1,2,4,5)], pc2_dfc[,c(4,5)]) #combine the gene signatures

#--------------------------------------------------------------------#
#Calculate disease centers
ad_centroid = data.frame( mean(subset(re_int, dis_updated=="AD")$pc1), mean(subset(re_int, dis_updated=="AD")$pc2) )
pv_centroid = data.frame( mean(subset(re_int, dis_updated=="PV")$pc1), mean(subset(re_int, dis_updated=="PV")$pc2) )
names(ad_centroid) = c("xpt","ypt")
names(pv_centroid) = c("xpt","ypt")

#--------------------------------------------------------------------#
#Plot everything together: cells, sample centroids, disease centroids
allcells_allgenes_fig = ggplot() +
  geom_point(data=allcell_allgene, alpha=0.05,aes(x=pc1, y=pc2, color=dis_updated), size=0.4) + #individual cells
  geom_point(data=re_int, alpha=1,aes(x=pc1, y=pc2, color=dis_updated, shape=dis_updated), size=4) + #sample centroids
  #geom_errorbarh(data=re_int, aes(x=pc1,xmax = pc1 + pc1_se, xmin = pc1 - pc1_se, color=dis_updated)) +
  #geom_errorbar(data=re_int, aes(y=pc2,ymax = pc2 + pc2_se, ymin = pc2 - pc2_se, color=dis_updated)) +
  stat_ellipse(data=allcell_allgene, aes(x=pc1, y=pc2, color=dis_updated)) +
  theme_classic() +
  ggtitle("All cells, full transcriptome") + 
  xlab("Dim 1") +
  ylab("Dim 2") +
  theme(legend.position="right") +
  scale_color_manual(name="",values = c("AD" = "darkgoldenrod3",
                                        "PV" = "olivedrab4",
                                        "Atypical rash case" = "darkblue")) +
  scale_fill_manual(name="",values = c("AD" = "lightgoldenrod2",
                                       "PV" = "darkolivegreen4",
                                       "Atypical rash case" = "darkolivegreen4")) +
  geom_point(data=ad_centroid, aes(x=xpt,y=ypt), color="darkgoldenrod3", shape=10, size=7, stroke = 2) +
  geom_point(data=pv_centroid, aes(x=xpt,y=ypt), color="darkolivegreen4", shape=10, size=7, stroke =2) +
  theme(axis.text=element_text(size=20,color='black'),
        axis.title=element_text(size=20,color='black'),
        plot.title=element_text(size=20,color='black')) 


png("/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/results/cluster_summary_fig/allcells_allgenes_v1.png",
    width=8,height=8,units="in",res=1000)
  print(allcells_allgenes_fig)
dev.off()  

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#Now the same thing but with cluster 2 cells only

clust2_cds = human_human_big_cds[,colData(human_human_big_cds)$ID %in% "2"] #subset to all normal patients and a diseased patient

clust2_cds <- preprocess_cds(clust2_cds, 
                                      method="PCA",
                                      num_dim = 100,
                                      norm_method = "none") #these data have already been normalized by Seurat

#optimal min distance is 0.45
clust2_cds <- reduce_dimension(
  clust2_cds,
  max_components = 2,
  reduction_method = c("UMAP"),
  preprocess_method = "PCA", #must now used ALIGNEMNT PRE-PROCESSING METHOD bc it replaces the other pre-process shit
  umap.metric = "cosine",
  umap.min_dist = 0.45, #update optimal min distance value here
  umap.n_neighbors = 30L,
  umap.fast_sgd = FALSE,
  umap.nn_method = "annoy",
  cores = 1,
  verbose = FALSE
)

#5. Cluster and partition ?cluster_cells cells into supergroups community detection (Phenograph algorithm) Community detection of cells as part of Phenograph algorithm https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4508757/
#5b.  ?cluster_cells
clust2_cds <- cluster_cells(clust2_cds,
                            reduction_method = c("UMAP"),
                            k = 10, #a bigger k will result in lower resolution 
                            cluster_method = c("leiden"),
                            num_iter = 2,
                            partition_qval = 0.05, #Numeric, the q-value cutoff to determine when to partition
                            weight = FALSE,
                            resolution = NULL,
                            random_seed = NULL,
                            verbose = F)

#Now make new dataframe of 
vars_of_interest = c("ident","donor","dis","treat","chem","final_clustering","clusters","ID","type","Ident1","Ident2","dis_updated")
allcell_allgene = data.frame(clust2_cds@int_colData@listData$reducedDims@listData$UMAP)
allcell_allgene = data.frame(allcell_allgene, colData(clust2_cds)[,match(vars_of_interest,names(colData(clust2_cds)))]) #combine UMAP coords with metadata factors of interest
names(allcell_allgene)[c(1,2)] = c("pc1","pc2")

#Get average across pc1 and pc2
#UNLOAD DPLYR PACKAGE#
pc1_dfc= summarySE(allcell_allgene, measurevar="pc1", groupvars=c("dis_updated","ident"))
names(pc1_dfc)[c(5)] = c("pc1_se")
pc2_dfc= summarySE(allcell_allgene, measurevar="pc2", groupvars=c("dis_updated","ident"))
names(pc2_dfc)[c(5)] = c("pc2_se")

#--------------------------------------------------------------------#
#Now make hull plot: https://www.r-bloggers.com/2021/07/ggforce-make-a-hull-plot-to-visualize-clusters-in-ggplot2/
re_int = cbind(pc1_dfc[,c(1,2,4,5)], pc2_dfc[,c(4,5)]) #combine the gene signatures

#--------------------------------------------------------------------#
#Calculate disease centers
ad_centroid = data.frame( mean(subset(re_int, dis_updated=="AD")$pc1), mean(subset(re_int, dis_updated=="AD")$pc2) )
pv_centroid = data.frame( mean(subset(re_int, dis_updated=="PV")$pc1), mean(subset(re_int, dis_updated=="PV")$pc2) )
names(ad_centroid) = c("xpt","ypt")
names(pv_centroid) = c("xpt","ypt")

#--------------------------------------------------------------------#
#Plot everything together: cells, sample centroids, disease centroids
clust2cells_allgenes_fig = ggplot() +
  geom_point(data=allcell_allgene, alpha=0.15,aes(x=pc1, y=pc2, color=dis_updated), size=0.4) + #individual cells
  geom_point(data=re_int, alpha=1,aes(x=pc1, y=pc2, color=dis_updated, shape=dis_updated), size=4) + #sample centroids
  #geom_errorbarh(data=re_int, aes(x=pc1,xmax = pc1 + pc1_se, xmin = pc1 - pc1_se, color=dis_updated)) +
  #geom_errorbar(data=re_int, aes(y=pc2,ymax = pc2 + pc2_se, ymin = pc2 - pc2_se, color=dis_updated)) +
  stat_ellipse(data=allcell_allgene, aes(x=pc1, y=pc2, color=dis_updated)) +
  theme_classic() +
  ggtitle("Trm cells, full transcriptome") + 
  xlab("Dim 1") +
  ylab("Dim 2") +
  theme(legend.position="right") +
  scale_color_manual(name="",values = c("AD" = "darkgoldenrod3",
                                        "PV" = "olivedrab4",
                                        "Atypical rash case" = "darkblue")) +
  scale_fill_manual(name="",values = c("AD" = "lightgoldenrod2",
                                       "PV" = "darkolivegreen4",
                                       "Atypical rash case" = "darkolivegreen4")) +
  geom_point(data=ad_centroid, aes(x=xpt,y=ypt), color="darkgoldenrod3", shape=10, size=7, stroke = 2) +
  geom_point(data=pv_centroid, aes(x=xpt,y=ypt), color="darkolivegreen4", shape=10, size=7, stroke =2) +
  theme(axis.text=element_text(size=20,color='black'),
        axis.title=element_text(size=20,color='black'),
        plot.title=element_text(size=20,color='black')) 


png("/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/results/cluster_summary_fig/clust2cells_allgenes_v2.png",
    width=8,height=8,units="in",res=1000)
print(clust2cells_allgenes_fig)
dev.off()  

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#Try UMAP decomposotion of just the genes from the RashX algo's genes

ad_genes = c("TWIST1","LGALS1",    "IL32",   "CAPG",  "ITM2C", "MFHAS1", "ANXA1","SOS1", "CSGALNACT1","LMO4",  "IFITM2","S100A10",  "MT-ND5",  "CYSLTR1","PLA2G16", "SYNE2", "THADA",    "NEAT1","IL17RB", "RPL36A","ARHGAP21",    "NBAS",  "ACTG1","PRKX", "TGFBR3",   "TIMP1","TNFSF10", "AHNAK",    "MT-ND2",  "ISG15",  "RPL17",  "LONRF2",   "CD99","TSHZ2", "MMP25",   "IFITM1","MT-ND1",  "BIRC3",  "FAM102A", "LPCAT2","NRIP3", "CRIP1",  "CLU",   "PLP2", "ZFP36",  "ZFP36L2",  "TUBA1B","GATA3","SLC5A3",    "SFXN1", "FANK1","TAGLN2") 
pv_genes = c("CXCL13", "MTRNR2L12",  "CD7",   "MGAT4A","FTH1", "LAYN", "IL17F", "KLRB1", "GNLY","CPM", "CTSH", "GBP5", "SOX4", "CLEC2B",    "GZMB", "CD2",   "CEBPD","ODF2L","LAG3", "LRRN3", "ARHGEF12",  "PTPN13",   "TNFAIP3", "TRPS1", "SNX9", "METRNL",  "BTG1", "JUN","SPOCK2",   "GABARAPL1",  "PMEPA1", "HIST1H1E","RBPJ", "LINC01871", "MAP3K4","H1FX", "UBC",  "GALNT1",  "PNRC1","GABPB1-AS1", "RPS26", "MUC20-OT1",  "CHN1",  "NAP1L4",   "PTMS",  "F2R",   "CTLA4", "DAPK2","RAP1B","CCR6", "B3GALT2", "YPEL2",  "FYN",   "PPDPF","SLA2", "CBLB", "ADGRG1","SARAF") 
clust2_cds@assays@data$counts[1:10,1:10]

#Make gene_metadata file of reduced gene set
genes_df = data.frame(c(ad_genes,pv_genes))
names(genes_df) = "gene_short_name"
row.names(genes_df) = genes_df$gene_short_name

#Create new cds object with reduced genes
clust2_cds_reduced_genes <- new_cell_data_set(clust2_cds@assays@data$counts[c(ad_genes,pv_genes),],
                                                cell_metadata = colData(clust2_cds),
                                                gene_metadata = genes_df)

# Try UMAP decomposition of just these genes
clust2_cds_reduced_genes <- preprocess_cds(clust2_cds_reduced_genes, 
                             method="PCA",
                             num_dim = 100,
                             norm_method = "none") #these data have already been normalized by Seurat


#optimal min distance is 0.45
clust2_cds_reduced_genes <- reduce_dimension(
  clust2_cds_reduced_genes,
  max_components = 2,
  reduction_method = c("UMAP"),
  preprocess_method = "PCA", #must now used ALIGNEMNT PRE-PROCESSING METHOD bc it replaces the other pre-process shit
  umap.metric = "cosine",
  umap.min_dist = 0.45, #update optimal min distance value here
  umap.n_neighbors = 30L,
  umap.fast_sgd = FALSE,
  umap.nn_method = "annoy",
  cores = 1,
  verbose = FALSE
)

#5. Cluster and partition ?cluster_cells cells into supergroups community detection (Phenograph algorithm) Community detection of cells as part of Phenograph algorithm https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4508757/
#5b.  ?cluster_cells
clust2_cds_reduced_genes <- cluster_cells(clust2_cds_reduced_genes,
                            reduction_method = c("UMAP"),
                            k = 10, #a bigger k will result in lower resolution 
                            cluster_method = c("leiden"),
                            num_iter = 2,
                            partition_qval = 0.05, #Numeric, the q-value cutoff to determine when to partition
                            weight = FALSE,
                            resolution = NULL,
                            random_seed = NULL,
                            verbose = F)

#Now make new dataframe of 
vars_of_interest = c("ident","donor","dis","treat","chem","final_clustering","clusters","ID","type","Ident1","Ident2","dis_updated")
allcell_allgene = data.frame(clust2_cds_reduced_genes@int_colData@listData$reducedDims@listData$UMAP)
allcell_allgene = data.frame(allcell_allgene, colData(clust2_cds_reduced_genes)[,match(vars_of_interest,names(colData(clust2_cds_reduced_genes)))]) #combine UMAP coords with metadata factors of interest
names(allcell_allgene)[c(1,2)] = c("pc1","pc2")

#Get average across pc1 and pc2
#UNLOAD DPLYR PACKAGE#
pc1_dfc= summarySE(allcell_allgene, measurevar="pc1", groupvars=c("dis_updated","ident"))
names(pc1_dfc)[c(5)] = c("pc1_se")
pc2_dfc= summarySE(allcell_allgene, measurevar="pc2", groupvars=c("dis_updated","ident"))
names(pc2_dfc)[c(5)] = c("pc2_se")

#--------------------------------------------------------------------#
#Now make hull plot: https://www.r-bloggers.com/2021/07/ggforce-make-a-hull-plot-to-visualize-clusters-in-ggplot2/
re_int = cbind(pc1_dfc[,c(1,2,4,5)], pc2_dfc[,c(4,5)]) #combine the gene signatures

#--------------------------------------------------------------------#
#Calculate disease centers
ad_centroid = data.frame( mean(subset(re_int, dis_updated=="AD")$pc1), mean(subset(re_int, dis_updated=="AD")$pc2) )
pv_centroid = data.frame( mean(subset(re_int, dis_updated=="PV")$pc1), mean(subset(re_int, dis_updated=="PV")$pc2) )
names(ad_centroid) = c("xpt","ypt")
names(pv_centroid) = c("xpt","ypt")

#--------------------------------------------------------------------#
#Plot everything together: cells, sample centroids, disease centroids
clust2cells_allgenes_fig = ggplot() +
  geom_point(data=allcell_allgene, alpha=0.15,aes(x=pc1, y=pc2, color=dis_updated), size=0.4) + #individual cells
  geom_point(data=re_int, alpha=1,aes(x=pc1, y=pc2, color=dis_updated, shape=dis_updated), size=4) + #sample centroids
  #geom_errorbarh(data=re_int, aes(x=pc1,xmax = pc1 + pc1_se, xmin = pc1 - pc1_se, color=dis_updated)) +
  #geom_errorbar(data=re_int, aes(y=pc2,ymax = pc2 + pc2_se, ymin = pc2 - pc2_se, color=dis_updated)) +
  stat_ellipse(data=allcell_allgene, aes(x=pc1, y=pc2, color=dis_updated)) +
  theme_classic() +
  ggtitle("Trm cells, RashX algorithm genes") + 
  xlab("Dim 1") +
  ylab("Dim 2") +
  theme(legend.position="right") +
  scale_color_manual(name="",values = c("AD" = "darkgoldenrod3",
                                        "PV" = "olivedrab4",
                                        "Atypical rash case" = "darkblue")) +
  scale_fill_manual(name="",values = c("AD" = "lightgoldenrod2",
                                       "PV" = "darkolivegreen4",
                                       "Atypical rash case" = "darkolivegreen4")) +
  geom_point(data=ad_centroid, aes(x=xpt,y=ypt), color="darkgoldenrod3", shape=10, size=7, stroke = 2) +
  geom_point(data=pv_centroid, aes(x=xpt,y=ypt), color="darkolivegreen4", shape=10, size=7, stroke =2) +
  theme(axis.text=element_text(size=20,color='black'),
        axis.title=element_text(size=20,color='black'),
        plot.title=element_text(size=20,color='black')) 
clust2cells_allgenes_fig

png("/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/results/cluster_summary_fig/clust2cells_RashX_genes_v1.png",
    width=8,height=8,units="in",res=1000)
print(clust2cells_allgenes_fig)
dev.off() 


#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#Plot 
#Get the genes of interest: from the original RashX algo
  ad_genes = c("TWIST1","LGALS1",    "IL32",   "CAPG",  "ITM2C", "MFHAS1", "ANXA1","SOS1", "CSGALNACT1","LMO4",  "IFITM2","S100A10",  "MT-ND5",  "CYSLTR1","PLA2G16", "SYNE2", "THADA",    "NEAT1","IL17RB", "RPL36A","ARHGAP21",    "NBAS",  "ACTG1","PRKX", "TGFBR3",   "TIMP1","TNFSF10", "AHNAK",    "MT-ND2",  "ISG15",  "RPL17",  "LONRF2",   "CD99","TSHZ2", "MMP25",   "IFITM1","MT-ND1",  "BIRC3",  "FAM102A", "LPCAT2","NRIP3", "CRIP1",  "CLU",   "PLP2", "ZFP36",  "ZFP36L2",  "TUBA1B","GATA3","SLC5A3",    "SFXN1", "FANK1","TAGLN2") 
  pv_genes = c("CXCL13", "MTRNR2L12",  "CD7",   "MGAT4A","FTH1", "LAYN", "IL17F", "KLRB1", "GNLY","CPM", "CTSH", "GBP5", "SOX4", "CLEC2B",    "GZMB", "CD2",   "CEBPD","ODF2L","LAG3", "LRRN3", "ARHGEF12",  "PTPN13",   "TNFAIP3", "TRPS1", "SNX9", "METRNL",  "BTG1", "JUN","SPOCK2",   "GABARAPL1",  "PMEPA1", "HIST1H1E","RBPJ", "LINC01871", "MAP3K4","H1FX", "UBC",  "GALNT1",  "PNRC1","GABPB1-AS1", "RPS26", "MUC20-OT1",  "CHN1",  "NAP1L4",   "PTMS",  "F2R",   "CTLA4", "DAPK2","RAP1B","CCR6", "B3GALT2", "YPEL2",  "FYN",   "PPDPF","SLA2", "CBLB", "ADGRG1","SARAF") 
clust2_cds@assays@data$counts[1:10,1:10]


#--------------------------------------------------------------------#
#AD gene data_matrix
ad_mat = t(clust2_cds@assays@data$counts[ad_genes,])
ad_mat[1:10,1:10]
ad_int = data.frame(colData(clust2_cds)$dis_updated, colData(clust2_cds)$ident,
                    rowSums(ad_mat) )
names(ad_int) = c("dis","sample","gene_sig")

#UNLOAD DPLYR PACKAGE#
ad_dfc = summarySE(ad_int, measurevar="gene_sig", groupvars=c("dis","sample"))
head(ad_dfc)
names(ad_dfc)[c(4,6)] = c("ad_gene_sig","ad_se")

#--------------------------------------------------------------------#
#PV gene data_matrix
pv_mat = t(clust2_cds@assays@data$counts[pv_genes,])
pv_mat[1:10,1:10]
pv_int = data.frame(colData(clust2_cds)$dis_updated, colData(clust2_cds)$ident,
                    rowSums(pv_mat) )
names(pv_int) = c("dis","sample","gene_sig")
pv_dfc = summarySE(pv_int, measurevar='gene_sig', groupvars=c("dis","sample"))
head(pv_dfc)
names(pv_dfc)[c(4,6)] = c("pv_gene_sig","pv_se")

#--------------------------------------------------------------------#
#Now make hull plot: https://www.r-bloggers.com/2021/07/ggforce-make-a-hull-plot-to-visualize-clusters-in-ggplot2/
re_int = cbind(ad_dfc[,c(1,2,4, 6,7)], pv_dfc[,c(4,6)]) #combine the gene signatures

#Calculate centers
ad_centroid = data.frame( mean(subset(re_int, dis=="AD")$ad_gene_sig), mean(subset(re_int, dis=="AD")$pv_gene_sig) )
pv_centroid = data.frame( mean(subset(re_int, dis=="PV")$ad_gene_sig), mean(subset(re_int, dis=="PV")$pv_gene_sig) )
names(ad_centroid) = c("xpt","ypt")
names(pv_centroid) = c("xpt","ypt")

#--------------------------------------------------------------------#
#Plot everything together: cells, sample centroids, disease centroids

int_dat = data.frame(ad_int,pv_int)

clust2cells_geneprogram_fig = ggplot() +
  geom_point(data=int_dat, alpha=0.25,aes(x=gene_sig, y=gene_sig.1, color=dis), size=0.4) + #individual cells
  geom_point(data=re_int, alpha=1,aes(x=ad_gene_sig, y=pv_gene_sig, color=dis, shape=dis), size=4) + #sample centroids
  #geom_errorbarh(data=re_int, aes(x=pc1,xmax = pc1 + pc1_se, xmin = pc1 - pc1_se, color=dis_updated)) +
  #geom_errorbar(data=re_int, aes(y=pc2,ymax = pc2 + pc2_se, ymin = pc2 - pc2_se, color=dis_updated)) +
  #stat_ellipse(data=re_int, aes(x=ad_gene_sig, y=pv_gene_sig, color=dis)) +
  stat_ellipse(data=int_dat, aes(x=gene_sig, y=gene_sig.1, color=dis)) +
  theme_classic() +
  ggtitle("Trm cells, condition gene program") + 
  xlab("AD program") +
  ylab("PV program") +
  theme(legend.position="right") +
  scale_color_manual(name="",values = c("AD" = "darkgoldenrod3",
                                        "PV" = "olivedrab4",
                                        "Atypical rash case" = "darkblue")) +
  scale_fill_manual(name="",values = c("AD" = "lightgoldenrod2",
                                       "PV" = "darkolivegreen4",
                                       "Atypical rash case" = "darkolivegreen4")) +
  geom_point(data=ad_centroid, aes(x=xpt,y=ypt), color="darkgoldenrod3", shape=10, size=7, stroke = 2) +
  geom_point(data=pv_centroid, aes(x=xpt,y=ypt), color="darkolivegreen4", shape=10, size=7, stroke =2) +
  theme(axis.text=element_text(size=20,color='black'),
        axis.title=element_text(size=20,color='black'),
        plot.title=element_text(size=20,color='black')) 


png("/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/results/cluster_summary_fig/clust2cells_geneprogram_v1.png",
    width=8,height=8,units="in",res=1000)
print(clust2cells_geneprogram_fig)
dev.off()  


#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#Optimize on distance between disease centroids relative to sample mean range!

#read in list of potential genes: From the smaller list
degs = read.xlsx(
  xlsxFile="/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/data/clust_2_Sci_immunology.xlsx",
  colNames=T
)

#read in expaned list of potential genes: From the epxanded list
degs = read.xlsx(
  xlsxFile="/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/data/cluster_2_expanded_degs.xlsx",
  colNames=T
)
degs = subset(degs, p_val_adj<0.05) #get only significant genes
head(degs)

#Remove sex-linked genes

ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
gene_attributes = getBM(attributes=c('hgnc_symbol','chromosome_name', 'start_position', 'end_position', 'strand'),
      filters=c('hgnc_symbol'),
      values=list(degs$gene),
      mart=ensembl)

#gene_attributers is in the wrong order, so merge with original df
df2 <- merge(degs,gene_attributes,
             by.x="gene",
             by.y="hgnc_symbol",
             sort=F)
nrow(degs)
nrow(df2)
table(df2$gene)
table(df2$chromosome_name)

#Get rid of these strange chromsome names that are duplicated for some genes and look like CHR_HSCHR17_1_CTT things
df3 = df2[-grep("CHR",df2$chromosome_name),]

#Get rid of sex-linked and mitochondrial genes
df3 = df3[-grep("X",df3$chromosome_name),]
df3 = df3[-grep("Y",df3$chromosome_name),]
df3 = df3[-grep("MT",df3$chromosome_name),]


#Get just pv and ad genes
head(degs)
pv_genes = subset(df3, dis=="PV")
ad_genes = subset(df3, dis=="AD")
ad_genes = ad_genes[order(-ad_genes$avg_log2FC),] #re-order from high to low FC column name in bigger DEGs

performance_met_df = do.call(rbind, lapply(1:nrow(pv_genes), function(i){
  tryCatch({ #suppress error
  #i=5
  print(i)
  
  #--------------------------------------------------------------------#
  #AD gene data_matrix
  ad_mat = t(clust2_cds@assays@data$counts[ad_genes$gene[1:i],])
  ad_int = data.frame(colData(clust2_cds)$dis_updated, colData(clust2_cds)$ident,
                      rowSums(ad_mat) )
  names(ad_int) = c("dis","sample","gene_sig")
  
  #UNLOAD DPLYR PACKAGE#
  ad_dfc = summarySE(ad_int, measurevar="gene_sig", groupvars=c("dis","sample"))
  names(ad_dfc)[c(4,6)] = c("ad_gene_sig","ad_se")
  
  #--------------------------------------------------------------------#
  #PV gene data_matrix
  pv_mat = t(clust2_cds@assays@data$counts[pv_genes$gene[1:i],])
  pv_int = data.frame(colData(clust2_cds)$dis_updated, colData(clust2_cds)$ident,
                      rowSums(pv_mat) )
  names(pv_int) = c("dis","sample","gene_sig")
  pv_dfc = summarySE(pv_int, measurevar='gene_sig', groupvars=c("dis","sample"))
  names(pv_dfc)[c(4,6)] = c("pv_gene_sig","pv_se")
  
  #--------------------------------------------------------------------#
  #combine the gene signatures
  re_int = cbind(ad_dfc[,c(1,2,4, 6,7)], pv_dfc[,c(4,6)])
  
  #Calculate centers
  ad_centroid = data.frame( mean(subset(re_int, dis=="AD")$ad_gene_sig), mean(subset(re_int, dis=="AD")$pv_gene_sig) )
  pv_centroid = data.frame( mean(subset(re_int, dis=="PV")$ad_gene_sig), mean(subset(re_int, dis=="PV")$pv_gene_sig) )
  names(ad_centroid) = c("xpt","ypt")
  names(pv_centroid) = c("xpt","ypt")
  
  #--------------------------------------------------------------------#
  #Calculate gene performance 
  
  #sample_level_range = hypotenus of AD range vs PV range
  ad_range = range(re_int$ad_gene_sig)[2] - range(re_int$ad_gene_sig)[1]
  pv_range = range(re_int$pv_gene_sig)[2] - range(re_int$pv_gene_sig)[1]
  sample_range = (ad_range^2 + pv_range^2)^(1/2)
  
  #calculate distance between centroids
  dist_centroids = pointDistance(ad_centroid, pv_centroid,allpairs=T, lonlat=F) #calculate all-vs-all distance
  
  #calculate performance metric = amount of sample range taken up by distance between disease centroids
  performance_metric = dist_centroids / sample_range
  
  #Put them togehter 
  res_df = data.frame(i, ad_range, pv_range, sample_range, dist_centroids, performance_metric)
  names(res_df)[1] = "genes_in_program"
  res_df
  
  }, error=function(e){})
}))

#Plot out performance across a range of program sizes
program_discrimination = ggplot(data=performance_met_df,  aes(x=genes_in_program, y=performance_metric)) +
  geom_point(size=1.5) + #individual cells
  geom_line() +
  #stat_smooth(se=F) +
  theme_classic() +
  ggtitle("Relative discriminatory power based on gene program size") + 
  xlab("Genes in AD and PV programs") +
  ylab("Relative power") +
  theme(axis.text=element_text(size=20,color='black'),
        axis.title=element_text(size=20,color='black'),
        plot.title=element_text(size=20,color='black')) +
  scale_x_continuous(breaks=seq(0, max(performance_met_df$genes_in_program), 15))

png("/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/results/cluster_summary_fig/performance_discriminatory_power_v2.png",
    width=15,height=6,units="in",res=1000)
  print(program_discrimination)
dev.off()  

#Order from highest performing gene number to lowest
performance_met_df_ordered = performance_met_df[order(-performance_met_df$performance_metric),]

#Plot with optimal gene program size colored
program_discrimination_opt = program_discrimination +
  geom_vline(aes(xintercept = performance_met_df_ordered$genes_in_program[1]),color="red", size=1)
  
png("/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/results/cluster_summary_fig/performance_discriminatory_power_v2.png",
    width=15,height=6,units="in",res=1000)
print(program_discrimination_opt)
dev.off()  

#Get the optimized gene programs
optimal_program = cbind(pv_genes[1:performance_met_df_ordered$genes_in_program[1],],
  ad_genes[1:performance_met_df_ordered$genes_in_program[1],])

#Write optimal program to excel sheet
write.xlsx(
  optimal_program,
  "/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/results/optimal_program_v2.xlsx",
  sheetName = "Optimal program",
  col.names = TRUE,
  row.names = F,
  append = FALSE,
  showNA = TRUE,
  password = NULL
)


#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#Clust 2 cells with optimized gene program

#Get the genes of interest: from the original RashX algo
ad_genes = ad_genes[1:performance_met_df_ordered$genes_in_program[1],]$gene
pv_genes = pv_genes[1:performance_met_df_ordered$genes_in_program[1],]$gene
clust2_cds@assays@data$counts[1:10,1:10]

#--------------------------------------------------------------------#
#AD gene data_matrix
ad_mat = t(clust2_cds@assays@data$counts[ad_genes,])
ad_mat[1:10,1:10]
ad_int = data.frame(colData(clust2_cds)$dis_updated, colData(clust2_cds)$ident,
                    rowSums(ad_mat) )
names(ad_int) = c("dis","sample","gene_sig")

#UNLOAD DPLYR PACKAGE#
ad_dfc = summarySE(ad_int, measurevar="gene_sig", groupvars=c("dis","sample"))
head(ad_dfc)
names(ad_dfc)[c(4,6)] = c("ad_gene_sig","ad_se")

#--------------------------------------------------------------------#
#PV gene data_matrix
pv_mat = t(clust2_cds@assays@data$counts[pv_genes,])
pv_mat[1:10,1:10]
pv_int = data.frame(colData(clust2_cds)$dis_updated, colData(clust2_cds)$ident,
                    rowSums(pv_mat) )
names(pv_int) = c("dis","sample","gene_sig")
pv_dfc = summarySE(pv_int, measurevar='gene_sig', groupvars=c("dis","sample"))
head(pv_dfc)
names(pv_dfc)[c(4,6)] = c("pv_gene_sig","pv_se")

#--------------------------------------------------------------------#
#Now make hull plot: https://www.r-bloggers.com/2021/07/ggforce-make-a-hull-plot-to-visualize-clusters-in-ggplot2/
re_int = cbind(ad_dfc[,c(1,2,4, 6,7)], pv_dfc[,c(4,6)]) #combine the gene signatures

#Calculate centers
ad_centroid = data.frame( mean(subset(re_int, dis=="AD")$ad_gene_sig), mean(subset(re_int, dis=="AD")$pv_gene_sig) )
pv_centroid = data.frame( mean(subset(re_int, dis=="PV")$ad_gene_sig), mean(subset(re_int, dis=="PV")$pv_gene_sig) )
names(ad_centroid) = c("xpt","ypt")
names(pv_centroid) = c("xpt","ypt")

#--------------------------------------------------------------------#
#Plot everything together: cells, sample centroids, disease centroids

int_dat = data.frame(ad_int,pv_int)

clust2cells_geneprogram_fig = ggplot() +
  geom_point(data=int_dat, alpha=0.25,aes(x=gene_sig, y=gene_sig.1, color=dis), size=0.4) + #individual cells
  geom_point(data=re_int, alpha=1,aes(x=ad_gene_sig, y=pv_gene_sig, color=dis, shape=dis), size=4) + #sample centroids
  #geom_errorbarh(data=re_int, aes(x=pc1,xmax = pc1 + pc1_se, xmin = pc1 - pc1_se, color=dis_updated)) +
  #geom_errorbar(data=re_int, aes(y=pc2,ymax = pc2 + pc2_se, ymin = pc2 - pc2_se, color=dis_updated)) +
  #stat_ellipse(data=re_int, aes(x=ad_gene_sig, y=pv_gene_sig, color=dis)) +
  stat_ellipse(data=int_dat, aes(x=gene_sig, y=gene_sig.1, color=dis)) +
  theme_classic() +
  ggtitle("Trm cells, optimized gene program") + 
  xlab("AD program") +
  ylab("PV program") +
  theme(legend.position="right") +
  scale_color_manual(name="",values = c("AD" = "darkgoldenrod3",
                                        "PV" = "olivedrab4",
                                        "Atypical rash case" = "darkblue")) +
  scale_fill_manual(name="",values = c("AD" = "lightgoldenrod2",
                                       "PV" = "darkolivegreen4",
                                       "Atypical rash case" = "darkolivegreen4")) +
  geom_point(data=ad_centroid, aes(x=xpt,y=ypt), color="darkgoldenrod3", shape=10, size=7, stroke = 2) +
  geom_point(data=pv_centroid, aes(x=xpt,y=ypt), color="darkolivegreen4", shape=10, size=7, stroke =2) +
  theme(axis.text=element_text(size=20,color='black'),
        axis.title=element_text(size=20,color='black'),
        plot.title=element_text(size=20,color='black')) 


png("/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/results/cluster_summary_fig/clust2cells_optimized_geneprogram_v1.png",
    width=8,height=8,units="in",res=1000)
print(clust2cells_geneprogram_fig)
dev.off()  


#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#

#Empiric null testing permutation random gene program sets

optimal_size = performance_met_df_ordered$genes_in_program[1]
num_permutations = 1000

all_genes = row.names(clust2_cds@assays@data$counts)

empiric_null_set = do.call(rbind, lapply(1:num_permutations, function(i){
  tryCatch({ #suppress error
    #i=5
    print(i)
    
    #get 2 random gene set same size as optimal
    random_genes = all_genes[sample(1:length(all_genes), optimal_size, replace=F)]
    random_genes_2 = all_genes[sample(1:length(all_genes), optimal_size, replace=F)]
    
    #--------------------------------------------------------------------#
    #AD gene data_matrix
    ad_mat = t(clust2_cds@assays@data$counts[random_genes,])
    ad_int = data.frame(colData(clust2_cds)$dis_updated, colData(clust2_cds)$ident,
                        rowSums(ad_mat) )
    names(ad_int) = c("dis","sample","gene_sig")
    
    #UNLOAD DPLYR PACKAGE#
    ad_dfc = summarySE(ad_int, measurevar="gene_sig", groupvars=c("dis","sample"))
    names(ad_dfc)[c(4,6)] = c("ad_gene_sig","ad_se")
    
    #--------------------------------------------------------------------#
    #PV gene data_matrix
    pv_mat = t(clust2_cds@assays@data$counts[random_genes_2,])
    pv_int = data.frame(colData(clust2_cds)$dis_updated, colData(clust2_cds)$ident,
                        rowSums(pv_mat) )
    names(pv_int) = c("dis","sample","gene_sig")
    pv_dfc = summarySE(pv_int, measurevar='gene_sig', groupvars=c("dis","sample"))
    names(pv_dfc)[c(4,6)] = c("pv_gene_sig","pv_se")
    
    #--------------------------------------------------------------------#
    #combine the gene signatures
    re_int = cbind(ad_dfc[,c(1,2,4, 6,7)], pv_dfc[,c(4,6)])
    
    #Calculate centers
    ad_centroid = data.frame( mean(subset(re_int, dis=="AD")$ad_gene_sig), mean(subset(re_int, dis=="AD")$pv_gene_sig) )
    pv_centroid = data.frame( mean(subset(re_int, dis=="PV")$ad_gene_sig), mean(subset(re_int, dis=="PV")$pv_gene_sig) )
    names(ad_centroid) = c("xpt","ypt")
    names(pv_centroid) = c("xpt","ypt")
    
    #--------------------------------------------------------------------#
    #Calculate gene performance 
    
    #sample_level_range = hypotenus of AD range vs PV range
    ad_range = range(re_int$ad_gene_sig)[2] - range(re_int$ad_gene_sig)[1]
    pv_range = range(re_int$pv_gene_sig)[2] - range(re_int$pv_gene_sig)[1]
    sample_range = (ad_range^2 + pv_range^2)^(1/2)
    
    #calculate distance between centroids
    dist_centroids = pointDistance(ad_centroid, pv_centroid,allpairs=T, lonlat=F) #calculate all-vs-all distance
    
    #calculate performance metric = amount of sample range taken up by distance between disease centroids
    performance_metric = dist_centroids / sample_range
    
    #Put them togehter 
    res_df = data.frame(i, optimal_size,ad_range, pv_range, sample_range, dist_centroids, performance_metric)
    names(res_df)[1] = "permutation"
    res_df
    
  }, error=function(e){})
}))

sig_cut_off = quantile(empiric_null_set$performance_metric, probs = 0.99)

#Plot out performance across a range of program sizes
empiric_null_fig = ggplot(data=empiric_null_set,  aes(x=optimal_size, y=performance_metric)) +
  geom_violin() + #individual cells
  geom_boxplot(size=1.5) +
  geom_jitter(size=0.2,alpha=0.3, width=0.2) +
  geom_hline(yintercept=sig_cut_off, linetype="dashed", 
             color = "red", size=2) + #percentile cutoff
  geom_point(aes(x=empiric_null_set$optimal_size[1],y=performance_met_df_ordered$performance_metric[1]),colour="darkgreen",shape=18, size=7) + #actual performance metric
  theme_classic() +
  ggtitle("Permutation\nsignificance") + 
  xlab("") +
  ylab("Relative power") +
  theme(axis.text.x=element_text(size=0,color='black'),
        axis.text.y=element_text(size=20,color='black'),
        axis.title=element_text(size=20,color='black'),
        plot.title=element_text(size=20,color='black'))

png("/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/results/cluster_summary_fig/empiric_null_perm_v2.png",
    width=3,height=8,units="in",res=1000)
print(empiric_null_fig)
dev.off()  

#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#GSEA on optimal program
#----------------------------------------------------------------------------------------------------#
library(msigdbr)
library(fgsea)

#Get all gene sets for humans
all_gene_sets = msigdbr(species = "Homo sapiens")
head(all_gene_sets)

# fixing format to work with fgsea
pathwaysH = split(x = all_gene_sets$gene_symbol, f = all_gene_sets$gs_name)

#----------------------------------------------------------------------------#
#GSEA for PV program
pv_program = pv_genes[1:performance_met_df_ordered$genes_in_program[1],]
gene_ranks = as.vector(abs(pv_program$avg_log2FC)) #make ranked and named vector to put into fgsea
names(gene_ranks) = pv_program$gene

fgseaRes <- fgsea(pathways = pathwaysH, 
                  stats    = gene_ranks,
                  minSize  = 10,
                  maxSize  = 500,
                  scoreType="pos")

ordered_gsea = fgseaRes[order(fgseaRes$pval),]

#Write optimal program to excel sheet
?write.xlsx
write.xlsx2(
  ordered_gsea,
  "/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/results/optimal_program.xlsx",
  sheetName = "PV GSEA",
  col.names = TRUE,
  row.names = F,
  append = T
)
#--------------------------------------------------------------#
#----------------------------------------------------------------------------#
#GSEA for PV program
ad_program = ad_genes[1:performance_met_df_ordered$genes_in_program[1],]
gene_ranks = as.vector(abs(ad_program$avg_log2FC)) #make ranked and named vector to put into fgsea
names(gene_ranks) = ad_program$gene

fgseaRes <- fgsea(pathways = pathwaysH, 
                  stats    = gene_ranks,
                  minSize  = 10,
                  maxSize  = 500,
                  scoreType="pos")

ordered_gsea = fgseaRes[order(fgseaRes$pval),]

#Write optimal program to excel sheet
write.xlsx2(
  ordered_gsea,
  "/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/results/optimal_program.xlsx",
  sheetName = "AD GSEA",
  col.names = TRUE,
  row.names = F,
  append = T,
  showNA = TRUE,
  password = NULL
)

#--------------------------------------------------------------#

#----------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#Load new data with uknown pt samples!
#----------------------------------------------------------------------------------------------------#
#Read in External samples excel

#In unix command line on terminal unzip all the compressed data
#gunzip /Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/data/external\ samples/*/filtered_feature_bc_matrix/*.gz

humandat = readRDS("/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/data/Manuscript\ Object.rds")
data.frame(unique(humandat$donor))

#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#

BiocManager::install("scDblFinder")
library(scDblFinder)

#Read in external sample metadat
external_sample_metadat = read.xlsx(
  file="/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/data/external\ samples/external\ sample\ details\ MT\ 20230406.xlsx",
  colNames=T,
  sheetIndex=1
)

#Load cell-typed original data with core samples
skin.integrated <- readRDS("/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/data/Manuscript\ Object.rds")

args = commandArgs(trailingOnly=TRUE)
#name <- basename(args[1])
output_dir <- "/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/data/external\ samples/cell_typing_output/" # used it for both batch1 and batch2 - just change the name
samples_dir <- "/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/data/external\ samples/h5_files/" # used it for both batch1 and batch2 - just change the name
setwd(samples_dir)
all_samples <- list.files()
ext <- NULL

external_samp_sum = do.call(rbind, lapply(1:length(all_samples), function(d) {
  tryCatch({ #suppress error
  print(d)
  #d=1
  
  data <- Read10X_h5(all_samples[d])
  if(!is.null(names(data))){
    data[["Gene Expression"]] <- CollapseSpeciesExpressionMatrix(data[["Gene Expression"]], prefix = "GRCh38_", controls = "mm10___", ncontrols =0)
    skinX.query <- CreateSeuratObject(data[["Gene Expression"]])
    #saveRDS(skinX.query, paste0(output_dir,all_samples[d],".rds",sep=""))
  } else if(is.null(names(data))){
    skinX.query <- CreateSeuratObject(data)
    #saveRDS(skinX.query, paste0(output_dir,all_samples[d],".rds",sep=""))
  }

  DefaultAssay(skinX.query) <- "RNA"
  
  ############# this part to upload, user may upload more than one files, how to assign them automatically##########
  #1. Find doublets - only if the samples are sparse matrices or h5. If samples are in RDS format, skip this step
  if(is.null(ext)){
    #skinX.sce <- as.SingleCellExperiment(skinX.query)
    #skinX.sce <- scDblFinder(skinX.sce)
    #skinX.query <- as.Seurat(skinX.sce)
    skinX.query[["percent.mt"]] <- PercentageFeatureSet(skinX.query, pattern = "^MT-")
    skinX.query <- subset(skinX.query, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
    #skinX.query <- subset(skinX.query, scDblFinder.class== "singlet")
  }
  ####Normalization
  skinX.query <- NormalizeData(skinX.query)
  
  #6. Cell type classification using an integrated reference 
  #(i.e. Mapping and annotating external query datasets based on our reference datasets to define Trm1 cells).
  #find anchors and predict.id
  skinX.anchors <- FindTransferAnchors(reference = skin.integrated, query = skinX.query, reference.reduction = 'pca')
  predictionsX <- TransferData(anchorset = skinX.anchors, refdata = skin.integrated$ID)  
  skinX.query <- AddMetaData(skinX.query, metadata = predictionsX)
  
  #Add sample metadata from excel sheet
  skin_id = data.frame(strsplit(all_samples[d],"_"))[1,1]
  sample_metadat_match = external_sample_metadat[match(skin_id, external_sample_metadat$Skin..),]
  skinX.query$Sample.name = sample_metadat_match$Sample.name
  skinX.query$dis = sample_metadat_match$dis
  skinX.query$donor = skin_id
  skinX.query$skin_id = skin_id
  
  #Make predicted ID into ID for eventual merging
  skinX.query$ID = skinX.query$predicted.id
  
  #Get cluster 2 cell percentages
  cell_id = data.frame(skinX.query$predicted.id)
  names(cell_id) = "cell_ID"
  n_cells = nrow(cell_id)
  n_cells_clust2 = nrow(subset(cell_id, cell_ID == 2))
  percent_cell_clust2 = round((n_cells_clust2/n_cells)*100,2)
  
  #Write out seurat object that has predicted cell types
  saveRDS(skinX.query, paste0(output_dir,all_samples[d],".rds",sep=""))
  
  #Return percentage of clust 2 cells
  clust_2_percent_df = data.frame(skin_id, n_cells,n_cells_clust2,percent_cell_clust2)
  return(clust_2_percent_df)
  
  }, error=function(e){})
}))

#-------------------------------------------------------------------------------------------------#
#Merge query samples with original sample
#-------------------------------------------------------------------------------------------------#

#1. Load all query samples as rds objects
setwd("/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/data/external\ samples/cell_typing_output")
query_samples = list.files()

for(i in 1:length(query_samples)) {
    #i=1
    print(i)
    sum.name <- paste("seur.", i, sep = "")
    assign(sum.name, readRDS(query_samples[i])
  )
}

#merge.Seurat: can only merge 2 seurat objects at a time which is fucking annoying
skin.integrated2 <- subset(x = skin.integrated, subset = ID == "2")
rm(skin.integrated)
table(skin.integrated2$ID)

human_dat <- merge(skin.integrated2, subset(seur.1, subset = ID == "2"))
human_dat <- merge(human_dat, subset(seur.2, subset = ID == "2"))
human_dat <- merge(human_dat, subset(seur.3, subset = ID == "2"))
human_dat <- merge(human_dat, subset(seur.4, subset = ID == "2"))

rm(seur.1)
rm(seur.2)
rm(seur.3)
rm(seur.4)

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#RashX with new samples
#------------------------------------------------------------------------------------------#

#Seurat object --> CDS
exp_mat <- human_dat@assays[["RNA"]]@data #pull out NORMALIZED counts from Seurat object
cell_metadat <- human_dat@meta.data #pull out cell meta-data from Seurat object
gene_annot = data.frame(human_dat@assays[["RNA"]]@counts@Dimnames[[1]])#pull out gene names from Seurat object
names(gene_annot) = "gene_short_name"
row.names(gene_annot) = gene_annot$gene_short_name #row.names of gene_metadata must be equal to row.names of expression_data

clust2_cds <- new_cell_data_set(exp_mat,
                                cell_metadata = cell_metadat,
                                gene_metadata = gene_annot)

#Subset to core samples + 4 external samples
table(colData(clust2_cds)$donor)
pts = c("170","198","230","231","232","233","236","165","173","194","199","211","222","234","235",
        "266","276","279","289")
clust2_cds = clust2_cds[,colData(clust2_cds)$donor %in% pts]

#------------------------------------------------------------------------------------------#
#Organize metadata factor levels
colData(clust2_cds)$donor
table(colData(clust2_cds)$dis)

colData(clust2_cds)$sample_group = colData(clust2_cds)$dis
colData(clust2_cds)$sample_group = gsub("AD1","core",colData(clust2_cds)$sample_group)
colData(clust2_cds)$sample_group = gsub("AE","core",colData(clust2_cds)$sample_group)
colData(clust2_cds)$sample_group = gsub("Pso","core",colData(clust2_cds)$sample_group)
colData(clust2_cds)$sample_group = gsub("Test AD 1","external",colData(clust2_cds)$sample_group)
colData(clust2_cds)$sample_group = gsub("Test AD 2","external",colData(clust2_cds)$sample_group)
colData(clust2_cds)$sample_group = gsub("Test PV 1","external",colData(clust2_cds)$sample_group)
colData(clust2_cds)$sample_group = gsub("Test PV 2","external",colData(clust2_cds)$sample_group)
table(colData(clust2_cds)$sample_group)

colData(clust2_cds)$dis_cat = colData(clust2_cds)$dis
colData(clust2_cds)$dis_cat = gsub("AD1","AD",colData(clust2_cds)$dis_cat)
colData(clust2_cds)$dis_cat = gsub("AE","AD",colData(clust2_cds)$dis_cat)
colData(clust2_cds)$dis_cat = gsub("Pso","PV",colData(clust2_cds)$dis_cat)
colData(clust2_cds)$dis_cat = gsub("Test AD 1","AD",colData(clust2_cds)$dis_cat)
colData(clust2_cds)$dis_cat = gsub("Test AD 2","AD",colData(clust2_cds)$dis_cat)
colData(clust2_cds)$dis_cat = gsub("Test PV 1","PV",colData(clust2_cds)$dis_cat)
colData(clust2_cds)$dis_cat = gsub("Test PV 2","PV",colData(clust2_cds)$dis_cat)
table(colData(clust2_cds)$dis_cat)

#------------------------------------------------------------------------------------------#
#Process CDS object with 4 external samples and cotrol across project

clust2_cds <- preprocess_cds(clust2_cds, 
                             method="PCA",
                             num_dim = 100,
                             norm_method = "none") #these data have already been normalized by Seurat

clust2_cds = align_cds(cds = clust2_cds,
                    num_dim = 100,
                    preprocess_method = "PCA",
                    alignment_group = "sample_group",
                    verbose=T)

clust2_cds <- reduce_dimension(
  clust2_cds,
  max_components = 2,
  reduction_method = c("UMAP"),
  preprocess_method = "PCA", #must now used ALIGNEMNT PRE-PROCESSING METHOD bc it replaces the other pre-process shit
  umap.metric = "cosine",
  umap.min_dist = 0.45, #update optimal min distance value here
  umap.n_neighbors = 30L,
  umap.fast_sgd = FALSE,
  umap.nn_method = "annoy",
  cores = 1,
  verbose = FALSE
)

#5. Cluster and partition ?cluster_cells cells into supergroups community detection (Phenograph algorithm) Community detection of cells as part of Phenograph algorithm https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4508757/
#5b.  ?cluster_cells
clust2_cds <- cluster_cells(clust2_cds,
                            reduction_method = c("UMAP"),
                            k = 10, #a bigger k will result in lower resolution 
                            cluster_method = c("leiden"),
                            num_iter = 2,
                            partition_qval = 0.05, #Numeric, the q-value cutoff to determine when to partition
                            weight = FALSE,
                            resolution = NULL,
                            random_seed = NULL,
                            verbose = F)


#------------------------------------------------------------------------------------------#
#Get the genes of interest
ad_genes = c("TWIST1","LGALS1",    "IL32",   "CAPG",  "ITM2C", "MFHAS1", "ANXA1","SOS1", "CSGALNACT1","LMO4",  "IFITM2","S100A10",  "MT-ND5",  "CYSLTR1","PLA2G16", "SYNE2", "THADA",    "NEAT1","IL17RB", "RPL36A","ARHGAP21",    "NBAS",  "ACTG1","PRKX", "TGFBR3",   "TIMP1","TNFSF10", "AHNAK",    "MT-ND2",  "ISG15",  "RPL17",  "LONRF2",   "CD99","TSHZ2", "MMP25",   "IFITM1","MT-ND1",  "BIRC3",  "FAM102A", "LPCAT2","NRIP3", "CRIP1",  "CLU",   "PLP2", "ZFP36",  "ZFP36L2",  "TUBA1B","GATA3","SLC5A3",    "SFXN1", "FANK1","TAGLN2") 
pv_genes = c("CXCL13", "MTRNR2L12",  "CD7",   "MGAT4A","FTH1", "LAYN", "IL17F", "KLRB1", "GNLY","CPM", "CTSH", "GBP5", "SOX4", "CLEC2B",    "GZMB", "CD2",   "CEBPD","ODF2L","LAG3", "LRRN3", "ARHGEF12",  "PTPN13",   "TNFAIP3", "TRPS1", "SNX9", "METRNL",  "BTG1", "JUN","SPOCK2",   "GABARAPL1",  "PMEPA1", "HIST1H1E","RBPJ", "LINC01871", "MAP3K4","H1FX", "UBC",  "GALNT1",  "PNRC1","GABPB1-AS1", "RPS26", "MUC20-OT1",  "CHN1",  "NAP1L4",   "PTMS",  "F2R",   "CTLA4", "DAPK2","RAP1B","CCR6", "B3GALT2", "YPEL2",  "FYN",   "PPDPF","SLA2", "CBLB", "ADGRG1","SARAF") 
clust2_cds@assays@data$counts[1:10,1:10]

#--------------------------------------------------------------------#
#AD gene data_matrix
ad_mat = t(clust2_cds@assays@data$counts[ad_genes,])
ad_mat[1:10,1:10]
ad_int = data.frame(colData(clust2_cds)$dis, colData(clust2_cds)$donor, colData(clust2_cds)$dis_cat, colData(clust2_cds)$sample_group,
                    rowSums(ad_mat) )
names(ad_int) = c("dis","sample","dis_cat","sample_group","gene_sig")

#UNLOAD DPLYR PACKAGE#
ad_dfc = summarySE(ad_int, measurevar="gene_sig", groupvars=c("dis","sample","dis_cat","sample_group"))
names(ad_dfc)[c(6,8)] = c("ad_gene_sig","ad_se")
head(ad_dfc)

#--------------------------------------------------------------------#
#PV gene data_matrix
pv_mat = t(clust2_cds@assays@data$counts[pv_genes,])
pv_mat[1:10,1:10]
pv_int = data.frame(colData(clust2_cds)$dis, colData(clust2_cds)$donor,colData(clust2_cds)$dis_cat, colData(clust2_cds)$sample_group,
                    rowSums(pv_mat) )
names(pv_int) = c("dis","sample","dis_cat","sample_group","gene_sig")
pv_dfc = summarySE(pv_int, measurevar='gene_sig', groupvars=c("dis","sample","dis_cat","sample_group"))
names(pv_dfc)[c(6,8)] = c("pv_gene_sig","pv_se")
head(pv_dfc)

#--------------------------------------------------------------------#
#Now make hull plot: https://www.r-bloggers.com/2021/07/ggforce-make-a-hull-plot-to-visualize-clusters-in-ggplot2/
re_int = cbind(ad_dfc[,c(1,2,3,4, 6,8)], pv_dfc[,c(6,8)]) #combine the gene signatures

#--------------------------------------------------------------------#
#Calculate centers
ad_centroid = data.frame( mean(subset(re_int, dis_cat=="AD" & sample_group=="core")$ad_gene_sig), mean(subset(re_int, dis_cat=="AD" & sample_group=="core")$pv_gene_sig) )
pv_centroid = data.frame( mean(subset(re_int, dis_cat=="PV" & sample_group=="core")$ad_gene_sig), mean(subset(re_int, dis_cat=="PV" & sample_group=="core")$pv_gene_sig) )
names(ad_centroid) = c("xpt","ypt")
names(pv_centroid) = c("xpt","ypt")

#--------------------------------------------------------------------#
#Make data frame to drawn lines between external and internal samples

external_samps_pv = subset(re_int, sample_group=="external" & dis_cat=="PV")[,c(3,5,7)]
pvcent = data.frame("PV",pv_centroid)
names(pvcent) = names(external_samps_pv)
pv_vectors = rbind(pvcent,pvcent,external_samps_pv)
pv_vectors$group = c("a","b","a","b")

external_samps_ad = subset(re_int, sample_group=="external" & dis_cat=="AD")[,c(3,5,7)]
adcent = data.frame("AD",ad_centroid)
names(adcent) = names(external_samps_ad)
ad_vectors = rbind(adcent,adcent,external_samps_ad)
ad_vectors$group = c("a","b","a","b")

#--------------------------------------------------------------------#
#Plot everything together: cells, sample centroids, disease centroids
int_dat = data.frame(ad_int,pv_int)



clust2cells_geneprogram_fig = ggplot() +
  geom_point(data=int_dat, alpha=0.25,aes(x=gene_sig, y=gene_sig.1, color=dis_cat), size=0.4) + #individual cells
  geom_point(data=re_int,aes(x=ad_gene_sig, y=pv_gene_sig, color=dis_cat, shape=dis_cat, alpha=sample_group, size=sample_group)) + #sample centroids
  #geom_errorbarh(data=re_int, aes(x=pc1,xmax = pc1 + pc1_se, xmin = pc1 - pc1_se, color=dis_updated)) +
  #geom_errorbar(data=re_int, aes(y=pc2,ymax = pc2 + pc2_se, ymin = pc2 - pc2_se, color=dis_updated)) +
  #stat_ellipse(data=re_int, aes(x=ad_gene_sig, y=pv_gene_sig, color=dis)) +
  #geom_label_repel(data=re_int,aes(x=ad_gene_sig, y=pv_gene_sig, label = sample), box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50') +
  stat_ellipse(data=int_dat, aes(x=gene_sig, y=gene_sig.1, color=dis_cat)) +
  geom_line(data=pv_vectors, aes(x=ad_gene_sig, y=pv_gene_sig, color=dis_cat,group=group), size=1.7) +
  geom_line(data=ad_vectors, aes(x=ad_gene_sig, y=pv_gene_sig, color=dis_cat,group=group), size=1.7) +
  theme_classic() +
  ggtitle("Trm cells, condition gene program, extenral samples") + 
  xlab("AD program") +
  ylab("PV program") +
  theme(legend.position="right") +
  scale_color_manual(name="",values = c("AD" = "darkgoldenrod3",
                                        "PV" = "olivedrab4")) +
  scale_fill_manual(name="",values = c("AD" = "lightgoldenrod2",
                                       "PV" = "darkolivegreen4")) +
  scale_alpha_manual(name="",values = c("core" = 0.9,
                                        "external" = 1)) +
  scale_size_manual(name="",values = c("core" = 3.5,
                                        "external" = 6)) +
  geom_point(data=ad_centroid, aes(x=xpt,y=ypt), color="darkgoldenrod3", shape=10, size=7, stroke = 2) +
  geom_point(data=pv_centroid, aes(x=xpt,y=ypt), color="darkolivegreen4", shape=10, size=7, stroke =2) +
  theme(axis.text=element_text(size=20,color='black'),
        axis.title=element_text(size=20,color='black'),
        plot.title=element_text(size=20,color='black')) +
  theme(legend.position = "none")


png("/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/results/cluster_summary_fig/clust2cells_geneprogram_externalsamples_v1.png",
    width=10,height=10,units="in",res=1000)
print(clust2cells_geneprogram_fig)
dev.off()  

#----------------------------------------------------------------------------------------------------#
#Zoom in to samples
clust2cells_geneprogram_fig = ggplot() +
  geom_point(data=int_dat, alpha=0.1,aes(x=gene_sig, y=gene_sig.1, color=dis_cat), size=0.4) + #individual cells
  geom_point(data=re_int,aes(x=ad_gene_sig, y=pv_gene_sig, color=dis_cat, shape=dis_cat, alpha=sample_group, size=sample_group)) + #sample centroids
  #geom_errorbarh(data=re_int, aes(x=pc1,xmax = pc1 + pc1_se, xmin = pc1 - pc1_se, color=dis_updated)) +
  #geom_errorbar(data=re_int, aes(y=pc2,ymax = pc2 + pc2_se, ymin = pc2 - pc2_se, color=dis_updated)) +
  #stat_ellipse(data=re_int, aes(x=ad_gene_sig, y=pv_gene_sig, color=dis)) +
  #geom_label_repel(data=re_int,aes(x=ad_gene_sig, y=pv_gene_sig, label = sample), box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50') +
  stat_ellipse(data=int_dat, aes(x=gene_sig, y=gene_sig.1, color=dis_cat)) +
  geom_line(data=pv_vectors, aes(x=ad_gene_sig, y=pv_gene_sig, color=dis_cat,group=group), size=1.7) +
  geom_line(data=ad_vectors, aes(x=ad_gene_sig, y=pv_gene_sig, color=dis_cat,group=group), size=1.7) +
  theme_classic() +
  ggtitle("Trm cells, condition gene program, external samples") + 
  xlab("AD program") +
  ylab("PV program") +
  theme(legend.position="right") +
  scale_color_manual(name="",values = c("AD" = "darkgoldenrod3",
                                        "PV" = "olivedrab4")) +
  scale_fill_manual(name="",values = c("AD" = "lightgoldenrod2",
                                       "PV" = "darkolivegreen4")) +
  scale_alpha_manual(name="",values = c("core" = 0.5,
                                        "external" = 1)) +
  scale_size_manual(name="",values = c("core" = 3.5,
                                       "external" = 8)) +
  geom_point(data=ad_centroid, aes(x=xpt,y=ypt), color="darkgoldenrod3", shape=10, size=7, stroke = 2) +
  geom_point(data=pv_centroid, aes(x=xpt,y=ypt), color="darkolivegreen4", shape=10, size=7, stroke =2) +
  theme(axis.text=element_text(size=22,color='black'),
        axis.title=element_text(size=22,color='black'),
        plot.title=element_text(size=0,color='black')) +
  theme(legend.position = "none") +
  ylim(30,70) +
  xlim(40,75)

png("/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/results/cluster_summary_fig/clust2cells_geneprogram_externalsamples_v1_zoomedin.png",
    width=7,height=7,units="in",res=500)
print(clust2cells_geneprogram_fig)
dev.off()  




#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------##----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------##----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------##----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------##----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------##----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------##----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------##----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------##----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------##----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------##----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------##----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------##----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------##----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------##----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------##----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#




#Seurat object --> CDS
exp_mat <- humandat@assays[["RNA"]]@data #pull out NORMALIZED counts from Seurat object
cell_metadat <- humandat@meta.data #pull out cell meta-data from Seurat object
gene_annot = data.frame(humandat@assays[["RNA"]]@counts@Dimnames[[1]])#pull out gene names from Seurat object
names(gene_annot) = "gene_short_name"
row.names(gene_annot) = gene_annot$gene_short_name #row.names of gene_metadata must be equal to row.names of expression_data

#Make CDS object
human_human_big_cds <- new_cell_data_set(exp_mat,
                                         cell_metadata = cell_metadat,
                                         gene_metadata = gene_annot)

#Load new data with uknown pt samples!
rm(humandat,exp_mat)

length(unique(colData(human_human_big_cds)$donor3))

#Subset to patients of interest
pts = c("170","198","230","231","232","233","236","165","173","194","199","211","222","234","235","202","289","279","266","276")
human_human_big_cds = human_human_big_cds[,colData(human_human_big_cds)$donor3 %in% pts] #subset to all normal patients and a diseased patient
colData(human_human_big_cds)$dis_updated = colData(human_human_big_cds)$donor3
colData(human_human_big_cds)$dis_updated = gsub("170","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("198","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("230","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("231","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("232","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("233","AD",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("236","AD",colData(human_human_big_cds)$dis_updated)

colData(human_human_big_cds)$dis_updated = gsub("165","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("173","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("194","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("199","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("211","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("222","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("234","PV",colData(human_human_big_cds)$dis_updated)
colData(human_human_big_cds)$dis_updated = gsub("235","PV",colData(human_human_big_cds)$dis_updated)

colData(human_human_big_cds)$dis_updated = gsub("202","PV",colData(human_human_big_cds)$dis_updated) 

colData(human_human_big_cds)$dis_updated = gsub("289","test PV",colData(human_human_big_cds)$dis_updated) 
colData(human_human_big_cds)$dis_updated = gsub("279","test PV",colData(human_human_big_cds)$dis_updated) 
colData(human_human_big_cds)$dis_updated = gsub("266","test AD",colData(human_human_big_cds)$dis_updated) 
colData(human_human_big_cds)$dis_updated = gsub("276","test AD",colData(human_human_big_cds)$dis_updated) 

#Now subset to clust 2 T-cells
clust2_cds = human_human_big_cds[,colData(human_human_big_cds)$ID %in% "2"] #subset to all normal patients and a diseased patient

#-------------------------------------------------------------------------------------------------------------------------------------#
external_sample_metadat = read.xlsx(
  file="/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/data/external\ samples/external\ sample\ details\ MT\ 20230406.xlsx",
  colNames=T,
  sheetIndex=1
)
#-------------------------------------------------------------------------------------------------------------------------------------#
#WITHOUT QC!!
#-------------------------------------------------------------------------------------------------------------------------------------#
#Loop through the extenral samples without QCing them and then combine them with the existing human data
for(i in 1:nrow(external_sample_metadat)) {
  #i=1
  print(i)
  #1. load data from 10X output: 
  matrix_dir = paste("/Users/marktaylor/Desktop/Projects/Cho derm/JID review/data/external samples/",external_sample_metadat$Box.folder.Name[i],"/filtered_feature_bc_matrix/",sep="")
  barcode.path <- paste0(matrix_dir, "barcodes.tsv")
  features.path <- paste0(matrix_dir, "features.tsv")
  matrix.path <- paste0(matrix_dir, "matrix.mtx")
  
  mat <- readMM(file = matrix.path) #read in the expression matrix
  
  feature.names = read.delim(features.path,  #read in gene names
                             header = FALSE,
                             stringsAsFactors = FALSE)
  names(feature.names)[1:2] = c("ensemble_id","gene_short_name")
  
  barcode.names = read.delim(barcode.path, #read in barcode names
                             header = FALSE,
                             stringsAsFactors = FALSE)
  
  #create cell_metadata dataframe
  cell_metadat = external_sample_metadat[i,][rep(seq_len(nrow(external_sample_metadat[i,])), each = nrow(barcode.names)), ]
  row.names(cell_metadat) = barcode.names[,1]
  cell_metadat$ID = "unknown" #This is the cell type column that matches the previous columns
  
  #create gene_annotation_file
  gene_metadat = feature.names
  row.names(gene_metadat) = gene_metadat$ensemble_id
  
  #create expression_data matrix
  colnames(mat) = row.names(cell_metadat)
  row.names(mat) = row.names(gene_metadat)
  
  cds_externalsample <- new_cell_data_set(expression_data = mat,
                                    cell_metadata = cell_metadat,
                                    gene_metadata = gene_metadat)

  #Log normalize to combine with previous data that has already been normalized
  cds_externalsample <- preprocess_cds(cds_externalsample, 
                            method="PCA",
                            num_dim = 100,
                            norm_method = "log")
  
  #assign each cds object to a new name
  name <- paste("cds.", i, sep = "")
  assign(name, cds_externalsample)    #create a bunch of cds objects that we are going to combine
} 

#-------------------------------------------------------------------------------------------------------------------------------------#
#WITH QC!!
#-------------------------------------------------------------------------------------------------------------------------------------#
#Loop through the extenral samples without QCing them and then combine them with the existing human data

pdf(file="/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/data/external\ samples/external_sample_QC/cell_ranger_output_QC_without_cutoffs.pdf",height=12,width=10)
for(i in 1:nrow(external_sample_metadat)) {
  #i=1
  print(i)
  #1. load data from 10X output: 
  matrix_dir = paste("/Users/marktaylor/Desktop/Projects/Cho derm/JID review/data/external samples/",external_sample_metadat$Box.folder.Name[i],"/filtered_feature_bc_matrix/",sep="")
  barcode.path <- paste0(matrix_dir, "barcodes.tsv")
  features.path <- paste0(matrix_dir, "features.tsv")
  matrix.path <- paste0(matrix_dir, "matrix.mtx")
  
  mat <- readMM(file = matrix.path) #read in the expression matrix
  
  feature.names = read.delim(features.path,  #read in gene names
                             header = FALSE,
                             stringsAsFactors = FALSE)
  names(feature.names)[1:2] = c("ensemble_id","gene_short_name")
  
  barcode.names = read.delim(barcode.path, #read in barcode names
                             header = FALSE,
                             stringsAsFactors = FALSE)
  
  #create cell_metadata dataframe
  cell_metadat = external_sample_metadat[i,][rep(seq_len(nrow(external_sample_metadat[i,])), each = nrow(barcode.names)), ]
  row.names(cell_metadat) = barcode.names[,1]
  cell_metadat$ID = "unknown" #This is the cell type column that matches the previous columns
  
  #create gene_annotation_file
  gene_metadat = feature.names
  row.names(gene_metadat) = make.unique(gene_metadat$gene_short_name)
  
  #create expression_data matrix
  colnames(mat) = row.names(cell_metadat)
  row.names(mat) = row.names(gene_metadat)
  
  cds.balmain1 <- new_cell_data_set(expression_data = mat,
                                          cell_metadata = cell_metadat,
                                          gene_metadata = gene_metadat)
  
  #----------------------------------------------------------------------------------------------#
  #Start Seurat QC --> cds.balmain2
  #----------------------------------------------------------------------------------------------#
  #Quality Control after: https://broadinstitute.github.io/2019_scWorkshop/data-wrangling-scrnaseq.html#filtering-low-quality-cells
  
  counts_per_cell <- Matrix::colSums(mat)
  counts_per_gene <- Matrix::rowSums(mat)
  genes_per_cell <- Matrix::colSums(mat>0) # count gene only if it has non-zero reads mapped.
  cells_per_gene <- Matrix::rowSums(mat>0) # only count cells where the gene is expressed
  
  #plot counts
  par(mfrow=c(2,2))
  hist((counts_per_cell+1),xlab='counts per cell',main=paste(external_sample_metadat$dis[i],'counts per cell'),col='wheat')
  plot(counts_per_cell, xlab='counts per cell',ylab='genes per cell',genes_per_cell, log='xy', main=paste(external_sample_metadat$dis[i],'counts vs genes per cell'),col='wheat')
  
  #plots genes
  hist((genes_per_cell+1), xlab='genes per cell', main=paste(external_sample_metadat$dis[i],'genes per cell'), col='wheat')
  plot(sort(genes_per_cell), xlab='cell', ylab='genes per cell', main=paste(external_sample_metadat$dis[i],'genes per cell (ordered)'))#distribution of library complexity --> libraries with hardly any genes should probably be dropped
  
  #create Seurat object
  #make row names into gene symbol, not ENSEMBL ?CreatSeuratObject
  row.names(mat) = feature.names$gene_short_name
  seurat<-CreateSeuratObject(counts = mat, min.cells = 3, min.features = 100, project = "EK-3544")
  
  #QC: add mitochondrial percent
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "MT-")
  print(
    VlnPlot(seurat, features = c("percent.mt"), ncol = 1) + geom_violin(fill="wheat") + theme(legend.position = "none") + ggtitle(paste(external_sample_metadat$dis[i],"% mitochondrial")) + xlab("") + theme(axis.text.x = element_text(size=0))
  )
  
  #extract: median IQR
  
  #pre-processing matrix : mat
  n_cells = dim(mat)[2]
  counts_per_cell <- format_summary(data.frame(as.vector(summary(Matrix::colSums(mat))))[c(3,2,5),])
  counts_per_gene <- format_summary(data.frame(as.vector(summary(Matrix::rowSums(mat))))[c(3,2,5),])
  genes_per_cell <- format_summary(data.frame(as.vector(summary(Matrix::colSums(mat>0))))[c(3,2,5),]) # count gene only if it has non-zero reads mapped.
  cells_per_gene <- format_summary(data.frame(as.vector(summary(Matrix::rowSums(mat>0))))[c(3,2,5),]) # only count cells where the gene is expressed
  
  sum.name <- paste("sum.", i, sep = "")
  assign(sum.name, data.frame(external_sample_metadat$dis[i],
                              n_cells, counts_per_cell, counts_per_gene, genes_per_cell, cells_per_gene
  )
  )
  
} 

dev.off()

#write out QC summary
totcombsum = rbind(sum.1, sum.2,sum.3,sum.4) #combine the things assigned above
write.table(totcombsum, file="/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/data/external\ samples/external_sample_QC/cell_ranger_output_QC_without_cutoffs.txt",row.names=F,quote=F,sep="\t")

#----------------------------------------------------------------------------------------------------------------------------------------#
#Apply QC cutoffs
#----------------------------------------------------------------------------------------------------------------------------------------#

external_sample_metadat = read.xlsx(
  file="/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/data/external\ samples/external\ sample\ details\ MT\ 20230406.xlsx",
  colNames=T,
  sheetIndex=1
)

pdf(file="/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/data/external\ samples/external_sample_QC/cell_ranger_output_QC_with_cutoffs.pdf",height=12,width=10)
for(i in 1:nrow(external_sample_metadat)) {
  #i=1
  print(i)
  
  samp_names = external_sample_metadat$dis[i]
  mingene = external_sample_metadat$Minimum.gene.threshold
  maxgene = external_sample_metadat$Maximum.gene.threshold
  mito_cutoff = external_sample_metadat$Mitochondrial.threshold
  
  #1. load data from 10X output: 
  matrix_dir = paste("/Users/marktaylor/Desktop/Projects/Cho derm/JID review/data/external samples/",external_sample_metadat$Box.folder.Name[i],"/filtered_feature_bc_matrix/",sep="")
  barcode.path <- paste0(matrix_dir, "barcodes.tsv")
  features.path <- paste0(matrix_dir, "features.tsv")
  matrix.path <- paste0(matrix_dir, "matrix.mtx")
  
  mat <- readMM(file = matrix.path) #read in the expression matrix
  
  feature.names = read.delim(features.path,  #read in gene names
                             header = FALSE,
                             stringsAsFactors = FALSE)
  names(feature.names)[1:2] = c("ensemble_id","gene_short_name")
  
  barcode.names = read.delim(barcode.path, #read in barcode names
                             header = FALSE,
                             stringsAsFactors = FALSE)
  
  #create cell_metadata dataframe
  cell_metadat = external_sample_metadat[i,][rep(seq_len(nrow(external_sample_metadat[i,])), each = nrow(barcode.names)), ]
  row.names(cell_metadat) = barcode.names[,1]
  cell_metadat$ID = "unknown" #This is the cell type column that matches the previous columns
  
  #create gene_annotation_file
  gene_metadat = feature.names
  row.names(gene_metadat) = make.unique(gene_metadat$gene_short_name)
  
  #create expression_data matrix
  colnames(mat) = row.names(cell_metadat)
  row.names(mat) = row.names(gene_metadat)
  
  cds.balmain1 <- new_cell_data_set(expression_data = mat,
                                    cell_metadata = cell_metadat,
                                    gene_metadata = gene_metadat)
  

#----------------------------------------------------------------------------------------------#
#Start Seurat QC --> cds.balmain2
#----------------------------------------------------------------------------------------------#
#Quality Control after: https://broadinstitute.github.io/2019_scWorkshop/data-wrangling-scrnaseq.html#filtering-low-quality-cells

counts_per_cell <- Matrix::colSums(mat)
counts_per_gene <- Matrix::rowSums(mat)
genes_per_cell <- Matrix::colSums(mat>0) # count gene only if it has non-zero reads mapped.
cells_per_gene <- Matrix::rowSums(mat>0) # only count cells where the gene is expressed

#plot counts
par(mfrow=c(2,2))
tryCatch({  hist((genes_per_cell+1), xlab='genes per cell', main=paste(external_sample_metadat$dis[i],'genes per cell'), col='wheat') + abline(v=c(mingene[i],maxgene[i]), col="red",lwd=3)}, error=function(e){})

#plots genes
plot(counts_per_cell, xlab='counts per cell',ylab='genes per cell',genes_per_cell, log='xy', main=paste(external_sample_metadat$dis[i],'counts vs genees per cell'),col='wheat')  + abline(h=c(mingene[i],maxgene[i]), col="red",lwd=3)
tryCatch({  hist((genes_per_cell+1), xlab='genes per cell', main=paste(external_sample_metadat$dis[i],'genes per cell'), col='wheat') + abline(v=c(mingene[i],maxgene[i]), col="red",lwd=3)}, error=function(e){})
tryCatch({ plot(sort(genes_per_cell), xlab='cell', ylab='genes per cell', main=paste(external_sample_metadat$dis[i],'genes per cell (ordered)'))  + abline(h=c(mingene[i],maxgene[i]), col="red",lwd=3) }, error=function(e){})#distribution of library complexity --> libraries with hardly any genes should probably be dropped

#create Seurat object
#make row names into gene symbol, not ENSEMBL ?CreatSeuratObject
row.names(mat) = feature.names$gene_short_name
seurat<-CreateSeuratObject(counts = mat, min.cells = 3, min.features = 200, project = "EK-3205")

#QC: add mitochondrial percent
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "MT-")
print(
  VlnPlot(seurat, features = c("percent.mt"), ncol = 1) + geom_violin(fill="wheat") + theme(legend.position = "none") + ggtitle(paste(external_sample_metadat$dis[i],"% mitochondrial")) + xlab("") + geom_hline(yintercept=mito_cutoff[i], color="red",size=1.25) + theme(axis.text.x = element_text(size=0))
)

#Getting rid of low-quality cells/empty droplets (<mingene genes) ; doublets (>maxgene) ; stresst/dying cells (<10% mito)
seurat <- subset(seurat, subset = nFeature_RNA > mingene[i] & nFeature_RNA < maxgene[i] & percent.mt < external_sample_metadat$dis[i])

#normalize
#Klein lab normalization https://www.nature.com/articles/nature25741#MOESM1 The gene expression counts of each cell were then normalized using a variant of total-count normalization that avoids distortion from very highly expressed genes. Specifically, we calculated , the normalized transcript counts for gene j in cell i, from the raw counts xi,j as follows: , in which  and  is the average of Xi over all cells. To prevent very highly expressed genes (for example, haemoglobin) from correspondingly decreasing the relative expression of other genes, we excluded genes comprising >10% of the total counts of any cell when calculating  and Xi.
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)

#Seurat object --> CDS
exp_mat <- seurat@assays[["RNA"]]@data #pull out NORMALIZED counts from Seurat object
cell_metadat <- data.frame(cell_metadat[(1:ncol(exp_mat)),], seurat@meta.data, sample_metadat[i,]) #pull out cell meta-data from Seurat object and combine with metadata made above
row.names(cell_metadat) = row.names(seurat@meta.data)
gene_annot = data.frame(seurat@assays[["RNA"]]@data@Dimnames[[1]])#pull out gene names from Seurat object
names(gene_annot) = "gene_short_name"
row.names(gene_annot) = gene_annot$gene_short_name #row.names of gene_metadata must be equal to row.names of expression_data

#assign each cds object to a new name
name <- paste("cds.", i, sep = "")
assign(name, new_cell_data_set(exp_mat,
                               cell_metadata = cell_metadat,
                               gene_metadata = gene_annot)    #create a bunch of cds objects that we are going to combine
)

#summary BEFORE & AFTER QC metrics
#extract: median IQR

#pre-processing matrix : mat
n_cells = dim(mat)[2]
counts_per_cell <- format_summary(data.frame(as.vector(summary(Matrix::colSums(mat))))[c(3,2,5),])
counts_per_gene <- format_summary(data.frame(as.vector(summary(Matrix::rowSums(mat))))[c(3,2,5),])
genes_per_cell <- format_summary(data.frame(as.vector(summary(Matrix::colSums(mat>0))))[c(3,2,5),]) # count gene only if it has non-zero reads mapped.
cells_per_gene <- format_summary(data.frame(as.vector(summary(Matrix::rowSums(mat>0))))[c(3,2,5),]) # only count cells where the gene is expressed
#post-processing matrix : exp_mat
n_cells.exp = dim(exp_mat)[2]
counts_per_cell.exp <- format_summary(data.frame(as.vector(summary(Matrix::colSums(exp_mat))))[c(3,2,5),])
counts_per_gene.exp <- format_summary(data.frame(as.vector(summary(Matrix::rowSums(exp_mat))))[c(3,2,5),])
genes_per_cell.exp <- format_summary(data.frame(as.vector(summary(Matrix::colSums(exp_mat>0))))[c(3,2,5),]) # count gene only if it has non-zero reads mapped.
cells_per_gene.exp <- format_summary(data.frame(as.vector(summary(Matrix::rowSums(exp_mat>0))))[c(3,2,5),]) # only count cells where the gene is expressed

sum.name <- paste("sum.", i, sep = "")
assign(sum.name, data.frame(external_sample_metadat$dis[i],
                            n_cells, counts_per_cell, counts_per_gene, genes_per_cell, cells_per_gene,
                            n_cells.exp,counts_per_cell.exp, counts_per_gene.exp, genes_per_cell.exp, cells_per_gene.exp)
)

} 

dev.off()

#write out QC summary
totcombsum = rbind(sum.1, sum.2,sum.3,sum.4) #combine the things assigned above
write.table(totcombsum, file="/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/data/external\ samples/external_sample_QC/cell_ranger_output_QC_WITH_CUTOFFS.txt",row.names=F,quote=F,sep="\t")

names(colData(human_human_big_cds))
unique(colData(human_human_big_cds)$donor): 
table(colData(human_human_big_cds)$dis)
  
########################################################################################################################
########################################################################################################################
########################################################################################################################

#Add column to correct on
colData(human_human_big_cds)$project = "original samples"
colData(clust2_cds)$project = "original samples"
colData(cds.1)$project = "external samples"
colData(cds.2)$project = "external samples"
colData(cds.3)$project = "external samples"
colData(cds.4)$project = "external samples"

#create 1 giant ass cds object from all the cds's created above
big_cds <- combine_cds(list(cds.4,cds.3,cds.2,cds.1,human_human_big_cds))

big_cds <- preprocess_cds(big_cds, 
                          method="PCA",
                          num_dim = 100,
                          norm_method = "none") #no normalization needed because it was already normalized by Seurat in QC

#3b. Pre-process with Seurat-processed object
big_cds = align_cds(cds = big_cds,
                    num_dim = 100,
                    preprocess_method = "PCA",
                    alignment_group = "project",
                    verbose=T)

#4b. ?reduce_dimension
big_cds <- reduce_dimension(
  big_cds,
  max_components = 2,
  reduction_method = c("UMAP"),
  preprocess_method = "Aligned", #must now used ALIGNEMNT PRE-PROCESSING METHOD bc it replaces the other pre-process shit
  umap.metric = "cosine",
  umap.min_dist = 0.1,
  umap.n_neighbors = 15L,
  umap.fast_sgd = FALSE,
  umap.nn_method = "annoy",
  cores = 1,
  verbose = FALSE
)

#5. Cluster and partition ?cluster_cells cells into supergroups community detection (Phenograph algorithm) Community detection of cells as part of Phenograph algorithm https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4508757/

#5b. 
?cluster_cells
big_cds <- cluster_cells(big_cds,
                         reduction_method = c("UMAP"),
                         k = 10, #a bigger k will result in lower resolution 
                         cluster_method = c("leiden"),
                         num_iter = 2,
                         partition_qval = 0.05, #Numeric, the q-value cutoff to determine when to partition
                         weight = FALSE,
                         resolution = NULL,
                         random_seed = NULL,
                         verbose = F)


#Now see if there is a set of nice colors
colourCount = length(unique(colData(big_cds)$ID))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

# Use this code if there is a duplicated sample_name column
# cell_metadata = colData(big_cds)
# cell_metadata = cell_metadata[ ,-grep("sample_name",names(cell_metadata))[2]] #delete duplicated sample_name column
# colData(big_cds) <- cell_metadata


png("/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/data/external\ samples/external_sample_figs/clusters_dims_12.png",
    width=15,height=15,units="in",res=400)

print(
plot_cells(big_cds, color_cells_by="ID", group_cells_by="cluster",
                 group_label_size = 2,
                 x=1,
                 y=2,
                 alpha=0.6) +
        ggtitle("") + xlab("") + ylab("") + 
        theme(axis.text=element_text(size=0,color='black'),
              axis.title=element_text(size=0),
              axis.ticks = element_line(size=0),
              legend.position="right") +
        scale_color_manual(values=getPalette(colourCount))
)

dev.off()


#----------------------------------------------------------------------------------------------------------------#

#Now go through various metadata columns and groups

metadatcols = c("ID", "donor","dis") #shared metadata

#Loop through metadata columns and print out the data itself 
for(i in 1:length(metadatcols)) {

  tryCatch({ #suppress error
    
    #i=3
    
    #--------------------------------------------------------------------------#
    #Plot all factor levels together
    match_index = match(metadatcols[i], names(colData(big_cds)), ) #get column index of metadata field I wanna look at
    colourCount = length(unique(colData(big_cds)[,match_index]))
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    
    png(paste("/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/data/external\ samples/external_sample_figs/metadata/",
              i,".",0,".",metadatcols[i],".png",sep=""),
        width=15,height=15,units="in",res=250)
    print(
      plot_cells(big_cds, color_cells_by=metadatcols[i], group_cells_by="cluster",
                 group_label_size = 0,
                 x=1,
                 y=2,
                 alpha=0.6) +
        ggtitle(metadatcols[i]) + xlab("") + ylab("") + 
        theme(axis.text=element_text(size=0,color='black'),
              axis.title=element_text(size=0),
              title=element_text(size=15),
              axis.ticks = element_line(size=0),
              legend.position="right") +
        scale_color_manual(values=getPalette(colourCount))
    )
    dev.off()
    
    #--------------------------------------------------------------------------#
    #Now plot each factor level independently
    factor_levels = na.omit(unique(colData(big_cds)[,match_index]))
    
    for(j in 1:length(factor_levels)){
      #j=1
      
      print(paste("i=",i," j=",j,sep=""))
      tryCatch({ #suppress error
        
        matching_cell_ids = row.names(colData(big_cds)[ colData(big_cds)[,match_index] %in% factor_levels[j] , ]) #get cell id's that match a factor level
        
        #Add dummy column
        colData(big_cds)$isolated_cell_population = "Non-match"
        
        #Now match isolated cell IDs with all cell IDs and replace non-match with match in dummy column
        sub_cell_matches = match(matching_cell_ids,colnames(big_cds)) #match extracted cells with all cell names
        colData(big_cds)$isolated_cell_population[sub_cell_matches] <- "Match" #now replace matching cells with new name
        
        #Exract low-D embeddings and make a graph object bitch
        graph_df = data.frame( big_cds@int_colData@listData$reducedDims$UMAP , colData(big_cds)$isolated_cell_population)
        names(graph_df) = c("UMAP1","UMAP2","cell_pop")
        
        #Now plot out
        png(paste("/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/data/external\ samples/external_sample_figs/metadata/",
                  i,".",j,".",metadatcols[i],"_",factor_levels[j],".png",sep=""),
            width=7,height=7,units="in",res=600)
        
        print(
          ggplot(data=graph_df, aes(x=UMAP1,y=UMAP2)) +
            theme_classic() +
            geom_point(size=0.01, aes(alpha=cell_pop,color=cell_pop)) +
            ggtitle( paste( metadatcols[i] , factor_levels[j], sep=": ") ) + xlab("") + ylab("") + 
            theme(axis.text=element_text(size=0,color='black'),
                  axis.title=element_text(size=0),
                  title=element_text(size=15),
                  axis.ticks = element_line(size=0),
                  legend.position="none") +
            scale_color_manual(values = c("Non-match" = "lightgrey", 
                                          "Match" = "darkred")) +
            scale_alpha_manual(values = c("Non-match" = 0.1, 
                                          "Match" = 0.8))
        )
        dev.off()
        
      }, error=function(e){})
      
    } #End factor level loop
    
  }, error=function(e){})
  
} #End metadata field loop

#-----------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------#
#Select cluster 2 cells

cluster_2_cells = choose_cells(big_cds)


#Get percentages of cluster 2 cells of priginal core group and now newly selected cluster 2 cells

metadat = colData(big_cds)
names(metadat)

#Get metadata of original core_group of samples
metadat_coregroup = subset(metadat, ID != "unknown")
n_clust_2_cells = nrow(subset(metadat_coregroup, ID == "2"))
round((n_clust_2_cells / nrow(metadat_coregroup) )*100,2)

#Get sample-wise percentages of cluster 2 cells in original core group
samples = unique(metadat_coregroup$donor)
clust_2_percentages = do.call(rbind, lapply(1:length(samples), function(i){
  #i=1
  sample = samples[i]
  n_clust_2_cells = nrow(subset(metadat_coregroup, ID == "2" & donor == samples[i]))
  total_cells_sample = nrow(subset(metadat_coregroup, donor == samples[i]))
  clust_2_sample_percent = round((n_clust_2_cells / total_cells_sample )*100,2)
  data.frame(sample,n_clust_2_cells,total_cells_sample,clust_2_sample_percent)
}))

#Get metadata from external samples
metadat_externalsamples = subset(metadat, ID == "unknown")

#Get percentage of cluster 2 cells in external samples
cluster_2_cells_metadat = colData(cluster_2_cells)
cluster_2_cells_metadat_externalsamples = subset(cluster_2_cells_metadat, ID == "unknown")
nrow(cluster_2_cells_metadat_externalsamples)
round((nrow(cluster_2_cells_metadat_externalsamples) / nrow(metadat_externalsamples))*100,2)

#Get sample-wise percentages of cluster 2 cells in EXTERNAL SAMPLES
samples = unique(metadat_externalsamples$Sample.name)
samples = unique(metadat_externalsamples$Skin.id..)
samples = unique(metadat_externalsamples$Skin..)

clust_2_percentages = do.call(rbind, lapply(1:length(samples), function(i){
  #i=1
  sample = samples[i]
  n_clust_2_cells = nrow(subset(cluster_2_cells_metadat_externalsamples, Sample.name == samples[i]))
  total_cells_sample = nrow(subset(metadat_externalsamples, Sample.name == samples[i]))
  clust_2_sample_percent = round((n_clust_2_cells / total_cells_sample )*100,2)
  data.frame(sample,n_clust_2_cells,total_cells_sample,clust_2_sample_percent)
}))



#-----------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------#
#Read in external sample data

file(i in 1:nrow(external_sample_metadat)) {
  #i=1
  external_file = paste("/Users/marktaylor/Desktop/Projects/Cho derm/JID review/data/external samples/",
                        external_sample_metadat$Box.folder.Name[i],"/",
                        "filtered_feature_bc_matrix.h5",sep="")
  external_samp = Read10X_h5(external_file, use.names = TRUE, unique.features = TRUE)
  
  
  
}


clust2_cds <- preprocess_cds(clust2_cds, 
                             method="PCA",
                             num_dim = 100,
                             norm_method = "none") #these data have already been normalized by Seurat

#optimal min distance is 0.45
clust2_cds <- reduce_dimension(
  clust2_cds,
  max_components = 2,
  reduction_method = c("UMAP"),
  preprocess_method = "PCA", #must now used ALIGNEMNT PRE-PROCESSING METHOD bc it replaces the other pre-process shit
  umap.metric = "cosine",
  umap.min_dist = 0.45, #update optimal min distance value here
  umap.n_neighbors = 30L,
  umap.fast_sgd = FALSE,
  umap.nn_method = "annoy",
  cores = 1,
  verbose = FALSE
)

#5. Cluster and partition ?cluster_cells cells into supergroups community detection (Phenograph algorithm) Community detection of cells as part of Phenograph algorithm https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4508757/
#5b.  ?cluster_cells
clust2_cds <- cluster_cells(clust2_cds,
                            reduction_method = c("UMAP"),
                            k = 10, #a bigger k will result in lower resolution 
                            cluster_method = c("leiden"),
                            num_iter = 2,
                            partition_qval = 0.05, #Numeric, the q-value cutoff to determine when to partition
                            weight = FALSE,
                            resolution = NULL,
                            random_seed = NULL,
                            verbose = F)

#Now make new dataframe of 
vars_of_interest = c("ident","donor","dis","treat","chem","final_clustering","clusters","ID","type","Ident1","Ident2","dis_updated")
allcell_allgene = data.frame(clust2_cds@int_colData@listData$reducedDims@listData$UMAP)
allcell_allgene = data.frame(allcell_allgene, colData(clust2_cds)[,match(vars_of_interest,names(colData(clust2_cds)))]) #combine UMAP coords with metadata factors of interest
names(allcell_allgene)[c(1,2)] = c("pc1","pc2")

#Get average across pc1 and pc2
#UNLOAD DPLYR PACKAGE#
pc1_dfc= summarySE(allcell_allgene, measurevar="pc1", groupvars=c("dis_updated","ident"))
names(pc1_dfc)[c(5)] = c("pc1_se")
pc2_dfc= summarySE(allcell_allgene, measurevar="pc2", groupvars=c("dis_updated","ident"))
names(pc2_dfc)[c(5)] = c("pc2_se")

#--------------------------------------------------------------------#
#Now make hull plot: https://www.r-bloggers.com/2021/07/ggforce-make-a-hull-plot-to-visualize-clusters-in-ggplot2/
re_int = cbind(pc1_dfc[,c(1,2,4,5)], pc2_dfc[,c(4,5)]) #combine the gene signatures

#--------------------------------------------------------------------#
#Calculate disease centers
ad_centroid = data.frame( mean(subset(re_int, dis_updated=="AD")$pc1), mean(subset(re_int, dis_updated=="AD")$pc2) )
pv_centroid = data.frame( mean(subset(re_int, dis_updated=="PV")$pc1), mean(subset(re_int, dis_updated=="PV")$pc2) )
names(ad_centroid) = c("xpt","ypt")
names(pv_centroid) = c("xpt","ypt")

#--------------------------------------------------------------------#
#Plot everything together: cells, sample centroids, disease centroids
clust2cells_allgenes_fig = ggplot() +
  geom_point(data=allcell_allgene, alpha=0.15,aes(x=pc1, y=pc2, color=dis_updated), size=0.4) + #individual cells
  geom_point(data=re_int, alpha=1,aes(x=pc1, y=pc2, color=dis_updated, shape=dis_updated), size=4) + #sample centroids
  #geom_errorbarh(data=re_int, aes(x=pc1,xmax = pc1 + pc1_se, xmin = pc1 - pc1_se, color=dis_updated)) +
  #geom_errorbar(data=re_int, aes(y=pc2,ymax = pc2 + pc2_se, ymin = pc2 - pc2_se, color=dis_updated)) +
  stat_ellipse(data=allcell_allgene, aes(x=pc1, y=pc2, color=dis_updated)) +
  theme_classic() +
  ggtitle("Trm cells, full transcriptome") + 
  xlab("Dim 1") +
  ylab("Dim 2") +
  theme(legend.position="right") +
  scale_color_manual(name="",values = c("AD" = "darkgoldenrod3",
                                        "PV" = "olivedrab4",
                                        "Atypical rash case" = "darkblue")) +
  scale_fill_manual(name="",values = c("AD" = "lightgoldenrod2",
                                       "PV" = "darkolivegreen4",
                                       "Atypical rash case" = "darkolivegreen4")) +
  geom_point(data=ad_centroid, aes(x=xpt,y=ypt), color="darkgoldenrod3", shape=10, size=7, stroke = 2) +
  geom_point(data=pv_centroid, aes(x=xpt,y=ypt), color="darkolivegreen4", shape=10, size=7, stroke =2) +
  theme(axis.text=element_text(size=20,color='black'),
        axis.title=element_text(size=20,color='black'),
        plot.title=element_text(size=20,color='black')) 


png("/Users/marktaylor/Desktop/Projects/Cho\ derm/JID\ review/results/cluster_summary_fig/clust2cells_allgenes_v1_testsamples.png",
    width=8,height=8,units="in",res=1000)
print(clust2cells_allgenes_fig)
dev.off()  



library(ggforce)
library(concaveman)

with_sample_label = ggplot(data=re_int, aes(x=ad_gene_sig, y=pv_gene_sig)) +
  geom_point(alpha=1,aes(color=dis), size=2) +
  geom_text_repel(data=re_int, size=8, aes(x=ad_gene_sig, y=pv_gene_sig,label=name_alphabet),box.padding = 0.7) + #add labels
  geom_errorbarh(aes(xmax = ad_gene_sig + ad_se, xmin = ad_gene_sig - ad_se, color=dis)) +
  geom_errorbar(aes(ymax = pv_gene_sig + pv_se, ymin = pv_gene_sig - pv_se, color=dis)) +
  #stat_ellipse(aes(color=dis)) +
  theme_classic() +
  ggtitle("Hyperdimensional proximity of Cluster 2 cells") + 
  xlab("AD gene signature") +
  ylab("PV gene signature") +
  theme(legend.position="right") +
  scale_color_manual(name="",values = c("AD" = "darkgoldenrod3",
                                        "PV" = "olivedrab4",
                                        "Atypical rash case" = "darkblue")) +
  scale_fill_manual(name="",values = c("AD" = "lightgoldenrod2",
                                       "PV" = "darkolivegreen4",
                                       "Atypical rash case" = "darkolivegreen4")) +
  geom_point(data=ad_centroid, aes(x=xpt,y=ypt), color="darkgoldenrod3", shape=10, size=5, stroke = 2) +
  geom_point(data=pv_centroid, aes(x=xpt,y=ypt), color="darkolivegreen4", shape=10, size=5, stroke =2) +
  theme(axis.text=element_text(size=20,color='black'),
        axis.title=element_text(size=20,color='black'),
        plot.title=element_text(size=20,color='black')) +
  geom_mark_hull(
    aes(fill=dis), concavity=5
  )

png("/Users/marktaylor/Desktop/Projects/Cho\ derm/Case\ Report/Figs/cell_clusters_integrated_alphabetnames_v2.png",
    width=8,height=8,units="in",res=1000)
print(with_sample_label)
dev.off()  



#-----------------------------------------------------------------------------------------#
# Mean gene heatmaps: doesn't look good
#-----------------------------------------------------------------------------------------#

#Get data to put into heatmap
gene_n = 10 #number of genes from AD and PV sigs to show in heatmap
getdata = data.frame(colData(human_human_big_cds)$donor3, 
                     t(human_human_big_cds@assays@data$counts[gene_sigs$ad_genes[1:gene_n],]),
                     t(human_human_big_cds@assays@data$counts[gene_sigs$pv_genes[1:gene_n],])
)
names(getdata)[1] = "sample"
head(getdata)

#Get sample-wise mean for each gene
extratcted_gene_dat = do.call(cbind, lapply(2:ncol(getdata), function(i){
  #i=2
  print(i)
  gene_means = summarySE(getdata, measurevar=colnames(getdata)[i], groupvars=c("sample")) #calculate sample-wise means
  gene_means.1 = data.frame(gene_means[,3]) #calculate sample-wise means
  names(gene_means.1) = colnames(getdata)[i]
  row.names(gene_means.1) = gene_means$sample
  gene_means.1
}))

head(extratcted_gene_dat)
extratcted_gene_dat.mat = as.matrix(extratcted_gene_dat)
head(extratcted_gene_dat.mat)

colnames(extratcted_gene_dat.mat)

#unscaled matrix
heatmap.2(extratcted_gene_dat.mat, col = bluered(100), 
          trace = "none", density.info = "none",
          Colv = F,
          Rowv = F,
          lhei=c(1, 10),
          #labCol = T,
          dendrogram="row",
          main = paste("Unscaled sample-wise means")
)
dev.off()

#scaled matrix
extratcted_gene_dat.mat.scaled = scale(extratcted_gene_dat.mat)
#unscaled matrix
heatmap.2(extratcted_gene_dat.mat.scaled, col = bluered(100), 
          trace = "none", density.info = "none",
          Colv = T,
          Rowv = T,
          lhei=c(1, 10),
          #labCol = T,
          dendrogram="row",
          main = paste("Scaled sample-wise means")
)
dev.off()

#-----------------------------------------------------------------------------------------#
# Summed gene heatmaps: LOOK BETTER
#-----------------------------------------------------------------------------------------#

#Get data to put into heatmap
gene_n = 50 #number of genes from AD and PV sigs to show in heatmap
getdata = data.frame(colData(human_human_big_cds)$name_alphabet, 
                     t(human_human_big_cds@assays@data$counts[gene_sigs$ad_genes[2:gene_n],]),
                     t(human_human_big_cds@assays@data$counts[gene_sigs$pv_genes[1:gene_n],])
)
names(getdata)[1] = "sample"
head(getdata)

#Get sample-wise mean for each gene
extratcted_gene_dat = do.call(cbind, lapply(2:ncol(getdata), function(i){
  #i=2
  print(i)
  gene_means = aggregate(getdata[,i], by=list(Group=getdata$sample), FUN=sum) #calculate sample-wise means
  gene_means.1 = data.frame(gene_means[,2]) #calculate sample-wise means
  names(gene_means.1) = colnames(getdata)[i]
  row.names(gene_means.1) = gene_means$Group
  gene_means.1
}))

head(extratcted_gene_dat)
extratcted_gene_dat.mat = as.matrix(extratcted_gene_dat)
head(extratcted_gene_dat.mat)

#unscaled matrix
heatmap.2(extratcted_gene_dat.mat, col = bluered(100), 
          trace = "none", density.info = "none",
          Colv = F,
          Rowv = F,
          lhei=c(1, 10),
          #labCol = T,
          dendrogram="row",
          main = paste("Unscaled sample-wise means")
)

##**************
#THIS IS THE ONE THAT WORKS BEST
#scaled matrix
extratcted_gene_dat.mat.scaled = scale(extratcted_gene_dat.mat)

png("/Users/marktaylor/Desktop/Projects/Cho\ derm/Case\ Report/Figs/heatmap_summed_scaled_alphabetnames.png",
    width=16,height=5.5,units="in",res=300)

#SCALED matrix
print( heatmap.2(extratcted_gene_dat.mat.scaled, col = redblue(100), 
                 trace = "none", density.info = "none",
                 Colv = T,
                 Rowv = T,
                 labCol=as.expression(lapply(colnames(extratcted_gene_dat.mat.scaled), function(a) bquote(italic(.(a))))),
                 lhei=c(0, 10),
                 #labCol = T,
                 dendrogram="row",
                 main = paste("Scaled sample-wise sums"),
                 margins =c(6.5,5),
                 lwid = c(1.5,15),
                 cexRow = 2,
                 cexCol = 1.1,
                 offsetCol=-0.2
)
)
dev.off()
graphics.off()

#----------------------------------------------#
#Test Now calculate distance between each indeterminate sample, and the PV/AD samples

p1 = as.matrix(subset(re_int))
library(raster)

all_coords_mat = as.matrix(data.frame(re_int$ad_gene_sig,re_int$pv_gene_sig))
row.names(all_coords_mat) = re_int$sample
dist_mat = pointDistance(all_coords_mat, allpairs=T, lonlat=F) #calculate all-vs-all distance
rownames(dist_mat) = re_int$sample
colnames(dist_mat) = re_int$sample

ind_samps = subset(re_int, dis=="Indeterminate")$sample
pv_sample = 
  #Loop through the inderminate samples and get those sweet averages betch
  
  all_res = do.call(rbind, lapply(1:length(ind_samps), function(i){
    #i=1
    dist_df = data.frame(re_int$dis,re_int$sample, dist_mat[ind_samps[i],]) #get distance of a sample from the matrix and combine with metadata
    names(dist_df) = c("dis","sample","dist")
    
    ad_dist_df = subset(dist_df, dis=="AD") #get distances form 1 indeterminate sample to all AD samples
    pv_dist_df = subset(dist_df, dis=="PV") #get distances form 1 indeterminate sample to all PV samples
    
    #calculate means to see which sided test i wanna do
    ad_mean = mean(ad_dist_df$dist)
    pv_mean = mean(pv_dist_df$dist)
    
    #now test based on mean direction, find which is least
    if(ad_mean > pv_mean) {
      wiltest = wilcox.test(ad_dist_df$dist, pv_dist_df$dist, alternative="greater")
      wilrest = data.frame("PV",ad_mean,pv_mean,ad_mean-pv_mean,wiltest$statistic,wiltest$p.value)
    } else {wiltest = wilcox.test(ad_dist_df$dist, pv_dist_df$dist, alternative="less")
    wilrest = data.frame("AD",ad_mean,pv_mean,ad_mean-pv_mean,wiltest$statistic,wiltest$p.value)
    }
    
    #comine with data and return
    res = data.frame(ind_samps[i],wilrest)
    names(res) = c("sample","proximity","AD_dist_mean","PV_dist_mean","AD_PV","W","p")
    res
  }))


#What if I used the US and LS to discriminate among bulk samples:
#1. read in the data
malissadat = read.table(file="/Users/marktaylor/Desktop/Projects/Balmain\ lab/mouse_scRNAseq/data/bulk_network_data/Melissa\ data\ set\ including\ metastases\ from\ Nat\ Med\ paper/expr_all_probes_20140311.txt",header=T)
malissadat[1:10,1:10]
dim(malissadat)
malissaphenodat = read.delim(file="/Users/marktaylor/Desktop/Projects/Balmain\ lab/mouse_scRNAseq/data/bulk_network_data/Melissa\ data\ set\ including\ metastases\ from\ Nat\ Med\ paper/sample_attributes_updated_20140912.txt",header=T)
malissaphenodat[1:10,1:10]
dim(malissaphenodat)
malissaprobedat = read.delim(file="/Users/marktaylor/Desktop/Projects/Balmain\ lab/mouse_scRNAseq/data/bulk_network_data/Melissa\ data\ set\ including\ metastases\ from\ Nat\ Med\ paper/gene_attributes_MM10_all_probes.txt",header=T)
malissaprobedat[1:10,]

ls_markers = read.delim(file="/Users/marktaylor/Desktop/Projects/Balmain\ lab/mouse_scRNAseq/results/top_markers/cluster_33_vs_all_other_top_markers.txt", header=T)
us_markers = read.delim(file="/Users/marktaylor/Desktop/Projects/Balmain\ lab/mouse_scRNAseq/results/top_markers/cluster_19_vs_all_other_top_markers.txt", header=T)

#Get top 100 markers
ls_markers = ls_markers$gene_id[1:100]
us_markers = us_markers$gene_id[1:100]

#Get per sample signatures
ls_signature = t(na.omit(malissadat[match(ls_markers,malissaprobedat$SYMBOL),]))
ls_signature = ls_signature[-1,] #remove the first row which is the identifier
ls_signature_df = data.frame(rowSums(ls_signature))

us_signature = t(na.omit(malissadat[match(us_markers,malissaprobedat$SYMBOL),]))
us_signature = us_signature[-1,] #remove the first row which is the identifier
us_signature_df = data.frame(rowSums(us_signature))

#Add gene
malissaphenodat$ls_signature = ls_signature_df[,1]
malissaphenodat$us_signature = us_signature_df[,1]

#malissaphenodat = subset(malissaphenodat, TUMOR_TYPE!="SK")

ggplot(data=malissaphenodat, aes(x=ls_signature, y=us_signature)) +
  geom_point(alpha=1,aes(color=TUMOR_TYPE), size=2) +
  #geom_text_repel(data=re_int, aes(x=ad_gene_sig, y=pv_gene_sig,label=sample),box.padding = 0.4) + #add labels
  theme_classic() +
  stat_ellipse(geom = "polygon", alpha = 0.3, aes(fill = TUMOR_TYPE, color=TUMOR_TYPE)) +
  ggtitle("LS-US Gene Signature Space") + 
  xlab("LS gene signature") +
  ylab("US gene signature") +
  theme(legend.position="right") +
  scale_color_manual(name="",values = c("CA" = "orange",
                                        "MET" = "darkred",
                                        "PAP" = "purple",
                                        "SK" = "darkblue")) +
  scale_fill_manual(name="",values = c("CA" = "orange",
                                       "MET" = "darkred",
                                       "PAP" = "purple",
                                       "SK" = "darkblue")) +
  theme(axis.text=element_text(size=15,color='black'),
        axis.title=element_text(size=15,color='black'),
        plot.title=element_text(size=20,color='black'))  

#Infer relative amounts of certain cell types in each sample

#------------------------------------------------------------------------------------------------------#
#Repeat on the Tumor data
#------------------------------------------------------------------------------------------------------#


#What if I used the US and LS to discriminate among bulk samples:
#1. read in the data
aim2dat = read.table(file="/Users/marktaylor/Desktop/Projects/Balmain\ lab/mouse_scRNAseq/data/bulk_network_data/Aim2-4/AIM2_AIM4_BS_CARC_HFC_NOTED_ENTREZ_03_16_15_EXPRESSION.txt",header=T)
aim2dat[1:10,1:10]
dim(aim2dat)
aim2phenodat = read.table(file="/Users/marktaylor/Desktop/Projects/Balmain\ lab/mouse_scRNAseq/data/bulk_network_data/Aim2-4/AIM2_AIM4_BS_CARC_HFC_NOTED_ENTREZ_03_16_15_PHENOTYPE_DATA.txt",header=T)
aim2phenodat[1:10,1:10]
dim(aim2phenodat)
aim2probedat = read.table(file="/Users/marktaylor/Desktop/Projects/Balmain\ lab/mouse_scRNAseq/data/bulk_network_data/Aim2-4/AIM2_AIM4_BS_CARC_HFC_NOTED_ENTREZ_03_16_15_PROBE_DATA.txt",header=T)
aim2probedat[1:10,]

ls_markers = read.delim(file="/Users/marktaylor/Desktop/Projects/Balmain\ lab/mouse_scRNAseq/results/top_markers/cluster_33_vs_all_other_top_markers.txt", header=T)
us_markers = read.delim(file="/Users/marktaylor/Desktop/Projects/Balmain\ lab/mouse_scRNAseq/results/top_markers/cluster_19_vs_all_other_top_markers.txt", header=T)

#Get top 100 markers
ls_markers = ls_markers$gene_id[1:100]
us_markers = us_markers$gene_id[1:100]

#Get per sample signatures
ls_signature = t(na.omit(aim2dat[match(ls_markers,aim2probedat$SYMBOL),]))
ls_signature = ls_signature[-1,] #remove the first row which is the identifier
ls_signature_df = data.frame(rowSums(ls_signature))

us_signature = t(na.omit(aim2dat[match(us_markers,aim2probedat$SYMBOL),]))
us_signature = us_signature[-1,] #remove the first row which is the identifier
us_signature_df = data.frame(rowSums(us_signature))

#Add gene
aim2phenodat$ls_signature = ls_signature_df[,1]
aim2phenodat$us_signature = us_signature_df[,1]

#malissaphenodat = subset(malissaphenodat, TUMOR_TYPE!="SK")

unique(aim2phenodat$TISSUE)

ggplot(data=aim2phenodat, aes(x=ls_signature, y=us_signature)) +
  geom_point(alpha=1,aes(color=TISSUE), size=2) +
  #geom_text_repel(data=re_int, aes(x=ad_gene_sig, y=pv_gene_sig,label=sample),box.padding = 0.4) + #add labels
  theme_classic() +
  stat_ellipse(geom = "polygon", alpha = 0.3, aes(fill = TISSUE, color=TISSUE)) +
  ggtitle("LS-US Gene Signature Space") + 
  xlab("LS gene signature") +
  ylab("US gene signature") +
  theme(legend.position="right") +
  scale_color_manual(name="",values = c("CARC" = "orange",
                                        "BACK" = "darkblue")) +
  scale_fill_manual(name="",values = c("CARC" = "orange",
                                       "BACK" = "darkblue")) +
  theme(axis.text=element_text(size=15,color='black'),
        axis.title=element_text(size=15,color='black'),
        plot.title=element_text(size=20,color='black'))  


#p16 status
unique(aim2phenodat$GENOTYPE_INFERRED)
ggplot(data=aim2phenodat, aes(x=ls_signature, y=us_signature)) +
  geom_point(alpha=1,aes(color=GENOTYPE_INFERRED), size=2) +
  #geom_text_repel(data=re_int, aes(x=ad_gene_sig, y=pv_gene_sig,label=sample),box.padding = 0.4) + #add labels
  theme_classic() +
  stat_ellipse(geom = "polygon", alpha = 0.3, aes(fill = GENOTYPE_INFERRED, color=GENOTYPE_INFERRED)) +
  ggtitle("LS-US Gene Signature Space") + 
  xlab("LS gene signature") +
  ylab("US gene signature") +
  theme(legend.position="right") +
  scale_color_manual(name="",values = c("KO" = "brown",
                                        "WT" = "darkgreen")) +
  scale_fill_manual(name="",values = c("KO" = "brown",
                                       "WT" = "darkgreen")) +
  theme(axis.text=element_text(size=15,color='black'),
        axis.title=element_text(size=15,color='black'),
        plot.title=element_text(size=20,color='black'))  




