### Embryonic Region Convergence Testing ###

#### Create mean shapes for Non-gavialid and Gavialid embryos ####
## list of species with embryos
Extant_Maturity$species$Embryo
Unique_Embryo_Species_list <- unique(Extant_Maturity$species$Embryo)
##
## subset embryonic GPA coords
Embryo_coords <- Extant_Maturity$coords$Embryo
dimnames(Embryo_coords)[[3]] <- Extant_Maturity$species$Embryo
##

## create an array for the mean shapes of each species
Embryo_mean_coords <- array(dim = c(14,2,length(Unique_Embryo_Species_list)), dimnames = list(c(1:14),c("X","Y"),paste(gsub(" ", "_", Unique_Embryo_Species_list),"Embryo",sep="_")))
##
## iteratively collect coordinates for each species and save the mean embryonic shape per species
for (i in 1:length(Unique_Embryo_Species_list)){
  current_sp <- Unique_Embryo_Species_list[[i]]
  current_specimen <- grep(current_sp, dimnames(Embryo_coords)[[3]])
  current_meanshape <- mshape(Embryo_coords[,,current_specimen])

  Embryo_mean_coords[,,i] <- current_meanshape
}
##
#Make Mean shapes for Gavialid and Non-Gavialid Embryos#
toMatch <- c("Gavialis gangeticus", "Tomistoma schlegelii")

Mean_NonGavialid_Embryo <- mshape(Embryo_mean_coords[,,-grep(paste(toMatch, collapse = "|"),Unique_Embryo_Species_list)])

Mean_Gavialid_Embryo <- mshape(Embryo_mean_coords[,,grep(paste(toMatch, collapse = "|"),Unique_Embryo_Species_list)])
##



#### Create mean shapes for all species Adults/Subadults ####

#### Create Version of Phylogeny with Mean Embryonic Tips ####
## get crown node##
toMatch <- c("Gavialis_gangeticus","Alligator_mississippiensis","Crocodylus_porosus")
Crown_tip_list <- grep(paste(toMatch, collapse = "|"),Pseudosuchia_calibrated.phy$tip.label)
Crown_node <- mrca.phylo(Pseudosuchia_calibrated.phy,Crown_tip_list)
##
## add two tips as sister taxa at the base of crocodylia
MeanEmbryo.phy <- bind.tip(Pseudosuchia_calibrated.phy,
                                 tip.label = "Mean_Embryo",
                                 edge.length = 1,
                                 where = Crown_node,
                                 position = 0)

MeanEmbryo.phy <- bind.tip(MeanEmbryo.phy,
                                 tip.label = "Mean_Gavialid_Embryo",
                                 edge.length = 1,
                                 where = grep("Mean_NonGavialid_Embryo", MeanEmbryo.phy$tip.label),
                                 position = 0)
##
## resolve dichotomy and give meaningless branch lengths to the newly created tips
MeanEmbryo.phy <- multi2di(MeanEmbryo.phy, random = FALSE)
MeanEmbryo.phy$edge.length[0.000 == MeanEmbryo.phy$edge.length] <- 0.01
##
## checking phylogeny was generated correctly
# plot(MeanEmbryo.phy, show.tip.label=FALSE)
# nodelabels("Crown_node", Crown_node)
# tiplabels("Mean_NonGavialid_Embryo",grep("Mean_NonGavialid_Embryo",MeanEmbryo.phy$tip.label), cex = 0.5, offset = 50)
# tiplabels("Mean_Gavialid_Embryo",grep("Mean_Gavialid_Embryo",MeanEmbryo.phy$tip.label), cex = 0.5, offset = 0)
##






#### Create mean PC score matrix for convergence testing ####
## make and fill array containing mean GPA cords ##
MeanEmbryo_Array <- array(dim = c(14,2,length(MeanEmbryo.phy$tip.label)), dimnames = list(c(1:14),c("X","Y"),MeanEmbryo.phy$tip.label))

Adult_means <- na.omit(match(dimnames(Mean_coords_phy_order)[[3]], dimnames(MeanEmbryo_Array)[[3]]))

MeanEmbryo_Array[,,Adult_means] <- Mean_coords_phy_order
MeanEmbryo_Array[,,"Mean_Embryo"] <- Mean_NonGavialid_Embryo
GavialisEmbryo_Array <- MeanEmbryo_Array
GavialisEmbryo_Array[,,"Mean_Embryo"] <- Mean_Gavialid_Embryo
##
## turn gpa coords into a matrix of PC scores
MeanEmbryo_PC_Matrix <- Coords2PC(MeanEmbryo_Array,Combined_CrocDorsal.pca,Combined_CrocDorsal.gpa)
GavialisEmbryo_PC_Matrix <- Coords2PC(GavialisEmbryo_Array,Combined_CrocDorsal.pca,Combined_CrocDorsal.gpa)

##


#### Perform convergence analysis ####
## create vectors of extinct species to be compared
fossils_near_embryos <- c("Simosuchus_clarki",
                         "Araripesuchus_patagonicus",
                         "Gobiosuchus_kielanae",
                         "Knoetschkesuchus_langenbergensis",
                         "Acynodon_iberoccitanus",
                         "Procaimanoidea_kayi",
                         "Protosuchus_richardsoni",
                         "Protosuchus_sp.",
                         "Navajosuchus_mooki",
                         "Acynodon_adriaticus",
                         "Mariliasuchus_amarali",
                         "Notosuchus_terrestris")

fossils_near_gavialid_embryos <- c("Metriorhynchus_brachyrhynchus",
                                   "Suchodus_durobrivensis",
                                   "Dakosaurus_andiniensis",
                                   "Machimosaurus_mosae",
                                   "Metriorhynchus_casamiquelai",
                                   "Dakosaurus_maximus",
                                   "Effigia_okeeffeae",
                                   "Shuvosaurus_inexpectatus")
##

##Do ancestral state reconstruction of PC data for Dataset with Gavialid and Non-Gavialid mean embryo shapes
MeanEmbryo_ALL_PCDATA = multianc(MeanEmbryo.phy,MeanEmbryo_PC_Matrix[,1:4]) %>%
  # Give the resulting marix proper rownames
  "rownames<-"(c(MeanEmbryo.phy$tip.label,seq(from  = length(MeanEmbryo.phy$tip.label)+1, to = length(MeanEmbryo.phy$tip.label)+MeanEmbryo.phy$Nnode, by = 1)))

#Calculate REALLY big ditance matrix (distances between all nodes, both tips and internal)
MeanEmbryo_dist = as.matrix(dist((MeanEmbryo_ALL_PCDATA))) %>%
  # Melt the distance matrix into a more workable for format for the revised C-metric code
  melt(as.matrix(MeanEmbryo_dist), varnames = c("n1", "n2"))

#Simulate PC evolution on Single embryo mean shape phylogeny
nsim=1000
MeanEmbryo_PCs.sim = make.simdata(nsim=nsim,phy=MeanEmbryo.phy,trait=MeanEmbryo_PC_Matrix[,1:4])

#Calculate C-metrics scores for the single embryo mean shape to each convergent taxon
CER_conv.df = c()
for (i in 1:length(fossils_near_embryos)){
  t1 = as.character("Mean_Embryo")
  t2 = as.character(fossils_near_embryos[i])
  tx = data.frame(calc.conv.rates.sig(MeanEmbryo_dist,
                                      MeanEmbryo.phy,
                                      t1,t2,
                                      MeanEmbryo_ALL_PCDATA,
                                      MeanEmbryo_PCs.sim))
  ty = cbind(t1,t2,tx)
  CER_conv.df = rbind(CER_conv.df,ty)
}

p.adjust(0.0041, n=length(fossils_near_embryos))
CER_conv_SigOnly.df = subset(CER_conv.df , C1sig  < 0.0041)


##Do ancestral state reconstruction of PC data for Dataset with Gavialid and Non-Gavialid mean embryo shapes
GavialisEmbryo_ALL_PCDATA = multianc(MeanEmbryo.phy,GavialisEmbryo_PC_Matrix[,1:4]) %>%
  # Give the resulting marix proper rownames
  "rownames<-"(c(MeanEmbryo.phy$tip.label,seq(from  = length(MeanEmbryo.phy$tip.label)+1, to = length(MeanEmbryo.phy$tip.label)+MeanEmbryo.phy$Nnode, by = 1)))

#Calculate REALLY big ditance matrix (distances between all nodes, both tips and internal)
GavialisEmbryo_dist = as.matrix(dist((GavialisEmbryo_ALL_PCDATA))) %>%
  # Melt the distance matrix into a more workable for format for the revised C-metric code
  melt(as.matrix(GavialisEmbryo_dist), varnames = c("n1", "n2"))

#Simulate PC evolution on Single embryo mean shape phylogeny
nsim=1000
GavialisEmbryo_PCs.sim = make.simdata(nsim=nsim,phy=MeanEmbryo.phy,trait=GavialisEmbryo_PC_Matrix[,1:4])


GT_Embryo_conv.df = c()
for (i in 1:length(fossils_near_gavialid_embryos)){
  t1 = as.character("Mean_Embryo")
  t2 = as.character(fossils_near_gavialid_embryos[i])
  tx = data.frame(calc.conv.rates.sig(GavialisEmbryo_dist,
                                      MeanEmbryo.phy,
                                      t1,t2,
                                      GavialisEmbryo_ALL_PCDATA,
                                      GavialisEmbryo_PCs.sim))
  ty = cbind(t1,t2,tx)
  GT_Embryo_conv.df = rbind(GT_Embryo_conv.df,ty)
}

p.adjust(0.00625, n=length(fossils_near_gavialid_embryos))
GT_Embryo_conv_SigOnly.df = subset(GT_Embryo_conv.df , C1sig  < 0.00625)
##


