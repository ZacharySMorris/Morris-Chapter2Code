###Convergence Testing###

##Set up and formatting of phylogenies##
##Need to run code from Total Morphospace to get the subsetted object##

##Set up and formatting of mean shape data##
Embryo_Species_list <- as.character(Maturity$species$`Mid-Stage Embryo`)
Embryo_Species_list <- append(Embryo_Species_list,as.character(Maturity_combined$species$`Late-Stage Embryo`))
Unique_Embryo_Species_list <- unique(Embryo_Species_list)

Adult_Species_list <- as.character(Maturity_combined$species$Adult)
Adult_Species_list <- append(Adult_Species_list,as.character(Maturity_combined$species$Subadult))
Unique_Adult_Species_list <- unique(Adult_Species_list)

Adult_coords <- abind(Maturity_combined$coords$Adult, Maturity_combined$coords$Subadult)
dimnames(Adult_coords)[[3]] <- Adult_Species_list

Embryo_coords <- abind(Maturity_combined$coords$`Mid-Stage Embryo`, Maturity_combined$coords$`Late-Stage Embryo`)
dimnames(Embryo_coords)[[3]] <- Embryo_Species_list

Adult_mean_coords <- array(dim = c(14,2,length(Unique_Adult_Species_list)), dimnames = list(c(1:14),c("X","Y"),gsub(" ", "_", Unique_Adult_Species_list)))
Embryo_mean_coords <- array(dim = c(14,2,length(Unique_Embryo_Species_list)), dimnames = list(c(1:14),c("X","Y"),paste(gsub(" ", "_", Unique_Embryo_Species_list),"Embryo",sep="_")))

for (i in 1:length(Unique_Embryo_Species_list)){
  current_sp <- Unique_Embryo_Species_list[[i]]
  current_specimen <- grep(current_sp, dimnames(Embryo_coords)[[3]])
  current_meanshape <- mshape(Embryo_coords[,,current_specimen])

  Embryo_mean_coords[,,i] <- current_meanshape
}

for (i in 1:length(Unique_Adult_Species_list)){
  current_sp <- Unique_Adult_Species_list[[i]]
  current_specimen <- grep(current_sp, dimnames(Adult_coords)[[3]])
  current_meanshape <- mshape(Adult_coords[,,current_specimen])

  Adult_mean_coords[,,i] <- current_meanshape
}

#Make Mean shapes for Gavialid and Non-Gavialid Embryos#
toMatch <- c("Gavialis gangeticus", "Tomistoma schlegelii")

Mean_NonGavialid_Embryo <- mshape(Embryo_mean_coords[,,-grep(paste(toMatch, collapse = "|"),Unique_Embryo_Species_list)])

Mean_Gavialid_Embryo <- mshape(Embryo_mean_coords[,,grep(paste(toMatch, collapse = "|"),Unique_Embryo_Species_list)])


##Sort out which species are represented by which ontogentic stages##
Embryo_phy_species <- gsub(" ", "_", Unique_Embryo_Species_list)
Embryo_tip_names <- paste(gsub(" ", "_", Unique_Embryo_Species_list),"Embryo",sep="_")

EmbryoTaxa.phy <- Pseudosuchia_dichotomous.phy

for (i in 1:length(Embryo_phy_species)){
  EmbryoTaxa.phy <- bind.tip(EmbryoTaxa.phy,
                             tip.label = Embryo_tip_names[[i]],
                             edge.length = 1,
                             where = c(match(Embryo_phy_species, Pseudosuchia_dichotomous.phy$tip.label))[[i]],
                             position = 0)
}

EmbryoTaxa.phy$edge.length <- EmbryoTaxa.phy$edge.length + 0.01

plot(EmbryoTaxa.phy, show.tip.label=FALSE)
tiplabels(EmbryoTaxa.phy$tip.label[grep(paste(Embryo_phy_species, collapse = "|"),EmbryoTaxa.phy$tip.label)],
          grep(paste(Embryo_phy_species, collapse = "|"),EmbryoTaxa.phy$tip.label),
          cex = 0.5, offset = -50)


toMatch <- c("Gavialis_gangeticus","Alligator_mississippiensis","Crocodylus_porosus")
Crown_tip_list <- grep(paste(toMatch, collapse = "|"),Pseudosuchia_dichotomous.phy$tip.label)

SingleMeanEmbryo.phy <- bind.tip(Pseudosuchia_dichotomous.phy,
                                 tip.label = "Mean_NonGavialid_Embryo",
                                 edge.length = 1,
                                 where = mrca.phylo(Pseudosuchia_dichotomous.phy,Crown_tip_list),
                                 position = 0)

SingleMeanEmbryo.phy <- multi2di(collapse.singles(SingleMeanEmbryo.phy), random = FALSE)
SingleMeanEmbryo.phy$edge.length <- SingleMeanEmbryo.phy$edge.length + 0.01

Single_Embryo_Crown_tip_list <- grep(paste(toMatch, collapse = "|"),SingleMeanEmbryo.phy$tip.label)
mrca.phylo(SingleMeanEmbryo.phy,Single_Embryo_Crown_tip_list)

plot(SingleMeanEmbryo.phy, show.tip.label=FALSE)
nodelabels("247",247)
tiplabels("Mean_NonGavialid_Embryo",grep("Mean_NonGavialid_Embryo",SingleMeanEmbryo.phy$tip.label), cex = 0.5, offset = 50)


DoubleMeanEmbryo.phy <- bind.tip(SingleMeanEmbryo.phy,
                                 tip.label = "Mean_Gavialid_Embryo",
                                 edge.length = 1,
                                 where = grep("Mean_NonGavialid_Embryo", SingleMeanEmbryo.phy$tip.label),
                                 position = 0)

DoubleMeanEmbryo.phy <- multi2di(DoubleMeanEmbryo.phy, random = FALSE)
DoubleMeanEmbryo.phy$edge.length[0.000 == DoubleMeanEmbryo.phy$edge.length] <- 0.01

plot(DoubleMeanEmbryo.phy, show.tip.label=FALSE)
nodelabels("247",247)
tiplabels("Mean_NonGavialid_Embryo",grep("Mean_NonGavialid_Embryo",DoubleMeanEmbryo.phy$tip.label), cex = 0.5, offset = 50)
tiplabels("Mean_Gavialid_Embryo",grep("Mean_Gavialid_Embryo",DoubleMeanEmbryo.phy$tip.label), cex = 0.5, offset = 50)


##Build MeanShape Arrays##
MeanEmbryo_Array <- array(dim = c(14,2,length(EmbryoTaxa.phy$tip.label)), dimnames = list(c(1:14),c("X","Y"),EmbryoTaxa.phy$tip.label))

Adult_means <- na.omit(match(dimnames(Adult_mean_coords)[[3]], dimnames(MeanEmbryo_Array)[[3]]))
Embryo_means <- na.omit(match(dimnames(Embryo_mean_coords)[[3]], dimnames(MeanEmbryo_Array)[[3]]))

MeanEmbryo_Array[,,Adult_means] <- Adult_mean_coords
MeanEmbryo_Array[,,Embryo_means] <- Embryo_mean_coords

MeanEmbryo_PC_Matrix <- Coords2PC(MeanEmbryo_Array,Combined_CrocDorsal.pca,Combined_CrocDorsal.gpa)

SingleMeanEmbryo_Array <- array(dim = c(14,2,length(SingleMeanEmbryo.phy$tip.label)), dimnames = list(c(1:14),c("X","Y"),SingleMeanEmbryo.phy$tip.label))

Single_Adult_means <- na.omit(match(dimnames(Adult_mean_coords)[[3]], dimnames(SingleMeanEmbryo_Array)[[3]]))

SingleMeanEmbryo_Array[,,Single_Adult_means] <- Adult_mean_coords
SingleMeanEmbryo_Array[,,"Mean_NonGavialid_Embryo"] <- Mean_NonGavialid_Embryo

SingleMeanEmbryo_PC_Matrix <- Coords2PC(SingleMeanEmbryo_Array,Combined_CrocDorsal.pca,Combined_CrocDorsal.gpa)

DoubleMeanEmbryo_Array <- array(dim = c(14,2,length(DoubleMeanEmbryo.phy $tip.label)), dimnames = list(c(1:14),c("X","Y"),DoubleMeanEmbryo.phy$tip.label))

Double_Adult_means <- na.omit(match(dimnames(Adult_mean_coords)[[3]], dimnames(DoubleMeanEmbryo_Array)[[3]]))

DoubleMeanEmbryo_Array[,,Double_Adult_means] <- Adult_mean_coords
DoubleMeanEmbryo_Array[,,"Mean_NonGavialid_Embryo"] <- Mean_NonGavialid_Embryo
DoubleMeanEmbryo_Array[,,"Mean_Gavialid_Embryo"] <- Mean_Gavialid_Embryo

DoubleMeanEmbryo_PC_Matrix <- Coords2PC(DoubleMeanEmbryo_Array,Combined_CrocDorsal.pca,Combined_CrocDorsal.gpa)

####Run convevol functions finally####

dwarf_species <- c()

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
fossils_near_embryos
fossils_near_gavialid_embryos <- c("Metriorhynchus_brachyrhynchus",
                                   "Suchodus_durobrivensis",
                                   "Dakosaurus_andiniensis",
                                   "Machimosaurus_mosae",
                                   "Metriorhynchus_casamiquelai",
                                   "Dakosaurus_maximus",
                                   "Effigia_okeeffeae",
                                   "Shuvosaurus_inexpectatus")

#Do ancestral state reconstruction of PC data for Dataset with all embryo species


is.binary(EmbryoTaxa.phy)

MeanEmbryo_ALL_PCDATA = multianc(EmbryoTaxa.phy, MeanEmbryo_PC_Matrix[,1:6]) %>%
  # Give the resulting marix proper rownames
  "rownames<-"(c(EmbryoTaxa.phy$tip.label,seq(from  = length(EmbryoTaxa.phy$tip.label)+1, to = length(EmbryoTaxa.phy$tip.label)+EmbryoTaxa.phy$Nnode, by = 1)))

#Calculate REALLY big ditance matrix (distances between all nodes, both tips and internal)
MeanEmbryo_dist = as.matrix(dist((MeanEmbryo_ALL_PCDATA))) %>%
  # Melt the distance matrix into a more workable for format for the revised C-metric code
  melt(as.matrix(MeanEmbryo_dist), varnames = c("n1", "n2"))

#Simulate PC evolution on Single embryo mean shape phylogeny
nsim=1000
MeanEmbryo_PCs.sim = make.simdata(nsim=nsim,phy=EmbryoTaxa.phy,trait=MeanEmbryo_PC_Matrix[,1:6])

#Calculate C-metrics scores for the single embryo mean shape to each convergent taxon
MeanEmbryo_conv.df = c()

conv.df = c()
for (i in 1:length(Embryo_tip_names)){
  t1 = as.character(Embryo_tip_names[i])
  for (j in 1:length(fossils_near_embryos)){
    t2 = as.character(fossils_near_embryos[j])
    tx = data.frame(calc.conv.rates.sig(MeanEmbryo_dist,
                                        EmbryoTaxa.phy,
                                        t1,t2,
                                        MeanEmbryo_ALL_PCDATA,
                                        MeanEmbryo_PCs.sim))
    ty = cbind(t1,t2,tx)
    MeanEmbryo_conv.df = rbind(MeanEmbryo_conv.df,ty)
  }
}

p.adjust(DoubleMeanEmbryo_conv.df$C1sig, "holm")

MeanEmbryo_conv_SigOnly.df = subset(MeanEmbryo_conv.df , C1sig  < 0.0002777777778)
##


SingleMeanEmbryo.phy$tip.label
rownames(SingleMeanEmbryo_PC_Matrix)


##Do ancestral state reconstruction of PC data for Dataset with single mean embryo shape
SingleMeanEmbryo_ALL_PCDATA = multianc(SingleMeanEmbryo.phy,SingleMeanEmbryo_PC_Matrix[,1:6]) %>%
  # Give the resulting marix proper rownames
  "rownames<-"(c(SingleMeanEmbryo.phy$tip.label,seq(from  = length(SingleMeanEmbryo.phy$tip.label)+1, to = length(SingleMeanEmbryo.phy$tip.label)+SingleMeanEmbryo.phy$Nnode, by = 1)))

#Calculate REALLY big ditance matrix (distances between all nodes, both tips and internal)
SingleMeanEmbryo_dist = as.matrix(dist((SingleMeanEmbryo_ALL_PCDATA))) %>%
  # Melt the distance matrix into a more workable for format for the revised C-metric code
  melt(as.matrix(SingleMeanEmbryo_dist), varnames = c("n1", "n2"))

#Simulate PC evolution on Single embryo mean shape phylogeny
nsim=1000
SingleMeanEmbryo_PCs.sim = make.simdata(nsim=nsim,phy=SingleMeanEmbryo.phy,trait=SingleMeanEmbryo_PC_Matrix[,1:6])

#Calculate C-metrics scores for the single embryo mean shape to each convergent taxon
SingleMeanEmbryo_conv.df = c()
for (i in 1:length(fossils_near_embryos)){
  t1 = as.character("Mean_NonGavialid_Embryo")
  t2 = as.character(fossils_near_embryos[i])
  tx = data.frame(calc.conv.rates.sig(SingleMeanEmbryo_dist,
                                      SingleMeanEmbryo.phy,
                                      t1,t2,
                                      SingleMeanEmbryo_ALL_PCDATA,
                                      SingleMeanEmbryo_PCs.sim))
  ty = cbind(t1,t2,tx)
  SingleMeanEmbryo_conv.df = rbind(SingleMeanEmbryo_conv.df,ty)
}

SingleMeanEmbryo_conv_SigOnly.df = subset(SingleMeanEmbryo_conv.df , C1sig  < 0.005)
##

is.binary(DoubleMeanEmbryo.phy)

##Do ancestral state reconstruction of PC data for Dataset with Gavialid and Non-Gavialid mean embryo shapes
DoubleMeanEmbryo_ALL_PCDATA = multianc(DoubleMeanEmbryo.phy,DoubleMeanEmbryo_PC_Matrix[,1:6]) %>%
  # Give the resulting marix proper rownames
  "rownames<-"(c(DoubleMeanEmbryo.phy$tip.label,seq(from  = length(DoubleMeanEmbryo.phy$tip.label)+1, to = length(DoubleMeanEmbryo.phy$tip.label)+DoubleMeanEmbryo.phy$Nnode, by = 1)))

#Calculate REALLY big ditance matrix (distances between all nodes, both tips and internal)
DoubleMeanEmbryo_dist = as.matrix(dist((DoubleMeanEmbryo_ALL_PCDATA))) %>%
  # Melt the distance matrix into a more workable for format for the revised C-metric code
  melt(as.matrix(DoubleMeanEmbryo_dist), varnames = c("n1", "n2"))

#Simulate PC evolution on Single embryo mean shape phylogeny
nsim=1000
DoubleMeanEmbryo_PCs.sim = make.simdata(nsim=nsim,phy=DoubleMeanEmbryo.phy,trait=DoubleMeanEmbryo_PC_Matrix[,1:6])

#Calculate C-metrics scores for the single embryo mean shape to each convergent taxon
DoubleMeanEmbryo_conv.df = c()
for (i in 1:length(fossils_near_embryos)){
  t1 = as.character("Mean_NonGavialid_Embryo")
  t2 = as.character(fossils_near_embryos[i])
  tx = data.frame(calc.conv.rates.sig(DoubleMeanEmbryo_dist,
                                      DoubleMeanEmbryo.phy,
                                      t1,t2,
                                      DoubleMeanEmbryo_ALL_PCDATA,
                                      DoubleMeanEmbryo_PCs.sim))
  ty = cbind(t1,t2,tx)
  DoubleMeanEmbryo_conv.df = rbind(DoubleMeanEmbryo_conv.df,ty)
}

for (i in 1:length(fossils_near_gavialid_embryos)){
  t1 = as.character("Mean_Gavialid_Embryo")
  t2 = as.character(fossils_near_gavialid_embryos[i])
  tx = data.frame(calc.conv.rates.sig(DoubleMeanEmbryo_dist,
                                      DoubleMeanEmbryo.phy,
                                      t1,t2,
                                      DoubleMeanEmbryo_ALL_PCDATA,
                                      DoubleMeanEmbryo_PCs.sim))
  ty = cbind(t1,t2,tx)
  DoubleMeanEmbryo_conv.df = rbind(DoubleMeanEmbryo_conv.df,ty)
}

p.adjust(0.0056, n=8)
DoubleMeanEmbryo_conv_SigOnly.df = subset(DoubleMeanEmbryo_conv.df , C1sig  < 0.00625)
##



##RRphylo versions##
data("DataFelids")
DataFelids$PCscoresfel->PCscoresfel
DataFelids$treefel->treefel
DataFelids$statefel->statefel

RRphylo(treefel,PCscoresfel)->RRfel

search.conv(RR=RRfel, y=PCscoresfel, min.dim=5, min.dist="node9", foldername = tempdir())

plot(treefel, cex = 0.7)

Mean_Embryo_RRphylo <- RRphylo(EmbryoTaxa.phy, MeanEmbryo_PC_Matrix[,1:6])
 Mean_Embryo_RRphylo

 MeanEmbryo_Sub_PC_Matrix <- MeanEmbryo_PC_Matrix[,1:6]

search.conv(Mean_Embryo_RRphylo, tree=EmbryoTaxa.phy, y=MeanEmbryo_Sub_PC_Matrix,
            min.dim=5, min.dist="node9",
            foldername = tempdir())



Single_Mean_Embryo_RRphylo <- RRphylo(SingleMeanEmbryo.phy, SingleMeanEmbryo_PC_Matrix[,1:6])
Single_Mean_Embryo_RRphylo

search.conv(Mean_Embryo_RRphylo, SingleMeanEmbryo.phy, SingleMeanEmbryo_PC_Matrix[,c(1:6)],
            nodes = match(c("Mean_NonGavialid_Embryo", fossils_near_embryos),SingleMeanEmbryo.phy$tip.label))


##Wireframes for significantly interesting taxa
Paedomorphic_taxa <- c("Araripesuchus_patagonicus",
                       "Gobiosuchus_kielanae",
                       "Knoetschkesuchus_langenbergensis",
                       "Procaimanoidea_kayi",
                       "Protosuchus_sp.",
                       "Protosuchus_richardsoni",
                       "Navajosuchus_mooki",
                       "Metriorhynchus_brachyrhynchus",
                       "Shuvosaurus_inexpectatus")

Paedomorphic_taxa_coords <- Adult_mean_coords[,,match(Paedomorphic_taxa, dimnames(Adult_mean_coords)[[3]])]
Mean_NonGavialid_Embryo_rotated <- rotate.LM(Mean_NonGavialid_Embryo, 4, 1, FALSE)
Mean_Gavialid_Embryo_rotated <- rotate.LM(Mean_Gavialid_Embryo, 4, 1, FALSE)

ref <- mshape(Combined_CrocDorsal.gpa$coords)
pdf(file = "Paedomorphic_GlobalMean_Warps.pdf", useDingbats = FALSE)
plotRefToTarget(rotate.LM(ref, 4, 1, plotLM = FALSE), Mean_NonGavialid_Embryo_rotated, method = "vector", links=Dorsal_links)
mtext("Mean to Non-gavialoid Embryo")
plotRefToTarget(rotate.LM(ref, 4, 1, plotLM = FALSE), Mean_Gavialid_Embryo_rotated, method = "vector", links=Dorsal_links)
mtext("Mean to Gavialoid Embryo")
for (i in 1:length(Paedomorphic_taxa)){
  Paedomorphic_sp_Rotated <- rotate.LM(Paedomorphic_taxa_coords[,,i], LM1=4, LM2=1, plotLM=FALSE)
  plotRefToTarget(rotate.LM(ref, 4, 1, plotLM = FALSE), Paedomorphic_sp_Rotated, method = "vector", links=Dorsal_links)
  mtext(gsub("_"," ", Paedomorphic_taxa[[i]]))
  }
dev.off()

rotate.LM(Paedomorphic_taxa_coords[,,1], 1, 4, plotLM = TRUE)
