### Phylogeny import and manipulation ###

## import phylogeny
Pseudosuchia.phy <- read.tree("Pseudosuchia_Grades_Final.phy")
save(Pseudosuchia.phy, file = "Chapter2Code/data/Pseudosuchia.phy.Rdata")
##

## import data for species geologic persistence
Species_Durations <- read.csv("Species_Geologic_Ages.csv",header = T)
save(Species_Durations, file = "Chapter2Code/data/Species_Durations.Rdata")

Species_Duration_Table <- as.matrix(Species_Durations[,-1])
rownames(Species_Duration_Table) <- gsub(" ", "_", Species_Durations$Species)
Species_Duration_Table <- Species_Duration_Table[match(Pseudosuchia.phy$tip.label,rownames(Species_Duration_Table)),]
save(Species_Duration_Table, file = "Chapter2Code/data/Species_Duration_Table.Rdata")
##

## make list of colors and symbols for phytips
phylo.colors <- Combined_CrocDorsal.PlottingValues$color[match(Pseudosuchia.phy$tip.label,gsub(" ", "_", classifier$Species))]
names(phylo.colors) <- Pseudosuchia.phy$tip.label
phylo.shapes <- Combined_CrocDorsal.PlottingValues$shape[match(Pseudosuchia.phy$tip.label,gsub(" ", "_", classifier$Species))]
names(phylo.shapes) <- Pseudosuchia.phy$tip.label
##

## make time calibrated and dichotomous phylogenies
Pseudosuchia_dichotomous.phy <- multi2di(collapse.singles(Pseudosuchia.phy), random = FALSE)
Pseudosuchia_calibrated.phy <- timePaleoPhy(Pseudosuchia_dichotomous.phy, Species_Duration_Table, type="equal", vartime = 0.11)
Pseudosuchia_calibrated.phy$edge.length <- Pseudosuchia_calibrated.phy$edge.length + 0.001
is.binary(Pseudosuchia_dichotomous.phy)
is.binary(Pseudosuchia_calibrated.phy)
##
Pseudosuchia_calibrated.phy$edge.length

# plot(Pseudosuchia.phy, show.tip.label = FALSE, show.node.label = FALSE, edge.width = 2); nodelabels(bg="white",adj=1.5)
# plot(Pseudosuchia_calibrated.phy, type = "phylogram",
#      use.edge.length = TRUE, align.tip.label = TRUE,
#      show.tip.label = FALSE, edge.width = 2)

#### plot scaled phylogeny ####
geoscalePhylo(Pseudosuchia_calibrated.phy, units = c("Period"), boxes = "Period",
              tick.scale = "Period", show.tip.label = FALSE, root.edge = TRUE,
              use.edge.length = TRUE, align.tip.label = FALSE)

lastPP <- get("last_plot.phylo",envir=.PlotPhyloEnv)

bar_lengths <- Species_Duration_Table[,1] - Species_Duration_Table[,2]
bar_rgb <- col2rgb(phylo.colors, alpha = 0.75)
rect(lastPP$xx[c(1:161)] + 0.1,
     lastPP$yy[c(1:161)] - 0.1,
     lastPP$xx[c(1:161)] + bar_lengths,
     lastPP$yy[c(1:161)] + 0.1,
     border = alpha(phylo.colors, 0.75)
     )

text(lastPP$xx[c(1:161)] + bar_lengths + 0.5,
     lastPP$yy[c(1:161)],
     labels = gsub("_", " ", Pseudosuchia_calibrated.phy$tip.label),
     cex = 0.3,
     adj = 0
)


clade_names <- c("Pseudosuchia","Crocodylomorpha","Neosuchia","Eusuchia","Crocodylia")
clade_nodes <- c(171,182,217,238,246)
clade_offsets <- c(-20,-18,-16,-14,-12)

for (i in 1:length(clade_names)){
  cladelabels(text = clade_names[i],node = clade_nodes[i], wing.length=-50, offset=clade_offsets[i])
}

clade_names <- c("Notosuchia","Thalattosuchia","Tethysuchia","Goniopholidae")
clade_nodes <- c(207,186,218,233)
clade_offsets <- c(4,4,4,4)

for (i in 1:length(clade_names)){
  cladelabels(text = clade_names[i],node = clade_nodes[i], wing.length=-50, offset=clade_offsets[i],orientation = "horizontal")
}

Gradelabels(text = c("Basal Pseudosuchia"),
            node = c(171),
            node_exclude = c(182),
            wing.length=-50, offset=4,
            orientation = "horizontal")
###




#### plot phylomorphospace plots ####
Mean_Adult_PC_Matrix <- Coords2PC(Mean_coords_phy_order,Combined_CrocDorsal.pca,Combined_CrocDorsal.gpa)


phy <- Pseudosuchia_calibrated.phy
morphospace <- Combined_CrocDorsal.pca

Mean_Adult_classifier <- classifier[match(Unique_Adult_Species_list,classifier$Species),]
Mean_Adult_classifier <- Mean_Adult_classifier[match(rownames(tip.data),gsub(" ", "_", Mean_Adult_classifier$Species)),]

tip.data <- Mean_Adult_PC_Matrix

  ##Do ancestral state reconstruction of PC data for Dataset with Gavialid and Non-Gavialid mean embryo shapes
  all.data <- multianc(phy,tip.data[,1:4]) %>%
  # Give the resulting marix proper rownames
  "rownames<-"(c(phy$tip.label,seq(from  = length(phy$tip.label)+1, to = length(phy$tip.label)+phy$Nnode, by = 1)))

PV$color <- Combined_CrocDorsal.PlottingValues$color[match(rownames(tip.data),gsub(" ", "_", classifier$Species))]
PV$shape <- Combined_CrocDorsal.PlottingValues$shape[match(rownames(tip.data),gsub(" ", "_", classifier$Species))]
PV$size <- Combined_CrocDorsal.PlottingValues$size[match(rownames(tip.data),gsub(" ", "_", classifier$Species))]

CladesNames <- c("Thalattosuchia",
                 "Notosuchia",
                 "Tethysuchia",
                 "Gavialoidea",
                 "Borealosuchidae",
                 "Alligatoroidea",
                 "Crocodyloidea"
                 )

#Set nodes for bases of clades for plotting
Key_Nodes <- c(186,207,218,247,255,261,290)

#Name groups for plotting
CladesNames <- c("Non-Pseudo",
                 "Basal Pseudosuchia",
                 "Basal Crocodylomorpha",
                 "Basal Neosuchia",
                 "Basal Eusuchia",
                 "Crocodylia"
                 )

#Set nodes for bases of clades for plotting
Key_Nodes <- c(162,171,182,217,238,246)

#Get nodes and edges for plotting
Nodes = list()
Edges = list()
for (i in 1:length(Key_Nodes)){
  if(i == length(Key_Nodes))
    X <- as.vector(c(Key_Nodes[i]:nrow(all.data)))
  else
    X <- as.vector(c(Key_Nodes[i]:(Key_Nodes[i+1]-1)))
  Nodes[[CladesNames[i]]] <- X
  Y <- as.vector(sapply(X, grep, x = phy$edge[,1]))
  Z <- as.vector(sapply(X, grep, x = phy$edge[,2]))
  Edges[[CladesNames[i]]] <- unique(c(unlist(Y),unlist(Z)))
}


Xlim<-c(floor(min(morphospace$pc.scores[,1])*10)/10,ceiling(max(morphospace$pc.scores[,1])*10)/10)
Ylim<-c(floor(min(morphospace$pc.scores[,2])*10)/10,ceiling(max(morphospace$pc.scores[,2])*10)/10)

pdf("Phylogeny_on_CombinedMorphospace_PC1_PC2.pdf",11,8.5,useDingbats = FALSE)

for (g in 1:length(Nodes)){

  # pdf(file=paste("PC1vPC2_Phylo_Names", names(Nodes)[[g]], ".pdf", sep = ""), width = 11, height = 8.5, useDingbats = FALSE)

  plot(0, 0, type = "n",
       xlim = Xlim,
       ylim = Ylim,
       xlab = paste("Principal Component 1 (", round(100*morphospace$pc.summary$importance[2,"PC1"], digits = 1), "%)", sep = ""),
       ylab = paste("Principal Component 2 (", round(100*morphospace$pc.summary$importance[2,"PC2"], digits = 1), "%)", sep = ""),
       axes = FALSE,
       frame.plot = FALSE,
       asp=F)

  axis(1, round(seq(Xlim[1],Xlim[2],by=0.1),1), pos=Ylim[1])
  axis(2, round(seq(Ylim[1],Ylim[2],by=0.1),1), pos=Xlim[1], las = 1)
  clip(Xlim[1],Xlim[2],Ylim[1],Ylim[2])
  abline(h=0, lty=3)
  abline(v=0, lty=3)

  clip(-0.4,0.4,-0.4,0.4)

  mtext(names(Nodes)[[g]], at = -0.25)

  for (i in c(Edges[[g]])){
    lines(all.data[(phy$edge[i,1:2]),1],
          all.data[(phy$edge[i,1:2]),2])
  }

  points(all.data[Nodes[[g]],1], all.data[Nodes[[g]],2], pch=21, bg = "grey", cex = 0.8)
  if(g == length(Nodes)) {
    points(all.data[min(Nodes[[g]])-1,1], all.data[min(Nodes[[g]])-1,2], pch=8, col = "light blue", cex = 1.5) #plot node for previous ancestral clade
    lines(all.data[c(min(Nodes[[g]]),min(Nodes[[g]])-1),1],all.data[c(min(Nodes[[g]]),min(Nodes[[g]])-1),2])
    # text(all.data[Nodes[[g]],1], all.data[Nodes[[g]],2], labels = Nodes[[g]], adj = 2)

    current_tips <- grep(CladesNames[g],Mean_Adult_classifier$No_Loricata)
    Mean_Adult_classifier$Species[current_tips]
    rownames(tip.data)[current_tips]
    points(tip.data[current_tips,1],
           tip.data[current_tips,2],
           pch=PV$shape[current_tips],
           cex = PV$size[current_tips],
           bg=alpha(PV$color[current_tips], 0.75))
  }else
  points(all.data[max(Nodes[[g]])+1,1], all.data[max(Nodes[[g]])+1,2], pch=8, col = "gold", cex = 1.5) #plot node for next decendent clade
  points(all.data[min(Nodes[[g]])-1,1], all.data[min(Nodes[[g]])-1,2], pch=8, col = "light blue", cex = 1.5) #plot node for previous ancestral clade
  lines(all.data[c(min(Nodes[[g]]),min(Nodes[[g]])-1),1],all.data[c(min(Nodes[[g]]),min(Nodes[[g]])-1),2])
  # text(all.data[Nodes[[g]],1], all.data[Nodes[[g]],2], labels = Nodes[[g]], adj = 2)

  current_tips <- grep(CladesNames[g],Mean_Adult_classifier$No_Loricata)
  Mean_Adult_classifier$Species[current_tips]
  rownames(tip.data)[current_tips]
  points(tip.data[current_tips,1],
         tip.data[current_tips,2],
         pch=PV$shape[current_tips],
         cex = PV$size[current_tips],
         bg=alpha(PV$color[current_tips], 0.75))
  # text(tip.data[current_tips,1],
  #      tip.data[current_tips,2],
  #      labels = rownames(tip.data)[current_tips], adj = 1)

  # dev.off()

}
dev.off()
##
