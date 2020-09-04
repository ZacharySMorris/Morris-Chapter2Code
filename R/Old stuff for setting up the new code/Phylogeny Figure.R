###Phylogeny Figure###

molecule.phy <- read.nexus("Pseudosuchia_Grades_Final2.nex")

Pseudosuchia.phy <- read.tree("Pseudosuchia_Grades_Final.phy")

Willberg.phy <- read.nexus("habitat scorings with time calibrated tree.nex")

# Pseudosuchia.phy <- root(compute.brlen(collapse.singles(Pseudosuchia.phy)), node = 162)

Pseudosuchia_grades.phy <- molecule.phy[[2]]

pdf(file="Pseudosuchia_Phylogeny.pdf",useDingbats = FALSE, width=8.5, height=33)
plot(Pseudosuchia.phy); nodelabels(bg="white",adj=1.5)
dev.off()

Pseudosuchia_phylo_Ecology <- classifier$Godoy_Lifestyle[match(gsub("_", " ", Pseudosuchia.phy$tip.label), classifier$Species)]
names(Pseudosuchia_phylo_Ecology) <- Pseudosuchia.phy$tip.label

Ecology_Model <- rerootingMethod(Pseudosuchia.phy,Pseudosuchia_phylo_Ecology, model = "SYM")

pdf(file = "Pseudosuchia_Phylo_ACE2.pdf", width = 10, height = 28, useDingbats = FALSE)
plotTree(Pseudosuchia.phy, offset = 2, show.tip.label = FALSE)
nodelabels(pie=Ecology_Model$marginal.anc, cex = 0.5, piecol = c("#386CB0","#1B9E77","#A6761D"))
tiplabels(pch = 21,
          bg = c("#386CB0","#1B9E77","#A6761D")[Pseudosuchia_phylo_Ecology],
          cex = 2,
          adj = c(0.55,0.5))
dev.off()



#Phylogeny for Figure 1
plot(Pseudosuchia.phy, use.edge.length = TRUE, align.tip.label = FALSE, edge.width = 4, label.offset = 15)
cladelabels(text = c("Pseudosuchia"), node = c(74), offset=15)

cladelabels(text = c("Alligatoroidea"), node = c(105), wing.length=-50, offset=10)
cladelabels(text = c("Crocodyloidea"), node = c(117), wing.length=-50, offset=10)

Gradelabels(text = c("Basal Pseudosuchia","Basal Loricata","Basal Crocodylomorpha","Basal Neosuchia","Basal Eusuchia"),
            node = c(78,83,86,94,99),
            node_exclude = c(83,86,94,99,104),
            wing.length=-50, offset=5,
            orientation = "horizontal")

Gradelabels(text = c("Basal Loricata"), node = c(83), node_exclude = c(83), wing.length=-50, offset=10)
Gradelabels(text = c("Basal Crocodylomorpha"), node = c(86), node_exclude = c(83), wing.length=-50, offset=10)
Gradelabels(text = c("Basal Neosuchia"), node = c(94),node_exclude = c(83), wing.length=-50, offset=10)
Gradelabels(text = c("Basal Eusuchia"), node = c(99), node_exclude = c(83), wing.length=-50, offset=10)
cladelabels(text = c("Crocodylia"), node = c(104), wing.length=-50, offset=10)


Pseudosuchia.phy$root.time <- max(nodeHeights(Pseudosuchia.phy))

TimeRanges <- read.csv("Pseudosuchia_time_range.csv", header=T)
rownames(TimeRanges) <- paste(TimeRanges[,1])
TimeRanges <- TimeRanges[,2:3]
TimeRanges <- TimeRanges[match(Pseudosuchia.phy$tip.label, rownames(TimeRanges)),]

pdf("Pseudosuchia_grades_phylo.pdf",height=11,width=8.5)

geoscalePhylo(Pseudosuchia.phy, units = c("Period"), boxes = "Period",
              tick.scale = "Period", show.tip.label = TRUE,
              use.edge.length = TRUE, align.tip.label = FALSE)

Gradelabels(text = c("Basal Pseudosuchia","Basal Loricata","Basal Crocodylomorpha","Basal Neosuchia","Basal Eusuchia"),
            node = c(76,83,87,94,99),
            node_exclude = c(83,87,94,99,104),
            wing.length=-50, offset=5,
            orientation = "horizontal")
cladelabels(text = c("Crocodylia"), node = c(104), wing.length=-50, offset=5, orientation = "horizontal")





lastPP <- get("last_plot.phylo",envir=.PlotPhyloEnv)
points(x = rep(max(lastPP$xx)+35,length(lastPP$yy[1:68])), y = lastPP$yy[1:68],
       cex = 1,
       pch = molecule.phylo.shapes,
       bg = alpha(molecule.phylo.colors, 0.75)
)

cladelabels(text = c("Basal Pseudosuchia", "Basal Loricata", "Basal Crocodylomorpha",
                     "Basal Neosuchia", "Basal Eusuchia", "Crocodylia"),
            node = c(66,9,76,81,86,90), offset=0)


#data(timescales, envir = environment())
#timescale.rescaled <- timescales
#units = "Age"
#segments(min(timescale.rescaled[timescale.rescaled[,"Type"] == units[1], "Start"]), par()$usr[3],
#         max(timescale.rescaled[timescale.rescaled[,"Type"] == units[1], "End"]), par()$usr[3])

dev.off()
