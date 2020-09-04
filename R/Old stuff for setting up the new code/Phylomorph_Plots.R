##Plot phylomorphospace groups##

Fossils_on_Ontogeny_species
Fossils_on_Ontogeny.pca$pc.scores

Aggregate.pca <- aggregate.data.frame(Fossils_on_Ontogeny.pca$pc.scores,
                                                    by = list(Fossils_on_Ontogeny_species),
                                                    FUN = mean)
Agg_Fossils_on_Ontogeny.pca <- data.frame(Aggregate.pca[,-1], row.names = Aggregate.pca[,1])

na.omit(rownames(Agg_Fossils_on_Ontogeny.pca)[match(gsub("_", " ", molecule.phy$tip.label),rownames(Agg_Fossils_on_Ontogeny.pca))])

Phylo_Fossils_on_Ontogney.pca <- na.omit(Agg_Fossils_on_Ontogeny.pca[match(gsub("_", " ", molecule.phy$tip.label),rownames(Agg_Fossils_on_Ontogeny.pca)),])

Adult_species_list <- classifier$Species[grep("Adult",classifier$Age)]


phy <- Pseudosuchia.phy
N <- length(phy$tip.label)
Nnode <- phy$Nnode

x <- Combined_CrocDorsal.pca$pc.scores
anc.states <- NULL
for (i in 1:ncol(x)) {
  x1 <- x[, i]
  tmp <- vector()
  for (j in 1:Nnode + N) {
    a <- multi2di(root(phy, node = j))
    tmp[j - N] <- ace(x1, a, method = "pic")$ace[1]
  }
  anc.states <- cbind(anc.states, tmp)
}

colnames(anc.states) <- colnames(x)
row.names(anc.states) <- 1:length(tmp)
all.data <- rbind(x, anc.states)

row.names(all.data)

pdf("Pseudo_phy.pdf",height = 11,width = 8.5)
plot(phy); nodelabels(bg="white")
dev.off()

#Name groups for plotting
CladesNames <- c("Non-Pseudo",
                 "Pseudosuchia",
                 "Loricata",
                 "Crocodylomorpha",
                 "Neosuchia",
                 "Eusuchia",
                 "Crocodylia")

#Set nodes for bases of clades for plotting
Key_Nodes <- c(71,78,83,87,94,99,104)

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

for (g in 1:length(Nodes)){
  Xlim<-1.1*c(min(Ontogeny.pca$pc.scores[,1]),max(Ontogeny.pca$pc.scores[,1]))
  Ylim<-1.1*c(min(Ontogeny.pca$pc.scores[,2]),max(Ontogeny.pca$pc.scores[,2]))

  pdf(file=paste("PC1vPC2_Phylo_Names", names(Nodes)[[g]], ".pdf", sep = ""), width = 11, height = 8.5, useDingbats = FALSE)

  plot(0, 0, type = "n",
       xlim = c(-0.3,max(Xlim)),
       ylim = c(-0.3,0.2),
       xlab = paste("Principal Componenet 1 (", round(100*Ontogeny.pca$pc.summary$importance[2,"PC1"], digits = 1), "%)", sep = ""),
       ylab = paste("Principal Componenet 2 (", round(100*Ontogeny.pca$pc.summary$importance[2,"PC2"], digits = 1), "%)", sep = ""),
       axes = FALSE,
       frame.plot = FALSE,
       asp=F)

  axis(1, c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3), pos=-0.3)
  axis(2, c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3), pos=-0.3)
  clip(-0.3,max(Xlim),-0.3,0.2)
  abline(h=0, lty=3)
  abline(v=0, lty=3)

  mtext(names(Nodes)[[g]], at = -0.25)

  for (i in c(Edges[[g]])){
    lines(all.data[(phy$edge[i,1:2]),1],
          all.data[(phy$edge[i,1:2]),2])
    }

  points(all.data[Nodes[[g]],1], all.data[Nodes[[g]],2], pch=21, bg = "grey", cex = 0.8)
  points(all.data[max(Nodes[[g]])+1,1], all.data[max(Nodes[[g]])+1,2], pch=8, col = "gold", cex = 1.5)
  points(all.data[min(Nodes[[g]])-1,1], all.data[min(Nodes[[g]])-1,2], pch=8, col = "light blue", cex = 1.5)
  text(all.data[Nodes[[g]],1], all.data[Nodes[[g]],2], labels = Nodes[[g]], adj = 2)

  points(Fossils_on_Ontogeny_Clades$PCvalues[[CladesNames[g]]][,1],
         Fossils_on_Ontogeny_Clades$PCvalues[[CladesNames[g]]][,2],
         pch=Fossils_on_Ontogeny_Clades$shape[[CladesNames[g]]],
         cex = Fossils_on_Ontogeny_Clades$size[[CladesNames[g]]],
         bg=alpha(Fossils_on_Ontogeny_Clades$color[[CladesNames[g]]], 0.75))
  text(Fossils_on_Ontogeny_Clades$PCvalues[[CladesNames[g]]][,1],
       Fossils_on_Ontogeny_Clades$PCvalues[[CladesNames[g]]][,2],
       labels = Fossils_on_Ontogeny_Clades$species[[CladesNames[g]]], adj = 1)

  dev.off()

}
