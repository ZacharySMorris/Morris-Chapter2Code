### Fossil Only Morphospace comparison ###

Extant_Only_CrocDorsal
Fossil_Only_CrocDorsal

Fossil_Only_CrocDorsal.gpa <- gpagen(Fossil_Only_CrocDorsal)
Fossil_Only_CrocDorsal.pca <- plotTangentSpace(Fossil_Only_CrocDorsal.gpa$coords, label = TRUE)

Extant_Only_CrocDorsal.gpa <- gpagen(Extant_Only_CrocDorsal)
Extant_on_Fossil_Only_CrocDorsal.pca <- Coords2PC(Extant_Only_CrocDorsal.gpa$coords,Fossil_Only_CrocDorsal.pca,Fossil_Only_CrocDorsal.gpa)

Extant_on_Fossil_Only_CrocDorsal.pca

PCA <- Fossil_Only_CrocDorsal.pca
PV <- Fossil_CrocDorsal.PlottingValues

Xlim<-c(floor(min(PCA$pc.scores[,1])*10)/10,ceiling(max(PCA$pc.scores[,1])*10)/10)
Ylim<-c(floor(min(PCA$pc.scores[,2])*10)/10,ceiling(max(PCA$pc.scores[,2])*10)/10)

# pdf("Crocodylia_Extant_Ontogeny_on_Combined_Morphospace.pdf",11,8.5,useDingbats = FALSE)

plot(0, 0, type = "n",
     xlim = Xlim,
     ylim = Ylim,
     xlab = paste("Principal Component 1 (", round(100*PCA$pc.summary$importance[2,"PC1"], digits = 1), "%)", sep = ""),
     ylab = paste("Principal Component 2 (", round(100*PCA$pc.summary$importance[2,"PC2"], digits = 1), "%)", sep = ""),
     axes = FALSE,
     frame.plot = FALSE,
     asp=F)

axis(1, round(seq(Xlim[1],Xlim[2],by=0.1),1), pos=Ylim[1])
axis(2, round(seq(Ylim[1],Ylim[2],by=0.1),1), pos=Xlim[1], las = 1)
clip(Xlim[1],Xlim[2],Ylim[1],Ylim[2])
abline(h=0, lty=3)
abline(v=0, lty=3)

clip(-0.4,0.4,-0.4,0.4)

points(PCA$pc.scores[,1], PCA$pc.scores[,2],
       pch=PV$shape,
       cex=PV$size,
       bg=alpha(PV$color, 0.75), asp=F)

points(Extant_on_Fossil_Only_CrocDorsal.pca[,1], Extant_on_Fossil_Only_CrocDorsal.pca[,2],
       pch=Extant_CrocDorsal.PlottingValues$shape,
       cex=Extant_CrocDorsal.PlottingValues$size,
       bg=alpha(Extant_CrocDorsal.PlottingValues$color, 0.75), asp=F)

x <- Extant_on_Fossil_Only_CrocDorsal.pca[,1]
y <- Extant_on_Fossil_Only_CrocDorsal.pca[,2]

z <- chull(x,y)
dfHull <- cbind(x[z],y[z])
linesHull <- rbind(dfHull, dfHull[1,])
polygon(linesHull, col = alpha("grey", 0.50))


Extant_CrocDorsal.PlottingValues


Test <- SubsettingGMM(Fossil_classifier,
                                 Fossil_Only_CrocDorsal.gpa,
                                 Fossil_Only_CrocDorsal.pca,
                                 Fossil_CrocDorsal.PlottingValues,
                                 "Clade2")

pdf("Fossil_Only_Morphospace_Clade2_PC3_PC4.pdf",11,8.5,useDingbats = FALSE)
PCA <- Fossil_Only_CrocDorsal.pca
PV <- Fossil_CrocDorsal.PlottingValues

for(i in 1:length(Test$taxa)){
  Xlim<-c(floor(min(PCA$pc.scores[,3])*10)/10,ceiling(max(PCA$pc.scores[,3])*10)/10)
  Ylim<-c(floor(min(PCA$pc.scores[,4])*10)/10,ceiling(max(PCA$pc.scores[,4])*10)/10)

  plot(0, 0, type = "n",
       xlim = Xlim,
       ylim = Ylim,
       xlab = paste("Principal Component 3 (", round(100*PCA$pc.summary$importance[2,"PC3"], digits = 1), "%)", sep = ""),
       ylab = paste("Principal Component 4 (", round(100*PCA$pc.summary$importance[2,"PC4"], digits = 1), "%)", sep = ""),
       axes = FALSE,
       frame.plot = FALSE,
       asp=F)

  axis(1, round(seq(Xlim[1],Xlim[2],by=0.1),1), pos=Ylim[1])
  axis(2, round(seq(Ylim[1],Ylim[2],by=0.1),1), pos=Xlim[1], las = 1)
  clip(Xlim[1],Xlim[2],Ylim[1],Ylim[2])
  abline(h=0, lty=3)
  abline(v=0, lty=3)

  clip(-0.4,0.4,-0.4,0.4)

  points(Test$PCvalues[[i]][,3], Test$PCvalues[[i]][,4],
         pch=Test$shape[[i]],
         cex=Test$size[[i]],
         bg=alpha(Test$color[[i]], 0.75), asp=F)

  mtext(Test$taxa[i])

}
dev.off()
