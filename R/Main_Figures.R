###Main Figures for Chapter 2###

###Figure 1###
##Panel A - Extant only on total morphospace
#Generate Plot#
PCA <- Combined_CrocDorsal.pca
PV <- Combined_CrocDorsal.PlottingValues

Xlim<-c(floor(min(PCA$pc.scores[,1])*10)/10,ceiling(max(PCA$pc.scores[,1])*10)/10)
Ylim<-c(floor(min(PCA$pc.scores[,2])*10)/10,ceiling(max(PCA$pc.scores[,2])*10)/10)

pdf("Crocodylia_Extant_Ontogeny_on_Combined_Morphospace.pdf",11,8.5,useDingbats = FALSE)

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

points(PCA$pc.scores[grep("Extant", classifier$Extant),1], PCA$pc.scores[grep("Extant", classifier$Extant),2],
       pch=PV$shape[grep("Extant", classifier$Extant)],
       cex=PV$size[grep("Extant", classifier$Extant)],
       bg=alpha(PV$color[grep("Extant", classifier$Extant)], 0.75), asp=F)

legend(-0.275, 0.2, Crocodylia_Ecomorphs$taxa,
       cex = 1,
       fill = alpha(Crocodylia_Ecomorphs$legendcolor, 0.75),
       bty = "n",
       horiz = TRUE)

legend(-0.275, 0.18, c("Alligatoridae", "Crocodylidae", "Gavialis", "Tomistoma"),
       cex = 1,
       pch = Crocodylia_Ecomorphs$legendshape,
       pt.bg = "black",
       pt.cex = 1.5,
       bty = "n",
       horiz = TRUE)

dev.off()

plot.major.axis(Clades2.MajorAxis.Slopes,Fossil_Morphospace_Clades2$PCvalues,Combined_CrocDorsal.pca,c(1,3))




##Panel B - Major Axes of ontogeny and Adult Ecology
Subset_data <- Adult_Ecomorphs
i='Crocodylia'
PCs <- c(1,2)
ConfIntList <- Adult_Total_Morphospace.MajorAxis_1
MA_number <- Adult_Total_Morphospace.MajorAxis_1$MA_number
LineColor <- unique(Adult_Ecomorphs$color$Crocodylia)

Xlim<-c(floor(min(PCA$pc.scores[,1])*10)/10,ceiling(max(PCA$pc.scores[,1])*10)/10)
Ylim<-c(floor(min(PCA$pc.scores[,2])*10)/10,ceiling(max(PCA$pc.scores[,2])*10)/10)

pdf("Crocodylia_Ontogny_and_Ecology_MajorAxes.pdf",11,8.5,useDingbats = FALSE)
CI_Poly_Plot(Crocodylia_Ecomorphs, Combined_CrocDorsal.pca, Crocodylia_Ecomorphs_Total_Morphospace.MajorAxis_1, Crocodylia_Ecomorphs$legendcolor, OnePlot=TRUE)

clip(-0.4,0.4,-0.4,0.4)
polygon(c(rev(ConfIntList$NEW_Xs[[i]]), ConfIntList$NEW_Xs[[i]]),
        c(rev(ConfIntList$preds_CIupper[[i]]), ConfIntList$preds_CIlower[[i]]),
        col = alpha('grey',0.6), border = NA)

clip(
  min(ConfIntList$NEW_Xs[[i]]),
  max(ConfIntList$NEW_Xs[[i]]),
  min(ConfIntList$preds_CIlower[[i]]),
  max(ConfIntList$preds_CIupper[[i]])
)

lines(ConfIntList$NEW_Xs[[i]], ConfIntList$preds_CIlower[[i]], lty = 'dashed', col = LineColor)
lines(ConfIntList$NEW_Xs[[i]], ConfIntList$preds_CIupper[[i]], lty = 'dashed', col = LineColor)
abline(ConfIntList$Intercept_list[[i]][MA_number],ConfIntList$Slope_list[[i]][MA_number], col = LineColor, lwd=3, lty=1 )


dataEllipse(Subset_data$PCvalues[[i]][,PCs[1]],
            Subset_data$PCvalues[[i]][,PCs[2]],
            add = TRUE, plot.points = FALSE, levels = c(0.75),
            col = LineColor, fill = FALSE)
dev.off()

pdf("Crocodylia_Ontogny_and_Ecology_MajorAxes_G_T_Unique.pdf",11,8.5,useDingbats = FALSE)
CI_Poly_Plot(Crocodylia_Ecomorphs_G_T_Unique, Combined_CrocDorsal.pca, Crocodylia_Ecomorphs_G_T_Unique_Total_Morphospace.MajorAxis_1, Crocodylia_Ecomorphs_G_T_Unique$legendcolor, OnePlot=TRUE)
dev.off()

pdf("Crocodylia_Ontogny_and_Ecology_MajorAxes_No_ExtremeGroups.pdf",11,8.5,useDingbats = FALSE)
CI_Poly_Plot(Crocodylia_Ecomorphs_Non_Extreme, Combined_CrocDorsal.pca, Crocodylia_Ecomorphs_Non_Extreme_Total_Morphospace.MajorAxis_1, Crocodylia_Ecomorphs_Non_Extreme$legendcolor, OnePlot=TRUE)
dev.off()

###

###Figure 2###
##Panel A - All Crocodylomorph Fossils on Extant Major Axes

##Panel B - Sub-plots of distinct clades
  #Basal Pseudosuchians
  #"Spehnosuchians"
  #Thalattosuchians
  #Notosuchians
  #Protosuchians
  #Tethysuchians
  #Gonipholids


###Allometric Plots###
GPA <- Combined_CrocDorsal.gpa
PCA <- Combined_CrocDorsal.pca
Subset_Groups <- Combined_Morphospace_Clades3
PC_Y <- 2
for(i in 1:length(Subset_Groups$taxa)){
  CS_log <- log(Subset_Groups$CSize[[i]])
  PCScores <- Subset_Groups$PCvalues[[i]]
  PV_color <- Subset_Groups$color[[i]]
  PV_size <- Subset_Groups$size[[i]]
  PV_shape <- Subset_Groups$shape[[i]]

  Xlim<-c(0,ceiling(max(log(GPA$Csize))))
  Ylim<-c(floor(min(PCA$pc.scores[,PC_Y])*10)/10,ceiling(max(PCA$pc.scores[,PC_Y])*10)/10)

  plot(0, 0, type = "n",
       xlim = Xlim,
       ylim = Ylim,
       xlab = paste("Centroid Size (log scale)"),
       ylab = paste("Principal Component ", PC_Y, " (", round(100*PCA$pc.summary$importance[2, PC_Y], digits = 1), "%)", sep = ""),
       axes = FALSE,
       frame.plot = FALSE,
       asp=F)

  axis(1, round(seq(Xlim[1],Xlim[2]),1), pos=Ylim[1])
  axis(2, round(seq(Ylim[1],Ylim[2],by=0.1),1), pos=Xlim[1], las = 1)
  clip(Xlim[1],Xlim[2],Ylim[1],Ylim[2])
  abline(h=0, lty=3)

  clip(0,Xlim[2],-0.4,0.4)

  points(CS_log, PCScores[,1],
         pch=PV_shape,
         cex=PV_size,
         bg=alpha(PV_color, 0.75), asp=F)
  mtext(Subset_Groups$taxa[[i]])

}

CS_log <- log(Combined_CrocDorsal.gpa$Csize)

sapply(Combined_CrocDorsal.gpa$coords, "dist")

Combined_HeadLengh <- c()
for(i in 1:length(Combined_CrocDorsal.gpa$Csize)){
  temp_dist <- dist(combined_CrocDorsal[,,i][c(1,4),])
  Combined_HeadLengh <- append(Combined_HeadLengh, temp_dist)
}


PCA <- Combined_CrocDorsal.pca
PV <- Combined_CrocDorsal.PlottingValues
PC_Y <- 1

Xlim<-c(0,ceiling(max(CS_log)))
Ylim<-c(floor(min(PCA$pc.scores[,PC_Y])*10)/10,ceiling(max(PCA$pc.scores[,PC_Y])*10)/10)

# pdf("Crocodylia_Extant_Ontogeny_on_Combined_Morphospace.pdf",11,8.5,useDingbats = FALSE)

plot(0, 0, type = "n",
     xlim = Xlim,
     ylim = Ylim,
     xlab = paste("Centroid Size (log scale)"),
     ylab = paste("Principal Component ", PC_Y, " (", round(100*PCA$pc.summary$importance[2, PC_Y], digits = 1), "%)", sep = ""),
     axes = FALSE,
     frame.plot = FALSE,
     asp=F)

axis(1, round(seq(Xlim[1],Xlim[2]),1), pos=Ylim[1])
axis(2, round(seq(Ylim[1],Ylim[2],by=0.1),1), pos=Xlim[1], las = 1)
clip(Xlim[1],Xlim[2],Ylim[1],Ylim[2])
abline(h=0, lty=3)

clip(0,Xlim[2],-0.4,0.4)

points(CS_log[grep("Extant", classifier$Extant)], PCA$pc.scores[grep("Extant", classifier$Extant),1],
       pch=PV$shape[grep("Extant", classifier$Extant)],
       cex=PV$size[grep("Extant", classifier$Extant)],
       bg=alpha(PV$color[grep("Extant", classifier$Extant)], 0.75), asp=F)

points(CS_log[grep("Fossil", classifier$Extant)], PCA$pc.scores[grep("Fossil", classifier$Extant),1],
       pch=PV$shape[grep("Fossil", classifier$Extant)],
       cex=PV$size[grep("Fossil", classifier$Extant)],
       bg=alpha(PV$color[grep("Fossil", classifier$Extant)], 0.75), asp=F)

text(CS_log[grep("Extant", classifier$Extant)], PCA$pc.scores[grep("Extant", classifier$Extant),1],
     names(Combined_CrocDorsal.gpa$Csize)[grep("Extant", classifier$Extant)])

log(Combined_CrocDorsal.gpa$Csize)[grep("osteolaemus osborni", names(Combined_CrocDorsal.gpa$Csize))]
classifier$image[grep("osteolaemus osborni", names(Combined_CrocDorsal.gpa$Csize))]

points(CS_log, PCA[,1],
       pch=PV_shape,
       cex=PV_size,
       bg=alpha(PV_color, 0.75), asp=F)

plot(Combined_HeadLengh[grep("Fossil", classifier$Extant)], Combined_CrocDorsal.gpa$Csize[grep("Fossil", classifier$Extant)],
     pch=PV$shape[grep("Fossil", classifier$Extant)],
     cex=PV$size[grep("Fossil", classifier$Extant)],
     bg=alpha(PV$color[grep("Fossil", classifier$Extant)], 0.75), asp=F)

plot(Combined_HeadLengh, Combined_CrocDorsal.gpa$Csize,
     pch=PV$shape,
     cex=PV$size,
     bg=alpha(PV$color, 0.75), asp=F)

plot(log(Combined_HeadLengh), PCA$pc.scores[,1],
     pch=PV$shape,
     cex=PV$size,
     bg=alpha(PV$color, 0.75), asp=F)

