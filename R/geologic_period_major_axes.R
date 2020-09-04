### geologic period statistical comparisons ###

## load in species durations ##
load(file = "Chapter2Code/data/Species_Durations.Rdata")
##

Species_Durations

#Make a matrix of time boundaries for geologic periods
Geologic_Period_Ranges <- matrix(data =  c(252, 201, 145, 66, 23, 2.5, 201, 145, 66, 23, 2.5, 0),
                        nrow = length(c("Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "Recent")),
                        ncol = 2, byrow = FALSE)
rownames(Geologic_Period_Ranges) <- c("Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "Recent")
colnames(Geologic_Period_Ranges) <- c("max_age","min_age")
save(Geologic_Period_Ranges, file = "data/Geologic_Period_Ranges.Rdata")


##Create and then fill a list of species present in each geologic period
Stage_Bin_Taxa <- list()

for (i in 1:nrow(Geologic_Period_Ranges)){
  Current_oldest <- max(Geologic_Period_Ranges[i,])
  Current_youngest <- min(Geologic_Period_Ranges[i,])

  current_period_taxa <- c()
  for (j in 1:length(Species_Durations$Species)){
    Current_sp <- Species_Durations$Species[[j]]

    #create a vector which includes the sequence of values within the species duration range
    Sp_range_seq <- seq(Species_Durations$Maximum.Age..mya.[j],Species_Durations$Minimum.Age..mya.[j])

    #Test whether any of the numbers in the species duration rante fall within the period range
    Range_comparison <- Sp_range_seq >= Current_youngest & Sp_range_seq <= Current_oldest

    if (any(Range_comparison) == TRUE){
      current_period_taxa <- append(current_period_taxa,paste(Current_sp))

    } else next
  }
  Stage_Bin_Taxa[[rownames(Geologic_Period_Ranges)[[i]]]] <- current_period_taxa
}
##


##Make figures
PCA <- Combined_CrocDorsal.pca
PV <- Combined_CrocDorsal.PlottingValues

Xlim<-c(floor(min(PCA$pc.scores[,1])*10)/10,ceiling(max(PCA$pc.scores[,1])*10)/10)
Ylim<-c(floor(min(PCA$pc.scores[,2])*10)/10,ceiling(max(PCA$pc.scores[,2])*10)/10)
pdf("GeologicStages_Fossils_on_Combined_Morphospace_Fixed.pdf",11,8.5, useDingbats = FALSE)

for (i in 1:length(Stage_Bin_Taxa)){

  toMatch <- paste(Stage_Bin_Taxa[[i]],collapse = "|")

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

  points(PCA$pc.scores[grep(toMatch,classifier$Species),1], PCA$pc.scores[grep(toMatch,classifier$Species),2],
         pch=PV$shape[grep(toMatch,classifier$Species)],
         cex=PV$size[grep(toMatch,classifier$Species)],
         bg=alpha(PV$color[grep(toMatch,classifier$Species)], 0.75), asp=F)
}
dev.off()
##

##Add these as factors into the classifier objects
classifier$GeologicPeriods <- as.character(classifier$Species)

for (i in 1:length(levels(classifier$Species))){
  current_sp <- levels(classifier$Species)[i]
  current_specimens <- grep(current_sp, classifier$Species)

  current_sp_stage_bins <- melt(Stage_Bin_Taxa)[grep(current_sp,melt(Stage_Bin_Taxa)$value),2]
  classifier$GeologicPeriods[current_specimens] <- as.character(paste(unique(current_sp_stage_bins), sep = "", collapse = "|"))
}

classifier$GeologicPeriods <- as.character(classifier$GeologicPeriods)
levels(classifier$GeologicPeriods) <- rownames(Geologic_Period_Ranges)

Adult_CrocDorsal.classifier$GeologicPeriods <- classifier$GeologicPeriods[grep("Adult|Subadult",classifier$Age)]
levels(Adult_CrocDorsal.classifier$GeologicPeriods) <- rownames(Geologic_Period_Ranges)

Fossil_classifier$GeologicPeriods <- classifier$GeologicPeriods[grep("Fossil",classifier$Extant)]
levels(Fossil_classifier$GeologicPeriods) <- rownames(Geologic_Period_Ranges)
##

## create geologic period GMM subsets for adult/subadults and extinct only datasets
Adult_GeologicPeriods <- SubsettingGMM(Adult_CrocDorsal.classifier,
                                       Adult_CrocDorsal.gpa,
                                       Adult_CrocDorsal.pca,
                                       Adult_CrocDorsal.PlottingValues,
                                       "GeologicPeriods", print.plot=TRUE)

Fossil_GeologicPeriods <- SubsettingGMM(Fossil_classifier,
                                        Fossil_CrocDorsal.gpa,
                                        Fossil_CrocDorsal.pca,
                                        Fossil_CrocDorsal.PlottingValues,
                                        "GeologicPeriods")
##

## perform marjor axis calculation and comparisons
Adult_GeologicPeriods.MajorAxis <- resample.major.axis(Adult_GeologicPeriods, method = "bootstrap", iter = 9999)

Adult_GeologicPeriods.Slopes <- major.axis.slopes(Adult_GeologicPeriods.MajorAxis)

pdf("Adult_GeologicPeriods.MajorAxis_PC1_PC2.pdf",11,8.5,useDingbats = FALSE)
plot.major.axis(Adult_GeologicPeriods.Slopes,Adult_GeologicPeriods$PCvalues,Combined_CrocDorsal.pca,PCs = c(1,2))
dev.off()

Adult_GeologicPeriods.MA.comp <- major.axis.comparison(Adult_GeologicPeriods.MajorAxis$groups,
                                                       Adult_GeologicPeriods.MajorAxis$Transformed.MA,
                                                       Adult_GeologicPeriods.MajorAxis$resampled_transformed.MA,
                                                       MA_number = 1,  PCs = c(1:4), PC_comp = 1)

write.csv(Adult_GeologicPeriods.MA.comp$Results,"Adult_GeologicPeriods_MA_Results.csv")

Fossil_GeologicPeriods.MajorAxis <- resample.major.axis(Fossil_GeologicPeriods, method = "bootstrap", iter = 9999)
Fossil_GeologicPeriods_PC1_3.MajorAxis <- resample.major.axis(Fossil_GeologicPeriods, PCs=c(1:3), method = "bootstrap", iter = 9999)

Fossil_GeologicPeriods.MajorAxis.Slopes <- major.axis.slopes(Fossil_GeologicPeriods.MajorAxis)

pdf("Fossil_GeologicPeriods.MajorAxis_PC1_PC2.pdf",11,8.5,useDingbats = FALSE)
plot.major.axis(Fossil_GeologicPeriods.MajorAxis.Slopes,Fossil_GeologicPeriods$PCvalues,Combined_CrocDorsal.pca,PCs = c(1,2))
dev.off()

Fossil_GeologicPeriods.MA.comp <- major.axis.comparison(Fossil_GeologicPeriods.MajorAxis$groups,
                                                        Fossil_GeologicPeriods.MajorAxis$Transformed.MA,
                                                        Fossil_GeologicPeriods.MajorAxis$resampled_transformed.MA,
                                                        MA_number = 1,  PCs = c(1:4), PC_comp = 1)
write.csv(Fossil_GeologicPeriods.MA.comp$Results,"Fossil_GeologicPeriods_MA_Results.csv")
##

## perform post-hoc comparisons between geologic periods and extant ecology MA
Adult_GeologicPeriods_v_Ecology <- posthoc.major.axis.comparison(Adult_GeologicPeriods.MajorAxis,
                                                                 Maturity.MajorAxis,
                                                                 group_list = c("Adult_Subadult", Adult_GeologicPeriods.MajorAxis$groups),
                                                                 MA_number = 1,  PCs = c(1:4), PC_comp = 1)

temp_MA <- Adult_GeologicPeriods_v_Ecology
temp_results <- NULL
temp_results <- array(data=NA,
                      dim = dim(temp_MA$Adj_Pvalues)*c(1,2,1),
                      dimnames = list(dimnames(temp_MA$Slopes)[[1]],c(paste(dimnames(temp_MA$Slopes)[[2]],"Slopes"),paste(dimnames(temp_MA$Slopes)[[2]],"Pvalues")),dimnames(temp_MA$Slopes)[[3]]))

for (l in 1:dim(temp_results)[3]){
  temp_results[,,l] <- cbind(temp_MA$Slopes[,,l],temp_MA$Adj_Pvalues[,,l])
}

write.csv(temp_results,"Adult_GeologicPeriods_v_Ecology_MA_Results.csv")

Fossil_GeologicPeriods_v_Ecology <- posthoc.major.axis.comparison(Fossil_GeologicPeriods.MajorAxis,
                                                                  Maturity.MajorAxis,
                                                                  group_list = c("Adult_Subadult", Fossil_GeologicPeriods.MajorAxis$groups),
                                                                  MA_number = 1,  PCs = c(1:4), PC_comp = 1)

temp_MA <- Fossil_GeologicPeriods_v_Ecology
temp_results <- NULL
temp_results <- array(data=NA,
                      dim = dim(temp_MA$Adj_Pvalues)*c(1,2,1),
                      dimnames = list(dimnames(temp_MA$Slopes)[[1]],c(paste(dimnames(temp_MA$Slopes)[[2]],"Slopes"),paste(dimnames(temp_MA$Slopes)[[2]],"Pvalues")),dimnames(temp_MA$Slopes)[[3]]))

for (l in 1:dim(temp_results)[3]){
  temp_results[,,l] <- cbind(temp_MA$Slopes[,,l],temp_MA$Adj_Pvalues[,,l])
}

write.csv(temp_results,"Fossil_GeologicPeriods_v_Ecology_MA_Results.csv")
##

### Old code to resolve

PmaxGMM(Fossil_GeologicPeriods_Groups, Combined_CrocDorsal.pca, PCs = c(1,2), PrintPoints = TRUE)

Fossil_GeologicPeriods_Groups$legendcolor <- Fossil_CrocDorsal.PlottingValues$color[match(Fossil_GeologicPeriods_Total_Morphospace.MajorAxis_1$Groups,Fossil_classifier$GeologicPeriods)]
names(Fossil_GeologicPeriods_Groups$legendcolor) <- Fossil_GeologicPeriods_Total_Morphospace.MajorAxis_1$Groups

pdf("Fossil_GeologicPeriods_Groups_MA1_95CI.pdf",11,8.5,useDingbats = FALSE)
CI_Poly_Plot(Fossil_GeologicPeriods_Groups,
             Combined_CrocDorsal.pca,
             Fossil_GeologicPeriods_Total_Morphospace.MajorAxis_1_relaxed,
             Fossil_GeologicPeriods_Groups$legendcolor,
             OnePlot=FALSE)
dev.off()


#### ANGLE VERSIONS ####

##Function to statistically test pairwise angular comparisons between sets of group major axes##
Fossil_GeologicPeriods_Angle_MA<- MA_angular_comp(Fossil_GeologicPeriods_PC1_3.MajorAxis$groups,
                                               Fossil_GeologicPeriods_PC1_3.MajorAxis$Transformed.MA,
                                               Fossil_GeologicPeriods_PC1_3.MajorAxis$resampled_transformed.MA,
                                               MA_number = 1)

Fossil_GeologicPeriods_v_Ecology_Angle <- posthoc_MASA_comp(Fossil_GeologicPeriods_PC1_3.MajorAxis,
                                                            Maturity_PC1_3.MajorAxis,
                                                            group_list = c("Adult_Subadult", Fossil_GeologicPeriods_PC1_3.MajorAxis$groups),
                                                            MA_number = 1)
##

