###Calculation and comaprisons of the Major Axes of Varaition across time###

##Take PBDB data and create time and geologic period factors for performing time based Major Axis analysis

#Make a matrix of time boundaries for geologic periods
Geologic_time <- matrix(data =  c(252, 201, 145, 66, 23, 2.5, 201, 145, 66, 23, 2.5, 0),
                        nrow = length(c("Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "Recent")),
                        ncol = 2, byrow = FALSE)
rownames(Geologic_time) <- c("Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "Recent")
colnames(Geologic_time) <- c("max_age","min_age")

#Make a matrix of time boundaries for 25Myr timebins
Time_Bins <- matrix(data = c(as.numeric(seq(from = 250, to = 0, length.out = 11))[1:10],
                             as.numeric(seq(from = 250, to = 0, length.out = 11))[2:11]),
                    nrow = 10,
                    ncol = 2,
                    byrow = FALSE)

rownames(Time_Bins) <- paste("Timebin", LETTERS[1:10], sep = " ")
colnames(Time_Bins) <- c("max_age","min_age")

##Create and then fill a list of species present in each geologic period
Stage_Bin_Taxa <- list()

for (i in 1:nrow(Geologic_time)){
  Current_oldest <- max(Geologic_time[i,])
  Current_youngest <- min(Geologic_time[i,])

  current_period_taxa <- c()
  for (j in 1:nrow(Species_Duration_Table)){
    Current_sp <- rownames(Species_Duration_Table)[[j]]

    #create a vector which includes the sequence of values within the species duration range
    Sp_range_seq <- seq(max(Species_Duration_Table[j,]),min(Species_Duration_Table[j,]))

    #Test whether any of the numbers in the species duration rante fall within the period range
    Range_comparison <- Sp_range_seq >= Current_youngest & Sp_range_seq <= Current_oldest

    if (any(Range_comparison) == TRUE){
      current_period_taxa <- append(current_period_taxa,Current_sp)

    } else next
  }
  Stage_Bin_Taxa[[rownames(Geologic_time)[[i]]]] <- current_period_taxa
}
##

##Create and then fill a list of species present in each timebin
Time_Bin_Taxa <- list()

for (i in 1:nrow(Time_Bins)){
  Current_oldest <- max(Time_Bins[i,])
  Current_youngest <- min(Time_Bins[i,])

  current_period_taxa <- c()
  for (j in 1:nrow(Species_Duration_Table)){
    Current_sp <- rownames(Species_Duration_Table)[[j]]

    #create a vector which includes the sequence of values within the species duration range
    Sp_range_seq <- seq(max(Species_Duration_Table[j,]),min(Species_Duration_Table[j,]))

    #Test whether any of the numbers in the species duration rante fall within the period range
    Range_comparison <- Sp_range_seq >= Current_youngest & Sp_range_seq <= Current_oldest

    if (any(Range_comparison) == TRUE){
      current_period_taxa <- append(current_period_taxa,Current_sp)

    } else next
  }
  Time_Bin_Taxa[[rownames(Time_Bins)[[i]]]] <- current_period_taxa
}
##

##Make figure
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


##Add these as factors into the classifier object
classifier$GeologicPeriods <- as.character(classifier$Species)
# classifier$TimeBins <- as.character(classifier$Species)

for (i in 1:length(levels(classifier$Species))){
  current_sp <- levels(classifier$Species)[i]
  current_specimens <- grep(current_sp, classifier$Species)

  # current_sp_time_bins <- melt(Time_Bin_Taxa)[grep(gsub(pattern = " ", replacement = "_", current_sp),melt(Time_Bin_Taxa)$value),2]
  # classifier$TimeBins[current_specimens] <- paste(unique(current_sp_time_bins), sep = "", collapse = "|")

  current_sp_stage_bins <- melt(Stage_Bin_Taxa)[grep(current_sp,melt(Stage_Bin_Taxa)$value),2]
  classifier$GeologicPeriods[current_specimens] <- as.character(paste(unique(current_sp_stage_bins), sep = "", collapse = "|"))
}

classifier$GeologicPeriods <- as.character(classifier$GeologicPeriods)
levels(classifier$GeologicPeriods) <- rownames(Geologic_time)
# levels(classifier$TimeBins) <- rownames(Time_Bins)
##

###Perform subsetting and analysis of Adult/Subadults across geologic periods
Adult_classifier <- data.frame(classifier[grep("Adult|Subadult",classifier$Age),])
levels(Adult_classifier$GeologicPeriods) <- rownames(Geologic_time)


GeologicPeriods_Groups <- SubsettingGMM(Adult_classifier,
                                        Adult_CrocDorsal.gpa,
                                        Adult_CrocDorsal.pca,
                                        Adult_CrocDorsal.PlottingValues,
                                        "GeologicPeriods")

PmaxGMM(GeologicPeriods_Groups, Adult_CrocDorsal.pca, PCs = c(1,2), PrintPoints = TRUE)

GeologicPeriods_Total_Morphospace.MajorAxis_1 <- Major_Axis_Regression_ConfInt(GeologicPeriods_Groups, PCs=c(1,2), MA_number = 1, method="bootstrap", CI_values=c(0.975,0.025))
GeologicPeriods_Total_Morphospace.MajorAxis_2 <- Major_Axis_Regression_ConfInt(GeologicPeriods_Groups, PCs=c(1,2), MA_number = 2, method="bootstrap", CI_values=c(0.975,0.025))

GeologicPeriods_Groups$legendcolor <- Adult_CrocDorsal.PlottingValues$color[match(GeologicPeriods_Total_Morphospace.MajorAxis_1$Groups,Adult_classifier$GeologicPeriods)]
names(GeologicPeriods_Groups$legendcolor) <- GeologicPeriods_Total_Morphospace.MajorAxis_1$Groups

CI_Poly_Plot(GeologicPeriods_Groups,
             Combined_CrocDorsal.pca,
             GeologicPeriods_Total_Morphospace.MajorAxis_1,
             GeologicPeriods_Groups$legendcolor,
             OnePlot=FALSE)

MA_CI_calculation(GeologicPeriods_Total_Morphospace.MajorAxis_1, GeologicPeriods_Total_Morphospace.MajorAxis_1)
MA_P_calculation(GeologicPeriods_Total_Morphospace.MajorAxis_1, GeologicPeriods_Total_Morphospace.MajorAxis_1)

MA_CI_calculation(GeologicPeriods_Total_Morphospace.MajorAxis_1, Extant_Adult_Total_Morphospace.MajorAxis_1)
MA_P_calculation(GeologicPeriods_Total_Morphospace.MajorAxis_1, Extant_Adult_Total_Morphospace.MajorAxis_1)


###

TimeBins_Groups <- SubsettingGMM(Adult_classifier,
                                 Adult_CrocDorsal.gpa,
                                 Adult_CrocDorsal.pca,
                                 Adult_CrocDorsal.PlottingValues,
                                 "TimeBins")

PmaxGMM(TimeBins_Groups, Adult_CrocDorsal.pca, PCs = c(1,2), PrintPoints = TRUE)

TimeBins_Total_Morphospace.MajorAxis_1 <- Major_Axis_Regression_ConfInt(TimeBins_Groups, PCs=c(1,2), MA_number = 1, method="bootstrap", CI_values=c(0.975,0.025))
TimeBins_Total_Morphospace.MajorAxis_2 <- Major_Axis_Regression_ConfInt(TimeBins_Groups, PCs=c(1,2), MA_number = 2, method="bootstrap", CI_values=c(0.975,0.025))

CI_Poly_Plot(TimeBins_Groups, Adult_CrocDorsal.pca, TimeBins_Total_Morphospace.MajorAxis_1, Adult_CrocDorsal.PlottingValues$legendcolor, OnePlot=FALSE)
CI_Poly_Plot(TimeBins_Groups, Adult_CrocDorsal.pca, TimeBins_Total_Morphospace.MajorAxis_, Adult_CrocDorsal.PlottingValues$legendcolor, OnePlot=FALSE)
###

Fossil_classifier <- data.frame(classifier[grep("Fossil",classifier$Extant),])
levels(Fossil_classifier$GeologicPeriods) <- rownames(Geologic_time)
levels(Fossil_classifier$TimeBins) <- rownames(Time_Bins)


Fossil_GeologicPeriods_Groups <- SubsettingGMM(Fossil_classifier,
                                               Fossil_CrocDorsal.gpa,
                                               Fossil_CrocDorsal.pca,
                                               Fossil_CrocDorsal.PlottingValues,
                                               "GeologicPeriods")

as.matrix(unlist(lapply(Fossil_GeologicPeriods_Groups$CSize, "length")))

PmaxGMM(Fossil_GeologicPeriods_Groups, Combined_CrocDorsal.pca, PCs = c(1,2), PrintPoints = TRUE)

Fossil_GeologicPeriods_Total_Morphospace.MajorAxis_1 <- Major_Axis_Regression_ConfInt(Fossil_GeologicPeriods_Groups, PCs=c(1,2), MA_number = 1, method="bootstrap", CI_values=c(0.9993055555,0.0006944445))
Fossil_GeologicPeriods_Total_Morphospace.MajorAxis_1_relaxed <- Major_Axis_Regression_ConfInt(Fossil_GeologicPeriods_Groups, PCs=c(1,2), MA_number = 1, method="bootstrap", CI_values=c(0.975,0.025))
Fossil_GeologicPeriods_Total_Morphospace.MajorAxis_PC3_relaxed <- Major_Axis_Regression_ConfInt(Fossil_GeologicPeriods_Groups, PCs=c(1,3), MA_number = 1, method="bootstrap", CI_values=c(0.975,0.025))
Fossil_GeologicPeriods_Total_Morphospace.MajorAxis_PC4_relaxed <- Major_Axis_Regression_ConfInt(Fossil_GeologicPeriods_Groups, PCs=c(2,4), MA_number = 1, method="bootstrap", CI_values=c(0.975,0.025))


Fossil_GeologicPeriods_Groups$legendcolor <- Fossil_CrocDorsal.PlottingValues$color[match(Fossil_GeologicPeriods_Total_Morphospace.MajorAxis_1$Groups,Fossil_classifier$GeologicPeriods)]
names(Fossil_GeologicPeriods_Groups$legendcolor) <- Fossil_GeologicPeriods_Total_Morphospace.MajorAxis_1$Groups

pdf("Fossil_GeologicPeriods_Groups_MA1_95CI.pdf",11,8.5,useDingbats = FALSE)
CI_Poly_Plot(Fossil_GeologicPeriods_Groups,
             Combined_CrocDorsal.pca,
             Fossil_GeologicPeriods_Total_Morphospace.MajorAxis_1_relaxed,
             Fossil_GeologicPeriods_Groups$legendcolor,
             OnePlot=FALSE)
dev.off()

CI_Poly_Plot(Fossil_GeologicPeriods_Groups,
             Combined_CrocDorsal.pca,
             Fossil_GeologicPeriods_Total_Morphospace.MajorAxis_PC3_relaxed,
             Fossil_GeologicPeriods_Groups$legendcolor,
             OnePlot=FALSE)

CI_Poly_Plot(Fossil_GeologicPeriods_Groups,
             Combined_CrocDorsal.pca,
             Fossil_GeologicPeriods_Total_Morphospace.MajorAxis_PC4_relaxed,
             Fossil_GeologicPeriods_Groups$legendcolor,
             OnePlot=FALSE)

MA_CI_calculation(Fossil_GeologicPeriods_Total_Morphospace.MajorAxis_1_relaxed, Fossil_GeologicPeriods_Total_Morphospace.MajorAxis_1_relaxed)
MA_P_calculation(Fossil_GeologicPeriods_Total_Morphospace.MajorAxis_1_relaxed, Fossil_GeologicPeriods_Total_Morphospace.MajorAxis_1_relaxed)

MA_CI_calculation(Fossil_GeologicPeriods_Total_Morphospace.MajorAxis_1, Crocodylia_Ecomorphs_G_T_Unique_Total_Morphospace.MajorAxis_1)
MA_CI_calculation(Fossil_GeologicPeriods_Total_Morphospace.MajorAxis_1_relaxed, Extant_Adult_Total_Morphospace.MajorAxis_1)
MA_P_calculation(Fossil_GeologicPeriods_Total_Morphospace.MajorAxis_1_relaxed, Extant_Adult_Total_Morphospace.MajorAxis_1)


###

Fossil_TimeBins_Groups <- SubsettingGMM(Fossil_classifier,
                                        Fossil_CrocDorsal.gpa,
                                        Fossil_CrocDorsal.pca,
                                        Fossil_CrocDorsal.PlottingValues,
                                        "TimeBins")

as.matrix(unlist(lapply(Fossil_TimeBins_Groups$CSize, "length")))

PmaxGMM(Fossil_TimeBins_Groups, Fossil_CrocDorsal.pca, PCs = c(1,2), PrintPoints = TRUE)

Fossil_TimeBins_Total_Morphospace.MajorAxis_1 <- Major_Axis_Regression_ConfInt(Fossil_TimeBins_Groups, PCs=c(1,2), MA_number = 1, method="bootstrap", CI_values=c(0.9997395834,0.00026041665))
Fossil_TimeBins_Total_Morphospace.MajorAxis_2 <- Major_Axis_Regression_ConfInt(Fossil_TimeBins_Groups, PCs=c(1,2), MA_number = 2, method="bootstrap", CI_values=c(0.975,0.025))

CI_Poly_Plot(Fossil_TimeBins_Groups, Fossil_CrocDorsal.pca, Fossil_TimeBins_Total_Morphospace.MajorAxis_1, Fossil_CrocDorsal.PlottingValues$legendcolor, OnePlot=FALSE)
CI_Poly_Plot(Fossil_TimeBins_Groups, Fossil_CrocDorsal.pca, Fossil_TimeBins_Total_Morphospace.MajorAxis_2, Fossil_CrocDorsal.pca$legendcolor, OnePlot=FALSE)

MA_CI_calculation(Fossil_TimeBins_Total_Morphospace.MajorAxis_1, Fossil_TimeBins_Total_Morphospace.MajorAxis_1)

unlist(Fossil_TimeBins_Total_Morphospace.MajorAxis_1$Slope_list)[grep("MA1 .", paste(names(unlist(Fossil_TimeBins_Total_Morphospace.MajorAxis_1$Slope_list)),"."))]
unlist(Fossil_TimeBins_Total_Morphospace.MajorAxis_1$Intercept_list)[grep("MA1 .", paste(names(unlist(Fossil_TimeBins_Total_Morphospace.MajorAxis_1$Slope_list)),"."))]

Fossil_TimeBins_Total_Morphospace.MajorAxis_1$slope_CI
Fossil_TimeBins_Total_Morphospace.MajorAxis_1$intercept_CI


MA_CI_calculation(Fossil_TimeBins_Total_Morphospace.MajorAxis_1, Crocodylia_Ecomorphs_G_T_Unique_Total_Morphospace.MajorAxis_1)
MA_CI_calculation(Fossil_TimeBins_Total_Morphospace.MajorAxis_1, Extant_Adult_Total_Morphospace.MajorAxis_1)
MA_P_calculation(Fossil_TimeBins_Total_Morphospace.MajorAxis_1, Extant_Adult_Total_Morphospace.MajorAxis_1)
###

Fossil_TimeBins_Total_Morphospace.MajorAxis_1$Slope_list$`Timebin B`
Fossil_TimeBins_Total_Morphospace.MajorAxis_1$slope_CI$`Timebin B`

Extant_Adult_Total_Morphospace.MajorAxis_1$slope_CI$Crocodylia
Extant_Adult_Total_Morphospace.MajorAxis_1$Slope_list$Crocodylia




