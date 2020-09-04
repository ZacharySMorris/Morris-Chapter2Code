### stem-crocodylian ecomorph statistical comparisons

# need these three subsetGMM lists
Combined_Morphospace_Clades1
Combined_Morphospace_Clades2
Combined_Morphospace_Clades3
#

##Perform Major Axis calculation and resampling for each subset
Clades1.MajorAxis <- resample.major.axis(Fossil_Morphospace_Clades1, method = "bootstrap", iter = 9999)

Clades1.MajorAxis.Slopes <- major.axis.slopes(Clades1.MajorAxis)

pdf("Clades1.MajorAxis_PC1_PC2.pdf",11,8.5,useDingbats = FALSE)
plot.major.axis(Clades1.MajorAxis.Slopes,Fossil_Morphospace_Clades1$PCvalues,Combined_CrocDorsal.pca,PCs = c(1,2))
dev.off()


Clades1.MA.comp <- major.axis.comparison(Clades1.MajorAxis$groups,
                                         Clades1.MajorAxis$Transformed.MA,
                                         Clades1.MajorAxis$resampled_transformed.MA,
                                         MA_number = 1,  PCs = c(1:4), PC_comp = 1)

write.csv(Clades1.MA.comp$Results,"Fossil_Clades1_MA_Results.csv")

Clades2.MajorAxis <- resample.major.axis(Fossil_Morphospace_Clades2, method = "bootstrap", iter = 9999)

Clades2.MajorAxis.Slopes <- major.axis.slopes(Clades2.MajorAxis)

pdf("Clades2.MajorAxis_PC1_PC2.pdf",11,8.5,useDingbats = FALSE)
plot.major.axis(Clades2.MajorAxis.Slopes,Fossil_Morphospace_Clades2$PCvalues,Combined_CrocDorsal.pca,PCs = c(1,2))
dev.off()

Clades2.MA.comp <- major.axis.comparison(Clades2.MajorAxis$groups,
                                         Clades2.MajorAxis$Transformed.MA,
                                         Clades2.MajorAxis$resampled_transformed.MA,
                                         MA_number = 1,  PCs = c(1:4), PC_comp = 1)

write.csv(Clades2.MA.comp$Results,"Fossil_Clades2_MA_Results.csv")

Clades3.MajorAxis <- resample.major.axis(Fossil_Morphospace_Clades3, method = "bootstrap", iter = 9999)

Clades3.MajorAxis.Slopes <- major.axis.slopes(Clades3.MajorAxis)

pdf("Clades3.MajorAxis_PC1_PC2.pdf",11,8.5,useDingbats = FALSE)
plot.major.axis(Clades3.MajorAxis.Slopes,Fossil_Morphospace_Clades3$PCvalues,Combined_CrocDorsal.pca,PCs = c(1,2))
dev.off()

Clades3.MA.comp <- major.axis.comparison(Clades3.MajorAxis$groups,
                                         Clades3.MajorAxis$Transformed.MA,
                                         Clades3.MajorAxis$resampled_transformed.MA,
                                         MA_number = 1,  PCs = c(1:4), PC_comp = 1)

write.csv(Clades3.MA.comp$Results,"Fossil_Clades3_MA_Results.csv")

##

## Perform post-hoc comparisons between clades and extant ecology MA
Clades1_v_Ecology <- posthoc.major.axis.comparison(Clades1.MajorAxis,
                                                   Maturity.MajorAxis,
                                                   group_list = c("Adult_Subadult", Clades1.MajorAxis$groups),
                                                   MA_number = 1,  PCs = c(1:4), PC_comp = 1)

temp_MA <- Clades1_v_Ecology
temp_results <- array(data=NA,
                      dim = dim(temp_MA$Adj_Pvalues)+c(0,1,0),
                      dimnames = list(dimnames(temp_MA$Slopes)[[1]],c("Slopes","Pvalues"),dimnames(temp_MA$Slopes)[[3]]))

for (l in 1:dim(temp_results)[3]){
  temp_results[,,l] <- cbind(temp_MA$Slopes[,,l],temp_MA$Adj_Pvalues[,,l])
}

write.csv(temp_results,"Fossil_Clades1_v_Ecology_MA_Results.csv")

Clades2_v_Ecology <- posthoc.major.axis.comparison(Clades2.MajorAxis,
                                                   Maturity.MajorAxis,
                                                   group_list = c("Adult_Subadult", Clades2.MajorAxis$groups),
                                                   MA_number = 1,  PCs = c(1:4), PC_comp = 1)

temp_MA <- Clades2_v_Ecology
temp_results <- array(data=NA,
                      dim = dim(temp_MA$Adj_Pvalues)+c(0,1,0),
                      dimnames = list(dimnames(temp_MA$Slopes)[[1]],c("Slopes","Pvalues"),dimnames(temp_MA$Slopes)[[3]]))

for (l in 1:dim(temp_results)[3]){
  temp_results[,,l] <- cbind(temp_MA$Slopes[,,l],temp_MA$Adj_Pvalues[,,l])
}

write.csv(temp_results,"Fossil_Clades2_v_Ecology_MA_Results.csv")

Clades3_v_Ecology <- posthoc.major.axis.comparison(Clades3.MajorAxis,
                                                   Maturity.MajorAxis,
                                                   group_list = c("Adult_Subadult", "Dyrosauridae","Pholidosauridae","Hylaeochampsidae","Marine Thalattosuchia","Terrestrial Notosuchia"),
                                                   MA_number = 1,  PCs = c(1:4), PC_comp = 1)

temp_MA <- Clades3_v_Ecology
temp_results <- array(data=NA,
                      dim = dim(temp_MA$Adj_Pvalues)+c(0,1,0),
                      dimnames = list(dimnames(temp_MA$Slopes)[[1]],c("Slopes","Pvalues"),dimnames(temp_MA$Slopes)[[3]]))

for (l in 1:dim(temp_results)[3]){
  temp_results[,,l] <- cbind(temp_MA$Slopes[,,l],temp_MA$Adj_Pvalues[,,l])
}

write.csv(temp_results,"Fossil_Clades3_v_Ecology_MA_Results.csv")
##

## Perform post-hoc comparisons between clades and extant ontogeny MA
Clades1_v_Ontogenies <- posthoc.major.axis.comparison(Clades1.MajorAxis,
                                                      Ecomorph_Ontogenies.MajorAxis,
                                                      MA_number = 1,  PCs = c(1:4), PC_comp = 1)

temp_MA <- Clades1_v_Ontogenies
temp_results <- NULL
temp_results <- array(data=NA,
                      dim = dim(temp_MA$Adj_Pvalues)*c(1,2,1),
                      dimnames = list(dimnames(temp_MA$Slopes)[[1]],c(paste(dimnames(temp_MA$Slopes)[[2]],"Slopes"),paste(dimnames(temp_MA$Slopes)[[2]],"Pvalues")),dimnames(temp_MA$Slopes)[[3]]))

for (l in 1:dim(temp_results)[3]){
  temp_results[,,l] <- cbind(temp_MA$Slopes[,,l],temp_MA$Adj_Pvalues[,,l])
}

write.csv(temp_results,"Fossil_Clades1_v_Ontogenies_MA_Results.csv")
##

Clades2_v_Ontogenies <- posthoc.major.axis.comparison(Clades2.MajorAxis,
                                                      Ecomorph_Ontogenies.MajorAxis,
                                                      MA_number = 1,  PCs = c(1:4), PC_comp = 1)

temp_MA <- Clades2_v_Ontogenies
temp_results <- NULL
temp_results <- array(data=NA,
                      dim = dim(temp_MA$Adj_Pvalues)*c(1,2,1),
                      dimnames = list(dimnames(temp_MA$Slopes)[[1]],c(paste(dimnames(temp_MA$Slopes)[[2]],"Slopes"),paste(dimnames(temp_MA$Slopes)[[2]],"Pvalues")),dimnames(temp_MA$Slopes)[[3]]))

for (l in 1:dim(temp_results)[3]){
  temp_results[,,l] <- cbind(temp_MA$Slopes[,,l],temp_MA$Adj_Pvalues[,,l])
}

write.csv(temp_results,"Fossil_Clades2_v_Ontogenies_MA_Results.csv")

Clades3_v_Ontogenies <- posthoc.major.axis.comparison(Clades3.MajorAxis,
                                                      Ecomorph_Ontogenies.MajorAxis,
                                                      group_list = c(Ecomorph_Ontogenies.MajorAxis$groups, "Marine Thalattosuchia","Terrestrial Notosuchia"),
                                                      MA_number = 1,  PCs = c(1:4), PC_comp = 1)

temp_MA <- Clades3_v_Ontogenies
temp_results <- NULL
temp_results <- array(data=NA,
                      dim = dim(temp_MA$Adj_Pvalues)*c(1,2,1),
                      dimnames = list(dimnames(temp_MA$Slopes)[[1]],c(paste(dimnames(temp_MA$Slopes)[[2]],"Slopes"),paste(dimnames(temp_MA$Slopes)[[2]],"Pvalues")),dimnames(temp_MA$Slopes)[[3]]))

for (l in 1:dim(temp_results)[3]){
  temp_results[,,l] <- cbind(temp_MA$Slopes[,,l],temp_MA$Adj_Pvalues[,,l])
}

write.csv(temp_results,"Fossil_Clades3_v_Ontogenies_MA_Results.csv")

Clades1_v_nonGavialidOntogenies <- posthoc.major.axis.comparison(Clades1.MajorAxis,
                                                                 Non_Gavialid_Ontogenies.MajorAxis,
                                                                 group_list = c("Non-extreme Ecomorph", Clades1.MajorAxis$groups),
                                                                 MA_number = 1,  PCs = c(1:4), PC_comp = 1)

temp_MA <- Clades1_v_nonGavialidOntogenies
temp_results <- NULL
temp_results <- array(data=NA,
                      dim = dim(temp_MA$Adj_Pvalues)*c(1,2,1),
                      dimnames = list(dimnames(temp_MA$Slopes)[[1]],c(paste(dimnames(temp_MA$Slopes)[[2]],"Slopes"),paste(dimnames(temp_MA$Slopes)[[2]],"Pvalues")),dimnames(temp_MA$Slopes)[[3]]))

for (l in 1:dim(temp_results)[3]){
  temp_results[,,l] <- cbind(temp_MA$Slopes[,,l],temp_MA$Adj_Pvalues[,,l])
}

write.csv(temp_results,"Fossil_Clades1_v_nonGavialidOntogenies_MA_Results.csv")

Clades2_v_nonGavialidOntogenies <- posthoc.major.axis.comparison(Clades2.MajorAxis,
                                                                 Non_Gavialid_Ontogenies.MajorAxis,
                                                                 group_list = c("Non-extreme Ecomorph", Clades2.MajorAxis$groups),
                                                                 MA_number = 1,  PCs = c(1:4), PC_comp = 1)

temp_MA <- Clades2_v_nonGavialidOntogenies
temp_results <- NULL
temp_results <- array(data=NA,
                      dim = dim(temp_MA$Adj_Pvalues)*c(1,2,1),
                      dimnames = list(dimnames(temp_MA$Slopes)[[1]],c(paste(dimnames(temp_MA$Slopes)[[2]],"Slopes"),paste(dimnames(temp_MA$Slopes)[[2]],"Pvalues")),dimnames(temp_MA$Slopes)[[3]]))

for (l in 1:dim(temp_results)[3]){
  temp_results[,,l] <- cbind(temp_MA$Slopes[,,l],temp_MA$Adj_Pvalues[,,l])
}

write.csv(temp_results,"Fossil_Clades2_v_nonGavialidOntogenies_MA_Results.csv")

Clades3_v_nonGavialidOntogenies <- posthoc.major.axis.comparison(Clades3.MajorAxis,
                                                                 Non_Gavialid_Ontogenies.MajorAxis,
                                                                 group_list = c(Non_Gavialid_Ontogenies.MajorAxis$groups, "Marine Thalattosuchia", "Terrestrial Notosuchia"),
                                                                 MA_number = 1,  PCs = c(1:4), PC_comp = 1)
temp_MA <- Clades3_v_nonGavialidOntogenies
temp_results <- NULL
temp_results <- array(data=NA,
                      dim = dim(temp_MA$Adj_Pvalues)*c(1,2,1),
                      dimnames = list(dimnames(temp_MA$Slopes)[[1]],c(paste(dimnames(temp_MA$Slopes)[[2]],"Slopes"),paste(dimnames(temp_MA$Slopes)[[2]],"Pvalues")),dimnames(temp_MA$Slopes)[[3]]))

for (l in 1:dim(temp_results)[3]){
  temp_results[,,l] <- cbind(temp_MA$Slopes[,,l],temp_MA$Adj_Pvalues[,,l])
}

write.csv(temp_results,"Fossil_Clades3_v_nonGavialidOntogenies_MA_Results.csv")
##


####Function to statistically test pairwise angular comparisons between sets of group major axes####
Clades1_Angle_MA<- MA_angular_comp(Clades1.MajorAxis$groups,
                                   Clades1.MajorAxis$Transformed.MA,
                                   Clades1.MajorAxis$resampled_transformed.MA,
                                   MA_number = 1)

Clades1_v_Ecology_Angle <- posthoc_MASA_comp(Clades1.MajorAxis,
                                             Maturity.MajorAxis,
                                             group_list = c("Adult_Subadult", Clades1.MajorAxis$groups),
                                             MA_number = 1)

Clades1_v_Ecomorph_Ontogeny_Angle <- posthoc_MASA_comp(Clades1.MajorAxis,
                                             Ecomorph_Ontogenies.MajorAxis,
                                             group_list = c(Ecomorph_Ontogenies.MajorAxis$groups, Clades1.MajorAxis$groups),
                                             MA_number = 1)

Clades1_v_Non_Gavialid_Ontogeny_Angle <- posthoc_MASA_comp(Clades1.MajorAxis,
                                             Non_Gavialid_Ontogenies.MajorAxis,
                                             group_list = c("Non-extreme Ecomorph", Clades1.MajorAxis$groups),
                                             MA_number = 1)

Clades2_Angle_MA<- MA_angular_comp(Clades2.MajorAxis$groups,
                                   Clades2.MajorAxis$Transformed.MA,
                                   Clades2.MajorAxis$resampled_transformed.MA,
                                   MA_number = 1)

Clades2_v_Ecology_Angle <- posthoc_MASA_comp(Clades2.MajorAxis,
                                             Maturity.MajorAxis,
                                             group_list = c("Adult_Subadult", Clades2.MajorAxis$groups),
                                             MA_number = 1)

Clades2_v_Ecomorph_Ontogeny_Angle <- posthoc_MASA_comp(Clades2.MajorAxis,
                                                       Ecomorph_Ontogenies.MajorAxis,
                                                       group_list = c(Ecomorph_Ontogenies.MajorAxis$groups, Clades2.MajorAxis$groups),
                                                       MA_number = 1)

Clades2_v_Non_Gavialid_Ontogeny_Angle <- posthoc_MASA_comp(Clades2.MajorAxis,
                                                           Non_Gavialid_Ontogenies.MajorAxis,
                                                           group_list = c("Non-extreme Ecomorph", Clades2.MajorAxis$groups),
                                                           MA_number = 1)

Clades3_Angle_MA<- MA_angular_comp(Clades3.MajorAxis$groups,
                                   Clades3.MajorAxis$Transformed.MA,
                                   Clades3.MajorAxis$resampled_transformed.MA,
                                   MA_number = 1)

Clades3_v_Ecology_Angle <- posthoc_MASA_comp(Clades3.MajorAxis,
                                             Maturity.MajorAxis,
                                             group_list = c("Adult_Subadult", "Dyrosauridae","Pholidosauridae","Hylaeochampsidae","Marine Thalattosuchia","Terrestrial Notosuchia"),
                                             MA_number = 1)

Clades3_v_Ecomorph_Ontogeny_Angle <- posthoc_MASA_comp(Clades3.MajorAxis,
                                                       Ecomorph_Ontogenies.MajorAxis,
                                                       group_list = c(Ecomorph_Ontogenies.MajorAxis$groups, "Dyrosauridae","Pholidosauridae","Hylaeochampsidae","Marine Thalattosuchia","Terrestrial Notosuchia"),
                                                       MA_number = 1)

Clades3_v_Non_Gavialid_Ontogeny_Angle <- posthoc_MASA_comp(Clades3.MajorAxis,
                                                           Non_Gavialid_Ontogenies.MajorAxis,
                                                           group_list = c("Non-extreme Ecomorph", "Dyrosauridae","Pholidosauridae","Hylaeochampsidae","Marine Thalattosuchia","Terrestrial Notosuchia"),
                                                           MA_number = 1)
##



