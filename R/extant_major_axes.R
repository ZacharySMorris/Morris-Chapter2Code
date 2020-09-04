### crocodylian ecomorph statistical comparisons ###

# need these three subsetGMM lists
Crocodylia_Ecomorphs
Crocodylia_Ecomorphs_Non_Extreme
Extant_Maturity
#

#### Perform Major Axis calculation and resampling for each subset ####
Ecomorph_Ontogenies.MajorAxis <- resample.major.axis(Crocodylia_Ecomorphs, method = "bootstrap", iter = 9999)
Ecomorph_Ontogenies_PC1_3.MajorAxis <- resample.major.axis(Crocodylia_Ecomorphs, PCs=c(1:3), method = "bootstrap", iter = 9999)

Ecomorph_Ontogenies.MajorAxis.Slopes <- major.axis.slopes(Ecomorph_Ontogenies.MajorAxis)

Ecomorph_Ontogenies.MA.comp <- major.axis.comparison(Ecomorph_Ontogenies.MajorAxis$groups,
                                                     Ecomorph_Ontogenies.MajorAxis$Transformed.MA,
                                                     Ecomorph_Ontogenies.MajorAxis$resampled_transformed.MA,
                                                     MA_number = 1,  PCs = c(1:4), PC_comp = 1)

write.csv(Ecomorph_Ontogenies.MA.comp$Results,"Extant_Ecommorph_Ontogenies_MA_Results.csv")

Non_Gavialid_Ontogenies.MajorAxis <- resample.major.axis(Crocodylia_Ecomorphs_Non_Extreme, method = "bootstrap", iter = 9999)
Non_Gavialid_Ontogenies_PC1_3.MajorAxis <- resample.major.axis(Crocodylia_Ecomorphs_Non_Extreme, PCs=c(1:3), method = "bootstrap", iter = 9999)

Non_Gavialid_Ontogenies.MajorAxis.Slopes <- major.axis.slopes(Non_Gavialid_Ontogenies.MajorAxis)

pdf("Ontogenies.MajorAxis_PC1_PC2.pdf",11,8.5,useDingbats = FALSE)
plot.major.axis(Ecomorph_Ontogenies.MajorAxis.Slopes,Crocodylia_Ecomorphs$PCvalues,Combined_CrocDorsal.pca,PCs = c(1,2))
plot.major.axis(Non_Gavialid_Ontogenies.MajorAxis.Slopes,Crocodylia_Ecomorphs_Non_Extreme$PCvalues,Combined_CrocDorsal.pca,PCs = c(1,2))
dev.off()

Non_Gavialid_Ontogenies.MA.comp <- major.axis.comparison(Non_Gavialid_Ontogenies.MajorAxis$groups,
                                                         Non_Gavialid_Ontogenies.MajorAxis$Transformed.MA,
                                                         Non_Gavialid_Ontogenies.MajorAxis$resampled_transformed.MA,
                                                         MA_number = 1,  PCs = c(1:4), PC_comp = 1)

write.csv(Non_Gavialid_Ontogenies.MA.comp$Results,"Non_Gavialid_Ontogenies_MA_Results.csv")

Maturity.MajorAxis <- resample.major.axis(Extant_Maturity, method = "bootstrap", iter = 9999)
Maturity_PC1_3.MajorAxis <- resample.major.axis(Extant_Maturity, PCs=c(1:3), method = "bootstrap", iter = 9999)

Maturity.MajorAxis$Loadings

Maturity.MajorAxis.Slopes <- major.axis.slopes(Maturity.MajorAxis)

pdf("Maturity.MajorAxis_PC1_PC2.pdf",11,8.5,useDingbats = FALSE)
plot.major.axis(Maturity.MajorAxis.Slopes,Extant_Maturity$PCvalues,Combined_CrocDorsal.pca,PCs = c(1,2))
dev.off()

Maturity.MA.comp <- major.axis.comparison(Maturity.MajorAxis$groups,
                                          Maturity.MajorAxis$Transformed.MA,
                                          Maturity.MajorAxis$resampled_transformed.MA,
                                          MA_number = 1,  PCs = c(1:4), PC_comp = 1)

Maturity.MA.comp$Results


#### Perform post-hoc comparisons between extant ecology and ontogeny MAs ####
Ecology_v_Ontogenies <- posthoc.major.axis.comparison(Ecomorph_Ontogenies.MajorAxis,
                                                      Maturity.MajorAxis,
                                                      group_list = c("Adult_Subadult", Ecomorph_Ontogenies.MajorAxis$groups),
                                                      MA_number = 1,  PCs = c(1:4), PC_comp = 1)

Ecology_v_Ontogenies_Results <- array(data=NA,
                                      dim = dim(Ecology_v_Ontogenies$Adj_Pvalues)+c(0,1,0),
                                      dimnames = list(dimnames(Ecology_v_Ontogenies$Slopes)[[1]],c("Slopes","Pvalues"),dimnames(Ecology_v_Ontogenies$Slopes)[[3]]))

for (l in 1:dim(Ecology_v_Ontogenies$Adj_Pvalues)[3]){
  Ecology_v_Ontogenies_Results[,,l] <- cbind(Ecology_v_Ontogenies$Slopes[,,l],Ecology_v_Ontogenies$Adj_Pvalues[,,l])
  }

write.csv(Ecology_v_Ontogenies_Results, "Ecology_v_Ontogenies_MA_Results.csv")


Ecology_v_nonGavialidOntogenies <- posthoc.major.axis.comparison(Non_Gavialid_Ontogenies.MajorAxis,
                                                                 Maturity.MajorAxis,
                                                                 group_list = c("Adult_Subadult", Non_Gavialid_Ontogenies.MajorAxis$groups),
                                                                 MA_number = 1,  PCs = c(1:4), PC_comp = 1)
##
Ecology_v_nonGavialidOntogenies_Results <- array(data=NA,
                                      dim = dim(Ecology_v_nonGavialidOntogenies$Adj_Pvalues)+c(0,1,0),
                                      dimnames = list(dimnames(Ecology_v_nonGavialidOntogenies$Slopes)[[1]],c("Slopes","Pvalues"),dimnames(Ecology_v_nonGavialidOntogenies$Slopes)[[3]]))

for (l in 1:dim(Ecology_v_nonGavialidOntogenies_Results)[3]){
  Ecology_v_nonGavialidOntogenies_Results[,,l] <- cbind(Ecology_v_nonGavialidOntogenies$Slopes[,,l],Ecology_v_nonGavialidOntogenies$Adj_Pvalues[,,l])
}

write.csv(Ecology_v_nonGavialidOntogenies_Results,"Ecology_v_nonGavialidOntogenies_MA_Results.csv")

Ecomorph_v_nonGavialidOntogenies <- posthoc.major.axis.comparison(Ecomorph_Ontogenies.MajorAxis,
                                                                  Non_Gavialid_Ontogenies.MajorAxis,
                                                                  MA_number = 1,  PCs = c(1:4), PC_comp = 1)


#### ANGLE VERSIONS ####

##Function to statistically test pairwise angular comparisons between sets of group major axes##
Ecomorph_Ontogenies_Angle_MA<- MA_angular_comp(Ecomorph_Ontogenies_PC1_3.MajorAxis$groups,
                                               Ecomorph_Ontogenies_PC1_3.MajorAxis$Transformed.MA,
                                               Ecomorph_Ontogenies_PC1_3.MajorAxis$resampled_transformed.MA,
                                               MA_number = 1)

Non_Gavialid_Ontogenies_Angle_MA<- MA_angular_comp(Non_Gavialid_Ontogenies_PC1_3.MajorAxis$groups,
                                               Non_Gavialid_Ontogenies_PC1_3.MajorAxis$Transformed.MA,
                                               Non_Gavialid_Ontogenies_PC1_3.MajorAxis$resampled_transformed.MA,
                                               MA_number = 1)

Maturity_Angle_MA<- MA_angular_comp(Maturity_PC1_3.MajorAxis$groups,
                                               Maturity_PC1_3.MajorAxis$Transformed.MA,
                                               Maturity_PC1_3.MajorAxis$resampled_transformed.MA,
                                               MA_number = 1)

Ecology_v_Ontogenies_Angle <- posthoc_MASA_comp(Ecomorph_Ontogenies_PC1_3.MajorAxis,
                                                Maturity_PC1_3.MajorAxis,
                                                group_list = c("Adult_Subadult", Ecomorph_Ontogenies_PC1_3.MajorAxis$groups),
                                                MA_number = 1)


Ecology_v_nonGavialidOntogenies <- posthoc_MASA_comp(Non_Gavialid_Ontogenies_PC1_3.MajorAxis,
                                                     Maturity_PC1_3.MajorAxis,
                                                     group_list = c("Adult_Subadult", Non_Gavialid_Ontogenies_PC1_3.MajorAxis$groups),
                                                     MA_number = 1)
##


