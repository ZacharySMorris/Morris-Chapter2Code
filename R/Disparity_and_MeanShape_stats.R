####Disparity and Mean Shape Comparisons among groups

Extant_classifier$G_T_Unique <-  as.vector(Extant_classifier$Shape)
Extant_classifier[grep("Tomistoma", Extant_classifier$Genus),"G_T_Unique"] <- paste(Extant_classifier$Genus[grep("Tomistoma", Extant_classifier$Genus)])
Extant_classifier$G_T_Unique <- factor(Extant_classifier$G_T_Unique)

Extant_classifier$Age_Merged <-  as.vector(Extant_classifier$Age)
Extant_classifier[grep("Subadult", Extant_classifier$Age),"Age_Merged"] <- "Adult"
Extant_classifier$Age_Merged <- factor(Extant_classifier$Age_Merged)

Extant_classifier$EcomorphVSAGE <- as.vector(Extant_classifier$G_T_Unique)
Extant_classifier[grep("Adult", Extant_classifier$Age_Merged),"EcomorphVSAGE"] <- "Adult"
Extant_classifier$EcomorphVSAGE <- factor(Extant_classifier$EcomorphVSAGE)

#Meanshape comparisons among ecomorphs
Ontogeny_Ecology <- geomorph.data.frame(coords = Combined_CrocDorsal.gpa$coords[,,grep("Extant", classifier$Extant)],
                                        Ecomorph = Extant_classifier$G_T_Unique,
                                        Age = Extant_classifier$Age_Merged,
                                        EcomorphVSAge = Extant_classifier$EcomorphVSAGE)

Ontogeny_fit <- procD.lm(coords ~ 1 + Age, data = Ontogeny_Ecology)
anova(Ontogeny_fit)
summary(pairwise(Ontogeny_fit, groups = Ontogeny_Ecology$Age))

Ecology_fit <- procD.lm(coords ~ 1 + Ecomorph, data = Ontogeny_Ecology)
anova(Ecology_fit)
summary(pairwise(Ecology_fit, groups = Ontogeny_Ecology$Ecomorph))

Ontogeny_Ecology_fit <- procD.lm(coords ~ 1 + EcomorphVSAge, data = Ontogeny_Ecology)
anova(Ontogeny_Ecology_fit)
UGH2 <- interaction(Ontogeny_Ecology$Age, Ontogeny_Ecology$Ecomorph)
summary(pairwise(Ontogeny_Ecology_fit, groups = Ontogeny_Ecology$EcomorphVSAge))

#Disparity comparisons for various factors
morphol.disparity(Combined_CrocDorsal.gpa$coords ~ 1, groups = classifier$Extant)

morphol.disparity(Combined_CrocDorsal.gpa$coords ~ 1, groups = classifier$Clade1)

morphol.disparity(Combined_CrocDorsal.gpa$coords ~ 1, groups = classifier$No_Loricata)

morphol.disparity(Combined_CrocDorsal.gpa$coords ~ 1, groups = classifier$Clade2)

morphol.disparity(Combined_CrocDorsal.gpa$coords ~ 1, groups = classifier$Clade3)

##Lifestyle Differences
Lifestyle <- geomorph.data.frame(Combined_CrocDorsal.gpa, lifestyle = classifier$Godoy_Lifestyle)
Lifestyle_extinct <- geomorph.data.frame(Fossil_CrocDorsal.gpa, lifestyle =  Fossil_classifier$Godoy_Lifestyle)

Lifestyle_fit <- procD.lm(coords ~ 1 + lifestyle, data = Lifestyle)
Lifestyle_extinct_fit <- procD.lm(coords ~ 1 + lifestyle, data = Lifestyle_extinct)
anova(Lifestyle_fit)
anova(Lifestyle_extinct_fit)
#pairwise comparisons of lifestyle mean shapes
summary(pairwise(Lifestyle_fit, groups = Lifestyle$lifestyle))
summary(pairwise(Lifestyle_extinct_fit, groups = Lifestyle_extinct$lifestyle))
##Calculate Mean Shapes
mshape(Combined_CrocDorsal.gpa$coords[,,grep("Terrestrial", classifier$Godoy_Lifestyle)])
mshape(Combined_CrocDorsal.gpa$coords[,,grep("Semi-aquatic", classifier$Godoy_Lifestyle)])
mshape(Combined_CrocDorsal.gpa$coords[,,grep("Aquatic", classifier$Godoy_Lifestyle)])

##testing differences in disparity of "lifestyles"
morphol.disparity(Combined_CrocDorsal.gpa$coords ~ 1, groups = Lifestyle$lifestyle)
##testing differences in disparity of "lifestyles", excluding extant specimens from disparity measurements
morphol.disparity(Combined_CrocDorsal.gpa$coords[,,grep("Fossil", classifier$Extant)] ~ 1, groups = Lifestyle$lifestyle[grep("Fossil", classifier$Extant)])

Lifestyle_extinct <- procD.lm(coords[] ~ 1 + lifestyle, data = Lifestyle)


#####
Clades_GDF <- geomorph.data.frame(Combined_CrocDorsal.gpa,
                                  PCData = Combined_CrocDorsal.pca$pc.scores,
                                  clade1 = classifier$No_Loricata,
                                  clade2 = classifier$Clade2,
                                  clade3 = classifier$Clade3)
Clades_extinct_GDF <- geomorph.data.frame(coords = Combined_CrocDorsal.gpa$coords[,,grep("Fossil", classifier$Extant)],
                                          PCData = Combined_CrocDorsal.pca$pc.scores[grep("Fossil", classifier$Extant),],
                                          clade1 = classifier$No_Loricata[grep("Fossil", classifier$Extant)],
                                          clade2 = classifier$Clade2[grep("Fossil", classifier$Extant)],
                                          clade3 = classifier$Clade3[grep("Fossil", classifier$Extant)])

##Differences in clade1 mean shapes
Clade1_fit <- procD.lm(coords ~ 1 + clade1, data = Clades_GDF)
Clade1_extinct_fit <- procD.lm(coords ~ 1 + clade1, data = Clades_extinct_GDF)
anova(Clade1_fit)
anova(Clade1_extinct_fit)
#pairwise comparisons of lifestyle mean shapes
summary(pairwise(Clade1_fit, groups = Clades_GDF$clade1))
summary(pairwise(Clade1_extinct_fit, groups = Clades_extinct_GDF$clade1))

##Differences in clade2 mean shapes
Clade2_fit <- procD.lm(coords ~ 1 + clade2, data = Clades_GDF)
Clade2_extinct_fit <- procD.lm(coords ~ 1 + clade2, data = Clades_extinct_GDF)
anova(Clade2_fit)
anova(Clade2_extinct_fit)
#pairwise comparisons of lifestyle mean shapes
summary(pairwise(Clade2_fit, groups = Clades_GDF$clade2))
summary(pairwise(Clade2_extinct_fit, groups = Clades_extinct_GDF$clade2))

##Differences in clade3 mean shapes
Clade3_fit <- procD.lm(coords ~ 1 + clade3, data = Clades_GDF)
Clade3_extinct_fit <- procD.lm(coords ~ 1 + clade3, data = Clades_extinct_GDF)
anova(Clade3_fit)
anova(Clade3_extinct_fit)
#pairwise comparisons of lifestyle mean shapes
summary(pairwise(Clade3_fit, groups = Clades_GDF$clade3))
summary(pairwise(Clade3_extinct_fit, groups = Clades_extinct_GDF$clade3))


Clade3_lm_fit <- lm.rrpp(Clades_GDF$PCData[,2:4] ~ Clades_GDF$PCData[,1] * clade3, data = Clades_GDF)
test_lm_fit <- lm.rrpp(Clades_GDF$PCData[,2:4] ~ Clades_GDF$PCData[,1], data = Clades_GDF)
test_lm_fit <- lm(Clades_GDF$PCData[,2:4] ~ Clades_GDF$PCData[,1])

summary(test_lm_fit)$adj.r.squared

anova(Clade3_lm_fit)
summary(pairwise(Clade3_lm_fit, groups = Clades_GDF$clade3))$r.squared

