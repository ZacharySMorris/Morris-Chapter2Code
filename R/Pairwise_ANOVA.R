#### Pairwise lm.rrpp method

#### make covariates for all analyses ####
classifier$Age_Merged <-  as.vector(classifier$Age)
classifier$Age_Merged <- factor(gsub("Adult|Subadult", "Adult_Subadult", classifier$Age_Merged))
classifier$Age_Merged <- factor(gsub("Juvenile|Hatchling", "Juvenile_Hatchling", classifier$Age_Merged))
classifier$Age_Merged <- factor(gsub("Late-Stage Embryo|Mid-Stage Embryo", "Embryo", classifier$Age_Merged))

Age_v_Clade1 <- as.character(classifier$Clade1)
Age_v_Clade1[grep("Extant",classifier$Extant)] <- as.character(classifier$Age_Merged[grep("Extant",classifier$Extant)])
Age_v_Clade1 <- as.factor(Age_v_Clade1)

Age_v_Clade2 <- as.character(classifier$Clade2)
Age_v_Clade2[grep("Extant",classifier$Extant)] <- as.character(classifier$Age_Merged[grep("Extant",classifier$Extant)])
Age_v_Clade2 <- as.factor(Age_v_Clade2)

Age_v_Clade3 <- as.character(classifier$Clade3)
Age_v_Clade3[grep("Extant",classifier$Extant)] <- as.character(classifier$Age_Merged[grep("Extant",classifier$Extant)])
Age_v_Clade3 <- as.factor(Age_v_Clade3)

classifier$Clade1
classifier$Clade2
classifier$Clade3
##

#### make necessary datasets ####
Adult_Combined_CrocDorsal.pca <- Combined_CrocDorsal.pca$pc.scores[grep("Adult|Subadult",classifier$Age_Merged),]
dim(Adult_Combined_CrocDorsal.pca)

Adult_v_Clade1 <- Age_v_Clade1[grep("Adult|Subadult",classifier$Age_Merged)]
Adult_v_Clade2 <- Age_v_Clade2[grep("Adult|Subadult",classifier$Age_Merged)]
Adult_v_Clade3 <- Age_v_Clade3[grep("Adult|Subadult",classifier$Age_Merged)]
##

#### Fit the models ####

Adult_PC1_to_4.lm <- lm.rrpp(Adult_Combined_CrocDorsal.pca[,2:4] ~ Adult_Combined_CrocDorsal.pca[,1])

GPA_PC1.lm <- lm.rrpp(Combined_CrocDorsal.gpa$coords ~ Combined_CrocDorsal.pca$pc.scores[,1])

length(Age_v_Clade1)

##


#### Perform the pairwise comparisons ####
Adult_v_Clade1.pairwise <- pairwise(Adult_PC1_to_4.lm, groups = as.factor(Adult_v_Clade1), covariate = Adult_Combined_CrocDorsal.pca[,1])
summary(Adult_v_Clade1.pairwise,test.type ="VC")

Adult_v_Clade2.pairwise <- pairwise(Adult_PC1_to_4.lm, groups = as.factor(Adult_v_Clade2), covariate = Adult_Combined_CrocDorsal.pca[,1])
summary(Adult_v_Clade2.pairwise,test.type ="VC")

Adult_v_Clade3.pairwise <- pairwise(Adult_PC1_to_4.lm, groups = as.factor(Adult_v_Clade3), covariate = Adult_Combined_CrocDorsal.pca[,1])
summary(Adult_v_Clade3.pairwise,test.type ="VC")
##
