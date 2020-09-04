

setwd("~/Dropbox/Dissertation/Chapter_2")

## Import new landmark data (TPS files) ##
Fossil_CrocDorsal <- readland.tps("Fossil_Croc_Dorsal.tps", specID="imageID")
##
## Import Morris et al., 2019 extant crocodylian landmark data ##
Extant_CrocDorsal <- readland.tps("Extant_Croc_Dorsal.tps", specID="imageID")
##

## Make and save combined landmark dataset ##
Combined_CrocDorsal <- abind(Extant_CrocDorsal,Fossil_CrocDorsal)
save(Combined_CrocDorsal, file = "Chapter2Code/data/Combined_CrocDorsal.Rdata")
##

## Perform GPA and PCA for combined dataset ##
Combined_CrocDorsal.gpa <- gpagen(Combined_CrocDorsal)
Combined_CrocDorsal.pca <- plotTangentSpace(Combined_CrocDorsal.gpa$coords, label = TRUE)
# Combined_CrocDorsal.pca <- gm.prcomp(Combined_CrocDorsal.gpa$coords)
# plot(Combined_CrocDorsal.pca)
par_default <- par() #save plot parameters for later use
##

## Import and save links among landmarks to make wireframe
skull_wireframe <- read.csv("14LM_Wireframe.csv", header = FALSE)
save(skull_wireframe, file = "Chapter2Code/data/skull_wireframe.Rdata")
##

## Import and save classifier data (CSV files)#
classifier <- read.csv("Chapter2_Covariates.csv", header=TRUE)
save(classifier, file = "Chapter2Code/data/classifier.Rdata")
##

# #Perform extant only alignment and PCA#
# Extant_CrocDorsal.gpa <- gpagen(combined_CrocDorsal[,,grep("Extant",classifier$Extant)])
# Extant_CrocDorsal.pca <- plotTangentSpace(Extant_CrocDorsal.gpa$coords, label = TRUE)
#
# #Perform fossil only alignment and PCA#
# Fossil_CrocDorsal.gpa <- gpagen(combined_CrocDorsal[,,grep("Fossil",classifier$Extant)])
# Fossil_CrocDorsal.pca <- plotTangentSpace(Fossil_CrocDorsal.gpa$coords, label = TRUE)

#Adult & Subadult subsets
Adult_CrocDorsal.gpa <- Combined_CrocDorsal.gpa
Adult_CrocDorsal.gpa$coords <- Combined_CrocDorsal.gpa$coords[,,grep("Adult|Subadult",classifier$Age)]
Adult_CrocDorsal.gpa$Csize <- Combined_CrocDorsal.gpa$Csize[grep("Adult|Subadult",classifier$Age)]

Adult_CrocDorsal.pca <- Combined_CrocDorsal.pca
Adult_CrocDorsal.pca$pc.scores <- Combined_CrocDorsal.pca$pc.scores[grep("Adult|Subadult",classifier$Age),]

Adult_CrocDorsal.PlottingValues <- lapply(TotalCrocDorsal.PlottingValues, function(x) x[grep("Adult|Subadult",classifier$Age)])
Adult_CrocDorsal.PlottingValues$legendcolor <- TotalCrocDorsal.PlottingValues$legendcolor
Adult_CrocDorsal.PlottingValues$legendshape <- TotalCrocDorsal.PlottingValues$legendshape

Adult_CrocDorsal.classifier <- classifier[grep("Adult|Subadult",classifier$Age),]

###

#Fossil subsets
Fossil_CrocDorsal.gpa <- Combined_CrocDorsal.gpa
Fossil_CrocDorsal.gpa$coords <- Combined_CrocDorsal.gpa$coords[,,grep("Fossil",classifier$Extant)]
Fossil_CrocDorsal.gpa$Csize <- Combined_CrocDorsal.gpa$Csize[grep("Fossil",classifier$Extant)]

Fossil_CrocDorsal.pca <- Combined_CrocDorsal.pca
Fossil_CrocDorsal.pca$pc.scores <- Combined_CrocDorsal.pca$pc.scores[grep("Fossil",classifier$Extant),]

Fossil_CrocDorsal.PlottingValues <- lapply(TotalCrocDorsal.PlottingValues, function(x) x[grep("Fossil",classifier$Extant)])
Fossil_CrocDorsal.PlottingValues$legendcolor <- TotalCrocDorsal.PlottingValues$legendcolor
Fossil_CrocDorsal.PlottingValues$legendshape <- TotalCrocDorsal.PlottingValues$legendshape

#Extant subsets
Extant_CrocDorsal.gpa <- Combined_CrocDorsal.gpa
Extant_CrocDorsal.gpa$coords <- Combined_CrocDorsal.gpa$coords[,,grep("Extant",classifier$Extant)]
Extant_CrocDorsal.gpa$Csize <- Combined_CrocDorsal.gpa$Csize[grep("Extant",classifier$Extant)]

Extant_CrocDorsal.pca <- Combined_CrocDorsal.pca
Extant_CrocDorsal.pca$pc.scores <- Combined_CrocDorsal.pca$pc.scores[grep("Extant",classifier$Extant),]

Extant_CrocDorsal.PlottingValues <- lapply(TotalCrocDorsal.PlottingValues, function(x) x[grep("Extant",classifier$Extant)])
Extant_CrocDorsal.PlottingValues$legendcolor <- TotalCrocDorsal.PlottingValues$legendcolor
Extant_CrocDorsal.PlottingValues$legendshape <- TotalCrocDorsal.PlottingValues$legendshape

Extant_classifier <- classifier[grep("Extant",classifier$Extant),]


length(TotalCrocDorsal.PlottingValues)

###
##Create and edit plotting values for combined dataset##
TotalCrocDorsal.PlottingValues <- PlottingValues(classifier,
                                                 ColorGroup = "Godoy_Lifestyle",
                                                 ShapeGroup = "No_Loricata")

#Color Lifestyle for fossils
TotalCrocDorsal.PlottingValues$color[grep("Semi-aquatic", classifier$Godoy_Lifestyle)] <- "#1B9E77"
TotalCrocDorsal.PlottingValues$color[grep("Terrestrial", classifier$Godoy_Lifestyle)] <- "#A6761D"
TotalCrocDorsal.PlottingValues$color[grep("Aquatic", classifier$Godoy_Lifestyle)] <- "#386CB0"

#Color extant ecomorphs
TotalCrocDorsal.PlottingValues$color[grep("Blunt", classifier$Shape)] <- "#F4DAA2"
TotalCrocDorsal.PlottingValues$color[grep("Slender", classifier$Shape)] <- "#b3cde3"
TotalCrocDorsal.PlottingValues$color[grep("Moderate", classifier$Shape)] <- "#ccebc5"
TotalCrocDorsal.PlottingValues$color[grep("Tomistoma", classifier$Genus)] <- "#f4cae4"
TotalCrocDorsal.PlottingValues$color[grep("Gavialis", classifier$Genus)] <- "#f8b3ae"

#Make sure that shapes are different for extant groups
TotalCrocDorsal.PlottingValues$shape[grep("Tomistoma", classifier$Genus)] <- 24
toMatch <- c("Crocodylus","Mecistops","Osteolaemus")
TotalCrocDorsal.PlottingValues$shape[grep(paste(toMatch,collapse="|"), classifier$Genus)] <- 22
toMatch <- c("Alligator","Melanosuchus","Caiman","Paleosuchus")
TotalCrocDorsal.PlottingValues$shape[grep(paste(toMatch,collapse="|"), classifier$Genus)] <- 21
TotalCrocDorsal.PlottingValues$shape[grep("Gavialis", classifier$Genus)] <- 23
TotalCrocDorsal.PlottingValues$shape[grep("Fossil", classifier$Extant)] <- 25

#Make sure that shapes are different for extinct groups
TotalCrocDorsal.PlottingValues$shape[grep("Notosuchia", classifier$Clade2)] <- 8
TotalCrocDorsal.PlottingValues$shape[grep("Metriorhynchidae", classifier$Clade3)] <- 3
TotalCrocDorsal.PlottingValues$shape[grep("Teleosauridae", classifier$Clade3)] <- 4
TotalCrocDorsal.PlottingValues$shape[grep("Dyrosauridae", classifier$Clade3)] <- 7
TotalCrocDorsal.PlottingValues$shape[grep("Pholidosauridae", classifier$Clade3)] <- 12


#Color grades
TotalCrocDorsal.PlottingValues$color[grep("Non-Pseudo", classifier$Clade1)] <- "#ffffb3"
TotalCrocDorsal.PlottingValues$color[grep("Basal Pseudosuchia", classifier$Clade1)] <- "#fb8072"
TotalCrocDorsal.PlottingValues$color[grep("Basal Loricata", classifier$Clade1)] <- "#8dd3c7"
TotalCrocDorsal.PlottingValues$color[grep("Basal Crocodylomorpha", classifier$Clade1)] <- "#916da8"
TotalCrocDorsal.PlottingValues$color[grep("Basal Neosuchia", classifier$Clade1)] <- "#fdb462"
TotalCrocDorsal.PlottingValues$color[grep("Basal Eusuchia", classifier$Clade1)] <- "#58a7d2"

#Fix Legends
TotalCrocDorsal.PlottingValues$legendcolor <- TotalCrocDorsal.PlottingValues$color[match(Clades1.MajorAxis$Groups,classifier$Clade1)]
names(TotalCrocDorsal.PlottingValues$legendcolor) <- Clades1.MajorAxis$Groups


TotalCrocDorsal.PlottingValues$legendshape <- TotalCrocDorsal.PlottingValues$shape[names(TotalCrocDorsal.PlottingValues$legendshape)]
###

Value_names <- names(TotalCrocDorsal.PlottingValues)
ExtantCrocDorsal.PlottingValues <- list()
FossilCrocDorsal.PlottingValues <- list()

for (i in 1:length(Value_names)){
  ExtantCrocDorsal.PlottingValues[[Value_names[[i]]]] <- TotalCrocDorsal.PlottingValues[[i]][grep("Extant", classifier$Extant)]
  FossilCrocDorsal.PlottingValues[[Value_names[[i]]]] <- TotalCrocDorsal.PlottingValues[[i]][grep("Fossil", classifier$Extant)]
}

ExtantCrocDorsal.PlottingValues$legendcolor <- TotalCrocDorsal.PlottingValues$legendcolor
ExtantCrocDorsal.PlottingValues$legendshape <- TotalCrocDorsal.PlottingValues$legendshape

FossilCrocDorsal.PlottingValues$legendcolor <- TotalCrocDorsal.PlottingValues$legendcolor
FossilCrocDorsal.PlottingValues$legendshape <- TotalCrocDorsal.PlottingValues$legendshape



