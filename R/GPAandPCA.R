### load in data and set up values for plotting ###

## load datas
load("data/Combined_CrocDorsal.Rdata") # raw landmark data
load("data/skull_wireframe.Rdata") # landmark links for wireframe
load("data/classifier.Rdata") # speciman covariate data
# load("Chapter2Code/data/classifier.Rdata") # phylogeny
# load("Chapter2Code/data/classifier.Rdata") # species occurance datas

## perform GPA and PCA for combined dataset ##
Combined_CrocDorsal.gpa <- gpagen(Combined_CrocDorsal)
Combined_CrocDorsal.pca <- plotTangentSpace(Combined_CrocDorsal.gpa$coords, label = TRUE)
# Combined_CrocDorsal.pca <- gm.prcomp(Combined_CrocDorsal.gpa$coords)
# plot(Combined_CrocDorsal.pca)
par_default <- par() #save plot parameters for later use
##

## create and edit plotting variables for combined dataset ##
Combined_CrocDorsal.PlottingValues <- PlottingValues(classifier,
                                                 ColorGroup = "Godoy_Lifestyle",
                                                 ShapeGroup = "No_Loricata")

  # set colors for lifestyles
  Combined_CrocDorsal.PlottingValues$color[grep("Semi-aquatic", classifier$Godoy_Lifestyle)] <- "#1B9E77"
  Combined_CrocDorsal.PlottingValues$color[grep("Terrestrial", classifier$Godoy_Lifestyle)] <- "#A6761D"
  Combined_CrocDorsal.PlottingValues$color[grep("Aquatic", classifier$Godoy_Lifestyle)] <- "#386CB0"
  #
  # set colors for extant ecomorphs
  Combined_CrocDorsal.PlottingValues$color[grep("Blunt", classifier$Shape)] <- "#F4DAA2"
  Combined_CrocDorsal.PlottingValues$color[grep("Slender", classifier$Shape)] <- "#b3cde3"
  Combined_CrocDorsal.PlottingValues$color[grep("Moderate", classifier$Shape)] <- "#ccebc5"
  Combined_CrocDorsal.PlottingValues$color[grep("Tomistoma", classifier$Genus)] <- "#f4cae4"
  Combined_CrocDorsal.PlottingValues$color[grep("Gavialis", classifier$Genus)] <- "#f8b3ae"
  #
  # set shapes for extant groups
  Combined_CrocDorsal.PlottingValues$shape[grep("Tomistoma", classifier$Genus)] <- 24
    toMatch <- c("Crocodylus","Mecistops","Osteolaemus")
  Combined_CrocDorsal.PlottingValues$shape[grep(paste(toMatch,collapse="|"), classifier$Genus)] <- 22
    toMatch <- c("Alligator","Melanosuchus","Caiman","Paleosuchus")
  Combined_CrocDorsal.PlottingValues$shape[grep(paste(toMatch,collapse="|"), classifier$Genus)] <- 21
  Combined_CrocDorsal.PlottingValues$shape[grep("Gavialis", classifier$Genus)] <- 23
  Combined_CrocDorsal.PlottingValues$shape[grep("Fossil", classifier$Extant)] <- 25
  #
  # set shapes for extinct groups
  Combined_CrocDorsal.PlottingValues$shape[grep("Notosuchia", classifier$Clade2)] <- 8
  Combined_CrocDorsal.PlottingValues$shape[grep("Metriorhynchidae", classifier$Clade3)] <- 3
  Combined_CrocDorsal.PlottingValues$shape[grep("Teleosauridae", classifier$Clade3)] <- 4
  Combined_CrocDorsal.PlottingValues$shape[grep("Dyrosauridae", classifier$Clade3)] <- 7
  Combined_CrocDorsal.PlottingValues$shape[grep("Pholidosauridae", classifier$Clade3)] <- 12
  #
  # # set color for grades
  # Combined_CrocDorsal.PlottingValues$color[grep("Non-Pseudo", classifier$Clade1)] <- "#ffffb3"
  # Combined_CrocDorsal.PlottingValues$color[grep("Basal Pseudosuchia", classifier$Clade1)] <- "#fb8072"
  # Combined_CrocDorsal.PlottingValues$color[grep("Basal Loricata", classifier$Clade1)] <- "#8dd3c7"
  # Combined_CrocDorsal.PlottingValues$color[grep("Basal Crocodylomorpha", classifier$Clade1)] <- "#916da8"
  # Combined_CrocDorsal.PlottingValues$color[grep("Basal Neosuchia", classifier$Clade1)] <- "#fdb462"
  # Combined_CrocDorsal.PlottingValues$color[grep("Basal Eusuchia", classifier$Clade1)] <- "#58a7d2"
  # #

  # # fix legends
  # Combined_CrocDorsal.PlottingValues$legendcolor <- Combined_CrocDorsal.PlottingValues$color[match(Clades1.MajorAxis$Groups,classifier$Clade1)]
  # names(Combined_CrocDorsal.PlottingValues$legendcolor) <- Clades1.MajorAxis$Groups
  # Combined_CrocDorsal.PlottingValues$legendshape <- Combined_CrocDorsal.PlottingValues$shape[names(Combined_CrocDorsal.PlottingValues$legendshape)]
  # #
##

## make separate extant, 'fossil', adult, and embryo versions of the data
  # extant subset
  Extant_CrocDorsal.gpa <- Combined_CrocDorsal.gpa
  Extant_CrocDorsal.gpa$coords <- Combined_CrocDorsal.gpa$coords[,,grep("Extant",classifier$Extant)]
  Extant_CrocDorsal.gpa$Csize <- Combined_CrocDorsal.gpa$Csize[grep("Extant",classifier$Extant)]

  Extant_CrocDorsal.pca <- Combined_CrocDorsal.pca
  Extant_CrocDorsal.pca$pc.scores <- Combined_CrocDorsal.pca$pc.scores[grep("Extant",classifier$Extant),]

  Extant_CrocDorsal.PlottingValues <- lapply(Combined_CrocDorsal.PlottingValues, function(x) x[grep("Extant",classifier$Extant)])
  Extant_CrocDorsal.PlottingValues$legendcolor <- Combined_CrocDorsal.PlottingValues$legendcolor
  Extant_CrocDorsal.PlottingValues$legendshape <- Combined_CrocDorsal.PlottingValues$legendshape

  Extant_classifier <- classifier[grep("Extant",classifier$Extant),]
  Extant_classifier$Shape <- factor(Extant_classifier$Shape)
  Extant_classifier$Shape_G.T <- factor(Extant_classifier$Shape_G.T)
  Extant_classifier$Age <- factor(Extant_classifier$Age)
  #
  # fossil subset
  Fossil_CrocDorsal.gpa <- Combined_CrocDorsal.gpa
  Fossil_CrocDorsal.gpa$coords <- Combined_CrocDorsal.gpa$coords[,,grep("Fossil",classifier$Extant)]
  Fossil_CrocDorsal.gpa$Csize <- Combined_CrocDorsal.gpa$Csize[grep("Fossil",classifier$Extant)]

  Fossil_CrocDorsal.pca <- Combined_CrocDorsal.pca
  Fossil_CrocDorsal.pca$pc.scores <- Combined_CrocDorsal.pca$pc.scores[grep("Fossil",classifier$Extant),]

  Fossil_CrocDorsal.PlottingValues <- lapply(Combined_CrocDorsal.PlottingValues, function(x) x[grep("Fossil",classifier$Extant)])
  Fossil_CrocDorsal.PlottingValues$legendcolor <- Combined_CrocDorsal.PlottingValues$legendcolor
  Fossil_CrocDorsal.PlottingValues$legendshape <- Combined_CrocDorsal.PlottingValues$legendshape

  Fossil_classifier <- classifier[grep("Fossil",classifier$Extant),]
  #
  # adult & subadult subset
  Adult_CrocDorsal.gpa <- Combined_CrocDorsal.gpa
  Adult_CrocDorsal.gpa$coords <- Combined_CrocDorsal.gpa$coords[,,grep("Adult|Subadult",classifier$Age)]
  Adult_CrocDorsal.gpa$Csize <- Combined_CrocDorsal.gpa$Csize[grep("Adult|Subadult",classifier$Age)]

  Adult_CrocDorsal.pca <- Combined_CrocDorsal.pca
  Adult_CrocDorsal.pca$pc.scores <- Combined_CrocDorsal.pca$pc.scores[grep("Adult|Subadult",classifier$Age),]

  Adult_CrocDorsal.PlottingValues <- lapply(Combined_CrocDorsal.PlottingValues, function(x) x[grep("Adult|Subadult",classifier$Age)])
  Adult_CrocDorsal.PlottingValues$legendcolor <- Combined_CrocDorsal.PlottingValues$legendcolor
  Adult_CrocDorsal.PlottingValues$legendshape <- Combined_CrocDorsal.PlottingValues$legendshape

  Adult_CrocDorsal.classifier <- classifier[grep("Adult|Subadult",classifier$Age),]
  #
##

