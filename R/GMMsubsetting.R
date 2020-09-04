### make GMM subset lists ###

# make GMM subset list of data and variables for extant ecomorphstages of maturity
Extant_classifier$Age_Merged <-  as.vector(Extant_classifier$Age)
Extant_classifier$Age_Merged <- factor(gsub("Adult|Subadult", "Adult_Subadult", Extant_classifier$Age_Merged))
Extant_classifier$Age_Merged <- factor(gsub("Juvenile|Hatchling", "Juvenile_Hatchling", Extant_classifier$Age_Merged))
Extant_classifier$Age_Merged <- factor(gsub("Late-Stage Embryo|Mid-Stage Embryo", "Embryo", Extant_classifier$Age_Merged))

Extant_Maturity <- SubsettingGMM(Extant_classifier,
                                  Extant_CrocDorsal.gpa,
                                  Extant_CrocDorsal.pca,
                                  Extant_CrocDorsal.PlottingValues,
                                  "Age_Merged")
#

# make a GMM subset list of data and variables for extant ecomorphs, with Gavialis + Tomistoma as seperate groups
Extant_classifier$G_T_Unique <-  as.vector(Extant_classifier$Shape)
Extant_classifier[grep("Tomistoma", Extant_classifier$Genus),"G_T_Unique"] <- paste(Extant_classifier$Genus[grep("Tomistoma", Extant_classifier$Genus)])
Extant_classifier$G_T_Unique <- factor(Extant_classifier$G_T_Unique)

Crocodylia_Ecomorphs <- SubsettingGMM(Extant_classifier,
                                      Extant_CrocDorsal.gpa,
                                      Extant_CrocDorsal.pca,
                                      Extant_CrocDorsal.PlottingValues,
                                      "G_T_Unique")
#

# make a GMM subset list of data and variables for extant non-extreme ecomorph, Gavialis, + Tomistoma groups
Extant_classifier$Non_G_T_Ontogeny <- factor(gsub("Blunt|Moderate|Slender", "Non-extreme Ecomorph", Extant_classifier$G_T_Unique))

Crocodylia_Ecomorphs_Non_Extreme <- SubsettingGMM(Extant_classifier,
                                                  Extant_CrocDorsal.gpa,
                                                  Extant_CrocDorsal.pca,
                                                  Extant_CrocDorsal.PlottingValues,
                                                  "Non_G_T_Ontogeny")
#

# make a GMM subset list of data and variables for pseudosuchian grades
Fossil_Morphospace_Clades1 <- SubsettingGMM(Fossil_classifier,
                                            Fossil_CrocDorsal.gpa,
                                            Fossil_CrocDorsal.pca,
                                            Fossil_CrocDorsal.PlottingValues,
                                            "No_Loricata")
#

# make a GMM subset list of data and variables for major pseudosuchian subclades
Fossil_Morphospace_Clades2 <- SubsettingGMM(Fossil_classifier,
                                            Fossil_CrocDorsal.gpa,
                                            Fossil_CrocDorsal.pca,
                                            Fossil_CrocDorsal.PlottingValues,
                                            "Clade2")
#

# make a GMM subset list of data and variables for additional pseudosuchian subclades
Fossil_Morphospace_Clades3 <- SubsettingGMM(Fossil_classifier,
                                            Fossil_CrocDorsal.gpa,
                                            Fossil_CrocDorsal.pca,
                                            Fossil_CrocDorsal.PlottingValues,
                                            "Clade3")
#

# make a GMM subset list of data and variables for pseudosuchian grades
Combined_Morphospace_Clades1 <- SubsettingGMM(classifier,
                                              Combined_CrocDorsal.gpa,
                                              Combined_CrocDorsal.pca,
                                              Combined_CrocDorsal.PlottingValues,
                                              "No_Loricata")
#

# make a GMM subset list of data and variables for major pseudosuchian subclades
Combined_Morphospace_Clades2 <- SubsettingGMM(classifier,
                                              Combined_CrocDorsal.gpa,
                                              Combined_CrocDorsal.pca,
                                              Combined_CrocDorsal.PlottingValues,
                                              "Clade2")
#

# make a GMM subset list of data and variables for additional pseudosuchian subclades
Combined_Morphospace_Clades3 <- SubsettingGMM(classifier,
                                              Combined_CrocDorsal.gpa,
                                              Combined_CrocDorsal.pca,
                                              Combined_CrocDorsal.PlottingValues,
                                              "Clade3")
#
