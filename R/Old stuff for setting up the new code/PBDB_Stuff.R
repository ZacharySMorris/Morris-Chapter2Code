###PBDB stuff###

#Download all Pseudosuchia entries for closer inspection#
Pseudosuchia_total_pbdb_occurances <- pbdb_occurrences(limit= "all",
                                                       base_name = "Pseudosuchia",
                                                       vocab= "pbdb",
                                                       show=c("phylo", "ident"))

write.csv(Pseudosuchia_total_pbdb_occurances, file = "Pseudosuchia_total_pbdb_occurances.csv")

#Use just the taxa in this dataset#
Species_List <- paste(unique(classifier$Species))

Species_List_sorted <- Species_List[match(gsub("_"," ", Pseudosuchia.phy$tip.label),Species_List)]

Species_Duration_Table <- matrix(data = rep(c(15.000,0.000),length(Species_List_sorted)),
                                 byrow = TRUE,
                                 nrow = length(Species_List_sorted), ncol=2,
                                  dimnames=list(c(Species_List_sorted),c("max_age","min_age")))

PBDB_Data <- pbdb_occurrences(limit= "all",
                              taxon_name = c(Species_List_sorted),
                              # vocab= "pbdb",
                              show=c("phylo", "ident","ref"))

write.csv(PBDB_Data, file = "PBDB_Data.csv")

PBDB_Refs <- pbdb_ref_occurrences(limit= "all",
                              taxon_name = c(Species_List_sorted),
                              # vocab= "pbdb",
                              show=c("ref"))


Species_durations <- pbdb_temp_range(PBDB_Test, rank = "species", do.plot = FALSE)

rownames(Species_durations) %in% Species_List_sorted
Species_durations[is.na(match(rownames(Species_durations),Species_List_sorted)),]

rownames(Species_durations)[is.na(match(rownames(Species_durations),Species_List_sorted))] <- c("Desmatosuchus haplocerus",
                                                                                                "Metriorhynchus superciliosum",
                                                                                                "Uberabasuchus terrificus",
                                                                                                "Hyposaurus natator",
                                                                                                "Gavialosuchus eggenburgensis")

Species_Duration_Table[match(rownames(Species_durations),Species_List_sorted),c(1:2)] <- as.matrix(Species_durations[,c(1:2)])
Species_Duration_Table

###Don't actually need, as the dates were grabbed correctly, just the names were wrong
###So i've modified the code above to swap the names around for these 4 species
#swap for Desmatosuchus, since the pbdb is dumb sometimes
# Desmatosuchus_temp <- pbdb_temp_range(pbdb_occurrences(limit= "all",
#                                                        taxon_name = "Episcoposaurus haplocerus",
#                                                        vocab= "pbdb",
#                                                        show=c("phylo", "ident")),
#                                       rank="species",
#                                       do.plot = FALSE)
#
# Species_Duration_Table["Desmatosuchus haplocerus",] <- as.numeric(Desmatosuchus_temp)
#
# #swap for Uberabasuchus, since the pbdb is dumb sometimes
# Uberabasuchus_temp <- pbdb_temp_range(pbdb_occurrences(limit= "all",
#                                                        taxon_name = c("Uberabasuchus terrificus","Peirosaurus torminni"),
#                                                        vocab= "pbdb",
#                                                        show=c("phylo", "ident")),
#                                       rank="species",
#                                       do.plot = FALSE)
#
# Species_Duration_Table["Uberabasuchus terrificus",] <- as.numeric(Uberabasuchus_temp)
#
# #swap for Hyposaurus natator, since the pbdb is dumb sometimes
# Hyposaurus_temp <- pbdb_temp_range(pbdb_occurrences(limit= "all",
#                                                        taxon_name = c("Hyposaurus natator"),
#                                                        vocab= "pbdb",
#                                                        show=c("phylo", "ident")),
#                                    rank="species",
#                                    do.plot = FALSE)
#
# Species_Duration_Table["Hyposaurus natator",] <- as.numeric(Hyposaurus_temp)
#
# Metriorhynchus_temp <- pbdb_temp_range(pbdb_occurrences(limit= "all",
#                                                         taxon_name = c("Metriorhynchus superciliosum"),
#                                                         vocab= "pbdb",
#                                                         show=c("phylo", "ident")),
#                                        rank="species",
#                                        do.plot = FALSE)
#
# Species_Duration_Table["Metriorhynchus superciliosum",] <- as.numeric(Metriorhynchus_temp)
###
###

#swap for Diplocynodon hantoniensis, since the pbdb is dumb sometimes
Diplocynodon_temp <- pbdb_temp_range(pbdb_occurrences(limit= "all",
                                                    taxon_name = c("Diplocynodon hantoniensis","Crocodilus hastingsiae"),
                                                    vocab= "pbdb",
                                                    show=c("phylo", "ident")),
                                     rank="species",
                                     do.plot = FALSE)

Species_Duration_Table["Diplocynodon hantoniensis",] <- as.numeric(Diplocynodon_temp)



Species_Duration_Table["Prodiplocynodon Utah specimen",] <- Species_Duration_Table["Prodiplocynodon langi",]

Species_Duration_Table["Allognathosuchus sp.",] <- Species_Duration_Table["Allognathosuchus wartheni",]

Species_Duration_Table["Protosuchus sp.",] <- Species_Duration_Table["Protosuchus richardsoni",]

Species_Duration_Table["Borealosuchus tullock specimen",] <- Species_Duration_Table["Borealosuchus formidabilis",]

Species_Duration_Table["Thoracosaurus",] <- Species_Duration_Table["Thoracosaurus macrorhynchus",]

Species_Duration_Table["Mexicosuchus downiearum",] <- c(83.6,72.1)

Amphicotylus_temp <- pbdb_temp_range(pbdb_occurrences(limit= "all",
                                                      taxon_name = c("Goniopholis felix"),
                                                      vocab= "pbdb",
                                                      show=c("phylo", "ident")),
                                     rank="species",
                                     do.plot = FALSE)

Species_Duration_Table["Amphicotylus felix",] <- as.numeric(Amphicotylus_temp)

Extant_Species_List <- as.character(unique(classifier$Species[grep("Extant", classifier$Extant)]))

Species_Duration_Table[Extant_Species_List,2] <- as.numeric(0.000)
Species_Duration_Table[Species_Duration_Table[,1] == 15.000,1] <- as.numeric(5.000)

rownames(Species_Duration_Table) <- gsub(" ", "_", rownames(Species_Duration_Table))

Pseudosuchia_dichotomous.phy <- timePaleoPhy(multi2di(collapse.singles(Pseudosuchia.phy), random = FALSE), Species_Duration_Table)
is.binary(Pseudosuchia_dichotomous.phy)
is.binary(Pseudosuchia.phy)




Species_Duration_Table <= timebin_edges[1] & Species_Duration_Table >= timebin_edges[2]

grep("TRUE", Temp_duration_table[,1:2])

grep("TRUE", Temp_duration_table[,2])

class(Temp_duration_table)

pos <- c(1:dim(na.omit(Species_Duration_Table))[1] - 0.9)
t_range <- cbind(na.omit(Species_Duration_Table), pos)
par(mar = c(4, 0, 1, 15))


plot(c(min(t_range[,"max_age"]), max(t_range[,"max_age"])), c(0, dim(t_range)[1]),
     type = "n", axes = FALSE, xlab = "Time (Ma)", ylab = "",
     xlim = c(max(t_range[,"max_age"]), min(t_range[,"max_age"])))
segments(x0 = t_range[,"min_age"], y0 = t_range[,"pos"], x1 = t_range[,"max_age"],
         y1 = t_range[,"pos"], col = "gray30", lwd = 6, lend = 2)
axis(1, at = seq(from = round(max(t_range[,"max_age"])),
                 to = round(min(t_range[,"max_age"])),
                 length.out = 10),
     col = "gray30",
     cex.axis = 0.8)


rect(unit.depths[t],
     tscale_n[, "End"],
     unit.depths[t + 1],
     tscale_n[, "Start"],
     col = rgb(tscale_n[,"Col_R"],tscale_n[, "Col_G"], tscale_n[, "Col_B"],maxColorValue = 255)
     )



geoscalePlot(t_range[,1],t_range[,3], units = c("Period"), boxes = "Period",
              tick.scale = "Period", show.tip.label = TRUE,
              use.edge.length = TRUE, align.tip.label = FALSE)


  text(x = t_range[,"min_age"] - 0.3, y = t_range[,"pos"], labels = row.names(t_range),
       adj = c(0, 0), cex = 0.5, col = "gray30")



  Pseudosuchia.phy$tip.label




