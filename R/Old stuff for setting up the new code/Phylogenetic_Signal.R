##Calculate phylogenetic signal##

##Load in phylogeny##

Pseudosuchia_calibrated.phy$tip.label
##

##Create mean adult/subadult shapes for all species##
## list of species with adults/subadults
Adult_coords <- Adult_CrocDorsal.gpa$coords
dimnames(Adult_coords)[[3]] <- Adult_CrocDorsal.classifier$Species
Unique_Adult_Species_list <- unique(Adult_CrocDorsal.classifier$Species)
##

# create an array for the mean shapes of each species
Adult_mean_coords <- array(dim = c(14,2,length(Unique_Adult_Species_list)), dimnames = list(c(1:14),c("X","Y"),Unique_Adult_Species_list))

for (i in 1:length(Unique_Adult_Species_list)){
  current_sp <- Unique_Adult_Species_list[[i]]
  current_specimen <- grep(current_sp, dimnames(Adult_coords)[[3]])
  current_meanshape <- mshape(Adult_coords[,,current_specimen])

  Adult_mean_coords[,,i] <- current_meanshape
}
##

## Reorder coords to match phylogeny ##
Mean_coords_phy_order <- Adult_mean_coords[,,match(Pseudosuchia_calibrated.phy$tip.label,gsub(" ", "_", dimnames(Adult_mean_coords)[[3]]))]
dimnames(Mean_coords_phy_order)[[3]] <- Pseudosuchia_calibrated.phy$tip.label
##

## Peform phylogenetic signal calculation ##
physignal(A=Mean_coords_phy_order, phy=Pseudosuchia_calibrated.phy)
##

