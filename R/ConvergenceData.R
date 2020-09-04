## code to prepare `ConvergenceData` dataset goes here

usethis::use_data(ConvergenceData, overwrite = TRUE)

file.create("data-raw/MeanEmbryo_PC_Matrix.Rdata")
save(MeanEmbryo_PC_Matrix, file = "data-raw/MeanEmbryo_PC_Matrix.Rdata")

file.create("data-raw/Embryo.phy.Rdata")
save(Embryo.phy, file = "data-raw/Embryo.phy.Rdata")

file.create("data-raw/Embryo_means.Rdata")
save(Embryo_means, file = "data-raw/Embryo_means.Rdata")

file.create("data-raw/MeanJuvenile_PC_Matrix.Rdata")
save(MeanJuvenile_PC_Matrix, file = "data-raw/MeanJuvenile_PC_Matrix.Rdata")

file.create("data-raw/Juvenile.phy.Rdata")
save(Juvenile.phy, file = "data-raw/Juvenile.phy.Rdata")

file.create("data-raw/Juvenile_means.Rdata")
save(Juvenile_means, file = "data-raw/Juvenile_means.Rdata")

file.create("data-raw/MeanAdult_PC_Matrix.Rdata")
save(MeanAdult_PC_Matrix, file = "data-raw/MeanAdult_PC_Matrix.Rdata")

file.create("data-raw/Adult.phy.Rdata")
save(Adult.phy, file = "data-raw/Adult.phy.Rdata")

file.create("data-raw/Adult_means.Rdata")
save(Adult_means, file = "data-raw/Adult_means.Rdata")

