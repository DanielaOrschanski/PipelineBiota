
#patients_dir = "~/Daniela/Biota/Muestras/73m"
#source = "KRAKEN"
#de_host = "Bowtie"

#library(PipelineBiota)
extraerAR_Species <- function(patients_dir, source, de_host ) {
  tabla_otus <- generateOTUsTableGrupal(patients_dir = patients_dir, source = source, conEukaryota = FALSE, de_host = de_host)

  if(any(tabla_otus$Species == "Cutibacterium modestum")) {
    tabla_otus$Species[which(tabla_otus$Species == "Cutibacterium modestum")] <- "Propionibacterium sp. oral taxon 193"
    tabla_otus$Genus[which(tabla_otus$Species ==  "Propionibacterium sp. oral taxon 193")] <- "Propionibacterium"
  }

  out <- group_TaxonomicLevels(patients_dir = patients_dir, tabla_otus = tabla_otus, source = source, de_host = de_host, conE = FALSE)
  list_Conteos <- out[[6]]
  Conteos_species <- list_Conteos[[8]]
  #write.xlsx(Conteos_species, file = "~/Daniela/Biota/Muestras/Conteos_Species_94p.xlsx")

  list_AR <- out[[1]]
  AR_filos <- list_AR[[3]]
  AR_generos <- list_AR[[7]]
  AR_species <- list_AR[[8]]
  colnames(AR_species)[-1] <- gsub(sprintf("_%s", source), "", colnames(AR_species)[-1])

  otus_unique <- tabla_otus[!duplicated(tabla_otus$Species), c("Species", "Phylum")]
  AR_species_Filo <- merge(AR_species, otus_unique, by = "Species", all.x = TRUE, all.y = FALSE)
  #AR_species_Filo <- AR_species_fem[order(AR_species_fem$PromedioFem, decreasing = TRUE),]
  AR_species_Filo <- AR_species_Filo[-which(AR_species_Filo$Phylum == "-"),]

  #return(AR_species_Filo)
  return(AR_species)
}

