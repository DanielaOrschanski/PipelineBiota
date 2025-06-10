#' @title Run DRAGEN
#' @description kvjdfnkvjdf
#' @param patient_dir fskbfjsf
#' @param project_id it kvsbdkj
#' @param bs_path path to executable Base Space
#' @examples tabla_otus_kraken <- generateOTUsTable_Individual(patient_dir ="~/Biota/Nuevas13/Muestras/33", source= "KRAKEN" )
#' @export
#' @import tidyr
#' @import readr
#' @import stringr

ProcessFirst <- function(patients_dir, project_id, bs_path) {
  library(PipelineBiota)
  library(data.table)
  library(tidyr)
  library(dplyr)
  library(magrittr)
  library(readr)
  #patients_dir <- "~/Daniela/Biota/Muestras/73m"

  QCcontrol(patients_dir)

  list_dirs <- list.dirs(patients_dir, recursive = FALSE)
  #list_dirs <- list_dirs[1:10]
  #list_dirs <- list_dirs[-c(1,4)]
  #patient <- list_dirs[53]
  #list_dirs <- pacientes_sinDRAGENreport
  for (patient in list_dirs) {
    patient_dir <- patient
    print(patient)
    #primero las subo, tengo que esperar un rato que se carguen para poder mandarlas al dragen

    #uploadtoBS(patient_dir = patient, project_id = project_id, bs_path = bs_path, de_host = "")

    RunKRAKEN(patient_dir, de_host = "Bowtie")
    kraken <- generateOTUsTable_Individual(patient_dir, source= "KRAKEN", de_host = "Bowtie", conEukaryota = FALSE)
    counts_Tax(patient_dir, source= "KRAKEN", de_host = "Bowtie", conEukaryota = FALSE)

    dir_bs <- "~/Daniela/Biota"
    bs_path  <- sprintf("%s/bs", dir_bs)
    project_id <- "426787361"
    id <- basename(patient_dir)
    #RunDRAGEN(patient_dir = patient, project_id = project_id, bs_path = bs_path, de_host = "")
  }

  #Todo KRAKEN --------------
  RunKRAKEN(patients_dir)
  kraken <- generateOTUsTableGrupal(patients_dir, source= "KRAKEN" )
  counts_taxonomy_k <- counts_Tax(patients_dir, source= "KRAKEN")

  out <- group_TaxonomicLevels(patients_dir, kraken, source= "KRAKEN")
  list_AR <- out[[1]]
  AR_domain <- list_AR[[1]]
  AR_kingdom <- list_AR[[2]]
  AR_phylum <- list_AR[[3]]
  AR_class <- list_AR[[4]]
  AR_order <- list_AR[[5]]
  AR_family <- list_AR[[6]]
  AR_genus <- list_AR[[7]]
  AR_species <- list_AR[[8]]
  AR_subspecies <- list_AR[[9]]

  AR_esp_imp <- out[[2]]
  AR_sub_esp_imp <- out[[3]]
}


ProcessSecond <- function(patients_dir, project_id, bs_path) {
  list_dirs <- list.dirs(patients_dir, recursive = FALSE)
    #Tengo que esperar a que se terminen de hacer los analisis en dragen para hacer esto:
  list_dirs <- list_dirs[1:10]
  patient_dir <- list_dirs[68]
  for (patient_dir in list_dirs) {
    print(patient_dir)
    download_DRAGENReport(patient_dir, project_id = project_id, bs_path = bs_path)
    dragen <- generateOTUsTable_Individual(patient_dir = patient_dir, source= "DRAGEN")
  }
}

RunPipelineGrupal <- function(patients_dir, runDRAGEN = FALSE) {

  QCcontrol(patients_dir, de_host = "Bowtie")
  list.dirs(patients_dir, full.names = TRUE, recursive = FALSE)

  #KRAKEN ---------------------------------------------------------------------------
  RunKRAKEN(patients_dir, de_host = "Bowtie")
  #kraken <- generateOTUsTableGrupal(patients_dir, source= "KRAKEN", de_host = "Bowtie", conEukaryota = TRUE)
  kraken <- generateOTUsTableGrupal(patients_dir, source= "KRAKEN", de_host = "Bowtie", conEukaryota = FALSE)
  
  counts_taxonomy_k <- counts_Tax(patients_dir, source= "KRAKEN", de_host = "Bowtie", conEukaryota = FALSE)
  out <- group_TaxonomicLevels(patients_dir = patients_dir, tabla_otus = kraken, source = "KRAKEN", de_host = "Bowtie", conEukaryota= FALSE)
  list_AR <- out[[1]]
  AR_domain <- list_AR[[1]]
  AR_kingdom <- list_AR[[2]]
  AR_phylum <- list_AR[[3]]
  AR_class <- list_AR[[4]]
  AR_order <- list_AR[[5]]
  AR_family <- list_AR[[6]]
  AR_genus <- list_AR[[7]]
  AR_species <- list_AR[[8]]
  AR_subspecies <- list_AR[[9]]

  colSums(AR_phylum[,-1])
  #write.xlsx(AR_phylum, file = "/media/4tb2/Daniela/Biota/Muestras/AR_Filos_83p.xlsx")
  #write.xlsx(AR_genus, file = "/media/4tb2/Daniela/Biota/Muestras/AR_Generos_83p.xlsx")
  
  AR_esp_imp <- out[[2]]
  AR_sub_esp_imp <- out[[3]]

  #DRAGEN -----------------------------------------------------------------
  if(runDRAGEN == TRUE) {
    dir_bs <- "~/Daniela/Biota"
    bs_path  <- sprintf("%s/bs", dir_bs)

    project_id <- "426787361"

    RunDRAGEN(patients_dir,  project_id = project_id, bs_path = bs_path)
    #No se puede para muchas
    dragen <- generateOTUsTableGrupal(patients_dir, source= "DRAGEN" )
    counts_taxonomy_D <- counts_Tax(patients_dir, source= "DRAGEN" )

    #Abundancias relativas ----------------------------------------------------------
    outD <- group_TaxonomicLevels(patients_dir, dragen, source= "DRAGEN")
    list_ARD <- outD[[1]]

    AR_domainD <- list_ARD[[1]]
    AR_kingdomD <- list_ARD[[2]]
    AR_phylumD <- list_ARD[[3]]
    AR_classD <- list_ARD[[4]]
    AR_orderD <- list_ARD[[5]]
    AR_familyD <- list_ARD[[6]]
    AR_genusD <- list_ARD[[7]]
    AR_speciesD <- list_ARD[[8]]
    AR_subSpeciesD <- list_ARD[[9]]

    AR_esp_impD <- outD[[2]]
    AR_sub_esp_impD <- outD[[3]]

    #Filtro los tables para quedarme con las mismas especies y comprar las abundancias ------------------------------

    dragenykraken <- merge(dragen, kraken, by = c( "Domain", "Kingdom", "Phylum", "Class", "Order", "Family","Genus", "Species", "SubSpecies"))
    colnames(dragenykraken)[-c(1:9)] <- str_replace(colnames(dragenykraken)[-c(1:9)], ".x", "_DRAGEN")
    colnames(dragenykraken)[-c(1:9)] <- str_replace(colnames(dragenykraken)[-c(1:9)], ".y", "_KRAKEN")

    dragen_08F <- dragenykraken[, -which(grepl("KRAKEN", colnames(dragenykraken)))]
    kraken_08F <- dragenykraken[, -which(grepl("DRAGEN", colnames(dragenykraken)))]
    all(dragen_08F[,1:9] == kraken_08F[,1:9])

    #----------------------------------------------------------------------------------------

    outD <- group_TaxonomicLevels(patient_dir, dragen_08F, source= "DRAGEN")
    list_ARD <- outD[[1]]

    AR_domainD <- list_ARD[[1]]
    AR_kingdomD <- list_ARD[[2]]
    AR_phylumD <- list_ARD[[3]]
    AR_classD <- list_ARD[[4]]
    AR_orderD <- list_ARD[[5]]
    AR_familyD <- list_ARD[[6]]
    AR_genusD <- list_ARD[[7]]
    AR_speciesD <- list_ARD[[8]]
    AR_subspeciesD <- list_ARD[[9]]
    AR_esp_impD <- outD[[2]]
    AR_sub_esp_impD <- outD[[3]]

    out <- group_TaxonomicLevels(patient_dir, kraken_08F, source= "KRAKEN")
    list_AR <- out[[1]]
    AR_domain <- list_AR[[1]]
    AR_kingdom <- list_AR[[2]]
    AR_phylum <- list_AR[[3]]
    AR_class <- list_AR[[4]]
    AR_order <- list_AR[[5]]
    AR_family <- list_AR[[6]]
    AR_genus <- list_AR[[7]]
    AR_species <- list_AR[[8]]
    AR_subspecie <- list_AR[[9]]
    AR_esp_imp <- out[[2]]
    AR_sub_esp_imp <- out[[3]]

    #--------------------------------------------------------------

    ConteosKD <- data.frame(rbind(counts_taxonomy_k, counts_taxonomy_D))
    #ConteosKD$source <- c("KRAKEN", "DRAGEN")
    #ConteosKD <- as.data.frame(t(ConteosKD))
    #colnames(ConteosKD) <- ConteosKD[ nrow(ConteosKD),]

  }

}


# Visualizaciones ----------------------------
