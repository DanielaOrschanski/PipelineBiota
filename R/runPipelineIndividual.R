#' @title Run PipelineBiota Individual
#' @import reshape2
#' @import openxlsx
#' @description Process de patient given from the raw fastq file until 2 reports with the results.
#' @param patient_dir path to the directory of the patient you want to analyze. It must contain 2 files: ID_R1_fastq and ID_R2_fastq.
#' @param runDRAGEN is logical. If it set to TRUE the analyzes will be generated with both KRAKEN and DRAGEN.
#' @return path for 2 reports
#' @examples RunPipelineIndividual(patient_dir = "~/Daniela/Biota/Muestras/73m/114")
#' @export

RunPipelineIndividual <- function(patient_dir, runDRAGEN = FALSE) {
  id <- basename(patient_dir)
  print(id)
  QCcontrol(patient_dir, de_host = "Bowtie")

  #KRAKEN ---------------------------------------------------------------------------
  RunKRAKEN(patient_dir, de_host = "Bowtie")
  
  #kraken <- generateOTUsTable_Individual(patient_dir, source= "KRAKEN", de_host = "Bowtie", conEukaryota = TRUE)
  kraken <- generateOTUsTable_Individual(patient_dir, source= "KRAKEN", de_host = "Bowtie", conEukaryota = FALSE)
  library(readxl)
  #counts_taxonomy_k_bo <- counts_Tax(patients_dir = patient_dir, source= "KRAKEN", de_host = "Bowtie", conEukaryota = TRUE)
  library(dplyr)
  #group_TaxonomicLevels(patients_dir = patient_dir, tabla_otus = kraken, source = "KRAKEN", de_host = "Bowtie", conEukaryota= FALSE)

  simple_report_file <- generateSimpleReport(patient_dir = patient_dir, path_metadata = "/media/4tb2/Daniela/Biota/PipelineBiota-master/data/Metadata-Completa-SinLimpiar.xlsx")
  doctor_report_file <- generateCompleteReport(patient_dir = patient_dir)
  
  #DRAGEN -----------------------------------------------------------------
  if( runDRAGEN == TRUE) {

    RunDRAGEN(patient_dir, bs_path = "~/Daniela/Biota/bs", de_host = "", conEukaryota = TRUE)
    dragen <- generateOTUsTable_Individual(patient_dir, source= "DRAGEN", de_host = "", conEukaryota = TRUE)
    counts_taxonomy_D <- counts_Tax(patient_dir, source= "DRAGEN", de_host = "", conEukaryota = TRUE)

    #Bowtie ---
    uploadtoBS(patient_dir = patient_dir, bs_path = "~/Daniela/Biota/bs", de_host = "Bowtie")
    RunDRAGEN(patient_dir, bs_path = "~/Daniela/Biota/bs", de_host = "Bowtie", conEukaryota = TRUE)
    dragen <- generateOTUsTable_Individual(patient_dir, source= "DRAGEN", de_host = "Bowtie", conEukaryota = TRUE)
    counts_taxonomy_D <- counts_Tax(patient_dir, source= "DRAGEN", de_host = "Bowtie", conEukaryota = TRUE)

    #COMPARACION DRAGEN VS KRAKEN
    #Filtro los tables para quedarme con las mismas especies y comprar las abundancias ------------------------------

    dragenykraken <- merge(dragen, kraken, by = c( "Domain", "Kingdom", "Phylum", "Class", "Order", "Family","Genus", "Species", "SubSpecies"))
    colnames(dragenykraken)[10:11] <- c(sprintf("%s_DRAGEN", id), sprintf("%s_KRAKEN", id))
    dragen_08F <- dragenykraken[, -which(colnames(dragenykraken) == sprintf("%s_KRAKEN", id))]
    kraken_08F <- dragenykraken[, -which(colnames(dragenykraken) == sprintf("%s_DRAGEN", id))]
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
    AR_esp_imp <- out[[2]]
    AR_sub_esp_imp <- out[[3]]

    #--------------------------------------------------------------

    ConteosKD <- data.frame(rbind(counts_taxonomy_k, counts_taxonomy_D))
    ConteosKD$source <- c("KRAKEN", "DRAGEN")
    ConteosKD <- as.data.frame(t(ConteosKD))
    colnames(ConteosKD) <- ConteosKD[ nrow(ConteosKD),]
  }

  #return(report_file)
}
