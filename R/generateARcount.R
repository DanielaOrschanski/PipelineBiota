#' @title Generate AR counts
#' @description Generates one excel file for each taxonomic level containing the relative abundances of all patients in patient_dir.
#' @param patients_dir path to the folder that stores all the samples folders.
#' @return  list with: outD, out, ConteosKD.
#' @examples counts_taxonomy <- counts_Tax(otus_table_kraken04)
#' @export

generateARcounts <- function(patients_dir) {
  OTUs_11muestras_08K <- generateOTUsTableGrupal(patients_dir, source= "KRAKEN" )
  OTUs_11muestras_08D <- generateOTUsTableGrupal(patients_dir, source= "DRAGEN" )

  #Filtro los tables para quedarme con las mismas especies y comprar las abundancias ------------------------------

  dragenykraken <- merge(OTUs_11muestras_08D, OTUs_11muestras_08K, by = c( "Domain", "Kingdom", "Phylum", "Class", "Order", "Family","Genus", "Species"))
  colnames(dragenykraken) <- gsub("\\.y", "KRAKEN", colnames(dragenykraken))
  colnames(dragenykraken) <- gsub("\\.x", "DRAGEN", colnames(dragenykraken))

  dragen_08F <- dragenykraken[, -which(grepl("KRAKEN", colnames(dragenykraken)))]
  kraken_08F <- dragenykraken[, -which(grepl("DRAGEN", colnames(dragenykraken)))]
  all(dragen_08F[,1:8] == kraken_08F[,1:8])

  #-------------------------------------------------------------------------------------------------
  outD <- group_TaxonomicLevels(dragen_08F, source= "DRAGEN", patients_dir)
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

  out <- group_TaxonomicLevels(kraken_08F, source= "KRAKEN", patients_dir)
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

  #--------------------------------------------------------------

  conteosK <- counts_Tax(OTUs_11muestras_08K)
  conteosD <- counts_Tax(OTUs_11muestras_08D)

  ConteosKD <- data.frame(rbind(conteosK, conteosD))
  ConteosKD$source <- c("KRAKEN", "DRAGEN")
  ConteosKD <- as.data.frame(t(ConteosKD))
  colnames(ConteosKD) <- ConteosKD[ nrow(ConteosKD),]
  ConteosKD <- ConteosKD[-c(9,10),]
  #ConteosKD <- ConteosKD[-9,]

  return(list(outD, out, ConteosKD))
}
