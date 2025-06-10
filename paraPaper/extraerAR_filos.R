
extraerAR_Filos <- function(patients_dir, source, de_host ) {
  tabla_otus <- generateOTUsTableGrupal(patients_dir = patients_dir, source = source, conEukaryota = FALSE, de_host = de_host)
  out <- group_TaxonomicLevels(patients_dir = patients_dir, tabla_otus = tabla_otus, source = source, de_host = de_host, conE = FALSE)
  list_AR <- out[[1]]
  AR_filos <- list_AR[[3]]

  colnames(AR_filos)[-1] <- gsub(sprintf("_%s", source), "", colnames(AR_filos)[-1])

  return(AR_filos)
}
