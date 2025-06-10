
list_ids <- list.dirs(patients_dir, recursive = F, full.names = F)
i=1
test_Virus <- data.frame("ID" = character(length(list_ids)),
                         "Virus" = character(length(list_ids)),
                         stringsAsFactors = FALSE)

for (i in 1:length(list_ids)) {
  id <- list_ids[i]
  print(id)
  virus_presentes <- indicadorVirus(id = id, MetadataB = MetadataB)
  test_Virus$ID[i] <- id
  if(length(virus_presentes) == 0) {
    virus_presentes <- ""
  }
  test_Virus$Virus[i] <- virus_presentes
}

write.xlsx(test_Virus, file = "~/Daniela/Biota/Muestras/73m/Prueba_Virus_Mediax5.xlsx")



#Test indicador tipo de piel

list_ids <- list.dirs(patients_dir, recursive = F, full.names = F)
i=1
test_tipopiel <- data.frame("ID" = character(length(list_ids)),
                         "TipoPiel" = character(length(list_ids)),
                         stringsAsFactors = FALSE)

for (i in 1:length(list_ids)) {
  id <- list_ids[i]
  print(id)
  tipo_piel <- indicadorTipoPiel(id = id, MetadataB = MetadataB)
  tipo_piel <- tipo_piel[[1]]
  test_tipopiel$ID[i] <- id
  test_tipopiel$TipoPiel[i] <- paste0(tipo_piel, collapse = ",")
}

write.xlsx(test_tipopiel, file = "~/Daniela/Biota/Muestras/73m/Prueba_TipoPiel.xlsx")
