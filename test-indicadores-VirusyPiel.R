
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
patients_dir <- "/home/daniela/Daniela/Biota/Muestras/73m"
MetadataB <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx"))

list_ids <- list.dirs(patients_dir, recursive = F, full.names = T)
i=1
test_tipopiel <- data.frame("ID" = character(length(list_ids)),
                         "TipoPielPredicho" = character(length(list_ids)),
                         "TipoPielReal" = character(length(list_ids)),
                         stringsAsFactors = FALSE)

for (i in 1:length(list_ids)) {
  id <- list_ids[i]
  print(id)
  tipo_piel <- indicadorTipoPiel(patient_dir = id, MetadataB = MetadataB)
  tipo_piel <- tipo_piel[[1]]
  test_tipopiel$ID[i] <- basename(id)
  test_tipopiel$TipoPielPredicho[i] <- paste0(tipo_piel, collapse = ",")

  if(length(nchar( MetadataB$Tipodepiel[which(MetadataB$ID == basename(id))])) ==0) {
    test_tipopiel$TipoPielReal[i] <- NA
  } else {
    test_tipopiel$TipoPielReal[i] <- MetadataB$Tipodepiel[which(MetadataB$ID == basename(id))]
  }

}

#Ver que da las de acne:
indicadorTipoPiel(patient_dir = "~/Daniela/Biota/Muestras/acne/25", MetadataB = MetadataB)
indicadorTipoPiel(patient_dir = "~/Daniela/Biota/Muestras/acne/217", MetadataB = MetadataB)

test_tipopiel$match <- apply(test_tipopiel, 1, function(row) {
  grepl(row["TipoPielReal"], row["TipoPielPredicho"])
})

table(test_tipopiel$match)

write.xlsx(test_tipopiel, file = "~/Daniela/Biota/Muestras/73m/Prueba_TipoPiel.xlsx")



# Test indicador rango etario:
list_ids <- list.dirs(patients_dir, recursive = F, full.names = F)
i=1
test_Rango <- data.frame("ID" = character(length(list_ids)),
                         "RangoPredicho" = character(length(list_ids)),
                         "Edad" = character(length(list_ids)),
                         stringsAsFactors = FALSE)

for (i in 1:length(list_ids)) {
  id <- list_ids[i]
  print(id)
  resultado <- indicadorRangoEtario(id = id, MetadataB = MetadataB)
  resultado <- resultado[[2]]

  test_Rango$ID[i] <- id
  test_Rango$RangoPredicho[i] <- resultado
  test_Rango$Edad[i] <- MetadataB$Edad[which(MetadataB$ID == id)]
}

write.xlsx(test_Virus, file = "~/Daniela/Biota/Muestras/73m/Prueba_Virus_Mediax5.xlsx")

