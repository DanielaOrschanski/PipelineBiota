#' @title Run HUMAnN
#' @description Executes concat of R1 and R2 fastq files, executes HUMAnN 3.0 and generates CPM.
#' @param patient_dir path to the folder of the patient that will be analyzed.
#' @param de_host indicates which alignment-based de-hosting software will be used.
#' @export

RunHuman <- function(patient_dir, de_host) {
  if(de_host  == "Bowtie") {
    de_host_file <- "Bo"
  } else if( de_host == "BWA") {
    de_host_file <- "bwa"
  } else if(de_host == "RSubread") {
    de_host_file <- "Rs"
  } else if(de_host == "") {
    de_host_file <- "T"
  } else if(de_host == "sinDH_PD") {
    de_host_file <- "sinDH_PD"
  } else {
    stop("de_host must be Bowtie, BWA, RSubread, sinDH_PD or empty string")
  }

  id <- basename(patient_dir)
  print(id)
  #1. Concatenar R1 y R2:
  concat_file <- sprintf("%s/Vias/%s/%s%s_concatR1R2.fastq.gz", patient_dir, de_host_file, id, de_host_file)

  if(!file.exists(concat_file) | file.info(concat_file)$size == 0) {
    message("Concatenando ...")
    R1 <- sprintf("%s/trimmed/%sDH%s_S04_L001_R1_001.fastq.gz", patient_dir, id, de_host_file)
    R2 <- sprintf("%s/trimmed/%sDH%s_S04_L001_R2_001.fastq.gz", patient_dir, id, de_host_file)

    if(de_host == "") {
      R1 <- sprintf("%s/trimmed/%s%s_S04_L001_R1_001.fastq.gz", patient_dir, id, de_host_file)
      R2 <- sprintf("%s/trimmed/%s%s_S04_L001_R2_001.fastq.gz", patient_dir, id, de_host_file)
    }

    dir.create(sprintf("%s/Vias", patient_dir))
    dir.create(sprintf("%s/Vias/%s", patient_dir, de_host_file))

    start_time_cat <- Sys.time()
    system(sprintf("cat %s %s > %s", R1, R2, concat_file))
    end_time_cat <- Sys.time()
    cat_duration <- end_time_cat - start_time_cat

    #chmod u+r $output_file
  }


  #2. Ejecutar humann del concatenado:

  vias_dir <-  sprintf("%s/Vias/%s", patient_dir, de_host_file)

  path_human <- "/home/daniela/miniconda3/envs/biobakery3/bin/humann"

  command <- paste(
    "bash -c 'source activate biobakery3 &&",  # Activar el entorno
    path_human,                                # Ejecutar humann
    sprintf("--input %s", concat_file),
    sprintf("--output %s", vias_dir),
    "--threads 15'",
    sep = " "
  )

  tsv_file <- sprintf("%s/%s%s_concatR1R2_pathabundance.tsv", vias_dir, id, de_host_file)
  if(!file.exists(tsv_file)) {
    message("Ejecutando human ...")
    start_time_human <- Sys.time()
    system(command = command, intern = FALSE)
    end_time_human <- Sys.time()
    human_duration <- end_time_human - start_time_human

  } else {
    message("Human was already executed before!")
  }

  #Agregar el conteo en CPM: ---------------------------
  rpk_file <- sprintf("%s/Vias/%s/%s%s_concatR1R2_pathabundance.tsv", patient_dir, de_host_file, id, de_host_file)
  dir.create(sprintf("/home/daniela/Daniela/Biota/Muestras/SubsetPathways/%s", de_host_file))
  cpm_file <- sprintf("/home/daniela/Daniela/Biota/Muestras/SubsetPathways/%s/%s%s_CPM_pathabundance.tsv", de_host_file, id, de_host_file)

  #cpm_file <- sprintf("/home/daniela/Daniela/Biota/Muestras/SubsetPathways/%s/%s%s_CPM_pathabundance.tsv", de_host_file, id, de_host_file)

  command <- paste(
    "bash -c 'source activate biobakery3 &&",  # Activar el entorno
    "humann_renorm_table",                                # Ejecutar humann
    sprintf("--input %s", rpk_file),
    sprintf("--output %s", cpm_file),
    "--units cpm",
    "--update-snames'",
    sep = " "
  )
  if(!file.exists(cpm_file)) {
    start_time_cpm <- Sys.time()
    system(command = command, intern = FALSE)
    end_time_cpm <- Sys.time()
    cpm_duration <- end_time_cpm - start_time_cpm
  } else {
    message("CPM file has already been generated for this patient!")
  }

  #print(sprintf("Concat duró %s - Human duró %s - CPM duró %s", cat_duration, human_duration, cpm_duration))

  return(cpm_file)

}


#' @title Run HUMAnN
#' @description Executes concat of R1 and R2 fastq files, executes HUMAnN 3.0 and generates CPM.
#' @param patient_dir path to the folder of the patient that will be analyzed.
#' @param de_host indicates which alignment-based de-hosting software will be used.
#' @export
generatePathwayReport <- function(id, de_host) {

  if(de_host  == "Bowtie") {
    de_host_file <- "Bo"
  } else if( de_host == "BWA") {
    de_host_file <- "bwa"
  } else if(de_host == "RSubread") {
    de_host_file <- "Rs"
  } else if(de_host == "") {
    de_host_file <- "T"
  }

  report_ARvias_file <- sprintf("~/Daniela/Biota/Muestras/73m/%s/Vias/AR%s%s_vias.xlsx",id, id, de_host_file )
  report_Countvias_file <- sprintf("~/Daniela/Biota/Muestras/73m/%s/Vias/Count%s%s_vias.xlsx",id, id, de_host_file )

  #if(!file.exists(report_vias_file)) {
  message("Generating Pathway Report ...")

  path_abundance <- read_tsv(sprintf("~/Daniela/Biota/Muestras/73m/%s/Vias/%s/%s%s_concatR1R2_pathabundance.tsv", id, de_host_file, id, de_host_file))

  colnames(path_abundance)[1] <- "Pathway"

  # Separar la columna en tres columnas: pathway, description, organism
  df_separado <- path_abundance %>%
    separate(Pathway, into = c("Pathway", "description_organism"), sep = ":", extra = "merge") %>%
    separate(description_organism, into = c("Description", "Organism"), sep = "\\|")

  sum(df_separado[,ncol(df_separado)])


  df_sin_unmapped <- df_separado[-which(grepl("UNMAPPED", df_separado$Pathway) | grepl("UNINTEGRATED", df_separado$Pathway)),]
  write.xlsx(df_sin_unmapped, file = report_Countvias_file)

  AR_vias <- prop.table(df_sin_unmapped[,ncol(df_sin_unmapped)])*100
  AR_vias <- cbind(df_sin_unmapped, AR_vias)
  AR_vias <- AR_vias[,-3]
  colnames(AR_vias)[3] <- paste0("Conteo_", id, sep="")
  colnames(AR_vias)[4] <- paste0("AR_", id, sep="")

  write.xlsx(AR_vias, file = report_ARvias_file)
  #} else {
  #  message("The Pathway report has already been generated for this patient!")
  #}

}

#' @title Generate Pathways Table
#' @description Joins CPM tables, adds information of classes from Metacyc and calculates Relative Abundances.
#' @param patients_dir path to the folder of the patients that will be analyzed.
#' @param de_host indicates which alignment-based de-hosting software will be used.
#' @export

generatePathwaysTable <- function(patients_dir, de_host) {

  if(de_host  == "Bowtie") {
    de_host_file <- "Bo"
  } else if( de_host == "BWA") {
    de_host_file <- "bwa"
  } else if(de_host == "RSubread") {
    de_host_file <- "Rs"
  } else if(de_host == "") {
    de_host_file <- "T"
  }

  list_patients <- list.dirs(patients_dir, recursive = FALSE, full.names = TRUE)
  #p <- list_patients[1]
  list_cpm_files <- c()
  for (p in list_patients) {
    print(p)
    cpm_file <- RunHuman(patient_dir = p, de_host = de_host)
    list_cpm_files <- c(list_cpm_files, cpm_file)
  }

  pathway_dir <- path.expand(sprintf("~/Daniela/Biota/Muestras/SubsetPathways/%s", de_host_file))
  #patients <- list.files(pathway_dir, recursive = FALSE, full.names = TRUE)
  #patients <- patients[which(grepl(sprintf("%s_CPM_pathabundance.tsv", de_host_file), patients))]
  patients <- list_cpm_files

  #Unir tsv de CPM:
  cpm_table <- path.expand(sprintf("~/Daniela/Biota/Muestras/SubsetPathways/%s/subset_CPM_pathabundance.tsv", de_host_file))
  if(file.exists(cpm_table)) {
    file.remove(cpm_table)
  }

  #quedarme solo con los pacientes que quiero:
  ids_quiero <- basename(list_patients)

  tsv_tengo <- list.files(pathway_dir)
  tsv_tengo <- tsv_tengo[grepl("_CPM_pathabundance.tsv", tsv_tengo)]
  ids_tengo <- gsub(sprintf("%s_CPM_pathabundance.tsv", de_host_file), "", tsv_tengo)

  ids_a_eliminar <- setdiff(ids_tengo, ids_quiero)
  for(id in ids_a_eliminar){
    print(sprintf("The ID %s will be removed from SubsetPathway/%s folder", id, de_host_file))
    file.remove(sprintf("%s/%s%s_CPM_pathabundance.tsv", pathway_dir, id, de_host_file))
  }

  command <- paste(
    "bash -c 'source activate biobakery3 &&",  # Activar el entorno
    "humann_join_tables",                                # Ejecutar humann
    sprintf("--input %s", pathway_dir),
    sprintf("--output %s", cpm_table),
    "--file_name pathabundance'",
    sep = " "
  )
  system(command = command, intern = FALSE)

  library(readr)
  library(dplyr)
  library(tidyr)
  path_abundance <- as.data.frame(read_tsv(cpm_table))
  colnames(path_abundance)[1] <-"FEATURE"
  colnames(path_abundance)[-1] <- gsub(sprintf("%s_concatR1R2_Abundance-CPM", de_host_file), "", colnames(path_abundance)[-1] )



  # Separar la columna en tres columnas: pathway, description, organism
  path_ab <- path_abundance
  colnames(path_ab)[1] <-"Pathway"
  path_ab <- path_ab %>%
    separate(Pathway, into = c("Pathway", "description_organism"), sep = ":", extra = "merge") %>%
    separate(description_organism, into = c("Description", "Organism"), sep = "\\|")

  #Eliminar unmapped:
  df_sin_unmapped <- path_ab[-which(grepl("UNMAPPED", path_ab$Pathway) | grepl("UNINTEGRATED", path_ab$Pathway)),]
  str(df_sin_unmapped)
  length(unique(df_sin_unmapped$Pathway))

  ids_patients <- list.dirs(patients_dir, recursive = FALSE, full.names = FALSE)
  all(ids_patients %in% colnames(df_sin_unmapped))
  df_sin_unmapped <- df_sin_unmapped[, which(colnames(df_sin_unmapped) %in% c(ids_patients, "Pathway", "Description", "Organism") )]

  #Para conservar la info de los organismos que contribuyen en cada via: ------------------------------------------------------
  vias_con_org <- df_sin_unmapped
  library(dplyr)
  library(stringr)

  vias_con_org <- vias_con_org %>%
    mutate(
      Generos = str_extract(Organism, "g__[^.]+"),
      Especies = str_extract(Organism, "s__[^.]+")
    )
  vias_con_org <- vias_con_org[,c(1,2,3, (ncol(vias_con_org))-1, ncol(vias_con_org), 4:(ncol(vias_con_org)-2) )]
  vias_con_org <- vias_con_org %>%
    mutate(
      Generos = str_replace(Generos, "g__", ""),
      Especies = str_replace(Especies, "s__", "")
    )

  library(openxlsx)
  write.xlsx(vias_con_org, file = sprintf("~/Daniela/Biota/Muestras/SubsetPathways/%s/%s_CPM_Vias_Organismos_%sp.xlsx", de_host_file, de_host_file, length(patients)))


  #Me quedo con las vias que no tienen organismos identificados: conteos totales --------------------------------
  df_sin_unmapped <- df_sin_unmapped[which(is.na(df_sin_unmapped$Organism)),]
  df_sin_unmapped <- df_sin_unmapped[-1,-3]
  df_sin_unmapped[,-c(1,2)] <- lapply(df_sin_unmapped[,-c(1,2)], as.numeric)
  str(df_sin_unmapped)

  #cpm_final_table <- sprintf("~/Daniela/Biota/Muestras/SubsetPathways/%s/%s_CPM_Vias_%sp.xlsx", de_host_file, de_host_file, length(patients))

  # Agregar info de clases:
  clases_vias <- read.table("~/Daniela/Biota/Muestras/SubsetPathways/ontology_All_pathways_MetaCyc.txt", header = TRUE, sep = "\t")
  #CPM_vias_completo <- as.data.frame(read_excel("~/Daniela/Biota/Muestras/SubsetPathways/Bo/CPM_vias_completo.xlsx"))

  CPM_vias <- df_sin_unmapped
  all(CPM_vias$Pathway %in% clases_vias$Pathways)

  v=1
  for (v in 1:nrow(CPM_vias)) {
    via <- CPM_vias$Pathway[v]
    print(via)
    clase <- clases_vias$Ontology...parents.of.class[which(clases_vias$Pathways == via)]
    CPM_vias$Clase[v] <- clase
  }

  length(unique(CPM_vias$Clase))
  CPM_vias <- cbind("Pathway" = CPM_vias$Pathway, "Clase" = CPM_vias$Clase, CPM_vias[,-which(colnames(CPM_vias) %in% c("Clase", "Pathway"))])

  write.xlsx(CPM_vias, file = sprintf("~/Daniela/Biota/Muestras/SubsetPathways/%s/%s_CPM_Vias_Clases_%sp.xlsx", de_host_file, de_host_file, length(patients)))

  #str(CPM_vias)

  AR_vias <- prop.table(as.matrix(CPM_vias[, -c(1,2,3)]), margin = 2) * 100
  AR_vias <- as.data.frame(cbind("Pathway" = CPM_vias[,1], "Clase" = CPM_vias[,2], "Description" = CPM_vias[,3], AR_vias))
  str(AR_vias)
  #AR_vias[,-c(1,2,3)] <- lapply(AR_vias[,-c(1,2,3)], as.numeric)
  AR_vias[,-c(1,2,3)] <- lapply(AR_vias[,-c(1,2,3)], function(x) as.numeric(gsub(",", ".", x)))

  colSums(AR_vias[,-c(1,2,3)])

  if(any(colSums(is.na(AR_vias)))) {
    #cantidad de NA en una columna = cantidad de filas
    AR_vias[, colSums(is.na(AR_vias)) == nrow(AR_vias)] <- 0
  }

  AR_final_table <- sprintf("~/Daniela/Biota/Muestras/SubsetPathways/%s/%s_AR_Vias_%sp.xlsx", de_host_file, de_host_file, length(patients))
  write.xlsx(AR_vias, file = AR_final_table, overwrite = TRUE)

  return(list(CPM_vias, AR_vias, vias_con_org))

}


#' @title Count Vias
#' @description Generates a df that counts the number of identified pathways for each sample.
#' @param patients_dir path to the folder of the patients that will be analyzed.
#' @param de_host indicates which alignment-based de-hosting software will be used.
#' @export

count_Vias <- function(patients_dir, de_host) {

  list_dirs <- list.dirs(patients_dir, recursive = FALSE, full.names = TRUE)
  #list_dirs <- list_dirs[1:50]
  #patient_dir <- list_dirs[1]

  if(length(nchar(list_dirs)) == 0 |  length(nchar(list_dirs[endsWith(list_dirs, "/trimmed")])) !=0 ) {
    list_dirs <- patients_dir
  }

  if(de_host  == "Bowtie") {
    de_host_file <- "Bo"
  } else if( de_host == "BWA") {
    de_host_file <- "bwa"
  } else if(de_host == "RSubread") {
    de_host_file <- "Rs"
  } else if(de_host == "") {
    de_host_file <- "T"
  } else if(de_host == "crudas") {
    de_host_file <- "cruda"
  } else if(de_host == "Biota") {
    de_host_file <- "trimmed"
  } else if(de_host == "sinDH_PD") {
    de_host_file <- "sinDH_PD"
  } else {
    stop("de_host must be Bowtie, BWA, RSubread, sinDH_PD or empty string")
  }

  count_Vias <- data.frame("ID" = NA,
                           "Vias" = NA,
                           "de_host" = de_host_file)
  i=1
  #patient_dir <- list_dirs[1]
  for (patient_dir in list_dirs) {

    id <- basename(patient_dir)

    print(patient_dir)
    tabla_vias <-  read_excel(paste(patient_dir, "/Vias/AR", id, de_host_file,"_vias.xlsx", sep=""))

    if(any(tabla_vias[,ncol(tabla_vias)] == 0)) {
      tabla_vias <- tabla_vias[-which(tabla_vias[,ncol(tabla_vias)] == 0),]
    }

    cants_original <- data.frame("ID" = id,
                                 "Vias" = length(unique(tabla_vias$Pathway)),
                                 "de_host" = de_host_file)

    count_Vias[i,] <- cants_original[1,]
    i <- i+1

  }

  #count_Vias_report <- as.data.frame(t(count_Vias))
  #colnames(count_Vias_report) <- count_Vias_report[1,]
  #count_Vias_report <- count_Vias_report[-1,]

  return(count_Vias)
}
