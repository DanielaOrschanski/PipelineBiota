#' @title upload to BS (Base Space)
#' @description kvjdfnkvjdf
#' @param patient_dir fskbfjsf allows fastq and fastq.gz
#' @param project_id it kvsbdkj
#' @param bs_path path to base space executable
#' @examples tabla_otus_kraken <- generateOTUsTable_Individual(patient_dir ="~/Biota/Nuevas13/Muestras/33", source= "KRAKEN" )
#' @export

uploadtoBS <- function(patient_dir, bs_path, de_host = "") {

  id <- basename(patient_dir)

  if(de_host  == "Bowtie") {
    de_host_file <- "Bo"
    project_id <- "438230795"
    carpeta_bs <- "dehostBowtie"
    nombre_p <- sprintf("%sDHBo", id)

  } else if( de_host == "BWA") {
    de_host_file <- "bwa"
    project_id <- "436873438"
    carpeta_bs <- "dehostBWA"
    nombre_p <- sprintf("%sDHbwa", id)
  } else if(de_host == "RSubread") {
    de_host_file <- "Rs"
    project_id = "436780349"
    carpeta_bs ="dehostRsubread"
    nombre_p <- sprintf("%sDHRs", id)
  } else if(de_host == "") {
    de_host_file <- "sin"
    project_id <- "436799364"
    carpeta_bs <- "MuestrasTrimmed"
    nombre_p <- sprintf("%sT", id)
  } else if(de_host == "crudas") {
    de_host_file <- "cruda"
    project_id <- "438331896"
    carpeta_bs <- "MuestrasCrudas"
    nombre_p <-  paste0(id, "C", sep="")
  }

  patient <- basename(patient_dir)
  patient_dir <- paste(patient_dir, "/trimmed", sep="")

  #Chequea en el trimmed pero sube el fastq original
  nuevo_path_report <- sprintf("%s/Resultados_DRAGEN/%s%s.DRAGEN-report.tsv", patient_dir, patient, de_host_file)

  if(file.exists(nuevo_path_report)) {
    return("This patient has already been processed with DRAGEN")
  }

  #Identifico R1 y R2

  print(patient)
  #setwd(patient_dir)
  list_files <- list.files(patient_dir, full.names = TRUE)

  #fileR1 <- list_files[grepl("R1_001.fastq", list_files) & !(grepl("fastqc", list_files)) & !(grepl(".txt", list_files))]
  #fileR2 <- list_files[grepl("R2_001.fastq", list_files)  & !(grepl("fastqc", list_files)) & !(grepl(".txt", list_files))]

  fileR1 <- list_files[grepl(sprintf("%sDH%s_S04_L001_R1_001.fastq.gz", patient, de_host_file), list_files)]
  fileR2 <- list_files[grepl(sprintf("%sDH%s_S04_L001_R2_001.fastq.gz", patient, de_host_file), list_files)]

  #fileR1 <- "~/Daniela/Biota/Muestras/73m/108/108C_S04_L001_R1_001.fastq.gz"
  #fileR2 <- "~/Daniela/Biota/Muestras/73m/108/108C_S04_L001_R2_001.fastq.gz"

  if(de_host == "") {
    fileR1 <- list_files[grepl("T_S04_L001_R1_001.fastq", list_files)]
    fileR2 <- list_files[grepl("T_S04_L001_R2_001.fastq", list_files)]
  }
  if(de_host == "crudas") {
    fileR1 <- list_files[grepl("C_S04_L001_R1_001.fastq", list_files)]
    fileR2 <- list_files[grepl("C_S04_L001_R2_001.fastq", list_files)]
  }

  #Cambio los nombres para que cumplan el formato de BaseSapce
  #new_names <- c()
  #for (r in c(fileR1, fileR2)) {
  #  file_name <- r
  #  file_parts <- strsplit(basename(file_name), split= "_")
  #  new_file_name <- paste0(dirname(file_name), "/", file_parts[[1]][1], "_", file_parts[[1]][2], "_L001_", file_parts[[1]][3], "_", file_parts[[1]][4])
  #  file.rename(file_name, new_file_name)
  #  new_names <- c(new_names, new_file_name)
  #}

  #fileR1 <- new_names[1]
  #fileR2 <- new_names[2]


  #Subo a base space ----------------------------------------
  command <- sprintf("%s list biosample", bs_path)
  output <- capture.output(system(command, intern = TRUE))
  output <- output[!grepl("^\\+|^\\s*$", output)]
  output <- gsub("^\\|\\s*|\\s*\\|$", "", output)
  data_list <- strsplit(output, "\\s*\\|\\s*")
  df_output <- do.call(rbind, data_list)

  colnames(df_output) <- trimws(df_output[2,])
  df_output <- df_output[-c(1:3),]

  biosamples <- df_output[,2]


  #if(!(nombre_p %in% biosamples)) {
    command <- sprintf("%s upload dataset --name %s -p %s %s %s", bs_path, carpeta_bs, project_id, fileR1, fileR2)
    system(command)
  #} else {
  #  return(message(sprintf("This patient %s has already been uploaded to BS", nombre_p)))
  #}

}


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

RunDRAGEN <- function(patient_dir, bs_path, reference = "hg38", de_host, conEukaryota = TRUE) {
  #Identificar el id de basespace de la muestra
  #Veo la lista de los ids de las muestras del proyecto
  #command <- sprintf("%s list biosamples --project-id %s", bs_path, project_id)
  if(conEukaryota == TRUE) {
    conE <- "conEUKARYOTA"
  } else {
    conE <- ""
  }


  id <- basename(patient_dir)
  if(id == "trimmed") {
    id <- basename(dirname(patient_dir))
  }

  if(de_host  == "Bowtie") {
    de_host_file <- "Bo"
    project_id = "438230795"
    carpeta_bs = "dehostBowtie"
    nombre_p <- sprintf("%sDHBo", id)
    dsin_dehost = TRUE

  } else if( de_host == "BWA") {
    de_host_file <- "bwa"
    project_id = "436873438"
    carpeta_bs = "dehostBWA"
    nombre_p <- sprintf("%sDHbwa", id)
    dsin_dehost = TRUE

  } else if(de_host == "RSubread") {
    de_host_file <- "Rs"
    project_id = "436780349"
    carpeta_bs ="dehostRsubread"
    nombre_p <- sprintf("%sDHRs", id)
    dsin_dehost = TRUE

  } else if(de_host == "") {
    de_host_file <- "sin"
    project_id = "436799364"
    carpeta_bs = "MuestrasTrimmed"
    nombre_p <- sprintf("%sT", id)
    dsin_dehost = FALSE #ejecuta dehosting con dragen a trimmeadas

  } else if(de_host == "crudas") {
    de_host_file <- "cruda"
    project_id <- "438331896"
    carpeta_bs <- "MuestrasCrudas"
    nombre_p <-  paste0(id, "C", sep="")
    dsin_dehost = FALSE #ejecuta dehosting con dragen a crudas

  } else if(de_host == "sinDH_PD") { # sin de host previo ni dehost dragen
    de_host_file <- "sinDH_PD"
    project_id = "436799364"
    carpeta_bs = "MuestrasTrimmed"
    nombre_p <- sprintf("%sT", id)
    dsin_dehost = TRUE
  }

  patient <- basename(patient_dir)
  patient_dir <- paste(patient_dir, "/trimmed", sep="")

  if (!(file.exists(sprintf("%s/Resultados_DRAGEN/%s%s_DRAGEN-report.tsv", patient_dir, patient, de_host_file)))) {

    #controlar que no se haya lanzado el analysis anteriormente ----
    command <- sprintf("%s list appsession --project-id %s", bs_path, project_id)
    #print(system(command))

    output <- capture.output(system(command, intern = TRUE))
    df_output <- as.data.frame(do.call(rbind, strsplit(output, " +", perl = TRUE)))
    df_output <- df_output[which(df_output$V3 == sprintf("DRAGENHG38_%s", de_host_file) ),]
    df_output <- df_output[, c(3,4,5,6)]
    #colnames(df_output) <- c("DRAGEN", "SAMPLE", "ID")

    if( (patient %in% df_output[,3]) | (patient %in% df_output[,4])| (patient %in% df_output[,2]) ) {
      message(sprintf("The DRAGEN aplication has already been launched for this patient: %s", patient))
      return(message(sprintf("The DRAGEN aplication has already been launched for this patient: %s", patient)))
    }

    #---------------------

    command <- sprintf("%s list biosamples --project-id %s", bs_path, project_id)
    output <- capture.output(system(command, intern = TRUE))
    df_output <- as.data.frame(do.call(rbind, strsplit(gsub("[|]", "", output), " +")))
    df_output <- df_output[,-c(1,2)]
   #df_output <- as.data.frame(do.call(rbind, strsplit(output, " | ", perl = TRUE)))

    #
    colnames(df_output) <- df_output[2,]
    df_output <- df_output[-c(1,2,3, nrow(df_output)),]
    df_output <- df_output[, 1:3]
    colnames(df_output)[1] <- "BioSampleName1"

    id_bs <- df_output$Id[which(df_output$BioSampleName == nombre_p)]
    if(length(nchar(id_bs)) == 0 ){
      id_bs <- df_output$BioSampleName[which(df_output$BioSampleName1 == nombre_p)]
    }

    #Puede ser que salga en formatos raros la tabla
    if(length(nchar(id_bs))==0) {
      id_bs <- df_output$`|`[which(df_output$`"|` == patient)]
    }

    if( length(nchar(id_bs)) == 0) {
      output <- system(command, intern = TRUE)
      output_text <- paste(output, collapse = "\n")
      output_lines <- unlist(strsplit(output_text, "\n"))
      data_lines <- output_lines[-which(grepl("-", output_lines))]
      data_lines <- data.frame(data_lines)
      data_lines <- as.data.frame(lapply(data_lines, function(x) gsub("^\\||\\|$", "", x)))
      data_lines <- data_lines %>%
        separate(data_lines, into = c("BioSampleName", "Id", "ContainerName", "ContainerPosition", "Status"), sep = "\\|")
      id_bs <- data_lines$Id[which(grepl(patient, data_lines$BioSampleName))]
      id_bs <- gsub(" ", "", id_bs)
    }

    if(length(nchar(id_bs)) == 0) {
      print(system(command, intern = TRUE))
      print(patient)
      id_bs <- readline(prompt = "Enter the id of your sample: ")
    }

    ref <- ifelse(reference == "hg19", "hg19-altaware-cnv-anchor.v8", "hg38-altaware-cnv-anchor.v8")
    # Definir el comando bs con todos los argumentos
    #system(sprintf('%s launch application -n "DRAGEN Metagenomics Pipeline" --app-version 3.5.12 --list', bs_path))

    if(dsin_dehost == TRUE) { #DRAGEN SIN DEHOSTING
      command <- sprintf('%s launch application -n "DRAGEN Metagenomics Pipeline" --app-version 3.5.12 -o app-session-name:"DRAGENHG38_%s %s" -l "DRAGENHG38_%s %s" -o project-id:%s -o sample-id:%s -o ht-ref:%s -o dehost-checkbox:false -o basespace-labs-disclaimer:Accepted -o db:minikraken20200312',
                         bs_path,
                         de_host_file,
                         patient,
                         de_host_file,
                         patient,
                         project_id,
                         id_bs,
                         ref)


      system(command)
    } else { # DRAGEN CON DEHOSTING
      command <- sprintf('%s launch application -n "DRAGEN Metagenomics Pipeline" --app-version 3.5.12 -o app-session-name:"DRAGENHG38_%s %s" -l "DRAGENHG38_%s %s" -o project-id:%s -o sample-id:%s -o ht-ref:%s -o basespace-labs-disclaimer:Accepted -o db:minikraken20200312',
                         bs_path,
                         de_host_file,
                         patient,
                         de_host_file,
                         patient,
                         project_id,
                         id_bs,
                         ref)
      system(command)
    }

  } else {
    return(message(sprintf("The DRAGEN report has already been generated for this patient: %s", patient)))
  }

}

#' @title Download DRAGEN Report
#' @description kvjdfnkvjdf
#' @param patient_dir fskbfjsf
#' @param project_id it kvsbdkj
#' @param bs_path path to bs
#' @examples path_report <- download_DRAGENReport(patient_dir ="~/Biota/Nuevas13/Muestras/33" )
#' @export
#' @import tidyr
# bs_path = "~/Daniela/Biota/PipelineBiota/data/bs"
download_DRAGENReport <- function(patient_dir, bs_path, de_host) {

  patient <- basename(patient_dir)
  print(patient)
  id <- patient
  patient_dir <- paste(patient_dir, "/trimmed", sep="")

  if(de_host  == "Bowtie") {
    de_host_file <- "Bo"
    project_id = "438230795"
    carpeta_bs = "dehostBowtie"
    nombre_p <- sprintf("%sDHBo", id)
    dsin_dehost = TRUE
  } else if( de_host == "BWA") {
    de_host_file <- "bwa"
    project_id = "436873438"
    carpeta_bs = "dehostBWA"
    nombre_p <- sprintf("%sDHbwa", id)
    dsin_dehost = TRUE
  } else if(de_host == "RSubread") {
    de_host_file <- "Rs"
    project_id = "436780349"
    carpeta_bs ="dehostRsubread"
    nombre_p <- sprintf("%sDHRs", id)
    dsin_dehost = TRUE
  } else if(de_host == "") {
    de_host_file <- "sin"
    project_id = "436799364"
    carpeta_bs = "MuestrasTrimmed"
    nombre_p <- sprintf("%sT", id)
    dsin_dehost = FALSE
  } else if(de_host == "sinDH_PD") { # sin de host previo ni dehost dragen
    de_host_file <- "sinDH_PD"
    project_id = "436799364"
    carpeta_bs = "MuestrasTrimmed"
    nombre_p <- sprintf("%sT", id)
    dsin_dehost = TRUE
  }

  nuevo_path_report <- sprintf("%s/Resultados_DRAGEN/%s%s.DRAGEN-report.tsv", patient_dir, patient, de_host_file)

  if(file.exists(nuevo_path_report)) {
    return(message(sprintf("The DRAGEN report was already downloaded for patient %s", patient)))
  }
  command <- sprintf("%s list appsession --project-id %s", bs_path, project_id)
  #print(system(command))

  output <- capture.output(system(command, intern = TRUE))
  df_output <- as.data.frame(do.call(rbind, strsplit(output, " +", perl = TRUE)))

  df_output <- df_output[which(df_output$V3 == sprintf("DRAGENHG38_%s", de_host_file) | df_output$V4 == sprintf("DRAGENHG38_%s", de_host_file)),]
  df_output <- df_output[, c(3,4,5,6,7)]
  #colnames(df_output) <- c("DRAGEN", "SAMPLE", "ID")

  #appsession_id <- readline(prompt = "Insert the Appsession ID: ")
  appsession_id <- unique(df_output[which(df_output[,2] == patient), 4])
  if(length(nchar(appsession_id)) == 0) {
    appsession_id <- unique(df_output[which(df_output[, 3] == patient), 5])
  }
  #df_output <- df_output[which(df_output$V3 == "DRAGENHG38"),]

  if(length(nchar(appsession_id)) == 0) {
    df_output <- df_output[, c(1,2,4)]
    colnames(df_output) <- c(sprintf("DRAGENHG38_%s", de_host_file), "SAMPLE", "ID")
    #appsession_id <- readline(prompt = "Insert the Appsession ID: ")
    appsession_id <- df_output$ID[which(df_output$SAMPLE == patient)]
  }

  if(length(nchar(appsession_id)) == 0) {
    output <- system(command, intern = TRUE)
    output_text <- paste(output, collapse = "\n")
    output_lines <- unlist(strsplit(output_text, "\n"))
    data_lines <- output_lines[-which(grepl("-", output_lines))]
    data_lines <- data.frame(data_lines)
    data_lines <- as.data.frame(lapply(data_lines, function(x) gsub("^\\||\\|$", "", x)))
    data_lines <- data_lines %>%
      separate(data_lines, into = c("BioSampleName", "Id", "ContainerName", "ContainerPosition", "Status"), sep = "\\|")
    appsession_id <- data_lines$Id[which(grepl(patient, data_lines$BioSampleName))]
    appsession_id <- gsub(" ", "", appsession_id)

  }

  if(length(nchar(appsession_id)) == 0) {
    print(system(command, intern = TRUE))
    print(patient)
    appsession_id <- readline(prompt = "Enter the id of your appsession: ")
  }

  if(length(appsession_id)>1) {
    appsession_id <- appsession_id[1]
  }

  dir.create(sprintf("%s/DRAGEN_Reports_%s", patient_dir, de_host_file))
  command <- sprintf("%s download appsession -i %s -o %s",
                     bs_path,
                     appsession_id,
                     sprintf("%s/DRAGEN_Reports_%s", patient_dir, de_host_file))
  system(command)

  folder_report <- list.dirs( sprintf("%s/DRAGEN_Reports_%s", patient_dir, de_host_file), full.names = TRUE, recursive = FALSE)
  file_report <- sprintf("%s/%s.microbe-classification-report.tsv", folder_report, nombre_p)
  nuevo_path_report <- sprintf("%s/Resultados_DRAGEN/%s%s.DRAGEN-report.tsv", patient_dir, patient, de_host_file)
  file.rename(from = file_report, to = nuevo_path_report)

  message(sprintf("The DRAGEN report was downloaded for patient %s", patient))

  return(nuevo_path_report)
}

