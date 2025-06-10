
RunDRAGENsinDeHost <- function(patient_dir, bs_path, reference = "hg38", de_host = "") {

  id <- basename(patient_dir)
  if(id == "trimmed") {
    id <- basename(dirname(patient_dir))
  }

  if(reference == "hg19") {
    de_host_file <- "sin"
  } else{
    de_host_file <- "sinDH_PD"
  }

  project_id = "436799364"
  carpeta_bs = "MuestrasTrimmed"
  nombre_p <- sprintf("%sT", id)


  patient <- basename(patient_dir)
  patient_dir <- paste(patient_dir, "/trimmed", sep="")


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


      if(reference == "hg19") {
        #con hg19 pero CON dehosting de dragen
        command <- sprintf('%s launch application -n "DRAGEN Metagenomics Pipeline" --app-version 3.5.12 -o app-session-name:"DRAGENHG19_%s %s" -l "DRAGENHG19_%s %s" -o project-id:%s -o sample-id:%s -o ht-ref:%s -o basespace-labs-disclaimer:Accepted -o db:minikraken20200312',
                           bs_path,
                           de_host_file,
                           patient,
                           de_host_file,
                           patient,
                           project_id,
                           id_bs,
                           ref)
      } else {
        command <- sprintf('%s launch application -n "DRAGEN Metagenomics Pipeline" --app-version 3.5.12 -o app-session-name:"DRAGENHG38_%s %s" -l "DRAGENHG38_%s %s" -o project-id:%s -o sample-id:%s -o ht-ref:%s -o dehost-checkbox:false -o basespace-labs-disclaimer:Accepted -o db:minikraken20200312',
                           bs_path,
                           de_host_file,
                           patient,
                           de_host_file,
                           patient,
                           project_id,
                           id_bs,
                           ref)
      }

      system(command)

}
