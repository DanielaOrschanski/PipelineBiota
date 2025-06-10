#' @title Run KRAKEN 2
#' @description kvjdfnkvjdf
#' @param patients_dir fskbfjsf
#' @param confidence is a score
#' @param path_kraken2 path to the executable kraken2
#' @examples RunKRAKEN(patients_dir = "/media/16TBDisk/Daniela/Biota/Muestras2daTanda/5/concat5_04/yoconmuestraheritas")
#' @export

RunKRAKEN <- function(patients_dir, confidence= 0.001, de_host) {

  #de_host = "Bowtie"
  libPath <- dirname(system.file(package = "PipelineBiota"))
  pipe_sof <- sprintf("%s/PipelineBiota-Softwares", libPath)
  softwares <- readLines(sprintf("%s/path_to_soft.txt", pipe_sof))
  linea_software <- grep("(?i)KRAKEN2", softwares, ignore.case = TRUE, value = TRUE)

  #setwd(pipe_sof)
  
  if(de_host  == "Bowtie") {
    de_host_file <- "Bo"
  } else if( de_host == "BWA") {
    de_host_file <- "bwa"
  } else if(de_host == "RSubread") {
    de_host_file <- "Rs"
  } else if(de_host == "") {
    de_host_file <- "sin"
  }

  #if(length(nchar(linea_software))==0) {
  #  message("KRAKEN2 will be installed")
  #  path_kraken2 <- downloadKRAKEN2()
  #} else {
  #  path_kraken2 <<- strsplit(linea_software, " ")[[1]][[2]]
  #}
  
  out_kraken2 <- downloadKRAKEN2()
  path_kraken2 <- out_kraken2[[1]]
  db_kraken2 <- out_kraken2[[2]]
  
  list_dirs <- list.dirs(patients_dir, full.names = TRUE, recursive = FALSE)
  list_files <- list.files(patients_dir, full.names = TRUE, recursive = FALSE)

  #Procesamiento de una muestra ---------------------------------------
  #Para identificar que el patients_dir hace referecia a un path de una muestra veo que la unica carpeta que tiene es la de DRAGEN Reports o que no tiene nada.
  if( any(grepl("fastq.gz", list_files)) ) {
    message("Se ejecutará el procesamiento de una única muestra")
    patient <- basename(patients_dir)
    print(patient)
    patient_dir <- paste(patients_dir, "/trimmed", sep="")
    #setwd(patient_dir)
    list_files <- list.files(patient_dir, full.names = TRUE)

    fileR1 <- list_files[grepl(sprintf("%sDH%s_S04_L001_R1_001.fastq.gz", patient, de_host_file), list_files)]
    fileR2 <- list_files[grepl(sprintf("%sDH%s_S04_L001_R2_001.fastq.gz", patient, de_host_file), list_files)]

    if(de_host == "") {
      fileR1 <- list_files[grepl("T_S04_L001_R1_001.fastq.gz", list_files)]
      fileR2 <- list_files[grepl("T_S04_L001_R2_001.fastq.gz", list_files)]
    }

    dir.create(sprintf("%s/Resultados_KRAKEN", patient_dir))

    if(!(file.exists(sprintf("%s/Resultados_KRAKEN/report_%s.sequences", patient_dir, de_host_file)))) {
      # Definir los argumentos del comando en una lista
      args <- c(
        #sprintf("--db %s/minikraken2_v2_8GB_201904_UPDATE", dirname(dirname(path_kraken2))),
        sprintf("--db %s", db_kraken2),
        sprintf("--confidence %s", confidence),
        "--threads 15",
        "--use-names",
        #"--report-zero-counts",
        sprintf("--output %s/Resultados_KRAKEN/output_%s.txt", patient_dir, de_host_file),
        sprintf("--report %s/Resultados_KRAKEN/report_%s.sequences", patient_dir, de_host_file),
        sprintf("--paired %s %s", fileR1, fileR2)
      )

      # Ejecutar el comando
      start_time_kraken <- Sys.time()
      system2(path_kraken2, args = args)
      end_time_kraken <- Sys.time()
      time_kraken <- end_time_kraken - start_time_kraken
      print(paste("Tiempo de ejecución de KRAKEN:", time_kraken))


    } else {
      return(message(sprintf("KRAKEN report has already been generated for this patient: %s", patient)))
    }

    #Procesamiento de varias muestras -------------------------------
  } else {
    message("Se ejecutará el procesamiento de varias muestras")
    for (patient in list_dirs) {
      #patient <- list_dirs[[2]]
      id <- basename(patient)
      patient <- paste(patient, "/trimmed", sep= "")
      print(patient)
      #setwd(patient)
      list_files <- list.files(patient, full.names = TRUE)

      fileR1 <- list_files[grepl(sprintf("%sDH%s_S04_L001_R1_001.fastq", id, de_host_file), list_files)]
      fileR2 <- list_files[grepl(sprintf("%sDH%s_S04_L001_R2_001.fastq", id, de_host_file), list_files)]

      #si se quiere procesar sin dehostear se agarra directamente el fastq salido del trimm
      if(de_host == "") {
        fileR1 <- list_files[grepl("T_S04_L001_R1_001.fastq", list_files)]
        fileR2 <- list_files[grepl("T_S04_L001_R2_001.fastq", list_files)]
      }

      gz <- ifelse(grepl(".gz", fileR1), "--gzip-compressed", "")
      dir.create(sprintf("%s/Resultados_KRAKEN", patient))

      if(!(file.exists(sprintf("%s/Resultados_KRAKEN/report_%s.sequences", patient, de_host_file)))) {
        # Definir los argumentos del comando en una lista
        args <- c(
          #sprintf("--db %s/minikraken2_v2_8GB_201904_UPDATE", dirname(dirname(path_kraken2))),
          sprintf("--db %s", db_kraken2),
          sprintf("--confidence %s", confidence),
          "--threads 15",
          "--use-names",
          #"--report-zero-counts",
          sprintf("--output %s/Resultados_KRAKEN/output_%s.txt", patient, de_host_file),
          sprintf("--report %s/Resultados_KRAKEN/report_%s.sequences", patient, de_host_file),
          sprintf("%s --paired %s %s", gz, fileR1, fileR2)
          #sprintf("--gzip-compressed --paired %s/%s %s/%s", patient, fileR1, patient, fileR2)
        )

        # Ejecutar el comando
        #system2("kraken2", args = args)
        system2(path_kraken2, args = args)
      } else {
        print("Ya se generó el reporte para este paciente")
      }

    }
  }

}


downloadKRAKEN2 <- function() {
  libPath <- dirname(system.file(package = "PipelineBiota"))
  pipe_sof <- sprintf("%s/PipelineBiota-Softwares", libPath)
  setwd(pipe_sof)
  KRAKEN2 <- sprintf("%s/KRAKEN2/kraken2/kraken2", pipe_sof)
  
  if(!file.exists(KRAKEN2)) {
    message("Kraken2 installation will now begin.")
    
    kraken_dir <- dirname(dirname(KRAKEN2))
    dir.create(kraken_dir, showWarnings = TRUE, recursive = TRUE)
    setwd(kraken_dir)
    
    # Descargar Kraken2 desde GitHub
    system("git clone https://github.com/DerrickWood/kraken2.git")
    
    setwd(file.path(kraken_dir, "kraken2"))
    system("make")
    system(sprintf("./install_kraken2.sh %s", file.path(kraken_dir, "kraken2")))
    
    # Verificar que se instaló correctamente
    system(sprintf("%s --version", KRAKEN2))
    
  } else {
    message("Kraken2 was already installed.")
  }
  
  #Base de datos de Kraken2: ---------------------------------------
  
  #kraken_dir <- dirname(dirname(KRAKEN2))
  #DB_Kraken2 <- sprintf("%s/minikraken2_v2_8GB_201904.tgz", kraken_dir)
  DB_Kraken2 <- "/media/16TBDisk/Daniela/Biota/Kraken2-DB-extrac/minikraken2_v2_8GB_201904_UPDATE"
  
  if(!file.exists(DB_Kraken2)) {
    
    message("Kraken2 DataBase will be downloaded now.")
    setwd(kraken_dir)
    
    # Descargar la base comprimida (8GB) desde el servidor oficial de Kraken
      #https://benlangmead.github.io/aws-indexes/k2
    options(timeout = 10000000)
    minikraken_url <- "https://genome-idx.s3.amazonaws.com/kraken/minikraken2_v2_8GB_201904.tgz"
    minikraken_file <- sprintf("%s/minikraken2_v2_8GB_201904.tgz", kraken_dir)
    download.file(minikraken_url, destfile = minikraken_file, mode = "wb")
    system(paste("tar -xvzf", minikraken_file))
    system(paste("tar -xvzf", minikraken_file, "-C", "/media/16TBDisk/Daniela/Biota/Kraken2-DB-extrac"))
    
    
    #Construir base con kraken:
    path_kraken_db <- 
    setwd("/media/16TBDisk/Daniela/Biota/Kraken2-DB")
    system(sprintf("%s/kraken2-build --standard --threads 15 --db /media/16TBDisk/Daniela/Biota/Kraken2-DB", dirname(KRAKEN2)))
    #/home/juan/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/KRAKEN2/kraken2/kraken2-build --standard --threads 15 --db /media/16TBDisk/Daniela/Biota/Kraken2-DB
    
    #system(sprintf("%s/kraken2-build --standard --threads 20 --db %s", dirname(KRAKEN2), dirname(DB_Kraken2)))
    dir_DB_Kraken2 <- sprintf("%s/minikraken2_v2_8GB_201904_UPDATE", dirname(dirname(KRAKEN2)))
    dir_DB_Kraken2 <- "/media/16TBDisk/Daniela/Biota/Kraken2-DB-extrac/minikraken2_v2_8GB_201904_UPDATE"
    DB_Kraken2 <- dir_DB_Kraken2
  } else {
    message("Kraken2 DB has already been downloaded")
  }
  
  return(list(KRAKEN2, DB_Kraken2))
}
