#' @title QC control
#' @import reshape2
#' @import openxlsx
#' @description It generates the following tasks:
#' - Quality control with FastQC
#' - Filter and trimming with Trim Galore
#' - De-hosting
#' - QC report (individual and grupal).
#' @param patients_dir path that points to the directory where all the folders of the samples are stored.
#' @param de_host Indicated the software that will be implemented to perform the de-hosting. It can be "Bowtie", "BWA", "RSubread" or "" (none de-hosting).
#' @param generate_QCReport_Individual if it is set to TRUE a QC report will be generated for each of the samples in patients_dir.
#' @param generate_QCReport_Grupal if it is set to TRUE a QC report will be generated for all the samples in patients_dir.
#' @param FastQC_trimmed if it is set to TRUE the QC report will be generated using QC metrics from trimmed FASTQ files instead of the raw FASTQ files.
#' @return qc_metrics is  dataframe that records qc mean, %Q30, number of sequences, percentage of CG, sequnece length and mapped metrics.
#' @export
#' @example QCcontrol(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host ="Bowtie", generate_QCReport_Individual = FALSE, generate_QCReport_Grupal = FALSE, FastQC_trimmed = TRUE)

QCcontrol <- function(patients_dir, de_host, generate_QCReport_Individual = FALSE, generate_QCReport_Grupal = TRUE, FastQC_trimmed = FALSE ) {

  if(de_host  == "Bowtie") {
    de_host_file <- "Bo"
  } else if( de_host == "BWA") {
    de_host_file <- "bwa"
  } else if(de_host == "RSubread") {
    de_host_file <- "Rs"
  } else if(de_host == "") {
    de_host_file <- "T"
  }  else {
    stop("de_host must be Bowtie, BWA, RSubread or empty string")
  }

  muestras <-  list.dirs(path = path.expand(patients_dir), full.names = TRUE, recursive = FALSE)
  #m <- muestras[3]

  if(any(grepl("_R1_001_fastqc", muestras)) | length(nchar(muestras)) <= 1) {
    #es una sola muestra
    message("Se va a analizar una muestra individual")
    muestras <- patients_dir
  }

  for (m in muestras) {
    print(m)
    id <- basename(m)
    file_list <- list.files(path = m, full.names = TRUE, recursive = FALSE)

      if(length(nchar(file_list[endsWith(file_list, "R1_001_fastqc.html")])) == 0 | length(nchar(file_list[endsWith(file_list, "R2_001_fastqc.html")])) == 0) {
        runFastQC(m, de_host = de_host)
      }

    val1 <- sprintf("%s/trimmed/%s_val_1.fq.gz", m, id)
    val2 <- sprintf("%s/trimmed/%s_val_2.fq.gz", m, id)

    if(file.exists(val1)) {
      trim1 <- sprintf("%s/trimmed/%sT_S04_L001_R1_001.fastq.gz", m, id)
      file.rename(val1, trim1)
    }
    if(file.exists(val2)) {
      trim2 <- sprintf("%s/trimmed/%sT_S04_L001_R1_001.fastq.gz", m, id)
      file.rename(val2, trim2)
    }
    
      file_list_trimmed <- list.files(path = sprintf("%s/trimmed", m), full.names = TRUE, recursive = FALSE)
      #if(length(nchar(file_list_trimmed[endsWith(file_list_trimmed, sprintf("%sT_S04_L001_R1_001.fastq.gz", id))])) == 0 | length(nchar(file_list_trimmed[endsWith(file_list_trimmed, sprintf("%sT_S04_L001_R2_001.fastq.gz", id))])) == 0) {
        m_trimmed <- runTrimgalore(m)
      #} 
      
      file_list_trimmed <- list.files(path = sprintf("%s/trimmed", m), full.names = TRUE, recursive = FALSE)
      m_trimmed <- sprintf("%s/trimmed", m)

      if(FastQC_trimmed == TRUE) {
        if(length(nchar(file_list_trimmed[endsWith(file_list_trimmed, sprintf("%s%s_S04_L001_R1_001.fastqc.html", id, de_host_file))])) == 0 | length(nchar(file_list[endsWith(file_list, sprintf("%s%s_S04_L001_R1_001.fastqc.html", id, de_host_file))])) == 0) {
          runFastQC(m_trimmed, de_host = de_host)
        }
      }

      deHosting(m, de_host = de_host)
      if(generate_QCReport_Individual == TRUE) {
        generateQCreport(patient_dir = m)
      }

  }

  if(any(grepl("_R1_001_fastqc", muestras)) | length(nchar(muestras)) <= 1) {
    #es una sola muestra
    message(sprintf("Quality control ready for sample: %s", patients_dir))
  } else if(generate_QCReport_Grupal == TRUE) {
    #Para analsiis poblacional:
    qc_metrics <- generate_QCReport_plotPBSQ(patients_dir = patients_dir, de_host = de_host , trimmed = FastQC_trimmed)
    return(qc_metrics)
  }

}

#' @title De Hosting
#' @import reshape2
#' @import openxlsx
#' @description kvjdfnkvjdf
#' @param patients_dir fskbfjsf
#' @param de_host. It can be "Bowtie", "BWA" or "RSubread".
#' @return qc_scores is  dataframe that records qc mean and %Q30
#' @export
deHosting <- function(patient_dir, de_host, compatible_DRAGEN = FALSE) {
  #remove host dna ----------------------

  if(de_host  == "Bowtie") {
    de_host_file <- "Bo"
  } else if( de_host == "BWA") {
    de_host_file <- "bwa"
  } else if(de_host == "RSubread") {
    de_host_file <- "Rs"
  }

  id <- basename(patient_dir)
  patient_dir_trim <- paste0(patient_dir, "/trimmed", sep="")
  file_list_trimmed <- list.files(patient_dir_trim, full.names = TRUE, recursive = FALSE)

  log_file <- path.expand(sprintf("%s/%s_%s_Log.txt", patient_dir_trim, id, de_host_file))
  error_file <- path.expand(sprintf("%s/%s_%s_Error.txt", patient_dir_trim, id, de_host_file))


  #if(length(nchar(file_list_trimmed[endsWith(file_list_trimmed, sprintf("DH%s_S04_L001_R1_001.fastq.gz", de_host_file))])) == 0 | length(nchar(file_list_trimmed[endsWith(file_list_trimmed, sprintf("DH%s_S04_L001_R2_001.fastq.gz", de_host_file))])) == 0 | !file.exists(log_file) ) {
  if(length(nchar(file_list_trimmed[endsWith(file_list_trimmed, sprintf("DH%s_S04_L001_R1_001.fastq.gz", de_host_file))])) == 0 | length(nchar(file_list_trimmed[endsWith(file_list_trimmed, sprintf("DH%s_S04_L001_R2_001.fastq.gz", de_host_file))])) == 0 ) {
    message(sprintf("Removing host DNA with %s...", de_host))

    r1_trim <- sprintf("%s/%sT_S04_L001_R1_001.fastq.gz", patient_dir_trim, id)
    r2_trim <- sprintf("%s/%sT_S04_L001_R2_001.fastq.gz", patient_dir_trim, id)
    out_r1_hrmv <- sprintf("%s/%sDH%s_S04_L001_R1_001.fastq", patient_dir_trim, id, de_host_file)
    out_r2_hrmv <- sprintf("%s/%sDH%s_S04_L001_R2_001.fastq", patient_dir_trim, id, de_host_file)

    out <- downloadHG38()
    hg38Fasta <- out[[2]]
    
    # DEHOSTING CON BOWTIE
    if(de_host == "Bowtie") {

      #system2("/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/Bowtie2/bowtie2-2.4.5-linux-x86_64/bowtie2",
      #      args = c(
      #        "-x", "/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/HG38/index/Homo_sapiens.GRCh38.dna_sm.primary_assembly",
      #        "-1", r1_trim,
      #        "-2", r2_trim,
      #        "-q", "--phred33", "--sensitive", "--end-to-end", "-p", "10",
      #       sprintf("|/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/Samtools/samtools-1.16.1/samtools view -bh -h - |/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/Samtools/samtools-1.16.1/samtools view -bh -h -f 12 -F 256 - |/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/Samtools/samtools-1.16.1/samtools sort -n - | bedtools bamtofastq -i - -fq %s -fq2 %s",
      #        out_r1_hrmv, out_r2_hrmv)
      #      ),
      #      stdout = log_file,
      #      stderr = error_file
      #     )

      #Otra opcion:
      # Paso 1: Ejecutar Bowtie2 y capturar las métricas
      B <- downloadBowtie2()
      Bowtie2 <- B[[1]]
      hg38_indexBowtie2 <- B[[2]]
      
      bowtie_output_bam <- sprintf("%s/%sDH%s.bam", patient_dir_trim, id, de_host_file)
      #setwd(patient_dir_trim)
      start_time_bowtie <- Sys.time()
      system2(
        Bowtie2,
        #"bowtie2",
        args = c(
          "-x", hg38_indexBowtie2,
          "-1", r1_trim,
          "-2", r2_trim,
          "-q", "--phred33", "--sensitive", "--end-to-end", "-p", "10",
          "-S", bowtie_output_bam 
        ),
        stdout = log_file,
        stderr = error_file
      )
      end_time_bowtie <- Sys.time()
      bowtie_duration <- end_time_bowtie - start_time_bowtie

      #bowtie2 -x /home/juan/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/HG38/index/HG38_index -1 /media/4tb2/Daniela/Biota/Muestras/37/trimmed/37T_S04_L001_R1_001.fastq.gz -2 /media/4tb2/Daniela/Biota/Muestras/37/trimmed/37T_S04_L001_R2_001.fastq.gz -q --phred33 --sensitive --end-to-end -p 10 -S /media/4tb2/Daniela/Biota/Muestras/37/trimmed/37DHBo.bam
      
      # Paso 2: Procesar la salida de Bowtie2 con Samtools y Bedtools
      start_time_samtools <- Sys.time()
      Samtools <- downloadSamtools()
      system(sprintf("%s --version", Samtools))
      system(sprintf(
        "%s view -bh -h %s | %s view -bh -h -f 12 -F 256 - | %s sort -n - | bedtools bamtofastq -i - -fq %s -fq2 %s",
        Samtools, bowtie_output_bam, Samtools, Samtools, out_r1_hrmv, out_r2_hrmv
      ))
      end_time_samtools <- Sys.time()
      samtools_duration <- end_time_samtools - start_time_samtools
      time_bo <- samtools_duration + bowtie_duration
      print(paste("Tiempo de ejecución de Bowtie con Samtools y Bedtools:", time_bo))

    # DEHOSTING CON BWA
    } else if (de_host == "BWA") {
      
      
      #hg38Fasta <- "/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/HG38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
      
      BWA_out <- downloadBWA()
      BWA <- BWA_out[[1]]
      indexBWA <- BWA_out[[2]]
      #indexBWA <- "/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/HG38/indexBWA"

      if(!file.exists(indexBWA) | length(list.files(indexBWA)) == 0) {
        # Ejecutar BWA index
        system2(
          BWA,
          args = c(
            "index",
            "-p", file.path(indexBWA, "Homo_sapiens.GRCh38.dna_sm.primary_assembly"), # Prefijo del índice
            hg38Fasta))
      }

      start_time_bwa <- Sys.time()
      system2(
        BWA,
        args = c(
          "mem",                                   # Subcomando de BWA
          "-t", "10",                              # Número de hilos
          file.path(indexBWA, "Homo_sapiens.GRCh38.dna_sm.primary_assembly"),  # Índice de referencia generado
          r1_trim,                                 # Archivo de lecturas forward
          r2_trim,                                 # Archivo de lecturas reverse
          sprintf(
            "| /home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/Samtools/samtools-1.16.1/samtools view -bh -h - | /home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/Samtools/samtools-1.16.1/samtools view -bh -h -f 12 -F 256 - | /home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/Samtools/samtools-1.16.1/samtools sort -n - | bedtools bamtofastq -i - -fq %s -fq2 %s",
            out_r1_hrmv, out_r2_hrmv
          )  ))
      end_time_bwa <- Sys.time()
      time_bwa <- end_time_bwa - start_time_bwa
      print(paste("Tiempo de ejecución de BWA:", time_bwa))


      #DEHOSTING CON RSUBREAD
    } else if (de_host == "RSubread") {

      indexRsubread <- "/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/HG38/indexRsubread/Homo_sapiens.GRCh38.dna_sm.primary_assembly"
      library(Rsubread)

      if(length(list.files("/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/HG38/indexRsubread")) == 0) {
        # Ejecutar Rsubread index
        dir.create("/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/HG38/indexRsubread")
        buildindex(
          basename = "/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/HG38/indexRsubread/Homo_sapiens.GRCh38.dna_sm.primary_assembly",
          reference = hg38Fasta)
      }

      # Alineamiento con Rsubread
      start_time_rs <- Sys.time()
      out_bam <- paste0(patient_dir_trim, "/", id ,"_R1R2_Rsubread.bam")
      align(
        index = indexRsubread,
        readfile1 = r1_trim,
        readfile2 = r2_trim,
        input_format = "FASTQ",
        output_file = out_bam,
        phredOffset = 33, # Para phred33
        nthreads = 10,
        PE_orientation = "fr", # orientación por defecto "forward-reverse"
        type = "dna", # Alineación al genoma
        maxMismatches = 3, # Número permitido de desajustes
        useAnnotation = FALSE
      )

      # Filtrado y conversión a FASTQ
      system(sprintf("samtools view -bh -h %s | samtools view -bh -h -f 12 -F 256 - | samtools sort -n - | bedtools bamtofastq -i - -fq %s -fq2 %s",
                     out_bam, out_r1_hrmv, out_r2_hrmv))
      end_time_rs <- Sys.time()
      time_rs <- end_time_rs - start_time_rs
      print(paste("Tiempo de ejecución de Rsubread:", time_rs))

      file.remove(out_bam)
      file.remove(paste0(patient_dir_trim, "/", id ,"_R1R2_Rsubread.bam.indel.vcf"))
    }


    system2("gzip", args = c(out_r1_hrmv, out_r2_hrmv))

    out_r1_hrmv <- sprintf("%s/%sDH%s_S04_L001_R1_001.fastq.gz", patient_dir_trim, id, de_host_file)
    out_r2_hrmv <- sprintf("%s/%sDH%s_S04_L001_R2_001.fastq.gz", patient_dir_trim, id, de_host_file)

    # Modificar los headers para que sean compatibles para subirse a BS Illumina: ------
    if(compatible_DRAGEN == "TRUE") {
      cmd <- paste(
        "zcat", out_r1_hrmv,
        "| sed 's/\\/1$/ 1:N:0:CGAGGCTG+CTCCTTAC/'",
        "| gzip >",
        paste0(out_r1_hrmv, "_fixed")
      )
      system(cmd)
  
      cmd <- paste(
        "zcat", out_r2_hrmv,
        "| sed 's/\\/2$/ 2:N:0:CGAGGCTG+CTCCTTAC/'",
        "| gzip >",
        paste0(out_r2_hrmv, "_fixed")
      )
      system(cmd)
  
      file.remove(out_r1_hrmv)
      file.remove(out_r2_hrmv)
  
      old_r1_hrmv <- sprintf("%s/%sDH%s_S04_L001_R1_001.fastq.gz_fixed", patient_dir_trim, id, de_host_file)
      old_r2_hrmv <- sprintf("%s/%sDH%s_S04_L001_R2_001.fastq.gz_fixed", patient_dir_trim, id, de_host_file)
  
  
      file.rename(old_r1_hrmv, out_r1_hrmv)
      file.rename(old_r2_hrmv, out_r2_hrmv)

    }
    
  } else {
    message("De-hosting already done.")
  }
}

#' @title Run FastQC
#' @description Executes the FastQC for R1 and R2.
#' @param patent_dir Path of the directory that contains R1 and R2 fastq files. It can be either the "trimmed" folder or the original folder.
#' @import stringr
#' @import viridis
#' @import reshape2
#' @export
runFastQC <- function(patient_dir, de_host) {

  if(de_host  == "Bowtie") {
    de_host_file <- "DHBo"
  } else if( de_host == "BWA") {
    de_host_file <- "DHbwa"
  } else if(de_host == "RSubread") {
    de_host_file <- "DHRs"
  } else if(de_host == "") {
    de_host_file <- "T"
  }  else {
    stop("de_host must be Bowtie, BWA, RSubread or empty string")
  }

  #FastQC <- downloadFastQC()
  #FastQC <- "/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/FastQC/FastQC/fastqc"
  FastQC <- downloadFastQC()
  file_list <- list.files(patient_dir)

  #Para que se pueda poner como entrada la carpeta de los trimmeados o la carpeta original:

  # TRIMEADO -------------------------------------------
  if  (startsWith(basename(patient_dir), "trimmed")) {
    #gzip <- ifelse(length(nchar(file_list[endsWith(file_list, sprintf("%s_S04_L001_R1_001.fastq.gz", de_host_file))])) == 0, "", ".gz")
    #fileR1 <- paste0(patient_dir, "/", file_list[endsWith(file_list, sprintf("%s_S04_L001_R1_001.fastq%s", de_host_file, gzip))], sep="")
    #fileR2 <- paste0(patient_dir, "/", file_list[endsWith(file_list, sprintf("%s_S04_L001_R2_001.fastq%s", de_host_file, gzip))], sep="")

    id <- basename(dirname(patient_dir))
    gzip <- ifelse(length(nchar(sprintf("%s/%s%s_S04_L001_R1_001.fastq.gz", patient_dir, id, de_host_file))) == 0, "", ".gz")
    fileR1 <- paste0(patient_dir, "/", sprintf("%s%s_S04_L001_R1_001.fastq%s", id, de_host_file, gzip), sep="")
    fileR2 <- paste0(patient_dir, "/", sprintf("%s%s_S04_L001_R2_001.fastq%s", id, de_host_file, gzip), sep="")

    if ((length(nchar(fileR1)) == 0) | (length(nchar(fileR2)) == 0)) {
      stop("There are no fastq files in this directory")
    }

    #R1
    if (!(length(nchar(file_list[endsWith(file_list, sprintf("%s_S04_L001_R1_001_fastqc", de_host_file))])) == 0)) {
      message("The FastQC for this sample R1 has already been done.")
      #return(paste0(patient_dir, "/", file_list[endsWith(file_list, "T_S04_L001_R1_001.fastqc")], sep=""))
    } else {
      system2(FastQC, fileR1)
      file_fastqc_zip <- paste0(patient_dir, "/", sprintf("%s%s_S04_L001_R1_001_fastqc.zip", id, de_host_file), sep="")
      unzip(file_fastqc_zip, exdir = patient_dir)

    }
    #R2
    if (!(length(nchar(file_list[endsWith(file_list, sprintf("%s_S04_L001_R2_001_fastqc", de_host_file))])) == 0)) {
      message("The FastQC for this sample R2 has already been done.")
      #return(paste0(patient_dir, "/", file_list[endsWith(file_list, "T_S04_L001_R1_001.fastqc")], sep=""))
    } else {
      system2(FastQC, fileR2)
      file_fastqc_zip <- paste0(patient_dir, "/", sprintf("%s%s_S04_L001_R2_001_fastqc.zip", id, de_host_file), sep="")
      unzip(file_fastqc_zip, exdir = patient_dir)

    }



  } else { # Cuando es con los fastq normales, sin trimmear: ------------------------

    gzip <- ifelse(length(nchar(file_list[endsWith(file_list, "R1_001.fastq.gz")])) == 0, "", ".gz")
    fileR1 <- paste0(patient_dir, "/", file_list[endsWith(file_list, sprintf("R1_001.fastq%s", gzip)) & !(grepl("fastqc", file_list))], sep="")
    fileR2 <- paste0(patient_dir, "/", file_list[endsWith(file_list, sprintf("R2_001.fastq%s", gzip)) & !(grepl("fastqc", file_list))], sep="")

    if ((length(nchar(fileR1)) == 0) | (length(nchar(fileR2)) == 0)) {
      stop("There are no fastq files in this directory")
    }

    #R1
    if ((length(nchar(file_list[endsWith(file_list, "R1_001_fastqc")])) != 0) & (length(nchar(file_list[endsWith(file_list, "R2_001_fastqc")])) != 0)) {
      message("The FastQC for this sample R1 has already been done.")
      #return(paste0(patient_dir, "/", file_list[endsWith(file_list, "R1_fastqc")], sep=""))
    } else {
      system2(FastQC, fileR1)
      file_fastqc_zip <- paste0(patient_dir, "/", file_list[endsWith(file_list, "R1_001_fastqc.zip")], sep="")
      unzip(file_fastqc_zip, exdir = patient_dir)
    }
    #R2
    if ((length(nchar(file_list[endsWith(file_list, "R2_001_fastqc")])) != 0) & (length(nchar(file_list[endsWith(file_list, "R2_001_fastqc")])) != 0)) {
      message("The FastQC for this sample R2 has already been done.")
      #return(paste0(patient_dir, "/", file_list[endsWith(file_list, "R1_fastqc")], sep=""))
    } else {
      system2(FastQC, fileR2)
      file_fastqc_zip <- paste0(patient_dir, "/", file_list[endsWith(file_list, "R2_001_fastqc.zip")], sep="")
      unzip(file_fastqc_zip, exdir = patient_dir)
    }

  }

  message("FastQC's analysis has finished!")

}

#' @title plotFastQC_PBSQ
#' @description generates the principal plot ("Per Base Sequence Quality") which compares the quality for all the samples.
#' @param patients_dir Path of the directory that contains one folder with fastq files (R1 and R2) of each patient.
#' @param trimmed it is set to TRUE when the fastQC files that you want to plot came from a trimmed file.
#' @param R will indicate if the curve of the plot will be constructed by the mean of R1 (R= "R1"), R2 (R= "R2") or by both (R= "R1R2").
#' @export
#' @import readr
#' @import stringr
#' @import ggplot2
#' @import reshape2
#' @import png
generate_QCReport_plotPBSQ <- function(patients_dir, de_host, trimmed = FALSE, R= "R1R2") {

  #Separación de los modulos del fastqc para generar graficos
  dir_list <- list.dirs(path = patients_dir, full.names = FALSE, recursive = FALSE)
  cant_patients <- length(dir_list)

  list_modulo_PBSQ <- list() #Va a almacenar los dataframes para el grafico Per Base Sequence Quality
  Table_Basic_Stats <- data.frame("ID" = c(), "TotalSequencesMean" = c() )

  #Recorro cada paciente para guardar las tablas para el grafico y sacar el promedio de las medias.
  for (p in 1:cant_patients) {

    #Toma los fastqc y se asegura de que estén unzippeados
    if (trimmed == FALSE) {
      file_list <- list.files(sprintf("%s/%s", patients_dir, dir_list[p]))

      file_fastqc_zip <- paste0(patients_dir, "/", dir_list[p],"/", file_list[endsWith(file_list, "R1_001_fastqc.zip")], sep="")
      unzip(file_fastqc_zip, exdir = sprintf("%s/%s", patients_dir, dir_list[p]))
      file_fastqc_zip <- paste0(patients_dir, "/", dir_list[p],"/", file_list[endsWith(file_list, "R2_001_fastqc.zip")], sep="")
      unzip(file_fastqc_zip, exdir = sprintf("%s/%s", patients_dir, dir_list[p]))

      file_list <- list.files(sprintf("%s/%s", patients_dir, dir_list[p]))

      if (!(length(nchar(file_list[endsWith(file_list, "R1_fastqc")])) == 0 )) {
        dir_fastqc_R1 <- paste0(patients_dir, "/", dir_list[p],"/", file_list[endsWith(file_list, "R1_fastqc")], sep="")
      } else if (!(length(nchar(file_list[endsWith(file_list, "R1_001_fastqc")])) == 0 )) {
        dir_fastqc_R1 <- paste0(patients_dir, "/", dir_list[p],"/", file_list[endsWith(file_list, "R1_001_fastqc")], sep="")
      }

      if(!(length(nchar(file_list[endsWith(file_list, "R2_fastqc")])) == 0 )) {
        dir_fastqc_R2 <- paste0(patients_dir, "/", dir_list[p],"/", file_list[endsWith(file_list, "R2_fastqc")], sep="")
      } else if(!(length(nchar(file_list[endsWith(file_list, "R2_001_fastqc")])) == 0 )) {
        dir_fastqc_R2 <- paste0(patients_dir, "/", dir_list[p],"/", file_list[endsWith(file_list, "R2_001_fastqc")], sep="")
      }

      dir_fastq_R1 <- paste0(patients_dir, "/", dir_list[p],"/", file_list[endsWith(file_list, "R1.fastq.gz")], sep="")
      dir_fastq_R2 <- paste0(patients_dir, "/", dir_list[p],"/", file_list[endsWith(file_list, "R2.fastq.gz")], sep="")


    } else {
      file_list <- list.files(sprintf("%s/%s/trimmed", patients_dir, dir_list[p]))

      file_fastqc_zip <- paste0(patients_dir, "/", dir_list[p],"/trimmed/", file_list[endsWith(file_list, "val_1_fastqc.zip")], sep="")
      unzip(file_fastqc_zip, exdir = sprintf("%s/%s/trimmed", patients_dir, dir_list[p]))
      file_fastqc_zip <- paste0(patients_dir, "/", dir_list[p],"/trimmed/", file_list[endsWith(file_list, "val_2_fastqc.zip")], sep="")
      unzip(file_fastqc_zip, exdir = sprintf("%s/%s/trimmed", patients_dir, dir_list[p]))

      file_list <- list.files(sprintf("%s/%s/trimmed", patients_dir, dir_list[p]))
      dir_fastqc_R1 <- paste0(patients_dir, "/", dir_list[p],"/trimmed/", file_list[endsWith(file_list, "val_1_fastqc")], sep="")
      dir_fastqc_R2 <- paste0(patients_dir, "/", dir_list[p],"/trimmed/", file_list[endsWith(file_list, "val_2_fastqc")], sep="")

      dir_fastq_R1 <- paste0(patients_dir, "/", dir_list[p],"/trimmed/", file_list[endsWith(file_list, "val_1_fq.gz")], sep="")
      dir_fastq_R2 <- paste0(patients_dir, "/", dir_list[p],"/trimmed/", file_list[endsWith(file_list, "val_2_fq.gz")], sep="")
    }


    if( length(nchar(dir_fastqc_R1)) == 0 | length(nchar(dir_fastqc_R2)) == 0 ) {
      stop("There is no fastQC file in this folder. Try trimmed == FALSE.")
    }
    print(basename(dir_fastqc_R1))

    report_R1 <- read_file(paste0(dir_fastqc_R1, "/fastqc_data.txt", sep=""))
    library(stringr)
    module_R1 <- str_split(report_R1, ">>")


    print(basename(dir_fastqc_R2))
    report_R2 <- read_file(paste0(dir_fastqc_R2, "/fastqc_data.txt", sep=""))
    module_R2 <- str_split(report_R2, ">>")

    #Basic stats:
    Basic_stats_R1 <- create_FQCdata(module_R1[[1]][2])
    Basic_stats_R2 <- create_FQCdata(module_R2[[1]][2])

    id <- basename(dir_list[p])
    seqsR1 <- as.numeric(Basic_stats_R1$Value[which(Basic_stats_R1$Measure == "Total Sequences")])
    seqsR2 <- as.numeric(Basic_stats_R2$Value[which(Basic_stats_R2$Measure == "Total Sequences")])
    promedio_totalseqs <- (seqsR1 + seqsR2) /2

    basesR1 <- Basic_stats_R1$Value[which(Basic_stats_R1$Measure == "Total Bases")]
    basesR2 <- Basic_stats_R2$Value[which(Basic_stats_R2$Measure == "Total Bases")]
    nbasesR1 <- as.numeric(strsplit(basesR1, split = " ")[[1]][1])
    nbasesR2 <- as.numeric(strsplit(basesR2, split = " ")[[1]][1])
    promedio_bases <- (nbasesR1 + nbasesR2)/2

    lenR1 <- as.numeric(Basic_stats_R1$Value[which(Basic_stats_R1$Measure == "Sequence length")])
    lenR2 <- as.numeric(Basic_stats_R2$Value[which(Basic_stats_R2$Measure == "Sequence length")])
    promedio_length <- (lenR1 + lenR2) /2

    GCR1 <- as.numeric(Basic_stats_R1$Value[which(Basic_stats_R1$Measure == "%GC")])
    GCR2 <- as.numeric(Basic_stats_R2$Value[which(Basic_stats_R2$Measure == "%GC")])
    promedio_GC <- (GCR1 + GCR2) /2

    basic_stats <- data.frame("Sample" = id,
                              "TotalSequences_R1" = seqsR1,
                              "TotalSequences_R2" = seqsR2,
                              "TotalSequences_Mean" = promedio_totalseqs,

                              "TotalBases_R1" = basesR1,
                              "TotalBases_R2" = basesR2,
                              "TotalBasesBases_Mean" = promedio_bases,

                              "SequenceLength_R1" = lenR1,
                              "SequenceLength_R2" = lenR2,
                              "SequenceLength_Mean" = promedio_length,

                              "GCPercentage_R1" = GCR1,
                              "GCPercentage_R2" = GCR2,
                              "GCPercentage_Mean" = promedio_GC
                              )

    Table_Basic_Stats <- rbind(Table_Basic_Stats, basic_stats)

    #Plot principal: -----------------------------------------------------------------
    Per_base_sequence_quality_R1 <- create_FQCdata(module_R1[[1]][4])
    order_fact <- Per_base_sequence_quality_R1$Base
    Per_base_sequence_quality_R1$Base <- factor(Per_base_sequence_quality_R1$Base, levels = order_fact)

    Per_base_sequence_quality_R2 <- create_FQCdata(module_R2[[1]][4])
    order_fact <- Per_base_sequence_quality_R2$Base
    Per_base_sequence_quality_R2$Base <- factor(Per_base_sequence_quality_R2$Base, levels = order_fact)

    #promedio la mean de R1 y R2 para poder tener una sola línea por cada paciente
    Per_base_sequence_quality <- Per_base_sequence_quality_R1
    #Per_base_sequence_quality$Mean <- ((Per_base_sequence_quality_R1$Mean + Per_base_sequence_quality_R2$Mean) / 2)

    if (R == "R2") {
      Per_base_sequence_quality$Mean <- Per_base_sequence_quality_R2$Mean
    } else if ( R == "both" | R == "R1R2") {
      Per_base_sequence_quality$Mean <- ((Per_base_sequence_quality_R1$Mean + Per_base_sequence_quality_R2$Mean) / 2)
    } else if (R == "R1") {
      Per_base_sequence_quality$Mean <- Per_base_sequence_quality_R1$Mean
    }

    Per_base_sequence_quality$Median <- ((Per_base_sequence_quality_R1$Median + Per_base_sequence_quality_R2$Median) / 2)
    Per_base_sequence_quality$`10th Percentile` <- ((Per_base_sequence_quality_R1$`10th Percentile`  + Per_base_sequence_quality_R2$`10th Percentile` ) / 2)
    Per_base_sequence_quality$`90th Percentile` <- ((Per_base_sequence_quality_R1$`90th Percentile`  + Per_base_sequence_quality_R2$`90th Percentile` ) / 2)
    Per_base_sequence_quality$`Lower Quartile` <- ((Per_base_sequence_quality_R1$`Lower Quartile`  + Per_base_sequence_quality_R2$`Lower Quartile` ) / 2)
    Per_base_sequence_quality$`Upper Quartile` <- ((Per_base_sequence_quality_R1$`Upper Quartile`  + Per_base_sequence_quality_R2$`Upper Quartile` ) / 2)

    list_modulo_PBSQ[[p]] <- Per_base_sequence_quality
    message(sprintf("The module of patient number %s for PBSQ graph has been loaded", p))

  }

  #ARMO UN REPORT DE QC COMPLETO CON %Q30 Y METRICAS DE MAPEO:

  trim <- ifelse(trimmed == TRUE, "trimmed", "raw")
  file_qc_scores <- paste(patients_dir, "/QC_scores_", trim,".xlsx", sep ="")
  #Genero ql qc_scores si no existe. PEro si ya existe evita volvr a hacerlo porque demora mucho tiempo:
  if(!file.exists(file_qc_scores)) {
    qc_scores <- calculate_qc_scores(patients_dir = patients_dir , trimmed = trimmed)
  } else {
    qc_scores <- as.data.frame(read_excel(file_qc_scores))
    qc_scores <- qc_scores[,-1]
  }
  Complete_QC_metrics <- merge(Table_Basic_Stats, qc_scores, by = "Sample")

  mapeo <- calculate_mapeo(patients_dir = patients_dir, de_host = de_host)
  Complete_QC_metrics <- merge(Complete_QC_metrics, mapeo, by = "Sample")

  write.xlsx(Complete_QC_metrics, file = sprintf("%s/QC_MetricsReport-%sp.xlsx", patients_dir, length(list.dirs(patients_dir, recursive = FALSE))) )


  #Ploteo una linea por paciente superpuestas en un mismo grafico -----------------
  #Base del plot
  plot <- ggplot(Per_base_sequence_quality, aes(x = Base, y = Mean)) +
    geom_boxplot(width = 0.5, fill = "black", color = "black") +
    geom_point(aes(y = Median), color = "yellow", size = 1, position = position_dodge(width = 0.75)) +
    geom_errorbar(
      aes(ymin = `Lower Quartile`, ymax = `Upper Quartile`),
      width = 0.2,
      position = position_dodge(width = 0.75),
      color = "black"
    ) +
    geom_linerange(
      aes(ymin = `10th Percentile`, ymax = `90th Percentile`),
      position = position_dodge(width = 0.75),
      color = "black"
    ) +
    geom_line(aes(y = Mean, group = 1), position = position_dodge(width = 0.75), color = "black", size = 1) +
    geom_hline(yintercept = 28, linetype = "dashed", color = "green", size = 0.5) +  # Add horizontal lines
    geom_hline(yintercept = 20, linetype = "dashed", color = "red", size = 0.5) +  # Add horizontal lines
    labs(title = "Per Base Sequence Quality",
         x = "Position in read (Base)",
         y = "Values") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  # Rotate x-axis labels
    scale_y_continuous(breaks = seq(0, 40, by = 2))+
    coord_cartesian(ylim = c(0, 41))  # Set y-axis limits

  #Genera un dataframe para poder agregar una linea por cada paciente:
  combined_dataframe <- do.call(cbind, lapply(seq_along(list_modulo_PBSQ), function(i) {
    df <- list_modulo_PBSQ[[i]]
    names(df) <- paste(names(df), "_", dir_list[i], sep = "_")
    df
  }))

  selected_columns <- combined_dataframe[, grepl("^Mean", names(combined_dataframe))]
  selected_columns$Base <-combined_dataframe[,1]
  class(selected_columns)
  dt <- as.data.table(selected_columns)
  df.long <- melt(dt, id.vars = "Base")
  #df.long <- melt(selected_columns,id.vars="Base")

  # Crea el gráfico completo con una linea por paciente:
  plot2 <- plot +
    geom_line(data = df.long, aes(x = Base, y = value, color = variable, group = variable), size = 1) +
    labs(title = "Per Base Sequence Quality - 1 curva por paciente", x = "Base", y = "Quality Value")

  indicador <- ifelse(trimmed == TRUE, "_trimmeado", "_SINtrimmeado")
  png(filename = sprintf("%s/plotFastQC_PBSQ%s.png", patients_dir, indicador), width = 800, height = 600)  # Ajusta el tamaño según tus necesidades
  print(plot2)
  dev.off()

  message(sprintf("The plot from FastQC was saved in %s/plotFastQC_PerBaseSequenceQuality.png", patients_dir))
  return(Complete_QC_metrics)

}

#' @title create FastQC data
#' @description Prepare the data for plotting
#' @param dat data from module of FastQC output
create_FQCdata <- function(dat) {
  DBname <- read_lines(dat)
  len <- length(DBname)
  DBname <- DBname[2:len]
  DBname <- read.table(text = DBname, sep = "\t")
  colnames(DBname) <- strsplit(read_lines(dat)[2],"\t", fixed=TRUE)[[1]]
  colnames(DBname)[1] <- str_replace(colnames(DBname)[1], "#", "")
  return(DBname)
}


#' @title Calculate QC scores from sequencing
#' @import qckitfastq
#' @description Prepare the data for plotting
#' @param patients_dir data from module of FastQC output
#' @param trimmed it can be T or F
calculate_qc_scores <- function(patients_dir, trimmed) {
  library(qckitfastq)
  scores_qc <- data.frame("Sample"= c(), "Q_Mean_R1" = c(), "%_>=Q30_R1"= c(),
                          "Q_Mean_R2" = c(), "%_>=Q30_R2"= c())
  i=1
  dir_list <- list.dirs(path = patients_dir, full.names = TRUE, recursive = FALSE)
  cant_patients <- length(dir_list)

  for (p in dir_list) {
    p <- dir_list[[i]]
    print(p)
    if  (trimmed == TRUE) {
      file_list <- list.files(sprintf("%s/trimmed/", p))
      gzip <- ifelse(length(nchar(file_list[endsWith(file_list, "T_S04_L001_R1_001.fastq.gz")])) == 0, "", ".gz")
      fileR1 <- paste0(p, "/trimmed/", file_list[endsWith(file_list, sprintf("T_S04_L001_R1_001.fastq%s", gzip))], sep="")
      fileR2 <- paste0(p, "/trimmed/", file_list[endsWith(file_list, sprintf("T_S04_L001_R2_001.fastq.%s", gzip))], sep="")
    } else {
      file_list <- list.files(p)
      gzip <- ifelse(length(nchar(file_list[endsWith(file_list, "R1_001.fastq.gz")])) == 0, "", ".gz")
      fileR1 <- paste0(p, "/", file_list[endsWith(file_list, sprintf("R1_001.fastq%s", gzip))], sep="")
      fileR2 <- paste0(p, "/", file_list[endsWith(file_list, sprintf("R2_001.fastq%s", gzip))], sep="")
    }

    if ((length(nchar(fileR1)) == 0) | (length(nchar(fileR2)) == 0)) {
      stop("There are no fastq files in this directory")
    }

    QC <- qual_score_per_read(fileR1)
    mean <- mean(QC$mu_per_read)
    Q30 <- (sum(QC$mu_per_read>=30)/length(QC$mu_per_read))*100

    scores_qc[i, "Sample"] <- basename(p)
    scores_qc[i, "Q_Mean_R1"] <- mean
    scores_qc[i, "%_>=Q30_R1"] <- Q30

    QC <- qual_score_per_read(fileR2)
    mean <- mean(QC$mu_per_read)
    Q30 <- (sum(QC$mu_per_read>=30)/length(QC$mu_per_read))*100

    scores_qc[i, "Q_Mean_R2"] <- mean
    scores_qc[i, "%_>=Q30_R2"] <- Q30

    i= i+1
  }
  #saveRDS(scores_qc, file= sprintf("%s/scores_QC.rds", patients_dir))

  trim <- ifelse(trimmed == TRUE, "trimmed", "raw")
  write.xlsx(scores_qc, file = paste(patients_dir, "/QC_scores_", trim,".xlsx", sep =""), rowNames= TRUE)

  return(scores_qc)

}


#' @title RunTrimgalore
#' @description Corre la funcion TrimGalore para eliminar los adapters y las lecturas (y bases) de baja calidad.
#' @param patient_dir Path donde se encuentra el archivo R1 de formato fasta o fastq.
#' @param fastQC_after_trim if it is set to TRUE it will run tha FastQC for the trimmed files.
#' @param trim_quality is the minimum value of the quality of each base within the sequence that will pass the filter.
#' @return path of the trimmed folder which contained the trimmed files and was created inside the patient folder.
#' @export
runTrimgalore <- function(patient_dir, trim_quality = 20) {
  # Chequeamos que esta descargado TrimGalore. En caso de no estarlo, lo descarga
  TrimGalore <- downloadTrimGalore()

  patient_id <- basename(patient_dir)
  file_list <- list.files(patient_dir)
  gzip <- ifelse(length(nchar(file_list[endsWith(file_list, "R1_001.fastq.gz")])) == 0, "", ".gz")
  fileR1 <- paste0(patient_dir, "/", file_list[endsWith(file_list, sprintf("R1_001.fastq%s", gzip))], sep="")
  fileR2 <- paste0(patient_dir, "/", file_list[endsWith(file_list, sprintf("R2_001.fastq%s", gzip))], sep="")

  if ((length(nchar(file_list[endsWith(file_list, "R1_001.fastq.gz")])) == 0) | (length(nchar(file_list[endsWith(file_list, "R2_001.fastq.gz")])) == 0)) {
    stop("There are no fastq files in this directory")
  }


  # Se fija si los archivos de entrada son de formato .gz
  gziped <- ifelse(stringr::str_detect(fileR1,".gz"),"--gzip","--dont_gzip")


  #Se fija si está hecho el trimmeado antes:
  trim_dir <- list.files(paste0(patient_dir, "/trimmed"), full.names = T)
  if (!(length(nchar(trim_dir[endsWith(trim_dir, "T_S04_L001_R1_001.fastq.gz")])) == 0) & !(length(nchar(trim_dir[endsWith(trim_dir, "T_S04_L001_R2_001.fastq.gz")])) == 0)) {
    #ifelse(gziped == "--gzip",
    #       file_trim <- paste0(patient_dir, "/", file_list[endsWith(file_list, "val_1.fq.gz")], sep=""),
    #       file_trim <- paste0(patient_dir, "/", file_list[endsWith(file_list, "val_1.fq")], sep=""))

    message("You have already trimmed this sample")
    return(sprintf("%s/trimmed", patient_dir))
  }

  dir.create(sprintf("%s/trimmed", patient_dir))

  outdir <- sprintf("%s/trimmed", patient_dir)

  # Trimeado por consola. Guarda el tiempo que tardo en ejecutarse
  t1 <- system.time(system2(command = TrimGalore,
                            args =c("--paired",
                                    gziped,
                                    sprintf("-q %s", trim_quality),
                                    paste0("--output_dir ", outdir),
                                    fileR1,
                                    fileR2,
                                    paste0("--basename ", basename(patient_dir)))))
  #stderr = file.path(outdir, "ErrLog.txt")))

  # Calculo el tamaño de los archivos de entrada sin trimmear
  of.size1 <- file.info(fileR1)$size
  of.size2 <- file.info(fileR2)$size

  # Si el formato era .gz, el output se llamara tambien como .gz, porque TrimGalore asi lo genera
  if (gziped == "--gzip") {
    ofile <- paste0(outdir, "/", basename(patient_dir),"_val_1.fq.gz")
  } else {
    ofile <- paste0(outdir, "/", basename(patient_dir),"_val_1.fq")
  }

  # Tamanos del archivo de salida
  ot.size1 <- file.info(ofile)$size
  ot.size2 <- file.info(stringr::str_replace_all(ofile,"val_1.","val_2."))$size

  # No entiendo porque hace de return esta asignacion de variables
  #return(c(ifile = fileR1, tfile = ofile , Time = t1, isize1 = of.size1, isize2 = of.size2,
  #         tsize1 = ot.size1, tsize2 = ot.size2))

  #renombrar
  #id <- basename(patient_dir)
  val1 <- trim_dir[endsWith(trim_dir, "val_1.fq.gz")]
  trim1 <- sprintf("%s/%sT_S04_L001_R1_001.fastq.gz", dirname(val1), patient_id)
  val2 <- trim_dir[endsWith(trim_dir, "val_2.fq.gz")]
  trim2 <- sprintf("%s/%sT_S04_L001_R2_001.fastq.gz", dirname(val2), patient_id)

  file.rename(val1, trim1)
  file.rename(val2, trim2)


  message("TrimGalore's analysis has finished!")
  return(dirname(ofile))
}

#################################
#patients_dir = "~/Daniela/Biota/Muestras/73m"
#de_host = "BWA"
#mapeo_Rs <- calculate_mapeo(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "RSubread")
#mapeo_Bo <- calculate_mapeo(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "Bowtie")
#mapeo_BWA <- calculate_mapeo(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "BWA")
#mapeo_D <- calculate_mapeo(patients_dir = "~/Daniela/Biota/Muestras/73m", de_host = "DRAGEN")


calculate_mapeo <- function(patients_dir, de_host) {

    mapeo <- data.frame("Sample"= c(), "PorcentajeMapeo" = c(), "Method" = c())

    if(de_host  == "Bowtie") {
      de_host_file <- "DHBo"
    } else if( de_host == "BWA") {
      de_host_file <- "DHbwa"
    } else if(de_host == "RSubread") {
      de_host_file <- "DHRs"
    } else if(de_host == "") {
      de_host_file <- "T"
    } else if(de_host == "DRAGEN") {
      de_host_file <- "DRAGEN"
    }

    i=1
    dir_list <- list.dirs(path = patients_dir, full.names = TRUE, recursive = FALSE)
    cant_patients <- length(dir_list)

    for (p in dir_list) {
      p <- dir_list[[i]]
      print(p)
      id <- basename(p)

      if(de_host == "BWA"){
        porcentaje_mapeo <- calculate_BWA_mapeo(patient_dir = p)

      } else if (de_host == "RSubread") {
        Rs_summary <- as.data.frame(read_delim(sprintf("%s/trimmed/%s_R1R2_Rsubread.bam.summary", p, id)))
        Rs_summary <- rbind(colnames(Rs_summary) , Rs_summary)
        colnames(Rs_summary) <- c("Metrica", "Valor")
        total_fragments <- as.numeric(Rs_summary$Valor[Rs_summary$Metrica == "Total_fragments"])
        mapped_fragments <- as.numeric(Rs_summary$Valor[Rs_summary$Metrica == "Mapped_fragments"])
        porcentaje_mapeo <- round(mapped_fragments/total_fragments*100,2)

      } else if (de_host == "Bowtie") {
        error_file <- sprintf("%s/trimmed/%s_Bo_Error.txt", p, id)
        if(!file.exists(error_file)) {
          dhbo_file_R1 <- sprintf("%s/trimmed/%sDHBo_S04_L001_R1_001.fastq.gz", p, id)
          dhbo_file_R2 <- sprintf("%s/trimmed/%sDHBo_S04_L001_R2_001.fastq.gz", p, id)

          viejo_dhbo_file_R1 <- sprintf("%s/trimmed/viejo_%sDHBo_S04_L001_R1_001.fastq.gz", p, id)
          viejo_dhbo_file_R2 <- sprintf("%s/trimmed/viejo_%sDHBo_S04_L001_R2_001.fastq.gz", p, id)

          file.rename(from = dhbo_file_R1, to = viejo_dhbo_file_R1)
          file.rename(from = dhbo_file_R2, to = viejo_dhbo_file_R2)

          deHosting(patient_dir = p, de_host = "Bowtie")
        }

        Bo_error <- as.data.frame(read_delim(error_file))
        porcentaje_mapeo <- Bo_error[nrow(Bo_error), 1]
        porcentaje_mapeo <- round(as.numeric(gsub("%", "", porcentaje_mapeo)),2)

      } else if (de_host == "DRAGEN") {

        patient_dir <- paste(p, "/trimmed", sep="")
        nombre_p <- sprintf("%sT", id)

        folder_report <- list.dirs( sprintf("%s/DRAGEN_Reports", patient_dir), full.names = TRUE, recursive = FALSE)
        metrics_report <- sprintf("%s/%s.mapping_metrics.csv", folder_report, nombre_p)
        metrics_report <- read_csv(metrics_report)
        metrics_report <- rbind(colnames(metrics_report), metrics_report)
        colnames(metrics_report) <- c("Mapping", "Sample", "Metric", "Value", "Percentage")
        porcentaje_mapeo <- as.numeric(metrics_report$Percentage[metrics_report$Metric == "Mapped reads" & grepl(id, metrics_report$Sample)])
        #num_readsR1 <- round(as.numeric(metrics_report$Value[metrics_report$Metric == "Unmapped reads" & grepl(id, metrics_report$Sample)])/2)
        #num_readsR2 <- num_readsR1
      }

      mapeo[i, "Sample"] <- basename(p)
      mapeo[i, "PorcentajeMapeo"] <- porcentaje_mapeo
      mapeo[i, "Method"] <- de_host_file

      i= i+1
    }

    write.xlsx(mapeo, file = sprintf("~/Daniela/Biota/PipelineBiota/paraPaper/Porcentaje_mapeo_%s.xlsx", de_host_file))
    return(mapeo)
}

calculate_BWA_mapeo <- function(patient_dir) {

  id <- basename(patient_dir)
  patient_dir_trim <- paste0(patient_dir, "/trimmed", sep="")
  file_list_trimmed <- list.files(patient_dir_trim, full.names = TRUE, recursive = FALSE)
  aligned_bam <- path.expand(sprintf("%s/%s_bwa.bam", patient_dir_trim, id))
  #bwa_summary <- path.expand(sprintf("%s/%s_bwa.summary", patient_dir_trim, id))

  if(!file.exists(aligned_bam)) {

    r1_trim <- sprintf("%s/%sT_S04_L001_R1_001.fastq.gz", patient_dir_trim, id)
    r2_trim <- sprintf("%s/%sT_S04_L001_R2_001.fastq.gz", patient_dir_trim, id)
    BWA <- downloadBWA()
    indexBWA <- "/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/HG38/indexBWA"
    # Guardar el BAM intermedio antes de bedtools
    system2(
      BWA,
      args = c(
        "mem", "-t", "10",
        file.path(indexBWA, "Homo_sapiens.GRCh38.dna_sm.primary_assembly"),
        r1_trim, r2_trim,
        sprintf("| /home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/Samtools/samtools-1.16.1/samtools view -bh -o %s", aligned_bam)
      ))
  }

  # Función para ejecutar comandos y leer salida
  run_cmd <- function(cmd) {
    con <- pipe(cmd)
    output <- readLines(con)
    close(con)
    return(output)
  }
  flagstat_output <- run_cmd(sprintf("/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/Samtools/samtools-1.16.1/samtools flagstat %s", aligned_bam))
  # Extraer valores clave
  total_fragments <- as.numeric(strsplit(flagstat_output[1], " ")[[1]][1])
  mapped_fragments <- as.numeric(strsplit(flagstat_output[7], " ")[[1]][1])
  porcentaje_mapeo <- round(mapped_fragments/total_fragments*100,2)

  return(porcentaje_mapeo)
}

