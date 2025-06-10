#' @title Generate OTUs Table Individual
#' @description kvjdfnkvjdf
#' @param patient_dir fskbfjsf
#' @param source it must be "DRAGEN" or "KRAKEN"
#' @return tabla_otus dataframe that formats KRAKEN or DRAGEN reports in a comparable way
#' @examples tabla_otus_kraken <- generateOTUsTable_Individual(patient_dir ="~/Biota/Nuevas13/Muestras/33", source= "KRAKEN" )
#' @export
#' @import data.table
#' @import openxlsx

generateOTUsTable_Individual <- function(patient_dir, source, conEukaryota = FALSE, de_host, paraBiota) {
  library(data.table)

  if(de_host  == "Bowtie") {
    de_host_file <- "Bo"
  } else if( de_host == "BWA") {
    de_host_file <- "bwa"
  } else if(de_host == "RSubread") {
    de_host_file <- "Rs"
  } else if(de_host == "") {
    de_host_file <- "sin"
  } else if(de_host == "sinDH_PD") { # sin de host previo ni dehost dragen
    de_host_file <- "sinDH_PD"
  } else {
    stop("de_host must be Bowtie, BWA, RSubread, sinDH_PD or empty string")
  }

  id <- basename(patient_dir)
  if( id == "trimmed") {
    id <- paste(basename(dirname(patient_dir)), "_", id, sep ="")
  }
  #if(source == "KRAKEN") {
    patient_dir <- paste(patient_dir, "/trimmed", sep="")
  #}


  if(conEukaryota == TRUE) {
    conE <- "conEUKARYOTA"
  } else {
    conE <- ""
  }

  if(file.exists(paste(patient_dir, "/Resultados_", source,"/TablaOTUS", conE, "_", id,de_host_file,"_", source, ".xlsx", sep=""))) {
    message("This OTU table has already been generated")
    otu_table <- as.data.frame(read_excel(paste(patient_dir, "/Resultados_", source, "/TablaOTUS", conE, "_", id, de_host_file,"_", source, ".xlsx", sep="")))
    return(otu_table)
  }

  if(source == "KRAKEN") {
    #ReportSequences <- read.table(paste(patient_dir, "/report.sequences", sep=""), header = FALSE, sep = "\t")
    ReportSequences <- fread(paste(patient_dir, sprintf("/Resultados_KRAKEN/report_%s.sequences", de_host_file), sep=""),  header = FALSE,sep = "\t")
  } else {

    list_files <- list.files(paste0(patient_dir, "/Resultados_DRAGEN", sep=""), full.names = TRUE, recursive = TRUE)
    path_report <- list_files[which(grepl(sprintf("%s%s.DRAGEN-report.tsv", id, de_host_file), list_files))]

    #ReportSequences <- read.table(path_report, header = FALSE, sep = "\t")
    ReportSequences <- fread(path_report, header = FALSE, sep = "\t")
  }

  ReportSequences$V4 <- as.factor(ReportSequences$V4)
  levels(ReportSequences$V4)
  #ReportSequences <- ReportSequences[-which(ReportSequences[,c(1,2,3)] == c(0,0,0))]
  tabla_otus <- data.frame("Domain" = c(), "Kingdom" = c(), "Phylum" = c(), "Class" = c(), "Order" = c(), "Family" = c(), "Genus" =c(), "Species" = c(), "SubSpecies" = c(), "Unknown" = c(), "CONTEO" = c(), "Acumulativo" =c())

  j=1
  domain <- "-"
  kingdom <- "-"
  phylum <- "-"
  class <- "-"
  order <- "-"
  family <- "-"
  genus <- "-"
  species <- "-"
  subspecies <- "-"
  unknown <- "-"


  for (i in 4:nrow(ReportSequences)) {
    #cant_espacios <- sum(attr(gregexpr(" ", ReportSequences$V6[i])[[1]], "match.length"))
    tax  <- gsub("^\\s+", "", ReportSequences$V6[i])

    #if (grepl("D", ReportSequences$V4[i]) | cant_espacios == 4) {
    if (ReportSequences$V4[i] == "D") {
      domain <- tax
      kingdom <- "-"
      phylum <- "-"
      class <- "-"
      order <- "-"
      family <- "-"
      genus <- "-"
      species <- "-"
      subspecies <- "-"
      unknown <- "-"

      #} else if (grepl("K",ReportSequences$V4[i]) | cant_espacios == 7) {
    } else if (ReportSequences$V4[i] == "K" ) {
      kingdom <- tax
      phylum <- "-"
      class <- "-"
      order <- "-"
      family <- "-"
      genus <- "-"
      species <- "-"
      subspecies <- "-"
      unknown <- "-"

      #} else if (grepl("P", ReportSequences$V4[i]) | cant_espacios == 8) {
    } else if (ReportSequences$V4[i] == "P") {
      phylum <- tax
      class <- "-"
      order <- "-"
      family <- "-"
      genus <- "-"
      species <- "-"
      subspecies <- "-"
      unknown <- "-"

      #} else if (grepl("C", ReportSequences$V4[i]) | cant_espacios == 10) {
    } else if (ReportSequences$V4[i] == "C") {
      class <- tax
      order <- "-"
      family <- "-"
      genus <- "-"
      species <- "-"
      subspecies <- "-"
      unknown <- "-"

      #} else if (grepl("O", ReportSequences$V4[i]) | cant_espacios == 12) {
    } else if (ReportSequences$V4[i] == "O") {
      order <- tax
      family <- "-"
      genus <- "-"
      species <- "-"
      subspecies <- "-"
      unknown <- "-"

      #} else if (grepl("F", ReportSequences$V4[i]) | cant_espacios == 14) {
    } else if (ReportSequences$V4[i] == "F") {
      family <- tax
      genus <- "-"
      species <- "-"
      subspecies <- "-"
      unknown <- "-"

      #} else if (grepl("G", ReportSequences$V4[i]) | cant_espacios == 16) {
    } else if (ReportSequences$V4[i] == "G") {
      genus <- tax
      species <- "-"
      subspecies <- "-"
      unknown <- "-"

      #} else if (grepl("S", ReportSequences$V4[i]) | cant_espacios == 19) {
    } else if (ReportSequences$V4[i] == "S") {
      species <- tax
      subspecies <- "-"
      unknown <- "-"

    } else if (ReportSequences$V4[i] == "S1" | ReportSequences$V4[i] == "S2") {
      subspecies <- tax
      unknown <- "-"
    } else {
      unknown <- ReportSequences$V6[i]

    }

    tabla_otus[j, 1] <- domain
    tabla_otus[j, 2] <- kingdom
    tabla_otus[j, 3] <- phylum
    tabla_otus[j, 4] <- class
    tabla_otus[j, 5] <- order
    tabla_otus[j, 6] <- family
    tabla_otus[j, 7] <- genus
    tabla_otus[j, 8] <- species
    tabla_otus[j, 9] <- subspecies
    tabla_otus[j, 10] <- unknown
    tabla_otus[j, 11] <- ReportSequences$V3[i]
    tabla_otus[j, 12] <- ReportSequences$V2[i]
    j = j+1
  }

  colnames(tabla_otus) <- c("Domain" , "Kingdom", "Phylum" , "Class", "Order", "Family", "Genus", "Species", "SubSpecies", "Unknown" ,  paste(id, "_CONTEO", sep = ""), paste(id, "Acumulativo", sep=""))
  library(openxlsx)

  #Me quedo con el acumulativo
  tabla_otus <- unique(tabla_otus)
  tabla_otus <- tabla_otus[which(tabla_otus$Unknown == "-"),]
  tabla_otus <- tabla_otus[, -c(10, 11)]

  colnames(tabla_otus)[ncol(tabla_otus)] <- paste(id, source, sep="_")
  if (conEukaryota == FALSE) {
    tabla_otus <- tabla_otus[-which(tabla_otus$Domain == "Eukaryota"),]
  }

  #if(paraBiota == TRUE) {
    tabla_otus$Species[which(tabla_otus$Species == "Propionibacterium sp. oral taxon 193")] <- "Cutibacterium modestum"
    tabla_otus$Genus[which(tabla_otus$Species == "Cutibacterium modestum")] <- "Cutibacterium"
  #}

  resultados_folder <- ifelse(source == "DRAGEN", "Resultados_DRAGEN", "Resultados_KRAKEN")
  path_tabla_otus <- paste(patient_dir, "/", resultados_folder, "/TablaOTUS", conE, "_", id, de_host_file, "_", source, ".xlsx", sep="")

  write.xlsx(tabla_otus, file = path_tabla_otus, overwrite = TRUE)
  return(tabla_otus)
}


###################################################################################

#' @title Generate OTUs Table Grupal
#' @description kvjdfnkvjdf
#' @param patients_dir path to the folder that stores all the samples folders.
#' @param source it must be "DRAGEN" or "KRAKEN"
#' @param paraBiota it can be TRUE or FALSE. If set to TRUE transform Cutibacterium modestum
#' @return marged_df dataframe that formats KRAKEN or DRAGEN reports in a comparable way. The columns contains the samples and the rows the OTUs
#' @examples tabla_otus_kraken <- generateOTUsTable_Grupal(patient_dir ="~/Biota/Nuevas13/Muestras", source= "KRAKEN" )
#' @export
generateOTUsTableGrupal <- function(patients_dir, source, conEukaryota = F, de_host, paraBiota) {

  if(de_host  == "Bowtie") {
    de_host_file <- "Bo"
  } else if( de_host == "BWA") {
    de_host_file <- "bwa"
  } else if(de_host == "RSubread") {
    de_host_file <- "Rs"
  } else if(de_host == "") {
    de_host_file <- "sin"
  } else if(de_host == "sinDH_PD") {
    de_host_file <- "sinDH_PD"
  } else if (de_host == "Biota") {
    de_host_file <- "trimmed"
  }


  library(readxl)
  patients <- list.dirs(patients_dir, full.names = TRUE, recursive = FALSE)
  #patients <- patients[14]
  #patients <- pacientes_conDRAGENreport
  #patients <- list_dirs

  if(conEukaryota == TRUE) {
    conE <- "conEUKARYOTA"
  } else {
    conE <-""
  }

  for(p in patients) {
    #p <- patients[14]
    id <- basename(p)
    print(id)

    if(de_host == "Biota") {
      otus_file <- sprintf("%s/%s/trimmed/TablaOTUS_%s_trimmed_%s.xlsx", patients_dir, id, id, source)
    } else {
      otus_file <- sprintf("%s/%s/trimmed/Resultados_%s/TablaOTUS%s_%s%s_%s.xlsx", patients_dir, id, source, conE, id, de_host_file, source)
    }

    if( !(file.exists(otus_file))) {
      generateOTUsTable_Individual(patient_dir = p, source = source, conEukaryota = conEukaryota, de_host = de_host, paraBiota = paraBiota)
    } else {
      print("La tabla de otus ya fue generada para esta muestra")
    }

  }

  id <- basename(patients[1])
  if(source == "KRAKEN") {
    #if(de_host == "Biota") {
    #  merged_df <- read_excel(sprintf("%s/%s/trimmed/TablaOTUS_%s_trimmed_%s.xlsx", patients_dir, id, id, source))
    #} else {
      merged_df <- read_excel(sprintf("%s/%s/trimmed/Resultados_KRAKEN/TablaOTUS%s_%s%s_%s.xlsx", patients_dir, id, conE, id, de_host_file, source))
    #}

  } else { # DRAGEN
    #if(de_host == "Biota") {
    #  merged_df <- read_excel(sprintf("%s/%s/trimmed/TablaOTUS_%s_%s.xlsx", patients_dir, id, id, source))
    #} else {
      merged_df <- read_excel(sprintf("%s/%s/trimmed/Resultados_DRAGEN/TablaOTUS%s_%s%s_%s.xlsx", patients_dir, id, conE, id, de_host_file, source))
    #}

  }
  common_cols <- names(merged_df)[1:9]

  #Me tengo que quedar con los que son>0 para poder comparar con los DRAGEN:
  merged_df <- merged_df[which(merged_df[,ncol(merged_df)]>0),]

  for (patient in patients[2:length(patients)]) {
    id <- basename(patient)
    print(id)



    if(source == "KRAKEN") {
      if(de_host == "Biota") {
        otusk_file <- sprintf("%s/%s/trimmed/TablaOTUS_%s_trimmed_%s.xlsx", patients_dir, id, id, source)
      } else {
        otusk_file <- sprintf("%s/%s/trimmed/Resultados_KRAKEN/TablaOTUS%s_%s%s_%s.xlsx", patients_dir, id, conE, id, de_host_file, source)
      }

      if (file.exists(otusk_file)) {
        tabla_otus <- read_excel(otusk_file)
      } else {
        tabla_otus <- generateOTUsTable_Individual(patient_dir = patient, source= source, conEukaryota = conEukaryota, de_host = de_host, paraBiota = paraBiota)
      }

    } else { #DRAGEN
      if(de_host == "Biota") {
        otusd_file <- sprintf("%s/%s/trimmed/TablaOTUS_%s_trimmed_%s.xlsx", patients_dir, id, id, source)
      } else {
        otusd_file <- sprintf("%s/%s/trimmed/Resultados_DRAGEN/TablaOTUS%s_%s%s_%s.xlsx", patients_dir, id, conE, id, de_host_file, source)
      }

      if (file.exists(otusd_file)) {
        tabla_otus <- read_excel(otusd_file)
      } else {
        tabla_otus <- generateOTUsTable_Individual(patient_dir = patient, source= source, conEukaryota = conEukaryota, de_host = de_host, paraBiota = paraBiota)
      }
    }

    tabla_otus <- tabla_otus[which(tabla_otus[,ncol(tabla_otus)]>0),]
    #merged <- merge(merged_df, tabla_otus, by = common_cols)
    merged_df <- merge(merged_df, tabla_otus, by = common_cols, all = TRUE)
  }

  if(any(is.na(merged_df))) {
    merged_df[is.na(merged_df)] <- 0
  }


  if(conEukaryota == F) {
    if(any(merged_df$Domain == "Eukaryota")) {
      merged_df <- merged_df[-which(merged_df$Domain == "Eukaryota"),]
    }
  }

  write.xlsx(merged_df, file = sprintf("%s/NEW-OTUs%s_%s_%s_%s.xlsx", patients_dir, source, de_host_file, conE, length(patients)), overwrite = TRUE)
  return(merged_df)
}



