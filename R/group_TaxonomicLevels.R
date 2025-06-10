#' @title Group Taxonomic Levels
#' @description Generates one excel file for each taxonomic level containing the relative abundances of all patients in patient_dir.
#' @param tabla_otus is a dataframe fskbfjsf
#' @param source it must be "DRAGEN" or "KRAKEN"
#' @param patients_dir path to the directory were are stored all the patients analyzed
#' @return 2 lists: one with 7 dataframes, one for each taxonomic level, that contains the relative abundances.
#' @export
#'

group_TaxonomicLevels <- function(patients_dir, tabla_otus, source, de_host, especies_seleccionadas = c(), subespecies_seleccionadas = c(), nombre_extra = "", conEukaryota) {

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
  }  else {
    stop("de_host must be Bowtie, BWA, RSubread, sinDH_PD or empty string")
  }


  Conteos_Cancer <- tabla_otus
  niveles <- colnames(Conteos_Cancer)[1:9]

  list_AR <- list()
  list_Conteos <- list()
  i=1
  #nivel <- niveles[8]
  for (nivel in niveles) {
    print(nivel)
    #Conteos_Cancer <- Conteos_Cancer[!is.na(Conteos_Cancer[, nivel]), ]
    Conteos_Cancer[is.na(Conteos_Cancer)] <- 0

    #separo solo los conteos y la columna genero (la convierto en categorica)
    generos <- as.data.frame(subset(Conteos_Cancer, select = -c(1:9)))
    ind <- which(colnames(Conteos_Cancer) == nivel)
    generos <- cbind(Conteos_Cancer[, ind:9], generos)
    l <- length(ind:9)

    #Cambiar esto para resolver el problema con los geneeros!!!!!!!!
    if(nivel == "Species") {
      filas_seleccionadas <- which(generos[, 2] == "-")
      generos <- generos[filas_seleccionadas, ]
      generos <- generos[,-c(2:l)]
    } else if(nivel == "SubSpecies") {
      generos <- generos[- which(generos[,1] == "-"), ]
      colnames(generos)[1] <- nivel
    } else {
      filas_seleccionadas <- apply(generos[, 2:l], 1, function(fila) all(fila == "-"))
      generos <- generos[filas_seleccionadas, ]
      generos <- generos[,-c(2:l)]
    }
    #Me quedo con las filas en los cuales tengo NOMBRE DE FAMILIA y todas las columnas que siguen hacia la derecha son -.

    if(conEukaryota == FALSE) {
      if(any(generos[, nivel] == "Eukaryota")) {
        generos <- generos[ -which(generos[, nivel] == "Eukaryota"),]
      }
      if(any(generos[, nivel] == "Homo sapiens")) {
        generos <- generos[ -which(generos[, nivel] == "Homo sapiens"),]
      }
    }

    if(any(generos[, nivel] == "-")) {
      generos <- generos[ -which(generos[, nivel] == "-"),]
    }


    length(unique(generos[, nivel])) #2781 especies diferentes
    # setdiff(generos[, nivel], unique(generos[,nivel]))
    rownames(generos) <- NULL
    conteos_por_nivel <- generos %>%
      # Agrupar por nivel
      group_by(!!sym(nivel)) %>%
      # Sumar los valores de las columnas para cada especie
      summarise(across(everything(), sum, na.rm = TRUE)) %>%
      ungroup()
    length(unique(generos$Species))

    list_Conteos[[i]] <- conteos_por_nivel
    write.xlsx(conteos_por_nivel, file = sprintf("%s/trimmed/Conteos_%s_%s_%s.xlsx", patients_dir, de_host_file, nivel, source), overwrite = TRUE)

    # Abundancias Relativas todas --------------------------------------------------------------
    #Hace que la suma de todos los generos en cada muestra den 100
    conteos_por_especie <- conteos_por_nivel
    nombres_especie <- conteos_por_especie[, nivel]
    abundancias_relativas <- prop.table(as.matrix(conteos_por_especie[, -1]), margin = 2) * 100
    # calcula las proporciones de cada valor en relaci?n con la suma de la columna (muestra).
    abundancias_relativas_df <- cbind(Especie = nombres_especie, as.data.frame(abundancias_relativas))
    str(abundancias_relativas_df)
    #colSums(abundancias_relativas_df[,-1])

    #id <- strsplit(colnames(Conteos_Cancer)[ncol(Conteos_Cancer)], split = "_")[[1]][1]
    colnames(abundancias_relativas_df)[1] <- nivel
    write.xlsx(abundancias_relativas_df, file = sprintf("%s/trimmed/AR_%s_%s_%s_%s.xlsx", patients_dir, de_host_file, nombre_extra, nivel, source), overwrite = TRUE)

    list_AR[[i]] <- abundancias_relativas_df
    i <- i+1

    #Abundancias relativas de especies importantes
    if(nivel == "Species") {
      especies_imp <- c("Cutibacterium acnes", "Streptococcus oralis", "Corynebacterium striatum", "Staphylococcus epidermidis", "Micrococcus luteus", "Acinetobacter johnsonii", "Streptococcus mitis", "Staphylococcus aureus", "Cutibacterium granulosum")
      AR_esp_imp <- abundancias_relativas_df[which(abundancias_relativas_df$Species %in% especies_imp), ]

      #Dando lista de especies:
      if(length(especies_seleccionadas) != 0) {
        AR_esp_selec <- abundancias_relativas_df[which(abundancias_relativas_df$Species %in% especies_seleccionadas), ]
        AR_esp_selec_prop <- prop.table(as.matrix(AR_esp_selec[, -1]), margin = 2) * 100
        AR_esp_selec <- cbind(Species = AR_esp_selec$Species, as.data.frame(AR_esp_selec_prop))

      } else {
        AR_esp_selec <- "No se seleccionaron especies"
      }

    } else if(nivel == "SubSpecies") {
      #Abundancias relativas de SUB especies importantes
      sub_especies_imp <- c("Cutibacterium acnes", "Corynebacterium", "Staphylococcus epidermidis", "Staphylococcus aureus")
      mask <- rowSums(sapply(sub_especies_imp, function(x) grepl(x, abundancias_relativas_df$SubSpecies))) > 0
      AR_sub_esp_imp <- abundancias_relativas_df[mask, ]

      #Dando lista de sub especies:
      if(length(subespecies_seleccionadas) != 0) {
        mask <- rowSums(sapply(subespecies_seleccionadas, function(x) grepl(x, abundancias_relativas_df$SubSpecies))) > 0
        AR_sub_esp_selec <- abundancias_relativas_df[mask, ]

        AR_sub_esp_selec_prop <- prop.table(as.matrix(AR_sub_esp_selec[, -1]), margin = 2) * 100
        AR_sub_esp_selec <- cbind(SubSpecies = AR_sub_esp_selec$SubSpecies, as.data.frame(AR_sub_esp_selec_prop))

      } else {
        AR_sub_esp_selec <- "No se seleccionaron sub especies"
      }
    }


  }
  return(list(list_AR, AR_esp_imp, AR_sub_esp_imp, AR_esp_selec, AR_sub_esp_selec, list_Conteos))
}


