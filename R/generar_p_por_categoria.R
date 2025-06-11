generar_p_por_categoria <- function(CPM_vias, de_host_file, nivel) {

  CPM_viasT <- as.data.frame(t(CPM_vias))
  colnames(CPM_viasT) <- CPM_viasT[1,]
  CPM_viasT <- CPM_viasT[-1,]

  #Agrego metadata:
  CPM_viasT <- cbind("ID" = rownames(CPM_viasT), CPM_viasT)

  MetadataB <- as.data.frame(read_excel(sprintf("%s/Metadata-soloColumnasUsables.xlsx", pipe_data)))
  MetadataB <- MetadataB[which(MetadataB$ID %in% CPM_viasT$ID),]
  all(CPM_viasT$ID %in% MetadataB$ID)
  all(CPM_viasT$ID == MetadataB$ID)

  colnames(MetadataB)
  CPM_vias_completo <- merge(CPM_viasT, MetadataB[, c("ID", "Sexo", "Rango etario")], by = "ID")
  str(CPM_vias_completo)
  n <- ncol(CPM_viasT)
  colnames(CPM_vias_completo)[n]
  colnames(CPM_vias_completo)[n+1]
  CPM_vias_completo[,2:n] <- lapply(CPM_vias_completo[,2:n], as.numeric)
  CPM_vias_completo[,(n+1):ncol(CPM_vias_completo)] <- lapply(CPM_vias_completo[,(n+1):ncol(CPM_vias_completo)], as.factor)

  #write.xlsx(CPM_vias_completo, file =  sprintf("~/Daniela/Biota/Muestras/SubsetPathways/%s/%s_%svias_completo_AR_83p.xlsx", de_host_file, de_host_file, length(2:n)))


  #Calculo pvalues:
  p_por_categoria <- data.frame("Categoria" = c(), "Via" = c(), "Wilcoxon"= c(), "GrupoDominante" = c())

  categorias <- colnames(CPM_vias_completo)[(n+1):ncol(CPM_vias_completo)]
  vias <- colnames(CPM_vias_completo)[2:n]
  c= n+1
  j = 2
  i=1

  #c= n+1
  for (c in (n+1):ncol(CPM_vias_completo)) {
    categoria <- colnames(CPM_vias_completo)[c]
    print(categoria)

    for (j in 2:(n)) {

      via <- colnames(CPM_vias_completo)[j]
      print(via)
      df_sin_na <- CPM_vias_completo[complete.cases(CPM_vias_completo[, categoria]), ]

      if (any(df_sin_na[,categoria] == "Desconocido")) {
        df_sin_na <- df_sin_na[-which(df_sin_na[,categoria] == "Desconocido"),]
      }

      if (any(is.na(df_sin_na[,categoria]))) {
        df_sin_na <- df_sin_na[-which(is.na(df_sin_na[,categoria])),]
      }

      if( categoria == "Maquillaje_Base" | categoria == "CuándoMaquillaje") {
        df_sin_na <- df_sin_na[-which(df_sin_na$Sexo == "Masculino"),]
      }

      #if( categoria == "AfeccionesPiel") {
      #  df_sin_na <- df_sin_na[-which(df_sin_na$`Rango etario` != "18-35"),]
      #}

      # Calcular la media para cada grupo en la categoría para sacar el grupo dominante:
      mean_values <- aggregate(df_sin_na[, via], by = list(df_sin_na[, categoria]), FUN = mean)
      colnames(mean_values) <- c("Grupo", "Media")
      grupo_dominante <- mean_values$Grupo[which.max(mean_values$Media)]
      print(grupo_dominante)

      if( length(levels(df_sin_na[, categoria])) > 2 ) {
        k_test <- kruskal.test(df_sin_na[, via] ~ df_sin_na[, categoria], data = df_sin_na)
        p_v <- k_test$p.value
        p_valor <- p_v
        p_val <- p_v
        p_adj <- p_v
      } else if (length(levels(df_sin_na[, categoria])) == 2 ) {
        w_test <- wilcox.test(df_sin_na[,via] ~ df_sin_na[,categoria], df_sin_na)
        p_valor <- w_test$p.value
      } else {
        p_valor <- 1
      }

      p_por_categoria[i,"Categoria"] <- categoria
      p_por_categoria[i,"Via"] <- via
      p_por_categoria[i,"Wilcoxon"] <- p_valor

      p_por_categoria[i, "GrupoDominante"] <-  ifelse(length(grupo_dominante) == 0, NA, as.character(grupo_dominante))

      i <- i+1
    }
  }

  unique(p_por_categoria$Categoria)
  p_por_categoria <- p_por_categoria[p_por_categoria$Categoria %in% c("Sexo", "Rango etario",
                                                                      "Tipodepiel", "Maquillaje_Base",
                                                                      "FacilidadBroncearse", "ActividadFísica",
                                                                      "CuandoMaquillaje"), ]
  write.xlsx(p_por_categoria, file =  sprintf("~/Daniela/Biota/PipelineBiota/paraPaper/%s_%s_p_por_categoria_83p.xlsx", nivel, de_host_file))
  return(p_por_categoria)

}
