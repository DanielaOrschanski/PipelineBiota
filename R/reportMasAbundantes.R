#' @title Report Most Abundant Bacteria
#' @import readxl
#' @import data.table
#' @import openxlsx
#' @import ggplot2
#' @import tidyr
#' @import gridExtra
#' @description It describes the relative abundance of the most abundant microorganisms (levels = phylum, genus, species, subspecies).
#' @param id patient's id
#' @param MetadataB dataframe that contains the metadata information. id must be within Metadata$ID.
#' @return list_plots list with 2 pie charts for every taxonomic level (one that represents the id's composition and the expected composition for its age range)
#' @export
#' @examples reportMasAbundantes(id = "108", MetadataB = MetadataB)
reportMasAbundantes <- function(patient_dir, MetadataB) {
  source = "KRAKEN"
  #niveles <- c("SubSpecies", "Species", "Genus", "Phylum")
  niveles <- c("Species", "Genus", "Phylum")
  i=1
  list_plots <- list()

  for(nivel in niveles) {
    #nivel <- niveles[1]
    print(nivel)

    Biomarcadores_SubEspecies_abundantes <- as.data.frame(read_excel(sprintf("%s/%sAbundantes_RangoEtario_sinoutliers.xlsx", pipe_data, nivel)))
    colnames(Biomarcadores_SubEspecies_abundantes)[1] <- nivel
    Biomarcadores_SubEspecies_abundantes[,1] <- gsub("\\.", " ", Biomarcadores_SubEspecies_abundantes[,1])
    length(unique(Biomarcadores_SubEspecies_abundantes$Species))

    #pedir id
    #id <- readline(prompt = "Ingrese el ID del paciente que desea analizar: ")
    
    #list_dirs <- list.dirs(patients_dir, recursive = FALSE, full.names = TRUE)
    #patient_dir <- path.expand(sprintf("%s/%s", patients_dir, id))

    #if(!(patient_dir %in% list_dirs)) {
    #  return("This id is not on your patients folder")
    #}

    #otus <-  as.data.frame(read_excel(sprintf("%s/trimmed/TablaOTUS_%s_trimmed_KRAKEN.xlsx", patient_dir,
    #AR_todos <- as.data.frame(read_excel(sprintf("~/Daniela/Biota/Muestras/73m/AR__%s_KRAKEN.xlsx", nivel)))
    #otus <- AR_todos[, which(colnames(AR_todos) == id)]

    #Se fija si esta la tabla de OTUS y las tablas de AR en el nivel que se necesita: --------------------

    library(dplyr)
    id <- basename(patient_dir)
    if(!(file.exists(paste(patient_dir, "/trimmed/Resultados_KRAKEN/TablaOTUS_", id, "Bo_", "KRAKEN", ".xlsx", sep="")))) {
      otus <- generateOTUsTable_Individual(patient_dir = patient_dir, source = "KRAKEN", de_host = "Bowtie")
      group_TaxonomicLevels(patients_dir = patient_dir, tabla_otus = otus, source= "KRAKEN", de_host = "Bowtie", conEukaryota = FALSE)
    }
    
    if(!(file.exists(sprintf("%s/trimmed/AR_Bo__%s_%s.xlsx", patient_dir, nivel, source)))) {
      otus <- as.data.frame(read_excel(paste(patient_dir, "/trimmed/Resultados_KRAKEN/TablaOTUS_", id, "Bo_", source, ".xlsx", sep="")))
      group_TaxonomicLevels(patients_dir = patient_dir, tabla_otus = otus, source= "KRAKEN", de_host = "Bowtie", conEukaryota = FALSE)
    }
    #-----------------------------------------------------------------------------

    otus <- as.data.frame(read_excel(sprintf("%s/trimmed/AR_Bo__%s_%s.xlsx", patient_dir, nivel, source)))

    #if(!any(colnames(AR_todos) == id)) {
    #  otus <- AR_todos[, which(colnames(AR_todos) == sprintf("%s_KRAKEN", id))]
    #}
    #otus <- data.frame(AR_todos[,1], otus)
    #colnames(otus) <- c(nivel, id)


    # me quedo con las  subespecies de la tabla y calculo AR -------------------------
    subespecies_ausar <- Biomarcadores_SubEspecies_abundantes[,1]
    subespecies_ausar <- unique(subespecies_ausar)
    otus_subespecies <- otus[which(otus[,nivel] %in% subespecies_ausar),]

    #AR_subespecies <- otus_subespecies[,9:10]
    AR_subespecies <- otus_subespecies
    str(AR_subespecies)
    AR_subespecies[,2:ncol(AR_subespecies)] <- prop.table(as.matrix(AR_subespecies[,2:ncol(AR_subespecies)]), margin = 2) * 100
    sum(AR_subespecies[,2:ncol(AR_subespecies)])

    #Extraigo los biomarcadores ---------------
    #AR_biomarcadores <- AR_subespecies[which(AR_subespecies$SubSpecies %in% Biomarcadores_SubEspecies_abundantes$Subespecie),]
    #colnames(AR_biomarcadores)[1] <- nivel
    AR_biomarcadores <- AR_subespecies


    #Veo a qué rango etario se aproxima más ----------------
    df_merged <- merge(AR_biomarcadores, Biomarcadores_SubEspecies_abundantes, by = nivel)

    edad <- MetadataB$Edad[which(MetadataB$ID == id)]
    rango_real <- MetadataB$`Rango etario`[which(MetadataB$ID == id)]

    # Digo si es alto/Eq/bajo segun su edad
    df_rango_real <- df_merged[which(df_merged$Rango == rango_real),]
    df_rango_real <- df_rango_real[order(df_rango_real[,2], decreasing = TRUE),]

    if(any(grepl("virus", df_rango_real[,1]))) {
      df_rango_real <- df_rango_real[-which(grepl("virus", df_rango_real[,1])),]
    }
    if(any(grepl("phage", df_rango_real[,1]))) {
      df_rango_real <- df_rango_real[-which(grepl("phage", df_rango_real[,1])),]

    }

    if(nivel == "Genus") {
      df_rango_medico <- df_rango_real[1:5,]
      df_rango_real <- df_rango_real[1:3,]

    } else if(nivel == "Phylum") {
      df_rango_medico <- df_rango_real[1:4,]
      df_rango_real <- df_rango_real[1:3,]

    } else { #especie
      df_rango_medico <- df_rango_real[1:10,]
      df_rango_real <- df_rango_real[1:3,]

    }


    library(dplyr)
    df_rango_real[,2] <- as.numeric(df_rango_real[,2] )
    df_rango_real <- df_rango_real %>%
      mutate(Resultados = ifelse(df_rango_real[,2] > Q1 & df_rango_real[,2] < Q3 | df_rango_real[,2] == Q1 | df_rango_real[,2] == Q3, "Equilibrado",
                                 ifelse(df_rango_real[,2] < Q1, "Bajo", "Alto")))

    library(dplyr)
    df_rango_medico[,2] <- as.numeric(df_rango_medico[,2] )
    df_rango_medico <- df_rango_medico %>%
      mutate(Resultados = ifelse(df_rango_medico[,2] > Q1 & df_rango_medico[,2] < Q3 | df_rango_medico[,2] == Q1 | df_rango_medico[,2] == Q3, "Equilibrado",
                                 ifelse(df_rango_medico[,2] < Q1, "Bajo", "Alto")))


    #Grafico de torta
    colnames(df_rango_real)[2] <- "AR"

    # Asegúrate de que 'nivel' es el nombre de la columna que quieres usar en 'fill'
    library(ggplot2)
    df_rango_real$AR <- as.numeric(df_rango_real$AR)

    if(nivel == "Species") {
      margen = 20
      units = 0.3
    } else if( nivel == "Phylum") {
      margen = 1
      units = 0.8
    } else {
      margen = 20
      units = 0.3
    }



    #Recomendacion y Funcion solo de especies para el reporte completo: -------------------------------------------------------------

    if(nivel == "Species") {
      Especies_Recomendaciones <- as.data.frame(read_excel(sprintf("%s/topTax-Recomendaciones.xlsx", pipe_data), sheet=1))
      paraCompleto <- Especies_Recomendaciones[which(Especies_Recomendaciones$Species %in% df_rango_medico$Species), c(1,2,3,5)]
      mensaje_paraCompleto <- c()
      v=1
      for(v in 1:nrow(paraCompleto)) {
        virus <- paraCompleto$Species[v]
        funcion <- paraCompleto$`Impacto en la salud de la piel`[v]
        recomendacion <- paraCompleto$`Vivi: Recomendación: si la Especie esta elevada`[v]

        #recomendacion solo si está elevada:
        estado <- df_rango_medico$Resultados[which(df_rango_medico$Species == virus)]

        if(!(is.na(funcion) & is.na(recomendacion)) & estado == "Alto") {
          mensaje_paraCompleto <- c(mensaje_paraCompleto, sprintf("La especie %s cumple con las siguientes funciones: %s Por lo tanto, se sugieren las siguientes recomendaciones: %s", virus, funcion, recomendacion))
        }
      }
    }

    #-------------------------------------------------------------------------
    # Identificar las 5 especies más abundantes
    top_5_species <- df_rango_real %>%
      arrange(desc(AR)) %>%
      head(5) %>%
      pull(!!sym(nivel))

    # Crear una nueva columna 'Species_group' que agrupe las demás especies como 'Otras'
    df_rango_real <- df_rango_real %>%
      mutate(Species_group = ifelse(!!sym(nivel) %in% top_5_species, !!sym(nivel), "Otras"))

    # Agrupar las especies bajo 'Otras' y sumar sus valores para AR y Media
    df_agrupado <- df_rango_real %>%
      group_by(Species_group) %>%
      summarise(AR = sum(AR),
                Media = sum(Media),
                Q1 = min(Q1),  # Puedes cambiar esto dependiendo de cómo quieras manejar Q1/Q3
                Q3 = max(Q3))

    # Crear un dataframe largo para graficar ambas barras
    library(tidyverse)
    df_largo <- df_agrupado %>%
      mutate(Rango_Text = paste(round(Q1), "-", round(Q3))) %>%
      pivot_longer(cols = c("AR", "Media"), names_to = "Metric", values_to = "Value")

    # Crear el factor para Species_group basado en el valor de AR
    df_largo$Species_group <- factor(df_largo$Species_group,
                                     levels = df_agrupado %>%
                                       arrange(desc(AR)) %>%  # Ordenar por AR en orden descendente
                                       pull(Species_group))  # Obtener los niveles de Species_group

    df_largo$Rango_Text <- ifelse(df_largo$Metric == "AR", round(df_largo$Value), df_largo$Rango_Text)

    # Graficar barras apiladas con la suma de las especies agrupadas bajo 'Otras'
    nombre_nivel <- ifelse(nivel == "Species", "Especies", ifelse(nivel == "SubSpecies", "Subespecies", ifelse(nivel == "Phylum", "Filos", "Géneros")))


    # Crear el gráfico con mejoras estéticas
    library(ggplot2)
    library(RColorBrewer)

    df_l <- df_largo
    df_l <- df_l[order(df_l$Value,decreasing = TRUE),]
    df_l <- df_l[which(df_l$Metric == "AR"),]

    #Llevar el elemento "Otras" al final de la lista. Los demas quedan de mayor AR a menor AR.
    etiquetas_personalizadas <- paste(unique(df_l$Species_group[which(df_l$Metric == "AR")]), "-", round(df_l$Value[which(df_l$Metric == "AR")], 2))
    if(any(grepl("Otras", etiquetas_personalizadas))) {
      otras_index <- grep("Otras", etiquetas_personalizadas)
      etiquetas_reordenadas <- c(etiquetas_personalizadas[-otras_index], etiquetas_personalizadas[otras_index])
    } else {
      etiquetas_reordenadas <- etiquetas_personalizadas
    }


    # Crear el gráfico con mejoras estéticas en los ejes y el texto dentro de las barras
    library(RColorBrewer)
    library(ggplot2)

    colores_fil <- c("#D3A1FF", "#8968CD", "#9932CC")

    colores_gen <- c( "#9ACD32",
                      "#7CCD7C",
                      "#548B54"
    )
    colores_esp <- c("#7FFFD4", "#66CDAA", "#458B74")


    if(nivel == "Phylum") {
      colores_pastel <- colores_fil
    } else if(nivel == "Genus") {
      colores_pastel <- colores_gen
    } else {
      colores_pastel <- colores_esp
    }
    colores_asignados <- setNames(colores_pastel[1:length(unique(df_largo$Species_group))], unique(df_largo$Species_group))


    bar_plot <- ggplot(df_largo, aes(x = Metric, y = round(Value), fill = Species_group)) +
      geom_bar(stat = "identity", position = "stack", width = 0.70) +  # Hacer las barras más delgadas
      geom_text(aes(label = Rango_Text),  # Redondeo al etiquetar los valores
                position = position_stack(vjust = 0.5),
                color = "gray20",  # Texto gris para destacar sobre los colores de las barras
                size = 2.5,         # Tamaño del texto ligeramente más grande
                fontface = "bold")+  # Números en negrita para mayor visibilidad
      scale_fill_manual(values = colores_pastel) +
      labs(x = "",  y = "") +
      ggtitle(sprintf("%s", nombre_nivel)) +
      scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
      #ylim(-1, 100) +  # Ajustar el límite superior del eje y para barras más largas
      theme(
        panel.background = element_rect(fill = "white", color = NA),  # Fondo blanco del panel
        plot.background = element_rect(fill = "white", color = NA),   # Fondo blanco del gráfico
        axis.text.x = element_text(
          angle = 0,
          hjust = 0.5,
          face = "bold",
          size = 10,          # Tamaño más grande
          color = "gray30"   # Color gris más suave
        ),
        axis.text.y = element_text(
          size = 7,          # Tamaño más grande
          color = "gray30"   # Color gris
        ),
        plot.title = element_text(
          hjust = 0.5,
          size = 9,           # Título más grande
          face = "bold",       # Negrita para destacar
          color = "gray20"     # Color gris oscuro
        ),
        legend.position = "bottom",
        legend.text = element_text(size = 7, face = "italic"),  # Texto de la leyenda más grande y moderno
        legend.key.size = unit(0.6, "cm"),  # Tamaño de los íconos en la leyenda
        plot.margin = margin(0, 0, 0, 0)  # Ajustar márgenes
      ) +
      # Cambiar las etiquetas del eje x
      scale_x_discrete(labels = c("AR" = "Muestra", "Media" = "Referencia")) +
      theme(axis.text.x = element_text(size = 9)) +
      # Eliminar el título de las leyendas
      guides(fill = guide_legend(ncol = 1, title = NULL))  # Ajustar el número de columnas en la leyenda


    print(bar_plot)


    ##########################################################################################################

    colnames(df_rango_medico)[2] <- "AR"
    colnames(df_rango_medico)[6] <- "Min"
    colnames(df_rango_medico)[7] <- "Max"
    colnames(df_rango_medico)[8] <- "Resultados"
    df_rango_medico <- df_rango_medico[,-c(3,4,5)]
    df_rango_medico[,2] <- round(df_rango_medico[,2], 2)

    if( colnames(df_rango_medico)[1] == "Species") {
      colnames(df_rango_medico)[1] <- "Especies"
    } else if( colnames(df_rango_medico)[1] == "Phylum") {
      colnames(df_rango_medico)[1] <- "Filos"
    } else if( colnames(df_rango_medico)[1] == "Genus") {
      colnames(df_rango_medico)[1] <- "Géneros"
    } else if( colnames(df_rango_medico)[1] == "SubSpecies") {
      colnames(df_rango_medico)[1] <- "SubEspecies"
    }


    #df_rango_medico <- df_rango_medico[, -"Species_group"]
    df_rango_medico[, c("Min", "Max")] <- round(df_rango_medico[, c("Min", "Max")],2)
    write.xlsx(df_rango_medico, file = sprintf("%s/Tabla_%s_masAbundantes.xlsx", patient_dir, nivel))


    #list_plots[[i]] <- combined_plot
    list_plots[[i]] <- bar_plot
    i <- i+1

  }
  library(gridExtra)
  plot_grid <- grid.arrange(list_plots[[3]], list_plots[[2]], list_plots[[1]], ncol = 3)

  #ggsave("~/Daniela/Biota/Muestras/73m/37/PRUEBA.png", plot_grid, height = 8, width = 14)
  dev.off()
  return(list(list_plots, plot_grid, mensaje_paraCompleto))

}

