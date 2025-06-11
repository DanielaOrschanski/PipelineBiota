# Carga el paquete ggplot2
library(ggplot2)


# Definir los niveles, tamaños y colores
niveles <- c("Filos",  "Clase", "Orden", "Familia", "Géneros", "Especies")
anchos <- c(7, 6,5, 4, 3, 2)  # Anchos para cada nivel
#colores <- c("#1f78b4", "#6baed6", "#c6dbef") # Tonos de azul
colores <- c("#0d4c92","#1f78b4", "#4ab4de", "#76c7e8", "#99d6f0", "#c6dbef")


# Crear el dataframe con los niveles, sus posiciones y anchos
df <- data.frame(
  Nivel = factor(niveles, levels = niveles),
  Posicion = c(6,5,4,3,2,1),
  Ancho = anchos
)

# Crear el gráfico
g <- ggplot(df, aes(x = 0, y = Posicion, fill = Nivel, width = Ancho)) +
  geom_tile(color = "black", height = 1) +
  geom_text(aes(label = Nivel), color = "black", size = 3, vjust = 0.5) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(limits = c(-4, 4), expand = c(0, 0)) +
  scale_fill_manual(values = colores) +
  theme_void() +
  theme(legend.position = "none") +
  #coord_flip()
  theme(legend.position = "none")

ggsave(g, file = sprintf("%s/Niveles_tax2.png", pipe_data), height = 1.2, width = 2.2)

library(ggplot2)

g <- ggplot(df, aes(x = 0, y = Posicion, fill = Nivel, width = Ancho)) +
  geom_tile(color = "grey20", height = 1, size = 0.5) + # Ajusta color y tamaño del borde
  geom_text(aes(label = Nivel), color = "white", size = 4, vjust = 0.5, fontface = "bold") + # Estilo de texto
  scale_y_continuous(expand = c(0.05, 0.05)) + # Añade un pequeño margen en el eje Y
  scale_x_continuous(limits = c(-3, 3), expand = c(0, 0)) +
  scale_fill_manual(values = colores) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "grey95", color = NA), # Fondo claro para el gráfico
    plot.margin = margin(10, 10, 10, 10) # Ajusta los márgenes
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

g
