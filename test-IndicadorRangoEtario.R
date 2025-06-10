
AR_SubSpecies_KRAKEN <- read_excel("~/Daniela/Biota/Muestras/73m/AR_SubSpecies_KRAKEN.xlsx")
str(df_completo$Edad)
df_completo$Edad <- as.numeric(as.character(df_completo$Edad))

ggplot(df_completo, aes(x = Edad, y = Corynebacterium_kroppenstedtii_DSM_44385)) +
  geom_point(color = "blue") +
  geom_hline(yintercept = 4.73, linetype = "dashed", color = "green") +
  geom_smooth(method = "loess", color = "red") +
  coord_flip() +
  labs(title = "Dispersión de Corynebacterium según Edad",
       x = "Edad", y = "Corynebacterium") +
  theme_minimal()



e <- solo12[1]
colnames(df_completo) <- gsub(" ", "_", colnames(df_completo))
solo12 <- gsub(" ", "_", solo12)

for(e in solo12) {
  print(e)
  Edad <- "Edad"
  g <- ggplot(df_completo, aes_string(x = Edad, y = e)) +
    geom_point(color = "blue") +
    geom_hline(yintercept = df_completo[which(df_completo$ID == "4"), e], linetype = "dashed", color = "green") +
    geom_smooth(method = "loess", color = "red") +
    coord_flip() +
    labs(title = "Dispersión de Corynebacterium según Edad",
         x = "Edad", y = paste0(e)) +
    theme_minimal()
  print(g)

}

solo12 <- unique(Biomarcadores_RangoEtario$Subespecie)
df_completo_solo12 <- df_completo[, c("Edad", "ID", solo12)]
str(df_completo_solo12)
df_completo_solo12[,3:ncol(df_completo_solo12)] <- prop.table(as.matrix(df_completo_solo12[,3:ncol(df_completo_solo12)]), margin = 1) * 100
sum(AR_biomarcadores[,3:ncol(AR_biomarcadores)])

mod <- lm(Edad~., data = df_completo_solo12[,-2])
summary(mod)
plot(mod)
predict(mod, df_completo_solo12[which(df_completo_solo12$ID == "4"),])

e <- solo12[1]
for(e in solo12) {
  print(e)
  Edad <- "Edad"
  g <- ggplot(df_completo_solo12, aes_string(x = Edad, y = e)) +
    geom_point(color = "blue") +
    geom_hline(yintercept = df_completo_solo12[which(df_completo_solo12$ID == "4"), e], linetype = "dashed", color = "green") +
    geom_smooth(method = "loess", color = "red") +
    coord_flip() +
    labs(title = "Dispersión de Corynebacterium según Edad",
         x = "Edad", y = paste0(e)) +
    theme_minimal()
  print(g)

}

#####################################################################################
id = "108"
patients_dir <- "~/Daniela/Biota/Muestras/73m"
