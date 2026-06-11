###############################################################
## ACTIVIDAD 3. PROYECTO TRANSVERSAL EN R
## Análisis bioestadístico de expresión génica en Los Simpson
###############################################################

### ============================================================================
### 0. LIMPIEZA DEL ENTORNO, REPRODUCIBILIDAD Y PAQUETES
### ============================================================================

## Limpio el entorno para evitar que objetos creados anteriormente
## interfieran con los resultados de esta actividad.

rm(list = ls())

## Fijo una semilla para que los procedimientos aleatorios,
## como el clustering, puedan reproducirse.

set.seed(2026)

## Evito que R utilice notación científica cuando no sea necesaria.

options(scipen = 999)

## Defino los paquetes que utilizaré durante el análisis.

paquetes <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "janitor",
  "factoextra",
  "FactoMineR",
  "pheatmap",
  "gtsummary",
  "flextable",
  "officer",
  "broom",
  "car",
  "pROC"
)

## Compruebo qué paquetes todavía no están instalados.

paquetes_nuevos <- paquetes[
  !paquetes %in% rownames(installed.packages())
]

## Instalo solamente los paquetes que falten.

if (length(paquetes_nuevos) > 0) {
  install.packages(paquetes_nuevos, dependencies = TRUE)
}

## Cargo las librerías necesarias.

library(dplyr)
library(tidyr)
library(ggplot2)
library(janitor)
library(factoextra)
library(FactoMineR)
library(pheatmap)
library(gtsummary)
library(flextable)
library(officer)
library(broom)
library(car)
library(pROC)

### ============================================================================
### 1. DIRECTORIO DE TRABAJO
### ============================================================================

## Defino la carpeta donde se encuentra la base de datos
## y donde guardaré los resultados de la actividad.

setwd("C:/Users/alber/Desktop/UNIR 26 -27/Primer cuatrimestre/Estadística y R para las Ciencias de la Salud/Actividades/Actividad grupal/Trabajos")



## Compruebo que la ruta sea correcta.

getwd()

## Reviso los archivos disponibles en la carpeta.

list.files()
### ============================================================================
### 2. CARGA Y EXPLORACIÓN INICIAL DE LA BASE DE DATOS
### ============================================================================

## Cargo la base de datos en formato CSV.

simpson_df <- read.csv(
  "Base_de_datos_Los_Simpson_completa.csv",
  sep = ",",
  stringsAsFactors = FALSE
)

## Reviso las dimensiones de la base de datos.

dim(simpson_df)

## Reviso la estructura general de las variables.

str(simpson_df)

## Reviso los nombres de todas las variables.

names(simpson_df)

## Reviso las primeras filas del dataset.

head(simpson_df)

## Obtengo un resumen general de las variables.

summary(simpson_df)

### ============================================================================
### 3. CONTROL DE CALIDAD INICIAL
### ============================================================================

## Compruebo el número total de valores perdidos en la base de datos.

sum(is.na(simpson_df))

is_na <- (is.na(simpson_df))

is_na

## Calculo los valores perdidos por variable.

na_por_variable <- data.frame(
  Variable = names(simpson_df),
  NA_total = colSums(is.na(simpson_df)),
  Porcentaje_NA = round(colMeans(is.na(simpson_df)) * 100, 2)
)

na_por_variable


## Muestro únicamente las variables que presentan valores perdidos.

na_por_variable_con_na <- na_por_variable %>%
  filter(NA_total > 0)

na_por_variable_con_na

## Compruebo si existen registros duplicados.

sum(duplicated(simpson_df))

### ============================================================================
### 4. SELECCIÓN DE VARIABLES DE EXPRESIÓN GÉNICA
### ============================================================================

## Compruebo la posición de los genes ADCY3 y UHMK1 dentro del dataset.

which(names(simpson_df) == "ADCY3")
which(names(simpson_df) == "UHMK1")

## Selecciono todas las variables génicas comprendidas entre ADCY3 y UHMK1.

genes_simpson <- names(simpson_df)[
  which(names(simpson_df) == "ADCY3"):
    which(names(simpson_df) == "UHMK1")
]

## Compruebo cuántos genes se han seleccionado.

length(genes_simpson)

## Muestro los nombres de los genes seleccionados.

genes_simpson

## Compruebo los valores perdidos en las variables de expresión génica.

colSums(is.na(simpson_df[, genes_simpson]))

## Calculo el total de valores perdidos en el conjunto de genes.

sum(is.na(simpson_df[, genes_simpson]))

### ============================================================================
### 5. COMPROBACIÓN DEL FORMATO DE LOS GENES
### ============================================================================

## Compruebo el tipo de dato de las variables génicas seleccionadas.

sapply(simpson_df[, genes_simpson], class)

## Compruebo si todas las variables génicas son numéricas.

sapply(simpson_df[, genes_simpson], is.numeric)

## Convierto las variables génicas a formato numérico por seguridad.

for (gen in genes_simpson) {
  simpson_df[[gen]] <- as.numeric(simpson_df[[gen]])
}

## Vuelvo a comprobar que todas las variables sean numéricas.

sapply(simpson_df[, genes_simpson], is.numeric)

### ============================================================================
### 6. TEST DE NORMALIDAD DE LOS GENES
### ============================================================================

## Aplico el test de Shapiro-Wilk a cada variable génica.
## Guardo el estadístico W, el valor p y la interpretación.

normalidad_genes <- data.frame(
  Variable = genes_simpson,
  Estadistico_W = NA,
  Valor_p = NA,
  Interpretacion = NA
)

for (i in seq_along(genes_simpson)) {
  
  prueba <- shapiro.test(simpson_df[[genes_simpson[i]]])
  
  normalidad_genes$Estadistico_W[i] <- prueba$statistic
  
  normalidad_genes$Valor_p[i] <- prueba$p.value
  
  normalidad_genes$Interpretacion[i] <- ifelse(
    prueba$p.value < 0.05,
    "Distribución no normal",
    "Distribución compatible con normalidad"
  )
}

## Creo una versión formateada para visualizar los resultados.

normalidad_genes_resultado <- normalidad_genes %>%
  mutate(
    Estadistico_W = round(Estadistico_W, 3),
    Valor_p = format(
      Valor_p,
      scientific = TRUE,
      digits = 3
    )
  )

normalidad_genes_resultado

## Resumo cuántos genes presentan distribución normal y cuántos no.

table(normalidad_genes$Interpretacion)

### ============================================================================
### 7. PREPARACIÓN DE LOS DATOS PARA EL PCA
### ============================================================================

## Creo un dataframe únicamente con las variables de expresión génica.

datos_genes <- simpson_df[, genes_simpson]

## Compruebo las dimensiones y la estructura del conjunto de genes.

dim(datos_genes)
str(datos_genes)

## Confirmo que no existan valores perdidos antes de realizar el PCA.

sum(is.na(datos_genes))

### ============================================================================
### 8. ANÁLISIS DE COMPONENTES PRINCIPALES
### ============================================================================

## Realizo el PCA con las variables génicas.
## Centro y escalo las variables para que todas tengan el mismo peso
## independientemente de su media y variabilidad original.

pca_genes <- prcomp(
  datos_genes,
  center = TRUE,
  scale. = TRUE
)

## Reviso un resumen general del PCA.

summary(pca_genes)

## Reviso los elementos guardados dentro del objeto PCA.

names(pca_genes)

### ============================================================================
### 9. VARIANZA EXPLICADA POR LOS COMPONENTES PRINCIPALES
### ============================================================================

## Extraigo los autovalores, el porcentaje de varianza explicada
## y la varianza acumulada de cada componente.

varianza_pca <- get_eigenvalue(pca_genes)

varianza_pca

## Creo una tabla ordenada con los resultados.

tabla_varianza_pca <- data.frame(
  Componente = rownames(varianza_pca),
  Autovalor = varianza_pca$eigenvalue,
  Varianza_explicada = varianza_pca$variance.percent,
  Varianza_acumulada = varianza_pca$cumulative.variance.percent
)

## Redondeo los resultados para facilitar su lectura.

tabla_varianza_pca <- tabla_varianza_pca %>%
  mutate(
    Autovalor = round(Autovalor, 3),
    Varianza_explicada = round(Varianza_explicada, 2),
    Varianza_acumulada = round(Varianza_acumulada, 2)
  )

tabla_varianza_pca

### ============================================================================
### 10. EXTRACCIÓN DE SCORES DE LAS PRIMERAS SEIS COMPONENTES
### ============================================================================

## Extraigo los scores de las primeras seis componentes principales.

scores_pca <- as.data.frame(
  pca_genes$x[, 1:6]
)

## Reviso las primeras filas de los scores.

head(scores_pca)

## Añado los scores al dataframe original.

simpson_df <- cbind(
  simpson_df,
  scores_pca
)

## Compruebo que PC1 a PC6 se hayan añadido correctamente.

names(simpson_df)

### ============================================================================
### 11. GRÁFICO DE VARIANZA EXPLICADA
### ============================================================================

## Represento el porcentaje de varianza explicado por cada componente principal.
## Este gráfico me ayuda a identificar qué componentes concentran más información.

grafico_varianza <- fviz_eig(
  pca_genes,
  addlabels = TRUE,
  ylim = c(0, 100)
) +
  ggtitle("Varianza explicada por los componentes principales") +
  theme_minimal()

grafico_varianza

### ============================================================================
### 12. COORDENADAS, COS2 Y CONTRIBUCIONES DE LAS VARIABLES
### ============================================================================

## Extraigo las coordenadas, la calidad de representación y
## las contribuciones de las variables génicas en el PCA.

variables_pca <- get_pca_var(pca_genes)

## Reviso las coordenadas de los genes.

head(variables_pca$coord)

## Reviso la calidad de representación de los genes.

head(variables_pca$cos2)

## Reviso la contribución de los genes a los componentes.

head(variables_pca$contrib)

### ============================================================================
### 13. VISUALIZACIÓN DE VARIABLES SEGÚN COS2
### ============================================================================

## Represento las variables génicas en el espacio de PC1 y PC2.
## El color indica la calidad de representación de cada variable.

grafico_variables_cos2 <- fviz_pca_var(
  pca_genes,
  axes = c(1, 2),
  col.var = "cos2",
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE
) +
  ggtitle("Variables génicas según su calidad de representación") +
  theme_minimal()

grafico_variables_cos2

### ============================================================================
### 14. CONTRIBUCIÓN DE LOS GENES A LOS COMPONENTES
### ============================================================================

## Muestro los genes que más contribuyen a la primera componente.

contribucion_pc1 <- fviz_contrib(
  pca_genes,
  choice = "var",
  axes = 1,
  top = 15
) +
  ggtitle("Contribución de los genes a PC1") +
  theme_minimal()

contribucion_pc1

## Muestro los genes que más contribuyen a la segunda componente.

contribucion_pc2 <- fviz_contrib(
  pca_genes,
  choice = "var",
  axes = 2,
  top = 15
) +
  ggtitle("Contribución de los genes a PC2") +
  theme_minimal()

contribucion_pc2
