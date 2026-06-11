###############################################################
## CLUSTERING (3) PCA DE INDIVIDUOS
###############################################################

#============================================================
## EMPEZAMOS CON EL BLOQUE DE CREACIÓN DEL PCA: ELIMINAR LUEGO PARA 
## EL Rmkd
#=====================================================================

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
library(flextable)


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


#================================================================
#================================================================

###############################################################
## CLUSTERING, PROPIAMENTE
###############################################################

## El análisis de clustering basado en los tres primeros componentes principales identificó 
# tres grupos de individuos con perfiles de expresión diferenciados, mostrando separación 
# principalmente a lo largo del eje PC1



# Ahora vamos a agrupar individuos, no genes 

### ============================================================================
### 1. SCORES DEL PCA (INDIVIDUOS)
### ============================================================================

# Extraemos los scores de los individuos Simpson

scores_ind <- as.data.frame(pca_genes$x[, 1:3])


### ============================================================================
### 2. CLUSTERING DE INDIVIDUOS
### ============================================================================

# Agrupamos por el método k-means, indicando tres centros (grupos)

set.seed(123)

kmeans_ind <- kmeans(scores_ind, centers = 3)

## Añadimos cluster

scores_ind$cluster <- as.factor(kmeans_ind$cluster)


### ============================================================================
### 3. VISUALIZACIÓN
### ============================================================================


ggplot(scores_ind, aes(x = PC1, y = PC2, color = cluster)) +
  
  geom_point(size = 3) +
  
  stat_ellipse(level = 0.95, linewidth = 1) +   # 🔥 elipses
  
  theme_minimal() +
  
  labs(
    title = "Clustering de individuos basado en PCA",
    x = "PC1",
    y = "PC2",
    color = "Cluster"
  )