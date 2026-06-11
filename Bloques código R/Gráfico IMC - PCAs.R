##############################################################
## PCA DE EXPRESIÓN GÉNICA Y VISUALIZACIÓN SEGÚN IMC
##############################################################


# Este gráfico permite evaluar si existe una separación de pacientes en función
# de su IMC dentro del espacio reducido del PCA. Cada punto representa un individuo,
# y la proximidad entre puntos indica similitud en los perfiles de expresión génica.
# Si las categorías de IMC aparecen agrupadas, podría sugerir una relación entre
# expresión génica y estado nutricional.

## Se representaron elipses de dispersión en el plano definido por las dos primeras componentes principales 
# para visualizar la variabilidad intra-grupo según categorías de IMC. 
# ¿La superposición de dichas elipses sugiere la ausencia de una clara separación entre los perfiles de 
# expresión génica en función del estado de IMC?




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
### 9. CREACIÓN DE VARIABLE IMC CATEGÓRICA
### ============================================================================

# En este dataset, el IMC ya es numérico, no requiere limpieza

# Se crean categorías clínicas de IMC
simpson_df$IMC_cat <- cut(
  simpson_df$imc_kg_m2,
  breaks = c(-Inf, 24.9, 29.9, Inf),
  labels = c("Normal", "Sobrepeso", "Obesidad")
)

# Se comprueba la distribución de pacientes por categoría
table(simpson_df$IMC_cat)

### ============================================================================
### 10. PREPARACIÓN PARA VISUALIZACIÓN
### ============================================================================

# Se extraen las coordenadas de los individuos en el PCA
pca_df <- as.data.frame(pca_genes$x)

# Se añaden variables clínicas y de identificación
pca_df$IMC_cat <- simpson_df$IMC_cat
pca_df$nombre <- simpson_df$nombre

### ============================================================================
### 11. VISUALIZACIÓN DEL PCA SEGÚN IMC
### ============================================================================

# Se representa el plano formado por PC1 y PC2
# Los individuos se colorean según su categoría de IMC

ggplot(pca_df, aes(x = PC1, y = PC2, color = IMC_cat)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCA de expresión génica en pacientes",
    subtitle = "Distribución de individuos según categorías de IMC",
    x = "Componente Principal 1",
    y = "Componente Principal 2",
    color = "Categoría IMC"
  ) +
  scale_color_manual(values = c(
    "Normal" = "#4CAF50",
    "Sobrepeso" = "#FFC107",
    "Obesidad" = "#F44336"
  ))

### ============================================================================
### 11. VISUALIZACIÓN DEL PCA CON ELIPSES POR IMC
### ============================================================================

# En este gráfico añadimos elipses para representar la dispersión de cada grupo
# de IMC en el espacio definido por PC1 y PC2. Estas elipses ayudan a visualizar
# si los grupos están separados o solapados.

ggplot(pca_df, aes(x = PC1, y = PC2, color = IMC_cat)) +
  
  # Puntos individuales (cada paciente)
  geom_point(size = 3, alpha = 0.8) +
  
  # Elipses por grupo de IMC
  stat_ellipse(aes(group = IMC_cat), linetype = 2, size = 1) +
  
  # Tema limpio
  theme_minimal() +
  
  # Etiquetas
  labs(
    title = "PCA de expresión génica en pacientes",
    subtitle = "Distribución por IMC con elipses de dispersión",
    x = "Componente Principal 1",
    y = "Componente Principal 2",
    color = "Categoría IMC"
  ) +
  
  # Colores personalizados
  scale_color_manual(values = c(
    "Normal" = "#4CAF50",
    "Sobrepeso" = "#FFC107",
    "Obesidad" = "#F44336"
  ))

