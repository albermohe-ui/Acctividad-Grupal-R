###############################################################
## TABLA COMPARATIVA COMPONENTES PCA POR TERCILES

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
## TABLA, PROPIAMENTE
###############################################################

# Ahora se clasificarán los individuos en terciles según su puntuación en PC1, PC2 y PC3. 
# Posteriormente, se compararán las variables génicas entre terciles utilizando estadísticas 
# descriptivas (mediana e IQR) y test de Kruskal-Wallis. Esto permite identificar genes cuya 
# expresión se asocia a los patrones de variabilidad capturados por el PCA.


### ============================================================================
### 9. CREACIÓN DE TERCILES DEL PCA
### ============================================================================

## Aquí cñasificamos cada individuo en “bajo / medio / alto” según su posición en PC1, PC2 y PC3

## Extraemos los scores del PCA (coordenadas de cada individuo en los componentes)

scores_pca <- as.data.frame(pca_genes$x)

## Creamos terciles (T1, T2, T3) para cada componente principal
## Dividimos los datos en 3 grupos según percentiles 33% y 66%

scores_pca <- scores_pca %>%
  mutate(
    
    PC1_t = cut(
      PC1,
      breaks = quantile(PC1, probs = c(0, 1/3, 2/3, 1)),
      include.lowest = TRUE,
      labels = c("T1", "T2", "T3")
    ),
    
    PC2_t = cut(
      PC2,
      breaks = quantile(PC2, probs = c(0, 1/3, 2/3, 1)),
      include.lowest = TRUE,
      labels = c("T1", "T2", "T3")
    ),
    
    PC3_t = cut(
      PC3,
      breaks = quantile(PC3, probs = c(0, 1/3, 2/3, 1)),
      include.lowest = TRUE,
      labels = c("T1", "T2", "T3")
    )
  )


### ============================================================================
### 10. UNIÓN DE DATOS ORIGINALES + TERCILES DEL PCA
### ============================================================================

# Usamos expresión génica ORIGINAL (sin transformar) y agrupamos por terciles del PCA

## Nos quedamos solo con los genes + los grupos de PCA
df_final <- cbind(
  simpson_df[, genes_simpson],   # datos originales de genes
  scores_pca[, c("PC1_t", "PC2_t", "PC3_t")]  # terciles
)


# Conversión explícita de genes a numéricos
df_final <- df_final %>%
  mutate(across(all_of(genes_simpson), ~ as.numeric(as.character(.))))

str(df_final)

### ============================================================================
### 11. TABLA DESCRIPTIVA - PC1
### ============================================================================

tabla_PC1 <- df_final %>%
  
  ## Eliminamos otros grupos para que no interfieran
  
  select(-PC2_t, -PC3_t) %>%
  
  ## Creamos tabla
  
  tbl_summary(
    by = PC1_t,   # agrupamos por terciles
    
    type = all_of(genes_simpson) ~ "continuous",
    
    ## Usamos mediana + rango intercuartílico
    ## (porque la mayoría de genes no siguen distribución normal)
    
    statistic = all_continuous() ~ "{median} ({p25} - {p75})",
    
    ## Número de decimales
    
    digits = all_continuous() ~ 1
  ) %>%
  
  ## Añadimos comparación estadística (Kruskal-Wallis)
  
  add_p(test = all_continuous() ~ "kruskal.test")



### ============================================================================
### 12. TABLA DESCRIPTIVA - PC2
### ============================================================================

tabla_PC2 <- df_final %>%
  
  select(-PC1_t, -PC3_t) %>%
  
  tbl_summary(
    by = PC2_t,
    type = all_of(genes_simpson) ~ "continuous",
    statistic = all_continuous() ~ "{median} ({p25} - {p75})",
    digits = all_continuous() ~ 1
  ) %>%
  
  add_p(test = all_continuous() ~ "kruskal.test")


### ============================================================================
### 13. TABLA DESCRIPTIVA - PC3
### ============================================================================

tabla_PC3 <- df_final %>%
  
  select(-PC1_t, -PC2_t) %>%
  
  tbl_summary(
    by = PC3_t,
    type = all_of(genes_simpson) ~ "continuous",
    statistic = all_continuous() ~ "{median} ({p25} - {p75})",
    digits = all_continuous() ~ 1
  ) %>%
  
  add_p(test = all_continuous() ~ "kruskal.test")


### ============================================================================
### 14. COMBINAR TABLAS
### ============================================================================

tabla_combinada <- tbl_merge(
  tbls = list(tabla_PC1, tabla_PC2, tabla_PC3),
  tab_spanner = c("**PC1**", "**PC2**", "**PC3**")
)


### ============================================================================
### 15. MOSTRAR TABLA COMBINADA
### ============================================================================

tabla_combinada
