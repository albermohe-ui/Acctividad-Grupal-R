
# ===============================
# 1. LIBRERÍAS
# ===============================

# gtsummary para tablas descriptivas
# dplyr para manipulación de datos
# readr para lectura de CSV

library(gtsummary)
library(dplyr)
library(readr)


# ===============================
# 2. CARGA DE DATOS
# ===============================

setwd("C:/Users/alber/Desktop/UNIR 26 -27/Primer cuatrimestre/Estadística y R para las Ciencias de la Salud/Actividades/Actividad grupal/Trabajos")

datos <- read_csv("Base_de_datos_Los_Simpson_completa.csv")


# Exploración inicial para comprobar estructura
summary(datos)
str(datos)


# ===============================
# 3. ELIMINAR VARIABLES INNECESARIAS
# ===============================

datos <- datos %>%
  select(-...1, -nombre, -orden)


# ===============================
# 4. IDENTIFICACIÓN AUTOMÁTICA DE GENES
# ===============================

# Detecta variables con nombres en mayúsculas (genes)

genes <- names(datos) %>%
  grep("^[A-Z0-9]+$", ., value = TRUE)


# ===============================
# 5. VARIABLES CONTINUAS EXTRA
# ===============================

# Variables que deben tratarse como numéricas aunque parezcan categóricas

vars_continuas_extra <- c(
  "legumbres_g_d", "aceites_g_d", "hdl_mg_dl",
  "albumina_g_dl", "bilirrubina_mg_dl",
  "frutos_secos_g_d", "lacteos_g_d",
  "af_moderada_h_semana", "af_intensa_h_semana"
)


# ===============================
# 6. LIMPIEZA GLOBAL DEL DATASET
# ===============================

datos <- datos %>%
  
  # Limpieza general de texto
  mutate(across(where(is.character), ~trimws(.))) %>%      # elimina espacios
  mutate(across(where(is.character), ~na_if(., ""))) %>%   # "" a NA
  
  
  # Corrección de variables clave
  mutate(
    sexo = recode(sexo,
                  "M" = "Masculino",
                  "F" = "Femenino",
                  "male" = "Masculino",
                  "female" = "Femenino"
    ),
    
    estado_civil = recode(estado_civil,
                          "casado" = "Casado",
                          "soltero" = "Soltero"
    ),
    
    etnicidad = recode(etnicidad,
                       "Caucásica" = "Caucásico"
    )
  ) %>%
  
  
  # Unificación automática de categorías
  mutate(
    
    situacion_vivienda = case_when(
      grepl("prop", situacion_vivienda, ignore.case = TRUE) ~ "Propia",
      grepl("alq|apart", situacion_vivienda, ignore.case = TRUE) ~ "Alquiler",
      grepl("famil", situacion_vivienda, ignore.case = TRUE) ~ "Familiar",
      TRUE ~ situacion_vivienda
    ),
    
    enfermedades = recode(enfermedades,
                          "Enfermedades cardiovasculares" = "Enfermedad cardiovascular",
                          "Diabetes tipo 2" = "Diabetes"
    ),
    
    trastornos = recode(trastornos,
                        "No" = "Ninguno"
    )
  ) %>%
  
  
  # Conversión a numérico de genes y variables continuas especiales
  mutate(across(all_of(c(genes, vars_continuas_extra)), as.numeric)) %>%
  
  
  # Conversión correcta a factor de variables categóricas
  mutate(across(c(
    sexo, estado_civil, nivel_educativo, ocupacion,
    lugar_residencia, etnicidad, situacion_vivienda,
    composicion_familiar, enfermedades, medicamento,
    tabaco, alcohol, ansiedad, depresion, trastornos
  ), ~factor(.))) %>%
  
  
  # Eliminación de niveles vacíos en factores
  mutate(across(where(is.factor), droplevels)) %>%
  
  
  # Eliminación de filas problemáticas mínimas
  filter(!is.na(sexo))


# ===============================
# 7. IDENTIFICAR VARIABLES NUMÉRICAS CLÍNICAS
# ===============================

numeric_vars <- datos %>%
  select(where(is.numeric))

numeric_vars_no_genes <- numeric_vars %>%
  select(-all_of(genes))


# ===============================
# 8. TEST DE NORMALIDAD (Shapiro-Wilk)
# ===============================

check_normality <- function(x) {
  
  x <- na.omit(x)
  
  if (is.numeric(x) && length(x) >= 10 && length(unique(x)) > 5) {
    p <- shapiro.test(x)$p.value
    return(p > 0.05)
  } else {
    return(FALSE)
  }
}


# Clasificación de variables

normal_vars <- names(numeric_vars_no_genes)[
  sapply(numeric_vars_no_genes, check_normality)
]

non_normal_vars <- setdiff(names(numeric_vars_no_genes), normal_vars)


# ===============================
# 9. TABLA DESCRIPTIVA FINAL
# ===============================

tabla <- datos %>%
  tbl_summary(
    
    # Forzar tipo continuo en genes y variables mal interpretadas
    type = list(
      all_of(c(genes, vars_continuas_extra)) ~ "continuous"
    ),
    
    
    # Estadísticos personalizados
    statistic = list(
      all_of(genes) ~ "{median} ({p25}, {p75})",
      all_of(normal_vars) ~ "{mean} ({sd})",
      all_of(non_normal_vars) ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    
    digits = all_continuous() ~ 2
  ) %>%
  
  # Añadir número de observaciones
  add_n() %>%
  
  # Resaltar etiquetas
  bold_labels()


# ===============================
# 10. MOSTRAR RESULTADO
# ===============================

tabla
