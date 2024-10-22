# Manipulación de datos
library(tidyverse)
library(ggpubr)
library(paletteer)
library(ggrepel)

# Abrir archibos excel
library(readxl)

# Visualizar valores faltantes
library(visdat)

# One-hot encoding
library(fastDummies)

# Poner varias graficas juntas
library(gridExtra)

# Imputación de datos
library(mice)

# Para analisis microbiano
library(mia)
library(vegan)

# Para visualizar analisis microbianos
library(miaTime)
library(miaViz)
library(scater)
library(patchwork)

# Efectos de Batch
library(ConQuR)
library(doParallel)
library(coda.base)

# Graficas Upset
library(UpSetR)
library(grid)

# DESEq2
library(DESeq2)



num_cores <- detectCores()
registerDoParallel(cores=num_cores-2)