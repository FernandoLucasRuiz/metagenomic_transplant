
load("Rdatas/lectura_conteos.RData")
load("Rdatas/covariables_limpieza.RData")
source("RMDs/Librerias.R")
source("RMDs/utils.R")



# vamos a colapsar las columnas de taxonomia para mia
taxonomia_combinada <- todas_muestras %>%
    unite("Combined", Phylum:Species, sep = "_", remove = F)

# buscamos los unique para ponerle IDs aleatorios
taxonos <- unique(taxonomia_combinada$Combined)
set.seed(123)  # Para reproducibilidad
random_ids <- replicate(length(taxonos), paste(sample(c(letters, 0:9), 8, replace = TRUE), collapse = ""))

taxonos_ids <- data.frame(cbind(taxonos, random_ids))

# juntar taxonos_ids a todas_muestras 
todas_muestras_taxonomia <- dplyr::left_join(taxonomia_combinada, taxonos_ids, by = c("Combined" = "taxonos"))



# Crear tabla de abundancia según MIA
## voy a sumar aquellos counts que tengan random_ids y codigo iguales ya que provienen de los "slash calls"
abundance_table <- todas_muestras_taxonomia %>%
    dplyr::filter(tipo_muestra == "LP") %>%
    group_by(random_ids, codigo) %>%
    summarize(Count = sum(as.numeric(Count), na.rm = TRUE)) %>%
    ungroup() %>%
    pivot_wider(names_from = codigo, values_from = Count) %>%
    mutate_all(~replace_na(., 0))

abundance_table <- as.data.frame(abundance_table)
rownames(abundance_table) <- abundance_table$random_ids
abundance_table <- abundance_table %>% select(-random_ids)

# Crear rowdata(tax)
tax <- todas_muestras_taxonomia %>%
    dplyr::filter(tipo_muestra == "LP") %>%
    dplyr::select(Phylum:Species, random_ids) %>%
    as.data.frame()

tax <- tax[!duplicated(tax$random_ids), ]

rownames(tax) <- tax$random_ids
tax <- dplyr::select(tax, -random_ids)

# Creas colData (samples). Esto incluye la metadata
samples <- todas_muestras_taxonomia %>%
    dplyr::filter(tipo_muestra == "LP") %>%
    dplyr::select(codigo, nombre_muestra, tipo_muestra, Batch) %>%
    as.data.frame()

samples <- samples[!duplicated(samples$codigo), ]

samples_metadata <- dplyr::left_join(samples, covs, by= "nombre_muestra" )

rownames(samples_metadata) <- samples_metadata$codigo
samples_metadata <- dplyr::select(samples_metadata, -codigo)




# taxa tiene que tener la taxonomia en las columnas y los samples en las filas
taxa <- data.frame(t(abundance_table))

# Para sacar la lista de los batch primero hay que ordenar los samples en el samples_metadata

## Asegurarse de que ambos dataframes tengan los mismos rownames
common_rows <- intersect(rownames(taxa), rownames(samples_metadata))

taxa <- taxa[common_rows, ]
samples_metadata <- samples_metadata[match(rownames(taxa), rownames(samples_metadata)), ]

covar_LP <- samples_metadata 

taxa_LP <- taxa[rownames(taxa) %in% rownames(covar_LP), ]

batchid_LP <- as.factor(covar_LP$Batch)

covar_LP <- covar_LP %>%
    dplyr::select(-c("nombre_muestra", "Batch", "tipo_muestra"))




# Definir las variables
variables <- c("Supervivencia", "MEAF_score", "AR", "AHT", "Biliary_complications")

# Generar todas las combinaciones posibles de 2 o más variables
combinaciones <- list()

for (i in 2:length(variables)) {  # Empezamos desde combinaciones de 2
    combinaciones[[i-1]] <- combn(variables, i, simplify = FALSE)
}

# Convertir a una sola lista de combinaciones
combinaciones_totales <- unlist(combinaciones, recursive = FALSE)

# Filtrar combinaciones que incluyan "MEAF_score"
combinaciones_con_meaf <- Filter(function(x) "MEAF_score" %in% x, combinaciones_totales)

# Convertir cada combinación en una cadena separada por comas
combinaciones_cadenas <- sapply(combinaciones_con_meaf, function(x) paste(x, collapse = " y "))


permanovas <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(permanovas) <- c("Batch Bray-Curtis", "Key Bray-Curtis", "Batch Aitchinson", "Key Aitchinson")

plots <- list()

for (i in levels(batchid_T1)){
    print(paste0("Analizando run = ", i))
    
    for (j in 1:length(combinaciones_con_meaf)){
        print(paste0("        Combinación: ", j, " ", combinaciones_cadenas[[j]]))
        tryCatch({
            covariables <- covar_T1 %>% 
                dplyr::select(combinaciones_con_meaf[[j]])
            
            
            taxa_coregida <- ConQuR(tax_tab=taxa_T1, 
                                    batchid=batchid_T1, 
                                    covariates=covariables, 
                                    batch_ref=i)
            
            saveRDS(taxa_coregida, paste0("RDSs/taxa_corrected_", combinaciones_cadenas[j] ,"_", i, "_T1.RDS"))
            
            print(paste0("     Taxa corregida"))
            
            
        }, error = function(e) {
            # Manejo del error, imprime un mensaje y continúa
            print(paste("Error en el RUN = ", i, ":", e$message))
        })
    }
}



