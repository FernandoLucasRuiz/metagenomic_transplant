
load("Rdatas/lectura_conteos.RData")
load("Rdatas/covariables_limpieza.RData")
load("Rdatas/cargar_tse.RData")
source("RMDs/Librerias.R")
source("RMDs/utils.R")

# taxa tiene que tener la taxonomia en las columnas y los samples en las filas
taxa <- data.frame(t(abundance_table))

## Asegurarse de que ambos dataframes tengan los mismos rownames
common_rows <- intersect(rownames(taxa), rownames(samples_metadata))

taxa <- taxa[common_rows, ]
samples_metadata <- samples_metadata[match(rownames(taxa), rownames(samples_metadata)), ]

batchid <- as.factor(samples_metadata$Batch)

covar_lp <- samples_metadata %>%
    dplyr::filter(tipo_muestra == "LP")

taxa_lp <- taxa[rownames(taxa) %in% rownames(covar_lp), ]

batchid_lp <- as.factor(covar_lp$Batch)

covar_lp <- covar_lp %>%
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

for (i in levels(batchid_lp)){
    print(paste0("Analizando run = ", i))
    
    for (j in 1:length(combinaciones_con_meaf)){
        print(paste0("        Combinación: ", j, " ", combinaciones_cadenas[[j]]))
        tryCatch({
            covariables <- covar_lp %>% 
                dplyr::select(combinaciones_con_meaf[[j]])
            
            
            taxa_coregida <- ConQuR(tax_tab=taxa_lp, 
                                    batchid=batchid_lp, 
                                    covariates=covariables, 
                                    batch_ref=i)
            
            saveRDS(taxa_coregida, paste0("../RDSs/taxa_corrected_", combinaciones_cadenas[j] ,"_", i, "_lp.RDS"))
            
            print(paste0("     Taxa corregida"))
            
            permanova <- PERMANOVA_R2(TAX= taxa_coregida, 
                                      batchid=batchid_lp,
                                      covariates=covar_lp[, combinaciones_con_meaf[[j]]], 
                                      key_index=1)
            
            print(paste0("     Permanova realizada"))
            
            nombre_base <- paste0(i, " ", combinaciones_con_meaf[[j]])
            
            permanovas[nombre_base,1] = permanova$tab_count[1,1]
            permanovas[nombre_base,2] = permanova$tab_count[2,1]
            permanovas[nombre_base,3] = permanova$tab_rel[1,1]
            permanovas[nombre_base,4] = permanova$tab_rel[2,1]
            
            
        }, error = function(e) {
            # Manejo del error, imprime un mensaje y continúa
            print(paste("Error en el RUN = ", i, ":", e$message))
        })
    }
}

