
source("RMDs/Librerias.R")
source("RMDs//utils.R")
load("Rdatas/lectura_conteos.RData")
load("Rdatas/covariables_limpieza.RData")


# Colapsar taxonomia para random ids

taxonomia_combinada_lp_t1_t3 <- otus_lp_t1_t3 %>%
    unite("Combined", Phylum:Species, sep = "_", remove = F)
taxonomia_combinada_bilis <- otus_bilis %>%
    unite("Combined", Phylum:Species, sep = "_", remove = F)
taxonomia_combinada_negative_controls <- otus_negative_controls %>%
    unite("Combined", Phylum:Species, sep = "_", remove = F)


taxonomia_combinada <- rbind(taxonomia_combinada_lp_t1_t3[,1:6], taxonomia_combinada_bilis[,1:6], taxonomia_combinada_negative_controls[,1:6]) 

taxonos <- unique(taxonomia_combinada$Combined)
set.seed(123)  # Para reproducibilidad
random_ids <- replicate(length(taxonos), {
    first_letter <- sample(letters, 1)  # Seleccionar una letra para el primer carácter
    rest <- sample(c(letters, 0:9), 15, replace = TRUE)  # Generar el resto del ID
    paste(c(first_letter, rest), collapse = "")  # Combinar la letra con el resto
})

taxonos_ids <- data.frame(cbind(taxonos, random_ids))

# juntar taxonos_ids a todas_muestras 
otus_lp_t1_t3_combined <- dplyr::left_join(taxonomia_combinada_lp_t1_t3, taxonos_ids, by = c("Combined" = "taxonos"))
otus_bilis_combined <- dplyr::left_join(taxonomia_combinada_bilis, taxonos_ids, by = c("Combined" = "taxonos"))
otus_negative_controls_combined <- dplyr::left_join(taxonomia_combinada_negative_controls, taxonos_ids, by = c("Combined" = "taxonos"))

# Controles negativos
abundance_table_negative_control <- otus_negative_controls_combined %>%
    group_by(random_ids, codigo) %>%
    summarize(Count = sum(as.numeric(Count), na.rm = TRUE)) %>%
    ungroup() %>%
    pivot_wider(names_from = codigo, values_from = Count) %>%
    mutate_all(~replace_na(., 0))

abundance_table_negative_control <- as.data.frame(abundance_table_negative_control)
rownames(abundance_table_negative_control) <- abundance_table_negative_control$random_ids
abundance_table_negative_control <- abundance_table_negative_control %>% select(-random_ids)


# LP

## Preparar datos

samples_LP <- otus_lp_t1_t3_combined %>%
    dplyr::filter(tipo_muestra == "LP") 

pretse_LP <- prepare_TSE_data(samples_LP, covs)


## Quitar contaminación

otu_combined <- merge(as.data.frame(pretse_LP$abundance_table), as.data.frame(abundance_table_negative_control), by = "row.names", all = TRUE)
rownames(otu_combined) <- otu_combined$Row.names
otu_combined$Row.names <- NULL

# Reemplazar los valores NA por ceros en la matriz combinada
otu_combined[is.na(otu_combined)] <- 0

# Crear un vector que indique qué columnas son controles (TRUE = control, FALSE = muestra real)
is_control <- c(rep(FALSE, ncol(as.data.frame(pretse_LP$abundance_table))), rep(TRUE, ncol(as.data.frame(abundance_table_negative_control))))

# Usar decontam para identificar contaminantes basados en prevalencia
contaminants_prev <- isContaminant(as.matrix(otu_combined), method = "prevalence", neg = is_control) %>%
    filter(contaminant == TRUE & p < 0.05)

head(contaminants_prev)


if (nrow(contaminants_prev) != 0){
    pretse_LP$abundance_table <- pretse_LP$abundance_table %>% 
        as.data.frame() %>%
        select(-rownames(contaminants_prev))
    
    pretse_LP$samples_metadata <- pretse_LP$samples_metadata[rownames(pretse_LP$samples_metadata) != rownames(contaminants_prev),] %>% na.omit()
}



## Profundidad de muestreo

# Extraer la matriz de abundancias del objeto TSE
abundances <- data.frame(t(pretse_LP$abundance_table))

# Calcular la profundidad mínima
abundances$sample <- rownames(abundances)

abundances_long <- pivot_longer(
    abundances, 
    cols = -sample,  # Especificamos que las columnas de especies deben ser pivotadas
    names_to = "Taxa",           # La nueva columna que contendrá los nombres de las especies
    values_to = "Abundance"        # La nueva columna que contendrá los valores de abundancia
)


min_n_seqs <- abundances_long %>% 
    group_by(sample) %>%
    summarize(n_seqs = sum(Abundance),
              min = min(n_seqs)) %>%
    pull(min)

quantil1_n_seqs <- abundances_long %>% 
    group_by(sample) %>%
    summarize(n_seqs = sum(Abundance)) %>%
    pull(n_seqs) %>%            
    quantile(probs = 0.25)

abundances <- abundances %>% 
    dplyr::select(-sample)


# Aplicar rarefacción real a las abundancias
rarefied_abundances_LP <- rrarefy(abundances, sample=quantil1_n_seqs)


## Quitar batch effect

### Preparar datos

# taxa tiene que tener la taxonomia en las columnas y los samples en las filas
taxa_LP <- data.frame(t(pretse_LP$abundance_table))

# Para sacar la lista de los batch primero hay que ordenar los samples en el samples_metadata

## Asegurarse de que ambos dataframes tengan los mismos rownames
common_rows <- intersect(rownames(taxa_LP), rownames(pretse_LP$samples_metadata))

taxa_LP <- taxa_LP[common_rows, ]
samples_metadata_LP <- pretse_LP$samples_metadata[match(rownames(taxa_LP), rownames(pretse_LP$samples_metadata)), ]

covar_LP <- samples_metadata_LP 

taxa_LP <- taxa_LP[rownames(taxa_LP) %in% rownames(covar_LP), ]

batchid_LP <- as.factor(covar_LP$Batch)

covar_LP <- covar_LP %>%
    dplyr::select(-c("nombre_muestra", "Batch", "tipo_muestra"))

covar_salida <- covar_LP %>% 
    dplyr::select(MEAF_score, Supervivencia, AR, AHT, Biliary_complications)



result_tuned <- Tune_ConQuR(tax_tab=taxa_LP,
                            batchid=batchid_LP,
                            covariates=covar_salida,
                            batch_ref_pool=c("RUN 1", "RUN 2", "RUN 3"),
                            logistic_lasso_pool=c(T, F),
                            quantile_type_pool=c("standard", "lasso"),
                            simple_match_pool=c(T, F),
                            lambda_quantile_pool=c(NA, "2p/n", "2p/logn"),
                            interplt_pool=c(T, F),
                            frequencyL=0,
                            frequencyU=1,
                            num_core = num_cores-2)

saveRDS(result_tuned, "RDSs/result_tuned.RDS")

