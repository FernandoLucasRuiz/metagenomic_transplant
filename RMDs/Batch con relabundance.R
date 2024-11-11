
print(Sys.time())

source("RMDs/Librerias.R")
source("RMDs/utils.R")
otus_lp_t1_t3 <- readRDS("RDSs/otus_lp_t1_t3.RDS")
otus_bilis <- readRDS("RDSs/otus_bilis.RDS")
otus_negative_controls <- readRDS("RDSs/otus_negative_controls.RDS")
covs <- readRDS("RDSs/covs.RDS")


# Preparar datos

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


samples_LP <- otus_lp_t1_t3_combined %>%
    dplyr::filter(tipo_muestra == "LP") 

pretse_LP <- prepare_TSE_data(samples_LP, covs)


# Create a TreeSE
tse <- TreeSummarizedExperiment(
    assays =  SimpleList(counts = pretse_LP$abundance_table),
    colData = DataFrame(pretse_LP$samples_metadata),
    rowData = DataFrame(pretse_LP$tax))

tse <- transformAssay(tse, assay.type = "counts", method = "relabundance")


# taxa tiene que tener la taxonomia en las columnas y los samples en las filas
taxa_LP <- data.frame(t(assay(tse, "relabundance")))

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


### Fine tunning automático

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

saveRDS(result_tuned, "RDSs/result_tuned_relabundance.RDS")

print(Sys.time())
