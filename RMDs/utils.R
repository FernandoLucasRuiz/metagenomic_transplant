# Convierte en un dataframe para graficar
dataframer_abundancias <- function(tse_object, abundancia_object){
    abundancia_agrupada_df <- data.frame(taxon = rownames(tse_object), abundancia = abundancia_object)
    
    rownames(abundancia_agrupada_df) <- NULL
    
    return(abundancia_agrupada_df)
}



# Ordenar taxa en orden
ordenar_taxa <- function(objecto_abundancia_df){
    abundancia_agrupada_df$taxon <- reorder(objecto_abundancia_df[["taxon"]], -objecto_abundancia_df[["abundancia"]])
    
    abundancia_agrupada_df$taxon <- factor(abundancia_agrupada_df[["taxon"]], levels = c(levels(abundancia_agrupada_df[["taxon"]])[levels(abundancia_agrupada_df[["taxon"]]) != "Other"], "Other"))
    
    return(abundancia_agrupada_df)
}




# Guardar conqur
guardar_conqur <- function(objeto_taxonomia, objeto_batchid, objeto_covar, objeto_batch, nombre_guardar){
    taxa_corrected <- ConQuR(tax_tab=objeto_taxonomia, 
                                   batchid=objeto_batchid, 
                                   covariates=objeto_covar, 
                                   batch_ref=objeto_batch, 
                                   num_core = num_cores-2)

    saveRDS(taxa_corrected, paste0("RDSs/taxa_corrected_", nombre_guardar, ".RDS"))
}





# Plotear conqur
plotear_conqur <- function(pathway_taxa_corrected, objeto_taxa, objeto_batchid){
    taxa_corrected <- readRDS(pathway_taxa_corrected)
    par(mfrow=c(2, 2))
    
    Plot_PCoA(TAX=objeto_taxa, factor=objeto_batchid, main="Before Correction, Bray-Curtis")
    Plot_PCoA(TAX=taxa_corrected, factor=objeto_batchid, main="ConQuR (Default), Bray-Curtis")
    
    Plot_PCoA(TAX=objeto_taxa, factor=objeto_batchid, dissimilarity="Aitch", main="Before Correction, Aitchison")
    Plot_PCoA(TAX=taxa_corrected, factor=objeto_batchid, dissimilarity="Aitch", main="ConQuR (Default), Aitchison")
}






# Prepare data for TSE

 prepare_TSE_data <- function(OTU_table, metadata){
    abundance_table <- OTU_table %>%
        dplyr::group_by(random_ids, codigo) %>%
        dplyr::summarize(Count = sum(as.numeric(Count), na.rm = TRUE)) %>%
        ungroup() %>%
        pivot_wider(names_from = codigo, values_from = Count) %>%
        mutate_all(~replace_na(., 0))
    
    abundance_table <- as.data.frame(abundance_table)
    rownames(abundance_table) <- abundance_table$random_ids
    abundance_table <- abundance_table %>% dplyr::select(-random_ids)
    
    # Crear rowdata(tax)
    tax <- OTU_table %>%
        dplyr::select(Phylum:Species, random_ids) %>%
        as.data.frame()
    
    tax <- tax[!duplicated(tax$random_ids), ]
    
    rownames(tax) <- tax$random_ids
    tax <- dplyr::select(tax, -random_ids)
    
    
    if ("Batch" %in% colnames(OTU_table)){
        samples <- OTU_table %>%
            dplyr::select(codigo, nombre_muestra, Batch, tipo_muestra) %>%
            as.data.frame()
    } else {
        # Creas colData (samples). Esto incluye la metadata
        samples <- OTU_table %>%
            dplyr::select(codigo, nombre_muestra, tipo_muestra) %>%
            as.data.frame()
    }
    
    
    samples <- samples[!duplicated(samples$codigo), ]
    
    samples_metadata <- dplyr::left_join(samples, metadata, by= "nombre_muestra" )
    
    rownames(samples_metadata) <- samples_metadata$codigo
    samples_metadata <- dplyr::select(samples_metadata, -codigo)
    
    
    # Match rows and columns
    abundance_table <- abundance_table[rownames(tax), rownames(samples_metadata)]
    
    # Let's ensure that the data is in correct (numeric matrix) format:
    abundance_table <- as.matrix(abundance_table)
    
    return(list(abundance_table = abundance_table, tax = tax, samples_metadata = samples_metadata))
    
}




# Chequeo preTSE
check_pre_TSE <- function(abundance_table, tax_table, samples_metadata_table){
    # Match rows and columns
    counts <- abundance_table[rownames(tax_table), rownames(samples_metadata_table)]
    
    # Let's ensure that the data is in correct (numeric matrix) format:
    counts <- as.matrix(counts)
    
    # coldata rownames match assay colnames
    samples_metadata_table <- samples_metadata_table[match(colnames(counts), rownames(samples_metadata_table)), ]
    
    
    has_printed <- FALSE
    
    if (all(rownames(samples_metadata_table) == colnames(counts)) != T) {
        print("rownames(samples_metadata_table) != colnames(counts)")
        has_printed <- TRUE 
    } 

    if (class(samples_metadata_table) != "data.frame"){
        print("samples_metadata_table to dataframe" )
        samples_metadata_table <- data.frame(samples_metadata_table)
        has_printed <- TRUE 
    }
    
    if (all(rownames(tax_table) == rownames(counts)) == F){
        print("rownames(tax_table) != rownames(counts)")
        has_printed <- TRUE 
    }
    
    if (class(tax_table) != "data.frame") {
        print("tax_table to dataframe")
        tax_table <- data.frame(tax_table)
        has_printed <- TRUE
    }
    
    if (!has_printed) {
        cat("Everything's OK")
        
        return(list(abundance = counts, tax = tax_table, samples_metadata = samples_metadata_table))
    }
}



# Abundancias relativas totales

abundancias_relativas_totales <- function(lista_objetos_tse, ranking, top = 10){
    
    if (is_null(names(lista_objetos_tse)) == T){
        stop("Pon nombres en la lista de TSEs")
        
    }
    
    abundancias_relativas <- data.frame(
        Taxon = character(),               
        Relative_Abundance = numeric(),      
        TSE = character())
    
    for (i in 1:length(lista_objetos_tse)){
        tse_ranking <- agglomerateByRank(lista_objetos_tse[[i]], rank = ranking)
        
        # Verificar si el assay "counts" está disponible
        if ("counts" %in% assayNames(tse_ranking)) {
            row_sums <- assay(tse_ranking, "counts") %>%
                as.data.frame() %>%
                rowSums()
            
        } else if ("relative_abundance" %in% assayNames(tse_ranking)) {
            row_sums <- assay(tse_ranking, "relative_abundance") %>%
                as.data.frame() %>%
                rowSums()
            
        } else {
            stop("Ni 'counts' ni 'relative_abundance' están disponibles en este objeto.")
            
        }
        
        total_sum <- sum(row_sums)
        
        relative_abundance <- row_sums / total_sum
        
        data <- data.frame(
            Taxon = rownames(assay(tse_ranking)),
            Relative_Abundance = relative_abundance
        )
        
        data$TSE <- names(lista_objetos_tse)[i]
        
        rownames(data) <- NULL
        
        if (nrow(abundancias_relativas) == 0){
            
            abundancias_relativas <- data
            
        } else {
            abundancias_relativas <- rbind(abundancias_relativas, data)
            
        }
        
    }
    
    top_taxones <- abundancias_relativas %>%
        group_by(TSE) %>%
        arrange(desc(Relative_Abundance)) %>%
        slice_head(n = top) %>%
        pull(Taxon)
    
    # Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
    taxa_renamed <- lapply(abundancias_relativas$Taxon, function(x){
        if (x %in% top_taxones) {x} else {"Other"}
    })
    
    abundancias_relativas$taxa_sub <- as.character(taxa_renamed)
    
    return(abundancias_relativas)
    
}



# Plotear abundancias relativas en pie chart

plotear_abundancia_relativa <- function(objeto_abundancias, titulo) {
    
    p <- ggplot(objeto_abundancias, aes(x = 1, y = Relative_Abundance, fill= taxa_sub)) +
        geom_bar(stat = "identity") +
        scale_fill_paletteer_d("ggthemes::Tableau_20") + 
        labs(title = paste0(titulo, " Relative Abundance"), x = "Taxons", y = "Relative Abundance") +
        theme_minimal() +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              legend.title = element_blank(),
              legend.position = "bottom",  # Posiciona la leyenda abajo
              legend.spacing.x = unit(0.5, 'cm'),  # Ajusta espacio entre elementos de la leyenda
              plot.title = element_text(hjust = 0.5, vjust = 10, size = 16),  # Sube el título y ajusta el tamaño
              #plot.margin = margin(t = 20)  # Añade margen superior para el título
        ) +
        coord_polar("y", start = 0) +
        facet_wrap(~TSE, nrow = 1) 
    
    return(p)
    
}


#plotear_diferencias_relabundance
plotear_diferencias_relabundance <- function(objeto_tse, taxa, covariable) {
    
    tse_filtrado <- agglomerateByRank(objeto_tse, taxa)
    
    relabundance <- as.data.frame(t(assay(tse_filtrado, "relabundance"))) %>% 
        rownames_to_column(var = "sample")
    
    covariates <- colData(tse_filtrado) %>%
        as.data.frame() %>%
        rownames_to_column(var = "sample")
    
    p <- dplyr::left_join(relabundance, covariates, by = "sample") %>%
        dplyr::select(intersect(names(.), names(relabundance)), all_of(covariable), -sample) %>%
        pivot_longer(-all_of(covariable), names_to = "Taxon", values_to = "Abundance") %>%
        ggplot(aes_string(x=covariable, y="Abundance", fill = "Taxon")) +
        geom_boxplot() +
        geom_jitter(aes_string(x=covariable, y="Abundance"), size = 0.5) +
        labs(title = paste("Relative Abundance at", taxa, "level"),
             x = covariable, 
             y = "Relative Abundance") +
        theme(
            axis.title.y = element_blank(),
            legend.title = element_blank(),
            legend.position = "none"
        ) +
        scale_y_log10() +
        facet_wrap(~ reorder(Taxon, -Abundance))
    
    return(p)
}






