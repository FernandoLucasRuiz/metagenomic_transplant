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


# DESEq de FLR

deseq_flr <- function(objeto_tse, rara = F, prevalence, taxa_relevante) {
    
    variables_salida <- c("Supervivencia", "AR", "AHT", "Biliary_complications")
    
    objetos_deseq2 <- list()
    
    for (i in variables_salida) {
        
        cat("\n\n DESeq2 de: ", i, "\n")
        
        counts <- as.data.frame(assay(objeto_tse))
        
        if (rara == T) {
            
            counts <- counts + 1
            
        }
        
        metadata <- data.frame(colData(objeto_tse)) %>% 
            rownames_to_column(var = "sample")
        
        metadata$Batch <- factor(metadata$Batch)
        
        if (i == "Supervivencia") {
            metadata[[i]] <- relevel(metadata[[i]], ref = "SI")
        } else {
            metadata[[i]] <- relevel(metadata[[i]], ref = "NO")
        }
        
        dds <- DESeqDataSetFromMatrix(countData=counts, 
                                      colData=metadata, 
                                      design=as.formula(paste0("~ ", "Batch + ", i)))
        
        dds <- DESeq(dds, quiet = T)
        res <- results(dds)
        summary(res)
        
        res <- res[order(res$padj),]
        res <- res %>%
            as.data.frame() %>%
            rownames_to_column(var = "code")
        
        taxa_names <- rowData(objeto_tse) %>%
            as.data.frame() %>%
            unite("Taxonomy", Phylum:Species, sep = "_", remove = T) %>%
            rownames_to_column(var = "code") 
        # # Separar la columna en múltiples columnas por "_"
        # separate(Taxonomy, into = c("Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "_", fill = "right", extra = "merge", remove = F) %>%
        # # Crear una nueva columna que contenga el último nombre y su prefijo correspondiente según la posición
        # mutate(Name_microorganism = case_when(
        #     !is.na(Species) & Species != "" ~ paste("Species", Species, sep = "_"),
        #     !is.na(Genus) & Genus != "" ~ paste("Genus", Genus, sep = "_"),
        #     !is.na(Family) & Family != "" ~ paste("Family", Family, sep = "_"),
        #     !is.na(Order) & Order != "" ~ paste("Order", Order, sep = "_"),
        #     !is.na(Class) & Class != "" ~ paste("Class", Class, sep = "_"),
        #     !is.na(Phylum) & Phylum != "" ~ paste("Phylum", Phylum, sep = "_")
        #     )) %>%
        # dplyr::select(code, Name_microorganism)
        
        res <- dplyr::left_join(res, taxa_names, by = "code") %>%
            as.data.frame()
        
        objetos_deseq2[[i]] <- res
        
    }
    
    for (i in names(objetos_deseq2)){
        
        i_filtrado <- objetos_deseq2[[i]] %>%
            dplyr::filter(padj < 0.05 & (log2FoldChange < -1 | log2FoldChange > 1)) %>%
            dplyr::select(code, log2FoldChange, lfcSE, padj, Taxonomy)
        
        i_filtrado$variable <- i
        
        i_filtrado$prevalencia <- prevalence
        
        for (j in 1:nrow(i_filtrado)) {
            
            if (i_filtrado[j, "log2FoldChange"] > 0) {
                
                i_filtrado[j, "biosis"] <- "Hiperbiosis"
                
            } else {
                
                i_filtrado[j, "biosis"] <- "Hipobiosis"
                
            }
            
        }
        
        if (nrow(taxa_relevante) == 0){
            
            taxa_relevante <- i_filtrado
            
        } else {
            
            taxa_relevante <- rbind(taxa_relevante, i_filtrado)
            
        }
        
    }
    
    return(list(objetos_deseq2 = objetos_deseq2, taxa_relevante = taxa_relevante))
    
}


# volcano plot deseq

volcanos_deseq_flr <- function(objeto_deseq){
    
    plots <- list()
    
    for (i in names(objeto_deseq)) {
        
        volcano_data <- as.data.frame(objeto_deseq[[i]])
        
        volcano_data[is.na(volcano_data)] <- 1
        
        # Crear un factor para diferenciar entre upregulados, downregulados y no significativos
        volcano_data$Regulation <- with(volcano_data, ifelse(padj < 0.05 & log2FoldChange > 1, "Hiperbiosis",
                                                             ifelse(padj < 0.05 & log2FoldChange < -1, "Hipobiosis", "Not Significant")))
        
        p <- ggplot(volcano_data, aes(x=log2FoldChange,
                                      y=-log10(padj), color=Regulation)) +
            geom_point(aes(size = Regulation, alpha = Regulation), show.legend = TRUE) +
            scale_size_manual(values = setNames(c(2, 4, 4), c("Not Significant", "Hiperbiosis" , "Hipobiosis"))) +
            scale_alpha_manual(values = setNames(c(0.1, 0.75, 0.75), c("Not Significant", "Hiperbiosis" , "Hipobiosis"))) +
            scale_color_manual(values = setNames(c("#26828EFF", "#FDE725FF", "#440154FF"), c("Not Significant", "Hiperbiosis" , "Hipobiosis"))) +
            geom_text_repel(aes(label=ifelse(Regulation != "Not Significant", code, "")),
                            color = "black",
                            max.overlaps = 20, # Reduce el número máximo de solapamientos
                            point.padding = unit(0.2, "lines"), # Menos padding alrededor de los puntos
                            size = 2.5, # Tamaño de fuente más pequeño
                            segment.size = 0.2, # Líneas de guía más finas
                            segment.color = 'grey50',
                            #max.segment.length = unit(0.5, "lines"), # Líneas de guía más cortas
                            arrow = arrow(length = unit(0.02, "npc"), type = "closed", ends = "last")) +
            labs(title = paste0("Volcano Plot: ", i), xlab = "Log2 Fold Change", ylab = "-Log10 P-value") + 
            theme_minimal() +
            geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey") +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
            theme(
                plot.title = element_text(size = 14, hjust = 0.5),  # Centrar el título
                legend.position = "bottom",                         # Posicionar la leyenda a la derecha
                legend.justification = "center",                   # Centrar la leyenda verticalmente
                legend.key.size = unit(0.5, "cm"),
                legend.text = element_text(size = 8),
                legend.title = element_blank(),
                plot.margin = margin(10, 30, 10, 10)               # Ajustar márgenes para evitar solapamientos
            )
        
        plots[[i]] <- p
        
    }
    
    return(plots)
    
}


# volcano plot maaslin2

volcanos_maaslin_flr <- function(objeto_maslin){
    
    plots <- list()
    
    for (i in names(objeto_maslin)) {
        
        volcano_data <- as.data.frame(objeto_maslin[[i]])
        
        volcano_data[is.na(volcano_data)] <- 1
        
        # Crear un factor para diferenciar entre upregulados, downregulados y no significativos
        volcano_data$Regulation <- with(volcano_data, ifelse(pval < 0.05 & coef > 0, "Hiperbiosis",
                                                             ifelse(pval < 0.05 & coef < 0, "Hipobiosis", "Not Significant")))
        
        p <- ggplot(volcano_data, aes(x=coef,
                                      y=-log10(pval), color=Regulation)) +
            geom_point(aes(size = Regulation, alpha = Regulation), show.legend = TRUE) +
            scale_size_manual(values = setNames(c(2, 4, 4), c("Not Significant", "Hiperbiosis" , "Hipobiosis"))) +
            scale_alpha_manual(values = setNames(c(0.1, 0.75, 0.75), c("Not Significant", "Hiperbiosis" , "Hipobiosis"))) +
            scale_color_manual(values = setNames(c("#26828EFF", "#FDE725FF", "#440154FF"), c("Not Significant", "Hiperbiosis" , "Hipobiosis"))) +
            geom_text_repel(aes(label=ifelse(Regulation != "Not Significant", feature, "")),
                            color = "black",
                            max.overlaps = 20, # Reduce el número máximo de solapamientos
                            point.padding = unit(0.2, "lines"), # Menos padding alrededor de los puntos
                            size = 2.5, # Tamaño de fuente más pequeño
                            segment.size = 0.2, # Líneas de guía más finas
                            segment.color = 'grey50',
                            max.segment.length = unit(10, "lines"), # Líneas de guía más cortas
                            arrow = arrow(length = unit(0.02, "npc"), type = "closed", ends = "last")) +
            labs(title = paste0("Volcano Plot: ", i), xlab = "Coef. (by MaAsLin2)", ylab = "-Log10 P-value") + 
            theme_minimal() +
            geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
            theme(
                plot.title = element_text(size = 14, hjust = 0.5),  # Centrar el título
                legend.position = "bottom",                         # Posicionar la leyenda a la derecha
                legend.justification = "center",                   # Centrar la leyenda verticalmente
                legend.key.size = unit(0.5, "cm"),
                legend.text = element_text(size = 8),
                legend.title = element_blank(),
                plot.margin = margin(10, 30, 10, 10)               # Ajustar márgenes para evitar solapamientos
            )
        
        plots[[i]] <- p
        
    }
    
    return(plots)
    
}


# ALDEX de FLR

aldex_flr <- function(objeto_tse, prevalence, taxa_relevante) {
    
    variables_salida <- c("Supervivencia", "AR", "AHT", "Biliary_complications")
    
    objetos_aldex <- list()
    
    for (i in variables_salida) {
        
        cat("\n\n ALDEX de: ", i, "\n")
        
        selex.sub <- assay(objeto_tse)
        
        conds <- as.character(colData(objeto_tse)[[i]])
        
        x <- aldex.clr(selex.sub, conds, denom="all", verbose=F)
        x.tt <- aldex.ttest(x, hist.plot=F, paired.test=FALSE, verbose=FALSE)
        x.effect <- aldex.effect(x, CI=T, verbose=F, include.sample.summary=F, 
                                 paired.test=FALSE, glm.conds=NULL, useMC=F)
        x.all <- data.frame(x.tt,x.effect) %>%
            rownames_to_column("code")
        
        if (i == "Supervivencia") {
            
            x.all$effect <- x.all$effect * -1
            
        }
        
        objetos_aldex[[i]] <- x.all
        
    }
    
    for (i in names(objetos_aldex)){
        
        i_filtrado <- objetos_aldex[[i]] %>%
            dplyr::filter(we.eBH < 0.05) %>%
            dplyr::select(code, effect, we.eBH)
        
        if (nrow(i_filtrado) == 0 ){
            next
            
        }
        
        i_filtrado$variable <- i
        
        i_filtrado$prevalencia <- prevalence
        
        for (j in 1:nrow(i_filtrado)) {
            
            if (i_filtrado[j, "effect"] > 0) {
                
                i_filtrado[j, "biosis"] <- "Hiperbiosis"
                
            } else {
                
                i_filtrado[j, "biosis"] <- "Hipobiosis"
                
            }
            
        }
        
        if (nrow(taxa_relevante) == 0){
            
            taxa_relevante <- i_filtrado
            
        } else {
            
            taxa_relevante <- rbind(taxa_relevante, i_filtrado)
            
        }
        
    }
    
    return(list(objetos_aldex = objetos_aldex, taxa_relevante = taxa_relevante))
    
}

# Volcanos ALdex2
volcanos_aldex_flr <- function(objetos_aldex){
    
    plots <- list()
    
    for (i in names(objetos_aldex)) {
        
        volcano_data <- as.data.frame(objetos_aldex[[i]])
        
        volcano_data[is.na(volcano_data)] <- 1
        
        # Crear un factor para diferenciar entre upregulados, downregulados y no significativos
        volcano_data$Regulation <- with(volcano_data, ifelse(we.eBH < 0.05 & effect > 0, "Hiperbiosis",
                                                             ifelse(we.eBH < 0.05 & effect < 0, "Hipobiosis", "Not Significant")))
        
        p <- ggplot(volcano_data, aes(x=effect,
                                      y=-log10(we.eBH), color=Regulation)) +
            geom_point(aes(size = Regulation, alpha = Regulation), show.legend = TRUE) +
            scale_size_manual(values = setNames(c(2, 4, 4), c("Not Significant", "Hiperbiosis" , "Hipobiosis"))) +
            scale_alpha_manual(values = setNames(c(0.1, 0.75, 0.75), c("Not Significant", "Hiperbiosis" , "Hipobiosis"))) +
            scale_color_manual(values = setNames(c("#26828EFF", "#FDE725FF", "#440154FF"), c("Not Significant", "Hiperbiosis" , "Hipobiosis"))) +
            geom_text_repel(aes(label=ifelse(Regulation != "Not Significant", code, "")),
                            color = "black",
                            max.overlaps = 20, # Reduce el número máximo de solapamientos
                            point.padding = unit(0.2, "lines"), # Menos padding alrededor de los puntos
                            size = 2.5, # Tamaño de fuente más pequeño
                            segment.size = 0.2, # Líneas de guía más finas
                            segment.color = 'grey50',
                            max.segment.length = unit(10, "lines"), # Líneas de guía más cortas
                            arrow = arrow(length = unit(0.02, "npc"), type = "closed", ends = "last")) +
            labs(title = paste0("Volcano Plot: ", i), xlab = "Coef. (by MaAsLin2)", ylab = "-Log10 P-value") + 
            theme_minimal() +
            geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
            theme(
                plot.title = element_text(size = 14, hjust = 0.5),  # Centrar el título
                legend.position = "bottom",                         # Posicionar la leyenda a la derecha
                legend.justification = "center",                   # Centrar la leyenda verticalmente
                legend.key.size = unit(0.5, "cm"),
                legend.text = element_text(size = 8),
                legend.title = element_blank(),
                plot.margin = margin(10, 30, 10, 10)               # Ajustar márgenes para evitar solapamientos
            )
        
        plots[[i]] <- p
        
    }
    
    return(plots)
    
}

# ancombc2_flr de FLR

ancombc2_flr <- function(objeto_tse, prevalence, taxa_relevante) {
    
    variables_salida <- c("Supervivencia", "AR", "AHT", "Biliary_complications")
    
    objetos_ancombc <- list()
    
    for (i in variables_salida) {
        
        cat("\n\n ANCOM-BC2 de: ", i, "\n")
        
        ancombc2_out <- ancombc2(
            objeto_tse,
            assay.type = "counts",
            fix_formula = paste0(i, " + Batch"),
            p_adj_method = "fdr",
            lib_cut = 0,
            group = i,
            struc_zero = TRUE,
            neg_lb = TRUE,
            alpha = 0.05,
            # multi-group comparison is deactivated automatically
            global = TRUE,
            n_cl = 4)
        
        objetos_ancombc[[i]] <- ancombc2_out
        
    }
    
    taxa_relevante <- data.frame()
    
    for (i in names(objetos_ancombc)){
        
        i_filtrado <- objetos_ancombc[[i]]$res %>%
            dplyr::select(taxon, 3, 19) 
        
        i_filtrado <- i_filtrado[i_filtrado[3] < 0.05, ]
        
        
        if (nrow(i_filtrado) == 0 ){
            next
            
        }
        
        i_filtrado$variable <- i
        
        i_filtrado$prevalencia <- prevalence
        
        for (j in 1:nrow(i_filtrado)) {
            
            if (i_filtrado[j, 2] > 0) {
                
                i_filtrado[j, "biosis"] <- "Hiperbiosis"
                
            } else {
                
                i_filtrado[j, "biosis"] <- "Hipobiosis"
                
            }
            
        }
        
        colnames(i_filtrado) <- c("taxon", "lfc", "FDR", "variable", "prevalencia", "biosis")
        
        if (nrow(taxa_relevante) == 0){
            
            taxa_relevante <- i_filtrado
            
        } else {
            
            taxa_relevante <- rbind(taxa_relevante, i_filtrado)
            
        }
        
    }
    
    return(list(objetos_ancombc = objetos_ancombc, taxa_relevante = taxa_relevante))
    
}


# volcano_ancombc2

volcanos_ancombc_flr <- function(objetos_ancombc){
    
    plots <- list()
    
    for (i in names(objetos_ancombc)) {
        
        volcano_data <- as.data.frame(objetos_ancombc[[i]]$res) %>%
            dplyr::select(taxon, 3, 19) 
        
        colnames(volcano_data) <- c("taxon", "lfc", "FDR")
        
        volcano_data[is.na(volcano_data)] <- 1
        
        # Crear un factor para diferenciar entre upregulados, downregulados y no significativos
        volcano_data$Regulation <- with(volcano_data, ifelse(FDR < 0.05 & lfc > 0, "Hiperbiosis",
                                                             ifelse(FDR < 0.05 & lfc < 0, "Hipobiosis", "Not Significant")))
        
        p <- ggplot(volcano_data, aes(x=lfc,
                                      y=-log10(FDR), color=Regulation)) +
            geom_point(aes(size = Regulation, alpha = Regulation), show.legend = TRUE) +
            scale_size_manual(values = setNames(c(2, 4, 4), c("Not Significant", "Hiperbiosis" , "Hipobiosis"))) +
            scale_alpha_manual(values = setNames(c(0.1, 0.75, 0.75), c("Not Significant", "Hiperbiosis" , "Hipobiosis"))) +
            scale_color_manual(values = setNames(c("#26828EFF", "#FDE725FF", "#440154FF"), c("Not Significant", "Hiperbiosis" , "Hipobiosis"))) +
            geom_text_repel(aes(label=ifelse(Regulation != "Not Significant", taxon, "")),
                            color = "black",
                            max.overlaps = 20, # Reduce el número máximo de solapamientos
                            point.padding = unit(0.2, "lines"), # Menos padding alrededor de los puntos
                            size = 2.5, # Tamaño de fuente más pequeño
                            segment.size = 0.2, # Líneas de guía más finas
                            segment.color = 'grey50',
                            max.segment.length = unit(10, "lines"), # Líneas de guía más cortas
                            arrow = arrow(length = unit(0.02, "npc"), type = "closed", ends = "last")) +
            labs(title = paste0("Volcano Plot: ", i), xlab = "Coef. (by MaAsLin2)", ylab = "-Log10 P-value") + 
            theme_minimal() +
            geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
            theme(
                plot.title = element_text(size = 14, hjust = 0.5),  # Centrar el título
                legend.position = "bottom",                         # Posicionar la leyenda a la derecha
                legend.justification = "center",                   # Centrar la leyenda verticalmente
                legend.key.size = unit(0.5, "cm"),
                legend.text = element_text(size = 8),
                legend.title = element_blank(),
                plot.margin = margin(10, 30, 10, 10)               # Ajustar márgenes para evitar solapamientos
            )
        
        plots[[i]] <- p
        
    }
    
    return(plots)
    
}



