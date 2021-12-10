
#' ADT-LR correlation pipeline Helper/Loader Function
#' @param subdir subdirectory in dir
#' @param dir appropriate formated citeseq directory
#' @inheritParams wrap_adt_corr
#'
#' @returns a list with two elements:
#' - a tibble with summarized ADT-LR correlations
#' - a tibble with summarized ADT
#'
#' @details Simply loads the needed seurat and liana objects and calls the
#' `wrap_adt_corr` function/pipeline.
#' `seurat_object_path` = seurat_object.RDS path
#' `liana_res_path` = liana_res.RDS path
run_adt_pipe <- function(subdir = subdir,
                         dir,
                         sobj_pattern = "_seurat.RDS",
                         liana_pattern,
                         op_resource,
                         cluster_key,
                         arbitrary_thresh,
                         organism = "human",
                         adt_pipe_type = "correlation",
                         .eval = "max",
                         ...){

    # read seurat
    seurat_object_path <- list.subfiles(subdir = subdir,
                                        dir = citeseq_dir,
                                        pattern = sobj_pattern)
    seurat_object <- readRDS(seurat_object_path)


    # read liana res
    liana_res_path <- list.subfiles(subdir = subdir,
                                    dir = citeseq_dir,
                                    pattern = liana_pattern)
    liana_res <- readRDS(liana_res_path) %>%
        liana_aggregate_enh(...,
                            .eval = .eval
                            )



    # Get Symbols of Receptors from OP
    receptor_syms <- str_to_symbol(op_resource$target_genesymbol, organism)
    # Get ADT Symbols
    adt_symbols <- rownames(seurat_object@assays$ADT)

    # Obtain adt genesymbols
    if(organism == "mouse"){
        # obtain both mouse and human but convert to title
        adt_aliases <- get_adt_aliases(adt_symbols = adt_symbols,
                                       organism = "mouse") %>%
            bind_rows(get_adt_aliases(adt_symbols = adt_symbols,
                                      organism = "human")) %>%
            mutate(across(everything(), str_to_title)) %>%
            distinct()
    } else{
        adt_aliases <- get_adt_aliases(adt_symbols = adt_symbols,
                                       organism = organism)
    }

    if(adt_pipe_type == "correlation"){
        message(str_glue("Calculating correlations on: {subdir}"))
        # call function that does the adt-lr correlation
        adt_lr_corr <- wrap_adt_corr(seurat_object = seurat_object,
                                     liana_res = liana_res,
                                     op_resource = op_resource,
                                     cluster_key = cluster_key,
                                     organism = organism,
                                     receptor_syms = receptor_syms,
                                     adt_symbols = adt_symbols,
                                     adt_aliases = adt_aliases)

    } else if(adt_pipe_type == "specificity"){
        message(str_glue("Calculating ROC summary for {subdir}"))

        adt_lr_roc <-
            generate_specificity_roc(
                seurat_object = readRDS(seurat_object_path),
                liana_res = liana_res,
                op_resource = op_resource,
                cluster_key = cluster_key,
                organism=organism,
                receptor_syms = receptor_syms,
                adt_symbols = adt_symbols,
                adt_aliases = adt_aliases,
                arbitrary_thresh = arbitrary_thresh)
    }
}



#' ADT-LR correlation pipeline wrapper function
#' @param seurat_object seurat_object with RNA and ADT assays
#' @param op_resource omnipath-formatted resource
#' @param cluster_key name of the vector in the seurat metadata
#' @param liana_res liana results obtained for the same dataset as above
#'
#' @returns a tibble with summarized ADT-LR correlations
wrap_adt_corr <- function(seurat_object,
                          liana_res,
                          op_resource,
                          cluster_key,
                          organism = "human",
                          receptor_syms,
                          adt_symbols,
                          adt_aliases){

    # Get ADT means
    adt_means <- get_adt_means(seurat_object = seurat_object,
                               liana_res = liana_res,
                               op_resource = op_resource,
                               cluster_key = cluster_key,
                               organism = organism,
                               receptor_syms = receptor_syms,
                               adt_symbols = adt_symbols,
                               adt_aliases = adt_aliases)

    # Join to LIANA and format
    liana_adt <- liana_format_adt(liana_res = liana_res,
                                  adt_means = adt_means)


    # Get correlations between RNA-LR scores and ADT means
    message("Calculating ADT-LR Correlations")
    adt_corr <- get_adt_correlations(liana_adt)

    # Get correlation between ADT-RNA assays (receptors alone)
    message("Calculating RNA-ADT Correlations")
    assays_corr <- get_assays_corr(seurat_object,
                                   adt_means = adt_means,
                                   adt_aliases = adt_aliases,
                                   cluster_key = cluster_key,
                                   receptor_syms = receptor_syms,
                                   organism = organism)

    # format and count number of ADTs used in correlation
    adt_corr %<>% bind_rows(assays_corr) %>%
        mutate(metric = gsub("\\_.*","", metric)) %>%
        mutate(n_genes = length(unique(adt_means$entity_symbol)))

    return(adt_corr)
}



#' Obtain ADT and RNA correlation
#' @param seurat_object seurat object with RNA assay
#' @param adt_means mean adt abundance per cluster, obtained via `get_adt_means`
#' @param adt_aliases as obtained from `get_adt_aliases`
#' @param cluster_key seurat cluster variable name
#'
#' @returns 2 row tibble, a row for correlation mean and scaled RNA-ADT
get_assays_corr <- function(seurat_object,
                            adt_means,
                            adt_aliases,
                            cluster_key = "seurat_clusters",
                            receptor_syms,
                            organism){

    sce_rna <- SingleCellExperiment::SingleCellExperiment(
        assays=list(
            counts = GetAssayData(seurat_object, assay = "RNA", slot = "counts"),
            data = GetAssayData(seurat_object, assay = "RNA", slot = "data")
        ),
        colData=DataFrame(label=seurat_object@meta.data[[cluster_key]])
    )

    # called adt_means but simply runs means on sce object and formats the output
    mean_summary <- scuttle::summarizeAssayByGroup(sce_rna,
                                                   ids = colLabels(sce_rna),
                                                   assay.type = "data",
                                                   statistics = c("mean"))
    rna_means <- mean_summary@assays@data$mean %>%
        # gene mean across cell types
        as_tibble(rownames = "entity_symbol") %>%
        mutate(entity_symbol = str_to_symbol(entity_symbol, organism)) %>%
        # keep only receptors present in OP
        filter(entity_symbol %in% receptor_syms) %>%
        pivot_longer(-entity_symbol,
                     names_to = "target",
                     values_to = "rna_mean") %>%
        group_by(entity_symbol) %>%
        # obtain across cluster scaled means
        mutate(rna_scale = scale(rna_mean)) %>%
        unnest(rna_scale) %>%
        ungroup()

    rna_adt_join <- adt_aliases %>%
        left_join(adt_means, by = c("alias_symbol"="entity_symbol")) %>%
        left_join(rna_means, by = c("alias_symbol" = "entity_symbol", "target")) %>%
        na.omit()

    mean_corr <- cor.test(rna_adt_join$adt_mean,
                          rna_adt_join$rna_mean,
                          method = "kendal") %>%
        broom::tidy() %>%
        mutate(metric = "mean")

    scale_corr <- cor.test(rna_adt_join$adt_scale,
                           rna_adt_join$rna_scale,
                           method = "kendal") %>%
        broom::tidy() %>%
        mutate(metric = "scale")

    adt_rna_corr <- bind_rows(mean_corr, scale_corr) %>%
        mutate(method = "RNA-ADT") %>%
        mutate(significant = if_else(p.value <= 0.05, TRUE, FALSE)) %>%
        dplyr::select(method, estimate, metric, significant)
    return(adt_rna_corr)
}

#' Correlation Helper Function
#' @param df nested df
#' @param var name of the variable
#'
#' @return cor.test result as tibble
corr_fun <- function(df, var, method = "kendall"){
    cor.test(df[[var]],
             df$value,
             method=method,
             exact = FALSE
    ) %>% broom::tidy()
}


#' Obtain ADT genes alias table
#' @returns a table with gene symbol alaiases
get_alias_table <- function(organism="human"){

    # first open the database connection
    if(organism=="human"){
        require(org.Hs.eg.db)
        dbCon <- org.Hs.eg_dbconn()
    } else if(organism=="mouse"){
        require(org.Mm.eg.db)
        dbCon <- org.Mm.eg_dbconn()
    }

    # SQL query to get alias table and gene_info table
    sql_query <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'

    # execute the query on the database
    alias_table <- DBI::dbGetQuery(dbCon, sql_query) %>%
        dplyr::select(-c("_id", "gene_name")) %>%
        dplyr::select(-c("_id")) %>% # present twice for some reason
        as_tibble()


    return(alias_table)
}


#' Helper function to load 10x citeseq matrix, cluster, and save to file
#' @param dir directory of all citeseq datasets
#' @param subdir subdirectories with h5 files
#' @param define the pattern to be loaded - i.e. .h5 for the 10x matrices
#'
#' @details saves the newly created seurat object in the appropriate subdir
load_and_cluster <- function(dir, subdir, pattern = ".h5"){
    message(str_glue("Loading: ",
                     list.subfiles(subdir = subdir,
                                   dir = dir,
                                   pattern = pattern)
    ))

    # Load matrix
    mats <- Seurat:::Read10X_h5(list.subfiles(subdir = subdir,
                                              dir = dir,
                                              pattern = pattern))
    # remove adt_prefix and filter controls
    rownames(mats$`Antibody Capture`) %<>%
        gsub("\\_.*","", .) %>%
        .[!grepl("control", .)]

    # Convert to Seurat object and run default analysis
    # (into function with silhuette score optimization)
    seurat_object <- SeuratObject::CreateSeuratObject(counts = mats$`Gene Expression`)

    set.seed(1234)
    seurat_object[["ADT"]] <- CreateAssayObject(counts = mats$`Antibody Capture`)
    seurat_object %<>%
        FindVariableFeatures(verbose = FALSE) %>%
        NormalizeData(verbose = FALSE) %>%
        ScaleData() %>%
        RunPCA(verbose = FALSE) %>%
        FindNeighbors(reduction = "pca") %>%
        # need to optimize silhuette scores, and save the 'best' to seurat_clusters
        FindClusters(resolution = 0.4, verbose = FALSE)

    # Change cluster names (CellChat does not allow 0s...)
    # Squidpy splits clusters at _
    seurat_object@meta.data %<>%
        mutate(seurat_clusters = str_glue("cluster.{seurat_clusters}")) %>%
        mutate(seurat_clusters = as.factor(seurat_clusters))
    Seurat::Idents(seurat_object) <- seurat_object@meta.data$seurat_clusters

    # Normalize ADT
    seurat_object <- NormalizeData(seurat_object,
                                   assay = "ADT",
                                   normalization.method = "CLR")

    # save the object
    saveRDS(seurat_object,
            file.path(dir, subdir, str_glue(subdir,
                                            "seurat.RDS",
                                            .sep="_")))
}

#' helper function to list files in a subdirectory
#' @param subdir subdirectory in which to look for the files
#' @param dir directory in which the subdir is located
#' @inheritParams list.files
#'
#' @returns returns the full path and the file name
list.subfiles <- function(subdir, dir, pattern = ".h5"){
    file_p <- list.files(file.path(dir, subdir), pattern = pattern)
    full_p <- file.path(dir, subdir, file_p)

    return(full_p)
}



#' run LIANA on newly created Seurat files
#' @param subdir subdirectory
#' @param dir directory of citeseq input and output
#' @param expr_prop proportion expressed
#' @param method liana wrap method
#'
#' @inheritDotParams liana_aggregate
wrap_liana_wrap <- function(subdir,
                            dir,
                            expr_prop,
                            organism = "human",
                            ...){
    message(str_glue("Loading: ",
                     list.subfiles(subdir = subdir,
                                   dir = dir,
                                   pattern = "_seurat.RDS")
    ))

    # load object
    seurat_object <- readRDS(list.subfiles(subdir = subdir,
                                           dir = dir,
                                           pattern = "_seurat.RDS"))

    # run liana
    liana_res <- liana_wrap(seurat_object,
                            expr_prop = expr_prop,
                            ...)
    # Format Method Results
    liana_res %<>%
        # standardize symbols
        map(function(res) res %>%
                mutate(across(c(ligand, receptor),
                              ~str_to_symbol(string = .x,
                                             organism=organism)))
            )

    message(str_glue("Saving liana_res to: ",
                     file.path(dir, subdir,
                               str_glue(subdir,
                                        "liana_res",
                                        expr_prop,
                                        ".RDS",
                                        .sep="_"))))

    # save results
    saveRDS(liana_res,
            file.path(dir, subdir, str_glue(subdir,
                                            "liana_res-{expr_prop}.RDS",
                                            .sep="-")))

}




#' Function to obtain all aliases for the ADT assay of the Seurat object
#' @param adt_symbols ADT names - i.e. rownames(sce or seurat_object)
#' @returns A tibble with ADT names and associated gene aliases
get_adt_aliases <- function(adt_symbols,
                            organism = "human"){

    # obtain aliases
    alias_table <- get_alias_table(organism)

    # check adts that are not in the alias_table
    missing_adts <- adt_symbols[!(str_to_symbol(adt_symbols, organism) %in%
                                      str_to_symbol(alias_table$alias_symbol, organism))]
    message(str_glue("Mismatched (to be taken from manual annotations): ",
                     glue::glue_collapse(missing_adts, sep = ", ")))

    # Manually check mismatches
    # i.e. check if any alias is present in OmniPath
    # some e.g. CD3 are proteins with multiple gene subunts

    if(organism=="human"){
        manual_list <- list(
            # aliases obtained from GeneCards
            "CD3" = c("CD3D", "CD3E", "CD3E", "CD3G", "CD3Z"),
            "CD45RA" = "PTPRC",
            "CD45RO" = "PTPRC",
            "HLA-DR" = "CD74"
        )
    } else if(organism=="mouse"){
        manual_list <- list(
            CD105 = c("Eng", "Endo"),
            CD107a = c("Lamp1", "Perk"),
            CD120b = c("Tnfrsf1b", "Tnfr2"),
            CD16.32 = "Fcgr3",
            CD198 = "Cxcr8",
            CD199 = c("Ccr9", "Cmkbr10"),
            CD201 = c("Procr", "Epcr"),
            CD21.CD35 = c("Cr2","Cr1"),
            CD278.1 = c("Icos, Ailim"),
            CD300c.d = c("Cd300c", "Clm6"),
            CD301a = c("Mgl", "Mgl1"),
            CD301b= c("Mgl2"),
            CD309.1 = c("Kdr","Flk1"),
            CD326 = c("Epcam", "Tacstd1"),
            CD34.1 = "CD34",
            CD45.1 = c("Ptprc"),
            CD45.2 = c("Ptprc"),
            CD45R.B220 = c("Ptprc"),
            CD49a = c("Itga1"),
            CD90.1 = c("Thy1", "Thy-1"),
            CD90.2 = c("Thy1", "Thy-1"),
            D62E = c("Sele", "Elam-1"),
            F4.80 = c("Adgre4", "Emr4"),
            FceRIa = c("Fcer1g", "Fce1g"),
            FolateReceptorb = c("Folr2", "Fbp2", "Folbp2"),
            IL.21Receptor = c("Il21r", "Nilr"),
            IL.33Ra = c("Il1rl1", "Ly84", "St2", "Ste2"),
            Ly.49A = c("Klra1", "Ly-49", "Ly-49a", "Ly49", "Ly49A"),
            Ly.6A.E = c("Sca-1", "Ly6e"),
            Ly.6G = "Ly6g6e",
            Ly.6G.Ly.6C = "Gr-1",
            MAdCAM.1 = "Madcam1",
            Mac.2 = "Lgals3",
            NK.1.1 = c("Klrb1c", "Ly55c", "Nkrp1c"),
            P2X7R = "P2rx7",
            # leave out T-cell Receptor variants
            # TCRVb13.1 = "Trbv13",
            # TCRVb5.1,
            # TCRVb8.1,
            # TCRVr1.1.Cr4,
            # TCRVr2,
            # TCRVr3,
            # TCRbchain,
            # TCRr.d,
            Tim.4 = c("Timd4", "Tim4"),
            anti.P2RY12 = c("P2ry12", "P2y12")#,
            # integrinb7,
            # B6,
            # B6.1
        )
    }

    adt_manual <- manual_list %>%
        enframe(name = "adt_symbol",
                value = "alias_symbol") %>%
        unnest(alias_symbol) %>%
        # only keep the ones that are in the object
        filter(!(alias_symbol %in% missing_adts)) %>%
        mutate(across(.cols = everything(), ~str_to_symbol(string = .x,organism=organism)))


    #' convert ADTs to aliases
    #' @param adt_names names to be queried (i.e. sce/seurat rownames)
    adt_match <- map(adt_symbols,
                     function(name){
                         alias_table %>%
                             filter(alias_symbol==str_to_symbol(name, organism)) %>%
                             dplyr::rename(adt_symbol = alias_symbol,
                                           alias_symbol = symbol)
                     }) %>% bind_rows()

    # expand to all aliases
    adt_aliases <- bind_rows(adt_manual,
                             adt_match)

    return(adt_aliases)
}




#' Obtain formatted across cluster means
#' @param seurat_object Seurat object with an active ADT assay
#' @param liana_res liana results as tibble
#' @param op_resource tibble with omnipath res
#' @param cluster_key name of the cluster key
#' @param organism organism - 'human' or 'mouse'
#' @param receptor_syms Receptor symbols obtained from OmniPath
#'
#' @returns a tibble with adt_means and adt_scale per cluster
get_adt_means <- function(seurat_object,
                          liana_res,
                          op_resource,
                          cluster_key,
                          organism,
                          receptor_syms,
                          adt_symbols,
                          adt_aliases){

    # convert to singlecell object
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays=list(counts = GetAssayData(seurat_object, assay = "ADT", slot = "counts"),
                    data = GetAssayData(seurat_object, assay = "ADT", slot = "data")),
        colData=DataFrame(label=seurat_object@meta.data[[cluster_key]])
    )

    # check if all matched
    mismatched_adts <- adt_symbols[!(str_to_symbol(adt_symbols, organism) %in%
                                         str_to_symbol(adt_aliases$adt_symbol, organism))]
    message("Missing (ideally none): ",
            glue::glue_collapse(mismatched_adts, sep = ", "))



    # Get ADT stats and bind to LIANA res
    message("Calculating ADT Means")
    mean_summary <- scuttle::summarizeAssayByGroup(sce,
                                                   ids = colLabels(sce),
                                                   assay.type = "data",
                                                   statistics = c("mean"))

    adt_means <- mean_summary@assays@data$mean %>%
        # gene mean across cell types
        as_tibble(rownames = "entity_symbol") %>%
        rowwise() %>%
        mutate(entity_symbol = str_to_symbol(entity_symbol, organism)) %>%
        # Join alias names and rename entity symbol to newly-obtained alias symbol
        # purpose being to be able to join to the gene symbols in the RNA assay
        left_join(adt_aliases, by = c("entity_symbol"="adt_symbol")) %>%
        dplyr::select(-c(entity_symbol)) %>%
        dplyr::select(entity_symbol = alias_symbol, everything()) %>%
        # keep only receptors present in OP
        filter(entity_symbol %in% receptor_syms) %>%
        pivot_longer(-entity_symbol,
                     names_to = "target",
                     values_to = "adt_mean") %>%
        group_by(entity_symbol) %>%
        # obtain across cluster scaled means
        mutate(adt_scale = scale(adt_mean)) %>%
        unnest(adt_scale) %>%
        ungroup()

    return(adt_means)

}


#' Format LIANA to ADT correlation pipeline
#' @param liana_res aggregated LIANA res table
#' @param adt_means adt means summary formatted as from `get_adt_means()`
#'
#' @returns a tibble with liana results for each receptor formatted according to
#' the ADT correlation pipeline
liana_format_adt <- function(liana_res, adt_means){
    liana_res %>%
        filter(receptor %in% adt_means$entity_symbol) %>%
        dplyr::select(source, ligand, target, receptor, ends_with("rank")) %>%
        pivot_longer(c(ligand,receptor),
                     names_to = "type",
                     values_to = "entity_symbol") %>%
        filter(type == "receptor") %>%
        distinct() %>%
        left_join(adt_means, by = c("target", "entity_symbol")) %>%
        dplyr::select(ends_with("rank"), starts_with("adt")) %>%
        pivot_longer(-c(adt_scale, adt_mean),
                     names_to = "method_name")
}

#' Obtain Correlations between LR-scores and ADT means
#'
#' @param liana_adt liana results joined to ADT summaries, as obtained from
#' `liana_format_adt()`
#'
#' @returns A tibble with method, metric, positive correlation estimates, and
#' whether the estimate is significant
#'
#' @details Note that since the ranks are ascending i.e. the most relevant
#' LR score for each method has a rank of 1, so I flip sign of the correlation
#' (i.e. highest rank number to highest expression)
get_adt_correlations <- function(liana_adt){
    liana_adt %>%
        # obtain correlations
        group_by(method_name) %>%
        nest() %>%
        mutate(mean_model = map(data, ~corr_fun(.x, var="adt_mean"))) %>%
        mutate(scale_model = map(data, ~corr_fun(.x, var="adt_scale"))) %>%
        ungroup() %>%
        # format correlations
        dplyr::select(-data) %>%
        pivot_longer(-method_name, names_to = "metric") %>%
        unnest(value) %>%
        mutate(significant = if_else(p.value <= 0.05, TRUE, FALSE)) %>%
        dplyr::select(method = method_name, metric, estimate, significant) %>%
        # reverse negative estimates (reverse rank - highest to first)
        mutate(estimate = -1*estimate)
}



#' Helper function to convert from str to symbol
#' @param
#' @inheritDotParams str_to_title
#'
str_to_symbol <- function(string, organism, ...){
    if(organism=="human"){
        str_to_upper(string, ...)
    } else if(organism=="mouse"){
        str_to_title(string, ...)
    }
}

#' Helper function used to prepare `adt_rank` elements
#'
#' @param df `get_rank_adt()` output.
#' @param arbitrary_thresh z-score threshold to calculate ROCs
#'
#' @return tidy data frame with meta information for each dataset and the
#'   response (1;0) and the predictor value which are required for ROC and
#'   curve analysis
#'
#' @details Yardstick ranks thresholds from lowest to highest, i.e. if we use
#' already ranked data then rank 1 would be the smallest value, rather than the
#' highest. Thus, we invert by multiplying by -1, and rank becomes the highest value
#' i.e. -1, while rank 50,000 becomes the lowest value i.e. -50,000
prepare_for_roc <- function(df, arbitrary_thresh){
    df %>%
        filter(!(method_name %in% c("mean_rank", "median_rank"))) %>%
        dplyr::rename(predictor = value) %>%
        # Reverse ranks, so that highest ranked interactions get the highest score - see details
        mutate(predictor = predictor * -1) %>%
        group_by(method_name) %>%
        dplyr::mutate(response = case_when(adt_scale >= arbitrary_thresh ~ 1,
                                           adt_scale < arbitrary_thresh ~ 0)) %>%
        mutate(response = factor(response, levels = c(1, 0))) %>% # 1 is truth
        dplyr::select(source.target.entity, method_name,
                      adt_scale, predictor, response)
}


#' Function to get Ranks for each method and ADT values
#' to reformat the `wrap_adt_corr` function
#'
#' @inheritParams get_adt_means
get_rank_adt <- function(seurat_object,
                         liana_res,
                         op_resource,
                         cluster_key,
                         organism,
                         receptor_syms,
                         adt_symbols,
                         adt_aliases){

    # Get ADT means
    adt_means <- get_adt_means(seurat_object = seurat_object,
                               liana_res = liana_res,
                               op_resource = op_resource,
                               cluster_key = cluster_key,
                               organism = organism,
                               receptor_syms = receptor_syms,
                               adt_symbols = adt_symbols,
                               adt_aliases = adt_aliases)

    # liana_format_adt
    liana_adt <- liana_res %>%
        filter(receptor %in% adt_means$entity_symbol) %>%
        dplyr::select(source, ligand, target, receptor, ends_with("rank")) %>%
        pivot_longer(c(ligand,receptor),
                     names_to = "type",
                     values_to = "entity_symbol") %>%
        filter(type == "receptor") %>%
        distinct() %>%
        left_join(adt_means, by = c("target", "entity_symbol")) %>%
        unite(col=source.target.entity, source, target, entity_symbol, sep = ".")  %>%
        dplyr::select(ends_with("rank"), starts_with("adt"), source.target.entity) %>%
        pivot_longer(-c(source.target.entity, adt_scale, adt_mean),
                     names_to = "method_name")

    return(liana_adt)
}




#' @title Generate AUROCs - required to run `calc_curve` function
#'
#' @inheritParams get_rank_adt
generate_specificity_roc <- function(seurat_object,
                                     liana_res,
                                     op_resource,
                                     cluster_key,
                                     organism,
                                     receptor_syms,
                                     adt_symbols,
                                     arbitrary_thresh,
                                     adt_aliases){

    # get Z-scaled ADT abundance and LR ranks DF
    adt_rank_results <- get_rank_adt(seurat_object = seurat_object,
                                     liana_res = liana_res,
                                     op_resource = op_resource,
                                     cluster_key = cluster_key,
                                     organism = organism,
                                     adt_symbols = adt_symbols,
                                     receptor_syms = receptor_syms,
                                     adt_aliases = adt_aliases)

    # To appropriate format (prepare_for_roc) and do ROC
    adt_rank_roc <- adt_rank_results %>%
        prepare_for_roc(arbitrary_thresh = arbitrary_thresh) %>%
        group_by(method_name) %>%
        group_nest(.key = "adt_rank") %>%
        # Calculate ROC
        mutate(roc = .data$adt_rank %>%
                   map(function(df) calc_curve(df,
                                               downsampling = FALSE,
                                               source = "source.target.entity",
                                               auc_only = TRUE))) %>%
        # Calculate PRROC
        mutate(prc = .data$adt_rank %>%
                   map(function(df) calc_curve(df,
                                               curve = "PR",
                                               downsampling = TRUE,
                                               times = 100,
                                               source = "source.target.entity",
                                               auc_only = TRUE)))
    return(adt_rank_roc)
}

