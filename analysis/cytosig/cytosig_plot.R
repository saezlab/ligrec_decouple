
# Read results
cytosig_eval <- readRDS("data/output/cytosig_out/brca_{}_cytosig_res.RDS") %>%
    select(dataset, cytosig_res) %>%
    unnest(cytosig_res)

aucs <- cytosig_eval %>%
    select(-c(cyto_liana, corr)) %>%
    unnest(prc) %>%
    select(dataset, method_name, roc, prc_auc = auc) %>%
    distinct() %>%
    unnest(roc) %>%
    select(dataset, method_name, roc_auc = auc, prc_auc) %>%
    distinct() %>%
    group_by(method_name) %>%
    mutate(roc_mean = mean(roc_auc),
           prc_mean = mean(prc_auc)) %>%
    ungroup() %>%
    mutate(method_name = gsub("\\..*","", method_name)) %>%
    mutate(method_name = recode_methods(method_name))


roc_min <- ifelse(min(aucs$roc_auc) > 0.5, 0.5, min(aucs$roc_auc))
prc_min <- ifelse(min(aucs$prc_auc) > 0.5, 0.5, min(aucs$prc_auc))

# min_lim <- floor(min(c(aucs$roc, aucs$prc)) * 100)/100
# max_lim <- ceiling(max(c(aucs$roc, aucs$prc)) * 100)/100

ggplot(aucs,
       aes(x=roc_mean,
           y=prc_mean,
           color=method_name)) +
    geom_point(shape = 9, size = 12, alpha=1) +
    geom_point(aes(x = roc_auc,
                   prc_auc,
                   shape=dataset),
               size = 6,
               alpha = 0.3) +
    theme(text = element_text(size=16)) +
    xlab('AUROC') +
    ylab('AUPRC') +
    xlim(roc_min, 0.7) +
    ylim(prc_min, 0.8) +
    geom_hline(yintercept = 0.5, colour = "pink",
               linetype = 2, size = 1.2) +
    geom_vline(xintercept = 0.5, colour = "pink",
               linetype = 2, size = 1.2) +
    theme_bw(base_size = 30) +
    guides(shape=guide_legend(title="Dataset"),
           color=guide_legend(title="Method"))
