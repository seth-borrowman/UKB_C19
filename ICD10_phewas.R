library(ggrepel)
library(PheWAS)
library(tidyverse)

### Import data ----
system("dx download APOE/apoe_genotypes.csv")
system("dx download APOE/covariates_new.txt")
system("dx download APOE/pheno_icd10_long.csv")

apoe <- read_csv("C:/Users/sborrowman/Downloads/phewas/apoe_genotypes.csv") %>%
    filter(epsilon != "e1/e3 or e2/e4") %>%
    select(eid, e1, e2, e3, e4) %>%
    rename(id = eid)
covar <- read_delim("C:/Users/sborrowman/Downloads/phewas/covariates_new.txt") %>%
    rename(id = IID) %>%
    select(-c(FID, Age2, Sex.Male)) %>%
    filter(id %in% apoe$id)
apoe <- apoe %>% filter(id %in% covar$id)

pheno <- read_csv('C:/Users/sborrowman/Downloads/phewas/pheno_icd10_long.csv')

# Make things easier for creating phenotypes
sex <- covar %>%
    select(id, Sex.Female) %>%
    mutate(sex = case_when(
        Sex.Female == 0 ~ "M",
        Sex.Female == 1 ~ "F"
    )) %>%
    select(-Sex.Female)
# Remove sex from further analysis - causes issues
covar <- covar %>% select(-Sex.Female)
### Create phenotypes ----
phenotypes <- createPhenotypes(
    pheno, min.code.count = 1, add.phecode.exclusions = T, translate = T,
    vocabulary.map = PheWAS::phecode_map_icd10,
    full.population.ids = covar$id,
    id.sex = sex)
write_csv(phenotypes, "phecode.csv")

# Set up folders ----
# Create table to aggregate significant results from all variants
sig_phewas = data.frame()

# Create directories to write PheWAS results and plots for each variant
ifelse(!dir.exists('png'), dir.create('png'), FALSE)
ifelse(!dir.exists('csv'), dir.create('csv'), FALSE)

# Run PheWAS for each variant ----
for (i in 2:4) {
    # Load geno data
    geno_data <- apoe[,c(1, 1 + i)]
    results <- phewas(
        phenotypes, geno_data, cores = detectCores() - 1,
        significance.threshold = c('p-value', 'bonferroni', 'fdr'),
        covariates = covar, additive.genotypes = T,
        alpha = 0.05)
    
    # Add PheWAS descriptions
    results_d <- addPhecodeInfo(results)
    
    # Get significant results
    res <- results_d[results_d$bonferroni & !is.na(results_d$p),]
    print("Results")
    print(res)
    sig_phewas <- rbind(sig_phewas, res)
    
    # Re-create the same threshold for plots as is used in bonferroni column
    sig_p <- (0.05)/(nrow(results_d[!is.na(results_d$p),]))
    print("sig_p")
    print(sig_p)
    print(nrow(results_d[!is.na(results_d$p),]))
    
    # Exctract significant result and save it as csv
    results_d <- results_d[!is.na(results_d$p),]
    results_d <- results_d[order(results_d$group),]
    # Add try() because this step sometimes fails and kills the loop
    try(results_d$order_num <- rep(1:nrow(results_d)), silent = T)
    try(results_d$order_num <- factor(results_d$order_num,
                                      levels = results_d$order_num), silent = T)
    write.csv(results_d, sprintf('csv/%s_phewas.csv', paste0("e", i)),
              row.names = FALSE)
    
    # Create Manhattan plot annotating significant phenotypes
    # Significant phenotypes are defined by passing Bonferroni correction
    png(filename = sprintf('png/%s_phewas.png', paste0("e", i)),
        width = 1400, height = 800)
    options(ggrepel.max.overlaps = Inf) 
    man_plot <- ggplot(
        results_d, 
        aes(x = order_num, y = -log(p))) +
        geom_point(aes(col = group, size = OR)) +
        theme_classic() +
        theme(
            axis.text.x = element_blank(),
            panel.grid.minor = element_line(colour = 'grey',
                                            linetype = 'dashed'), 
            axis.ticks = element_blank()) +
        labs(color = 'Category', size = 'Effect size',
             x = paste0("e", i), y = '-log(p-value)') +
        geom_text_repel(
            data = . %>%
                mutate(label = ifelse((p < sig_p) & (bonferroni == TRUE),
                                      as.character(description), '')),
            aes(label = label), size = 3, box.padding = unit(0.7, 'lines')) +
        geom_hline(yintercept = -log(sig_p), color = 'red',
                   linewidth = 1, alpha = 0.5) 
    print(man_plot)
    dev.off()
}

# Output aggregate results ----
write.csv(sig_phewas, 'significant_phewas.csv', row.names = FALSE)

sig_phewas_tb <- tibble::tibble(sig_phewas)
sig_phewas_agg <- dplyr::count(sig_phewas_tb, description, sort = TRUE)

write.csv(sig_phewas_agg, 'significant_phewas_agg.csv', row.names = FALSE)
