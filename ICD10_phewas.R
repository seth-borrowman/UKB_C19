library(ggrepel)
library(PheWAS)
library(tidyverse)

### Import data ----
apoe <- read_csv("Z:/Chandra Lab/002 Chandra Lab Projects/UKB Research/UKB APOE/apoe_genotypes.csv") %>%
  filter(epsilon != "e1/e3 or e2/e4") %>%
  select(eid, e1, e2, e3, e4) %>%
  rename(id = eid)
covar <- read_delim("Z:/Chandra Lab/002 Chandra Lab Projects/UKB Research/UKB APOE/covariates_new.txt") %>%
  rename(id = IID) %>%
  select(-c(FID, Age2, Sex.Male)) %>%
  filter(id %in% apoe$id)
apoe <- apoe %>% filter(id %in% covar$id)

pheno <- read_csv("Z:/Chandra Lab/002 Chandra Lab Projects/UKB Research/UKB APOE/pheno_icd10_long.csv")

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
  pheno,
  min.code.count = 1, add.phecode.exclusions = TRUE, translate = TRUE,
  vocabulary.map = PheWAS::phecode_map_icd10,
  full.population.ids = covar$id,
  id.sex = sex
)
write_csv(phenotypes, "phecode.csv")

# Set up folders ----
# Create table to aggregate significant results from all variants
sig_phewas <- data.frame()

# Create directories to write PheWAS results and plots for each variant
ifelse(!dir.exists("png"), dir.create("png"), FALSE)
ifelse(!dir.exists("csv"), dir.create("csv"), FALSE)

# Run PheWAS for each variant ----
for (i in 2:4) {
  # Load geno data
  geno_data <- apoe[, c(1, 1 + i)]
  results <- phewas(
    phenotypes, geno_data,
    cores = detectCores() - 1,
    significance.threshold = c("p-value", "bonferroni", "fdr"),
    covariates = covar, additive.genotypes = TRUE,
    alpha = 0.05
  )

  # Add PheWAS descriptions
  results_d <- addPhecodeInfo(results)

  # Get significant results
  res <- results_d[results_d$bonferroni & !is.na(results_d$p), ]
  print("Results")
  print(res)
  sig_phewas <- rbind(sig_phewas, res)

  # Re-create the same threshold for plots as is used in bonferroni column
  sig_p <- (0.05) / (nrow(results_d[!is.na(results_d$p), ]))
  print("sig_p")
  print(sig_p)
  print(nrow(results_d[!is.na(results_d$p), ]))

  # Exctract significant result and save it as csv
  results_d <- results_d[!is.na(results_d$p), ]
  results_d <- results_d[order(results_d$group), ]
  # Add try() because this step sometimes fails and kills the loop
  try(results_d$order_num <- rep(seq_len(nrow(results_d)), silent = TRUE))
  try(results_d$order_num <- factor(results_d$order_num,
        levels = results_d$order_num
      ), silent = TRUE)
  write.csv(results_d, sprintf("csv/%s_phewas.csv", paste0("e", i)),
    row.names = FALSE
  )

  # Create Manhattan plot annotating significant phenotypes
  # Significant phenotypes are defined by passing Bonferroni correction
  png(
    filename = sprintf("png/%s_phewas.png", paste0("e", i)),
    width = 1400, height = 800
  )
  options(ggrepel.max.overlaps = Inf)
  man_plot <- ggplot(
    results_d,
    aes(x = order_num, y = -log(p))
  ) +
    geom_point(aes(col = group, size = OR)) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      panel.grid.minor = element_line(
        colour = "grey",
        linetype = "dashed"
      ),
      axis.ticks = element_blank()
    ) +
    labs(
      color = "Category", size = "Effect size",
      x = paste0("e", i), y = "-log(p-value)"
    ) +
    geom_text_repel(
      data = . %>%
        mutate(label = ifelse((p < sig_p) & (bonferroni == TRUE),
          as.character(description), ""
        )),
      aes(label = label), size = 3, box.padding = unit(0.7, "lines")
    ) +
    geom_hline(
      yintercept = -log(sig_p), color = "red",
      linewidth = 1, alpha = 0.5
    )
  print(man_plot)
  dev.off()
}

# Output aggregate results ----
write.csv(sig_phewas, "significant_phewas.csv", row.names = FALSE)

sig_phewas_tb <- tibble::tibble(sig_phewas)
sig_phewas_agg <- dplyr::count(sig_phewas_tb, description, sort = TRUE)

write.csv(sig_phewas_agg, "significant_phewas_agg.csv", row.names = FALSE)

### Run with categorical genotype ----
apoe_epsilon <- read_csv("Z:/Chandra Lab/002 Chandra Lab Projects/UKB Research/UKB APOE/apoe_genotypes.csv") %>%
  select(eid, epsilon) %>%
  rename(id = eid) %>%
  mutate(epsilon = factor(epsilon,
    levels = c(
      "e3/e3",
      "e3/e4",
      "e2/e3",
      "e1/e3 or e2/e4",
      "e4/e4",
      "e2/e2"
    ),
    ordered = TRUE
  ))
covar_epsilon <- read_delim("Z:/Chandra Lab/002 Chandra Lab Projects/UKB Research/UKB APOE/covariates_new.txt") %>%
  rename(id = IID) %>%
  select(-c(FID, Age2, Sex.Male)) %>%
  filter(id %in% apoe_epsilon$id)

apoe_epsilon <- cbind(
  apoe_epsilon$id,
  model.matrix(~ epsilon - 1, data = apoe_epsilon)
) %>%
  as.data.frame() %>%
  rename(id = V1) %>%
  rename(`epsilone2/e4` = `epsilone1/e3 or e2/e4`) %>%
  filter(id %in% covar_epsilon$id)

# Make things easier for creating phenotypes
sex_epsilon <- covar_epsilon %>%
  select(id, Sex.Female) %>%
  mutate(sex = case_when(
    Sex.Female == 0 ~ "M",
    Sex.Female == 1 ~ "F"
  )) %>%
  select(-Sex.Female)
# Remove sex from further analysis - causes issues
covar_epsilon <- covar_epsilon %>% select(-Sex.Female)

phenotypes_epsilon <- createPhenotypes(
  pheno,
  min.code.count = 1, add.phecode.exclusions = TRUE, translate = TRUE,
  vocabulary.map = PheWAS::phecode_map_icd10,
  full.population.ids = covar_epsilon$id,
  id.sex = sex_epsilon
)

# Create directories to write PheWAS results and plots for each variant
ifelse(!dir.exists("epsilon"), dir.create("epsilon"), FALSE)

combined_results <- data.frame()
sig_phewas <- data.frame()

for (i in 3:7) {
  # Load geno data
  whichrows <- which(apoe_epsilon$`epsilone3/e3` == 1 |
                       apoe_epsilon[, i] == 1)
  geno_data <- apoe_epsilon[whichrows, c(1, i)]
  genotype <- colnames(geno_data)[2]
  results <- phewas(
    phenotypes_epsilon, geno_data,
    cores = detectCores() - 1,
    significance.threshold = c("p-value", "bonferroni", "fdr"),
    covariates = covar_epsilon, additive.genotypes = FALSE,
    alpha = 0.05, MASS.confint.level = 0.95
  )

  # Add PheWAS descriptions
  results_d <- addPhecodeInfo(results)

  # Get significant results
  res <- results_d[results_d$bonferroni & !is.na(results_d$p), ]
  print("Results")
  print(res)
  sig_phewas <- rbind(sig_phewas, res)

  # Re-create the same threshold for plots as is used in bonferroni column
  sig_p <- (0.05) / (nrow(results_d[!is.na(results_d$p), ]))
  print("sig_p")
  print(sig_p)
  print(nrow(results_d[!is.na(results_d$p), ]))

  # Extract results and save it as csv
  results_d <- results_d[!is.na(results_d$p), ]
  results_d <- results_d[order(results_d$group), ]
  # Add try() because this step sometimes fails and kills the loop
  try(results_d$order_num <- rep(seq_len(nrow(results_d)), silent = TRUE))
  try(results_d$order_num <- factor(results_d$order_num,
        levels = results_d$order_num
      ), silent = TRUE)
  write.csv(results_d, sprintf("epsilon/%s_phewas.csv", i),
    row.names = FALSE
  )
  combined_results <- rbind(combined_results, results_d)

  # Create Manhattan plot annotating significant phenotypes
  # Significant phenotypes are defined by passing Bonferroni correction
  png(
    filename = sprintf("epsilon/%s_phewas.png", i),
    width = 1400, height = 800
  )
  options(ggrepel.max.overlaps = Inf)
  man_plot <- ggplot(
    results_d,
    aes(x = order_num, y = -log(p))
  ) +
    geom_point(
      data = . %>%
        mutate(shape = case_when(
          OR > 1 ~ "OR > 1",
          OR < 1 ~ "OR < 1",
          .default = "OR = 1"
        )) %>%
        mutate(shape = factor(shape,
          levels = c("OR > 1", "OR < 1", "OR = 1"),
          ordered = TRUE
        )),
      aes(col = group, shape = shape, fill = group)
    ) +
    scale_shape_manual(values = c(24, 25, 21)) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      panel.grid.minor = element_line(
        colour = "grey",
        linetype = "dashed"
      ),
      axis.ticks = element_blank()
    ) +
    labs(
      color = "Category", fill = "Category", shape = "Odds Ratio",
      x = substr(genotype, 8, 12), y = "-log(p-value)"
    ) +
    geom_text_repel(
      data = . %>%
        mutate(label = ifelse((p < sig_p) & (bonferroni == TRUE),
          as.character(description), ""
        )),
      aes(label = label), size = 3, box.padding = unit(0.7, "lines")
    ) +
    geom_hline(
      yintercept = -log(sig_p), color = "red",
      linewidth = 1, alpha = 0.5
    )
  print(man_plot)
  dev.off()
}
write_csv(combined_results, "combined_results.csv")

### Plot combined ORs ----

# HBP
hbp_results <- combined_results %>%
  filter(phenotype == 401.1) %>%
  mutate(snp = substr(snp, 9, 13)) %>%
  rbind(., c(
    NA, "e3/e3", NA, NA, NA, NA, NA, 1, NA, NA,
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    NA, NA, NA
  )) %>%
  mutate(snp = factor(snp,
    levels = c(
      "e2/e2",
      "e2/e3",
      "e3/e3",
      "e2/e4",
      "e3/e4",
      "e4/e4"
    ),
    ordered = TRUE
  )) %>%
  mutate(OR = as.numeric(OR)) %>%
  mutate(lower = as.numeric(lower)) %>%
  mutate(upper = as.numeric(upper)) %>%
  mutate(summary = case_when(
    snp == "e3/e3" ~ "Reference",
    .default = paste0(
      round(OR, 2), " (", round(lower, 2), ", ",
      round(upper, 2), ")"
    )
  )) %>%
  arrange(snp) %>%
  mutate(summary = factor(summary, levels = row_number(), label = summary))


hbp_results %>%
  ggplot(aes(x = OR, y = as.numeric(snp))) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.25) +
  geom_vline(xintercept = 1, linetype = 2, alpha = 0.5) +
  labs(x = "OR (95% CI)", title = "Essential Hypertension") +
  scale_y_continuous(
    breaks = seq_along(levels(hbp_results$snp)),
    labels = levels(hbp_results$snp),
    sec.axis = sec_axis(~.,
      breaks = seq_along(levels(hbp_results$summary)),
      labels = levels(hbp_results$summary)
    )
  ) +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(),
    axis.line.y.left = element_line(),
    axis.line.y.right = element_blank(),
    axis.text.x = element_text(colour = "black", size = 11),
    axis.text.y = element_text(colour = "black", size = 11),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(),
    axis.title.y = element_blank(),
    plot.title = element_text(
      color = "black", size = 14,
      hjust = 0.5, face = "bold"
    )
  )
ggsave("Hypertension_OR.png", dpi = 600)

# Appendicitis
appendix_results <- combined_results %>%
  filter(phenotype == 540.1) %>%
  mutate(snp = substr(snp, 9, 13)) %>%
  rbind(., c(
    NA, "e3/e3", NA, NA, NA, NA, NA, 1, NA, NA,
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    NA, NA, NA
  )) %>%
  mutate(snp = factor(snp,
    levels = c(
      "e2/e2",
      "e2/e3",
      "e3/e3",
      "e2/e4",
      "e3/e4",
      "e4/e4"
    ),
    ordered = TRUE
  )) %>%
  mutate(OR = as.numeric(OR)) %>%
  mutate(lower = as.numeric(lower)) %>%
  mutate(upper = as.numeric(upper)) %>%
  mutate(summary = case_when(
    snp == "e3/e3" ~ "Reference",
    .default = paste0(
      round(OR, 2), " (", round(lower, 2), ", ",
      round(upper, 2), ")"
    )
  )) %>%
  arrange(snp) %>%
  mutate(summary = factor(summary, levels = row_number(), label = summary))


appendix_results %>%
  ggplot(aes(x = OR, y = as.numeric(snp))) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.25) +
  geom_vline(xintercept = 1, linetype = 2, alpha = 0.5) +
  labs(x = "OR (95% CI)", title = "Appendicitis") +
  scale_y_continuous(
    breaks = seq_along(levels(appendix_results$snp)),
    labels = levels(appendix_results$snp),
    sec.axis = sec_axis(~.,
      breaks = seq_along(levels(appendix_results$summary)),
      labels = levels(appendix_results$summary)
    )
  ) +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(),
    axis.line.y.left = element_line(),
    axis.line.y.right = element_blank(),
    axis.text.x = element_text(colour = "black", size = 11),
    axis.text.y = element_text(colour = "black", size = 11),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(),
    axis.title.y = element_blank(),
    plot.title = element_text(
      color = "black", size = 14,
      hjust = 0.5, face = "bold"
    )
  )
ggsave("Appendicitis_OR.png", dpi = 600)

# Hypovolemia
hypovolemia_results <- combined_results %>%
  filter(phenotype == 276.5) %>%
  mutate(snp = substr(snp, 9, 13)) %>%
  rbind(., c(
    NA, "e3/e3", NA, NA, NA, NA, NA, 1, NA, NA,
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    NA, NA, NA
  )) %>%
  mutate(snp = factor(snp,
    levels = c(
      "e2/e2",
      "e2/e3",
      "e3/e3",
      "e2/e4",
      "e3/e4",
      "e4/e4"
    ),
    ordered = TRUE
  )) %>%
  mutate(OR = as.numeric(OR)) %>%
  mutate(lower = as.numeric(lower)) %>%
  mutate(upper = as.numeric(upper)) %>%
  mutate(summary = case_when(
    snp == "e3/e3" ~ "Reference",
    .default = paste0(
      round(OR, 2), " (", round(lower, 2), ", ",
      round(upper, 2), ")"
    )
  )) %>%
  arrange(snp) %>%
  mutate(summary = factor(summary, levels = row_number(), label = summary))


hypovolemia_results %>%
  ggplot(aes(x = OR, y = as.numeric(snp))) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.25) +
  geom_vline(xintercept = 1, linetype = 2, alpha = 0.5) +
  labs(x = "OR (95% CI)", title = "Hypovolemia") +
  scale_y_continuous(
    breaks = seq_along(levels(hypovolemia_results$snp)),
    labels = levels(hypovolemia_results$snp),
    sec.axis = sec_axis(~.,
      breaks = seq_along(levels(hypovolemia_results$summary)),
      labels = levels(hypovolemia_results$summary)
    )
  ) +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(),
    axis.line.y.left = element_line(),
    axis.line.y.right = element_blank(),
    axis.text.x = element_text(colour = "black", size = 11),
    axis.text.y = element_text(colour = "black", size = 11),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(),
    axis.title.y = element_blank(),
    plot.title = element_text(
      color = "black", size = 14,
      hjust = 0.5, face = "bold"
    )
  )
ggsave("Hypovolemia_OR.png", dpi = 600)

# Retention of Urine
rou_results <- combined_results %>%
  filter(phenotype == 599.2) %>%
  mutate(snp = substr(snp, 9, 13)) %>%
  rbind(., c(
    NA, "e3/e3", NA, NA, NA, NA, NA, 1, NA, NA,
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    NA, NA, NA
  )) %>%
  mutate(snp = factor(snp,
    levels = c(
      "e2/e2",
      "e2/e3",
      "e3/e3",
      "e2/e4",
      "e3/e4",
      "e4/e4"
    ),
    ordered = TRUE
  )) %>%
  mutate(OR = as.numeric(OR)) %>%
  mutate(lower = as.numeric(lower)) %>%
  mutate(upper = as.numeric(upper)) %>%
  mutate(summary = case_when(
    snp == "e3/e3" ~ "Reference",
    .default = paste0(
      round(OR, 2), " (", round(lower, 2), ", ",
      round(upper, 2), ")"
    )
  )) %>%
  arrange(snp) %>%
  mutate(summary = factor(summary, levels = row_number(), label = summary))


rou_results %>%
  ggplot(aes(x = OR, y = as.numeric(snp))) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.25) +
  geom_vline(xintercept = 1, linetype = 2, alpha = 0.5) +
  labs(x = "OR (95% CI)", title = "Disordered Retention of Urine") +
  scale_y_continuous(
    breaks = seq_along(levels(rou_results$snp)),
    labels = levels(rou_results$snp),
    sec.axis = sec_axis(~.,
      breaks = seq_along(levels(rou_results$summary)),
      labels = levels(rou_results$summary)
    )
  ) +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(),
    axis.line.y.left = element_line(),
    axis.line.y.right = element_blank(),
    axis.text.x = element_text(colour = "black", size = 11),
    axis.text.y = element_text(colour = "black", size = 11),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(),
    axis.title.y = element_blank(),
    plot.title = element_text(
      color = "black", size = 14,
      hjust = 0.5, face = "bold"
    )
  )
ggsave("RoU_OR.png", dpi = 600)

# Incontinence
incont_results <- combined_results %>%
  filter(phenotype == 599.4) %>%
  mutate(snp = substr(snp, 9, 13)) %>%
  rbind(., c(
    NA, "e3/e3", NA, NA, NA, NA, NA, 1, NA, NA,
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    NA, NA, NA
  )) %>%
  mutate(snp = factor(snp,
    levels = c(
      "e2/e2",
      "e2/e3",
      "e3/e3",
      "e2/e4",
      "e3/e4",
      "e4/e4"
    ),
    ordered = TRUE
  )) %>%
  mutate(OR = as.numeric(OR)) %>%
  mutate(lower = as.numeric(lower)) %>%
  mutate(upper = as.numeric(upper)) %>%
  mutate(summary = case_when(
    snp == "e3/e3" ~ "Reference",
    .default = paste0(
      round(OR, 2), " (", round(lower, 2), ", ",
      round(upper, 2), ")"
    )
  )) %>%
  arrange(snp) %>%
  mutate(summary = factor(summary, levels = row_number(), label = summary))


incont_results %>%
  ggplot(aes(x = OR, y = as.numeric(snp))) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.25) +
  geom_vline(xintercept = 1, linetype = 2, alpha = 0.5) +
  labs(x = "OR (95% CI)", title = "Urinary Incontinence") +
  scale_y_continuous(
    breaks = seq_along(levels(incont_results$snp)),
    labels = levels(incont_results$snp),
    sec.axis = sec_axis(~.,
      breaks = seq_along(levels(incont_results$summary)),
      labels = levels(incont_results$summary)
    )
  ) +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(),
    axis.line.y.left = element_line(),
    axis.line.y.right = element_blank(),
    axis.text.x = element_text(colour = "black", size = 11),
    axis.text.y = element_text(colour = "black", size = 11),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(),
    axis.title.y = element_blank(),
    plot.title = element_text(
      color = "black", size = 14,
      hjust = 0.5, face = "bold"
    )
  )
ggsave("Incontinence_OR.png", dpi = 600)

# Anxiety
anxiety_results <- combined_results %>%
  filter(phenotype == 300) %>%
  mutate(snp = substr(snp, 9, 13)) %>%
  rbind(., c(
    NA, "e3/e3", NA, NA, NA, NA, NA, 1, NA, NA,
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    NA, NA, NA
  )) %>%
  mutate(snp = factor(snp,
    levels = c(
      "e2/e2",
      "e2/e3",
      "e3/e3",
      "e2/e4",
      "e3/e4",
      "e4/e4"
    ),
    ordered = TRUE
  )) %>%
  mutate(OR = as.numeric(OR)) %>%
  mutate(lower = as.numeric(lower)) %>%
  mutate(upper = as.numeric(upper)) %>%
  mutate(summary = case_when(
    snp == "e3/e3" ~ "Reference",
    .default = paste0(
      round(OR, 2), " (", round(lower, 2), ", ",
      round(upper, 2), ")"
    )
  )) %>%
  arrange(snp) %>%
  mutate(summary = factor(summary, levels = row_number(), label = summary))


anxiety_results %>%
  ggplot(aes(x = OR, y = as.numeric(snp))) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.25) +
  geom_vline(xintercept = 1, linetype = 2, alpha = 0.5) +
  labs(x = "OR (95% CI)", title = "Anxiety Disorders") +
  scale_y_continuous(
    breaks = seq_along(levels(anxiety_results$snp)),
    labels = levels(anxiety_results$snp),
    sec.axis = sec_axis(~.,
      breaks = seq_along(levels(anxiety_results$summary)),
      labels = levels(anxiety_results$summary)
    )
  ) +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(),
    axis.line.y.left = element_line(),
    axis.line.y.right = element_blank(),
    axis.text.x = element_text(colour = "black", size = 11),
    axis.text.y = element_text(colour = "black", size = 11),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(),
    axis.title.y = element_blank(),
    plot.title = element_text(
      color = "black", size = 14,
      hjust = 0.5, face = "bold"
    )
  )
ggsave("Anxiety_OR.png", dpi = 600)

# Dementias
dementia_results <- combined_results %>%
  filter(phenotype == 290.1) %>%
  mutate(snp = substr(snp, 9, 13)) %>%
  rbind(., c(
    NA, "e3/e3", NA, NA, NA, NA, NA, 1, NA, NA,
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    NA, NA, NA
  )) %>%
  mutate(snp = factor(snp,
    levels = c(
      "e2/e2",
      "e2/e3",
      "e3/e3",
      "e2/e4",
      "e3/e4",
      "e4/e4"
    ),
    ordered = TRUE
  )) %>%
  mutate(OR = as.numeric(OR)) %>%
  mutate(lower = as.numeric(lower)) %>%
  mutate(upper = as.numeric(upper)) %>%
  mutate(summary = case_when(
    snp == "e3/e3" ~ "Reference",
    .default = paste0(
      round(OR, 2), " (", round(lower, 2), ", ",
      round(upper, 2), ")"
    )
  )) %>%
  arrange(snp) %>%
  mutate(summary = factor(summary, levels = row_number(), label = summary))


dementia_results %>%
  ggplot(aes(x = OR, y = as.numeric(snp))) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.25) +
  geom_vline(xintercept = 1, linetype = 2, alpha = 0.5) +
  labs(x = "OR (95% CI)", title = "Dementias") +
  scale_y_continuous(
    breaks = seq_along(levels(dementia_results$snp)),
    labels = levels(dementia_results$snp),
    sec.axis = sec_axis(~.,
      breaks = seq_along(levels(dementia_results$summary)),
      labels = levels(dementia_results$summary)
    )
  ) +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(),
    axis.line.y.left = element_line(),
    axis.line.y.right = element_blank(),
    axis.text.x = element_text(colour = "black", size = 11),
    axis.text.y = element_text(colour = "black", size = 11),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(),
    axis.title.y = element_blank(),
    plot.title = element_text(
      color = "black", size = 14,
      hjust = 0.5, face = "bold"
    )
  )
ggsave("Dementia_OR.png", dpi = 600)

# PD
park_results <- combined_results %>%
  filter(phenotype == 332) %>%
  mutate(snp = substr(snp, 9, 13)) %>%
  rbind(., c(
    NA, "e3/e3", NA, NA, NA, NA, NA, 1, NA, NA,
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    NA, NA, NA
  )) %>%
  mutate(snp = factor(snp,
    levels = c(
      "e2/e2",
      "e2/e3",
      "e3/e3",
      "e2/e4",
      "e3/e4",
      "e4/e4"
    ),
    ordered = TRUE
  )) %>%
  mutate(OR = as.numeric(OR)) %>%
  mutate(lower = as.numeric(lower)) %>%
  mutate(upper = as.numeric(upper)) %>%
  mutate(summary = case_when(
    snp == "e3/e3" ~ "Reference",
    .default = paste0(
      round(OR, 2), " (", round(lower, 2), ", ",
      round(upper, 2), ")"
    )
  )) %>%
  arrange(snp) %>%
  mutate(summary = factor(summary, levels = row_number(), label = summary))


park_results %>%
  ggplot(aes(x = OR, y = as.numeric(snp))) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.25) +
  geom_vline(xintercept = 1, linetype = 2, alpha = 0.5) +
  labs(x = "OR (95% CI)", title = "Parkinson's Disease") +
  scale_y_continuous(
    breaks = seq_along(levels(park_results$snp)),
    labels = levels(park_results$snp),
    sec.axis = sec_axis(~.,
      breaks = seq_along(levels(park_results$summary)),
      labels = levels(park_results$summary)
    )
  ) +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(),
    axis.line.y.left = element_line(),
    axis.line.y.right = element_blank(),
    axis.text.x = element_text(colour = "black", size = 11),
    axis.text.y = element_text(colour = "black", size = 11),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(),
    axis.title.y = element_blank(),
    plot.title = element_text(
      color = "black", size = 14,
      hjust = 0.5, face = "bold"
    )
  )
ggsave("Parkinson_OR.png", dpi = 600)
