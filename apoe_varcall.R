library(tidyverse)

system("dx download APOE/APO_var.raw")
raw <- read_delim("APO_var.raw") %>%
    select(IID, `19:44908684:T:C_T`, `19:44908822:C:T_C`) %>%
    na.omit() %>%
    rename(eid = IID) %>%
    rename(rs429358 = `19:44908684:T:C_T`) %>%
    rename(rs7412 = `19:44908822:C:T_C`) %>%
    # Create named variant groups
    mutate(epsilon = case_when(
        rs429358 == 0 & rs7412 == 0 ~ "e1/e1",
        rs429358 == 1 & rs7412 == 0 ~ "e1/e2",
        rs429358 == 1 & rs7412 == 1 ~ "e1/e3 or e2/e4",
        rs429358 == 0 & rs7412 == 1 ~ "e1/e4",
        rs429358 == 2 & rs7412 == 0 ~ "e2/e2",
        rs429358 == 2 & rs7412 == 1 ~ "e2/e3",
        rs429358 == 2 & rs7412 == 2 ~ "e3/e3",
        rs429358 == 1 & rs7412 == 2 ~ "e3/e4",
        rs429358 == 0 & rs7412 == 2 ~ "e4/e4"
    )) %>%
    # Create additive allele variables
    mutate(e1 = case_when(
        epsilon == "e1/e1" ~ 2,
        epsilon == "e1/e2" |
            epsilon == "e1/e3 or e2/e4" |
            epsilon == "e1/e4" ~ 1,
        .default = 0
    )) %>%
    mutate(e2 = case_when(
        epsilon == "e2/e2" ~ 2,
        epsilon == "e1/e2" |
            epsilon == "e1/e3 or e2/e4" |
            epsilon == "e2/e3" ~ 1,
        .default = 0
    )) %>%
    mutate(e3 = case_when(
        epsilon == "e3/e3" ~ 2,
        epsilon == "e1/e3 or e2/e4" |
            epsilon == "e2/e3" |
            epsilon == "e3/e4" ~ 1,
        .default = 0
    )) %>%
    mutate(e4 = case_when(
        epsilon == "e4/e4" ~ 2,
        epsilon == "e1/e4" |
            epsilon == "e1/e3 or e2/e4" |
            epsilon == "e3/e4" ~ 1,
        .default = 0
    ))

# Get frequency of genotypes
raw %>%
    group_by(epsilon) %>%
    summarise(count = n()) %>%
    mutate(freq = round(count / sum(count), 3)) %>%
    arrange(desc(freq)) %>%
    write_csv(., "apoe_freq.csv")

write_csv(raw, "apoe_genotypes.csv")

system("dx upload apoe*.csv --path APOE/")