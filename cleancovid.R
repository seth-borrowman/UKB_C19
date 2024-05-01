library(tidyverse)

system("dx download APOE/covid.csv")

original <- read_csv("covid.csv")
# There seem to be issues with one row
original <- original[-2770,]

# Create dictionary for field id to name conversion
dictionary <- c(
    # Category 163
    "p31040" = "plasma_abeta40",
    "p31041" = "plasma_abeta42",
    "p31042" = "plasma_glialfap",
    "p31043" = "plasma_nfl",
    "p31044" = "plasma_ptau",
    "p31045" = "nuerobio_batch",
    "p31046" = "neurobio_plate",
    "p31047" = "neurobio_well",
    "p31048" = "ptau_comment",
    "p31049" = "plex_comment",
    # Category 1511
    "p29157" = "first_c19_date",
    "p29159" = "recent_c19_date",
    "p29158" = "first_c19_diag",
    "p29160" = "recent_c19_diag",
    "p29156" = "c19_count",
    "p29161" = "c19_recovery",
    "p29207" = "c19_qs_compl",
    "p29194" = "c19_qs_start",
    # Category 992
    "p28010" = "symptom_qs_method",
    "p28011" = "fever",
    "p28012" = "wheezing",
    "p28013" = "chills",
    "p28014" = "chest_pain",
    "p28015" = "fatigue",
    "p28016" = "headache",
    "p28017" = "aches",
    "p28018" = "nausea",
    "p28019" = "sorethroat",
    "p28020" = "abdpain",
    "p28021" = "cough_dry",
    "p28022" = "diarrhea",
    "p28023" = "runnynose",
    "p28024" = "smell_taste",
    "p28025" = "sob",
    "p28026" = "cough_wet",
    "p28027" = "c19_medattention",
    "p28028" = "c19_isolation",
    "p28029" = "c19_hospital",
    "p28030" = "symptom_date",
    "p28031" = "symptom_test_date",
    "p28032" = "post_rec_date",
    "p28033" = "web_rec_date",
    # Category 997
    "p27990" = "c19ab_result",
    "p27993" = "c19ab_despatched",
    "p27991" = "c19ab_analyzed",
    "p27992" = "c19ab_recd"
)

# Create new column name vector
new_colnames <- colnames(original)
for (i in 2:ncol(original)) {
    split_name <- str_split(colnames(original)[i], pattern = "_")
    
    new_colnames[i] <- dictionary[split_name[[1]][1]] %>%
        unname() %>%
        ifelse(is.na(split_name[[1]][2]),
               .,
               paste0(., "_", split_name[[1]][2]))
}

updated <- original
colnames(updated) <- new_colnames

write_csv(updated, "covid_updated.csv")
system("dx upload covid_updated.csv --path APOE/")
