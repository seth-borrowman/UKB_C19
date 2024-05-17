library(tidyverse)

system("dx download APOE/apoe_genotypes.csv")
system("dx download APOE/covid_updated.csv")
system("dx download APOE/covariates_new.txt")

apoe <- read_csv("apoe_genotypes.csv")
covid <- read_csv("covid_updated.csv")
# Some problems with unneccesary columns
covid <- covid[, -c(235, 236)]
covar <- read_delim("covariates_new.txt")

covar <- covar %>%
  rename(eid = IID) %>%
  select(-FID)


# Let's go exploring ----
# Associations with plasma abeta levels
abeta <- covid %>%
  select(
    eid,
    plasma_abeta40_i2, plasma_abeta40_i3,
    plasma_abeta42_i2, plasma_abeta42_i3,
    neurobio_plate_i2, neurobio_plate_i3,
    nuerobio_batch_i2, nuerobio_batch_i3
  ) %>%
  merge(., apoe, by = "eid") %>%
  merge(., covar, by = "eid")

hist(abeta$plasma_abeta40_i2)
ab40i2_data <- abeta %>%
  select(-c(
    eid, plasma_abeta40_i3,
    plasma_abeta42_i2, plasma_abeta42_i3,
    e1, e2, e3, e4, rs429358, rs7412, Age2,
    neurobio_plate_i3, nuerobio_batch_i3
  ))
ab40i2 <- lm(plasma_abeta40_i2 ~ ., data = ab40i2_data)
summary(ab40i2)

hist(abeta$plasma_abeta40_i3)
ab40i3_data <- abeta %>%
  select(-c(
    eid, plasma_abeta40_i2,
    plasma_abeta42_i2, plasma_abeta42_i3,
    e1, e2, e3, e4, rs429358, rs7412, Age2,
    neurobio_plate_i2, nuerobio_batch_i2
  ))
ab40i3 <- lm(plasma_abeta40_i3 ~ ., data = ab40i3_data)
summary(ab40i3)

hist(abeta$plasma_abeta42_i2)
ab42i2_data <- abeta %>%
  select(-c(
    eid, plasma_abeta40_i2,
    plasma_abeta40_i3, plasma_abeta42_i3,
    e1, e2, e3, e4, rs429358, rs7412, Age2,
    neurobio_plate_i3, nuerobio_batch_i3
  ))
ab42i2 <- lm(plasma_abeta42_i2 ~ ., data = ab42i2_data)
summary(ab42i2)

hist(abeta$plasma_abeta42_i3)
ab42i3_data <- abeta %>%
  select(-c(
    eid, plasma_abeta40_i2,
    plasma_abeta40_i3, plasma_abeta42_i2,
    e1, e2, e3, e4, rs429358, rs7412, Age2,
    neurobio_plate_i2, nuerobio_batch_i2
  ))
ab42i3 <- lm(plasma_abeta42_i3 ~ ., data = ab42i3_data)
summary(ab42i3)

# Any associations with change in plasma abeta levels?
delta_ab40_data <- abeta %>%
  select(-c(
    eid, plasma_abeta42_i2, plasma_abeta42_i3,
    e1, e2, e3, e4, rs429358, rs7412, Age2
  )) %>%
  mutate(delta = plasma_abeta40_i3 - plasma_abeta40_i2) %>%
  select(-c(plasma_abeta40_i2, plasma_abeta40_i3))
hist(delta_ab40_data$delta)
ab40_delta <- lm(delta ~ ., data = delta_ab40_data)
summary(ab40_delta)

delta_ab42_data <- abeta %>%
  select(-c(
    eid, plasma_abeta40_i2, plasma_abeta40_i3,
    e1, e2, e3, e4, rs429358, rs7412, Age2
  )) %>%
  mutate(delta = plasma_abeta42_i3 - plasma_abeta42_i2) %>%
  select(-c(plasma_abeta42_i2, plasma_abeta42_i3))
hist(delta_ab42_data$delta)
ab42_delta <- lm(delta ~ ., data = delta_ab42_data)
summary(ab42_delta)

# NFL and Glialfap
nflgli <- covid %>%
  select(
    eid, plasma_nfl_i2, plasma_nfl_i3,
    plasma_glialfap_i2, plasma_glialfap_i3,
    neurobio_plate_i2, neurobio_plate_i3,
    nuerobio_batch_i2, nuerobio_batch_i3
  ) %>%
  merge(., apoe, by = "eid") %>%
  merge(., covar, by = "eid")

# Plasma NeuroFilament Light
hist(nflgli$plasma_nfl_i2) # possible outliers?
nfli2_data <- nflgli %>%
  select(-c(
    eid, plasma_nfl_i3,
    plasma_glialfap_i2, plasma_glialfap_i3,
    e1, e2, e3, e4, rs429358, rs7412, Age2,
    neurobio_plate_i3, nuerobio_batch_i3
  ))
nfli2 <- lm(plasma_nfl_i2 ~ ., data = nfli2_data)
summary(nfli2)

hist(nflgli$plasma_nfl_i3) # possible outliers?
nfli3_data <- nflgli %>%
  select(-c(
    eid, plasma_nfl_i2,
    plasma_glialfap_i2, plasma_glialfap_i3,
    e1, e2, e3, e4, rs429358, rs7412, Age2,
    neurobio_plate_i2, nuerobio_batch_i2
  ))
nfli3 <- lm(plasma_nfl_i3 ~ ., data = nfli3_data)
summary(nfli3)

# Plasma Glial fibrillary acidic protein
hist(nflgli$plasma_glialfap_i2) # possible outliers?
glii2_data <- nflgli %>%
  select(-c(
    eid, plasma_nfl_i3,
    plasma_nfl_i2, plasma_glialfap_i3,
    e1, e2, e3, e4, rs429358, rs7412, Age2,
    neurobio_plate_i3, nuerobio_batch_i3
  ))
glii2 <- lm(plasma_glialfap_i2 ~ ., data = glii2_data)
summary(glii2)

hist(nflgli$plasma_glialfap_i3) # possible outliers?
glii3_data <- nflgli %>%
  select(-c(
    eid, plasma_nfl_i3,
    plasma_nfl_i2, plasma_glialfap_i2,
    e1, e2, e3, e4, rs429358, rs7412, Age2,
    neurobio_plate_i2, nuerobio_batch_i2
  ))
glii3 <- lm(plasma_glialfap_i3 ~ ., data = glii3_data)
summary(glii3)

# Any associations with change in plasma nfl/gli levels?
delta_nfl_data <- nflgli %>%
  select(-c(
    eid, plasma_glialfap_i2, plasma_glialfap_i3,
    e1, e2, e3, e4, rs429358, rs7412, Age2
  )) %>%
  mutate(delta = plasma_nfl_i3 - plasma_nfl_i2) %>%
  select(-c(plasma_nfl_i2, plasma_nfl_i3))
hist(delta_nfl_data$delta)
nfl_delta <- lm(delta ~ ., data = delta_nfl_data)
summary(nfl_delta)

delta_gli_data <- nflgli %>%
  select(-c(
    eid, plasma_nfl_i2, plasma_nfl_i3,
    e1, e2, e3, e4, rs429358, rs7412, Age2
  )) %>%
  mutate(delta = plasma_glialfap_i3 - plasma_glialfap_i2) %>%
  select(-c(plasma_glialfap_i2, plasma_glialfap_i3))
hist(delta_gli_data$delta)
gli_delta <- lm(delta ~ ., data = delta_gli_data)
summary(gli_delta)

# phosphotau
ptau <- covid %>%
  select(
    eid,
    plasma_ptau_i2, plasma_ptau_i3,
    neurobio_plate_i2, neurobio_plate_i3,
    nuerobio_batch_i2, nuerobio_batch_i3
  ) %>%
  merge(., apoe, by = "eid") %>%
  merge(., covar, by = "eid")

hist(ptau$plasma_ptau_i2) # possible outliers?
ptaui2_data <- ptau %>%
  select(-c(
    eid, plasma_ptau_i3,
    e1, e2, e3, e4, rs429358, rs7412, Age2,
    neurobio_plate_i3, nuerobio_batch_i3
  ))
ptaui2 <- lm(plasma_ptau_i2 ~ ., data = ptaui2_data)
summary(ptaui2)

hist(ptau$plasma_ptau_i3) # possible outliers?
ptaui3_data <- ptau %>%
  select(-c(
    eid, plasma_ptau_i2,
    e1, e2, e3, e4, rs429358, rs7412, Age2,
    neurobio_plate_i2, nuerobio_batch_i2
  ))
ptaui3 <- lm(plasma_ptau_i3 ~ ., data = ptaui3_data)
summary(ptaui3)

delta_ptau_data <- ptau %>%
  select(-c(eid, e1, e2, e3, e4, rs429358, rs7412, Age2)) %>%
  mutate(delta = plasma_ptau_i3 - plasma_ptau_i2) %>%
  select(-c(plasma_ptau_i2, plasma_ptau_i3))
hist(delta_ptau_data$delta)
ptau_delta <- lm(delta ~ ., data = delta_ptau_data)
summary(ptau_delta)


# COVID recovery ----
covid_recovery <- covid[, c(1, 24:30)] %>%
  mutate(first_c19_date = as.Date(first_c19_date, "%Y-%m-%d")) %>%
  mutate(recent_c19_date = as.Date(recent_c19_date, "%Y-%m-%d")) %>%
  mutate(most_recent_date = case_when(
    is.na(first_c19_date) & !is.na(recent_c19_date) ~ recent_c19_date,
    !is.na(first_c19_date) & is.na(recent_c19_date) ~ first_c19_date,
    first_c19_date > recent_c19_date ~ first_c19_date,
    recent_c19_date > first_c19_date ~ recent_c19_date,
    first_c19_date == recent_c19_date ~ first_c19_date
  )) %>%
  mutate(most_recent_method = case_when(
    is.na(first_c19_date) & !is.na(recent_c19_date) ~ recent_c19_diag,
    !is.na(first_c19_date) & is.na(recent_c19_date) ~ first_c19_diag,
    first_c19_date > recent_c19_date ~ first_c19_diag,
    recent_c19_date > first_c19_date ~ recent_c19_diag,
    first_c19_date == recent_c19_date ~ first_c19_diag
  )) %>%
  filter(
    c19_count != "I do not know if I have had COVID-19",
    c19_count != "I do not know how many times I have had COVID-19",
    c19_count != "Prefer not to answer"
  ) %>%
  mutate(c19_count = as.numeric(c19_count)) %>%
  mutate(time_passed = difftime(c19_qs_compl,
    most_recent_date,
    units = "days"
  )) %>%
  mutate(recovered = case_when(
    c19_recovery == "No, not at all" ~ "No",
    c19_recovery == "No, getting worse" ~ "No",
    c19_recovery == "Partially" ~ "No",
    c19_recovery == "Yes, mostly" ~ "Yes",
    c19_recovery == "Yes, completely" ~ "Yes"
  )) %>%
  mutate(long_covid = case_when(
    time_passed > 130 &
      recovered == "No" ~ 1,
    time_passed > 1 &
      recovered == "Yes" ~ 0
  ))

longcov <- merge(covid_recovery, apoe) %>%
  merge(., covar) %>%
  select(-c(
    eid, first_c19_date, first_c19_diag, recent_c19_date,
    recent_c19_diag, time_passed, recovered, c19_recovery,
    rs429358, rs7412, e1, e2, e3, e4, Age2, Sex.Male
  )) %>%
  filter(c19_count > 0)
longcov_mod <- lm(long_covid ~ ., data = longcov)
summary(longcov_mod)
table(longcov$epsilon, longcov$long_covid)
