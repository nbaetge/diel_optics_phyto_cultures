Final_Data
================
Nicholas Baetge
Last compiled on 10 November, 2023

- 

<!-- -->

    ```{r,load packages, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
library(tidyverse)
library(lubridate)
```

    ```

- 

############# 

# IMPORT & COMBINE DATA

############## 

- 

<!-- -->

    ```{r,bottle file, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
times <- read_csv(
    "~/Box Sync/Phyto_bbp/DATA/FINAL/diel_optics_phyto_cultures/analyses/fcm/Cultures_and_Optics/processed_fcm_bottle.csv"
  ) %>% select(exp, tp, datetime:plot_datetime) %>% distinct() 

poc <- read_csv(
    "~/Box Sync/Phyto_bbp/DATA/FINAL/diel_optics_phyto_cultures/analyses/poc/processed_chn.csv"
  ) %>% 
  mutate(poc_optics = case_when(exp == "TP22-14" ~ mean_c_mg_m3  * 0.013,
                                   exp == "SYN22-5" ~ mean_c_mg_m3  * 0.009,
                                   exp == "OL22-3" ~ mean_c_mg_m3  * 0.01),
         sd_poc_optics = case_when(exp == "TP22-14" ~ sd_c_mg_m3  * 0.013,
                                   exp == "SYN22-5" ~ sd_c_mg_m3  * 0.009,
                                   exp == "OL22-3" ~ sd_c_mg_m3  * 0.01),
         pon_optics = case_when(exp == "TP22-14" ~ mean_n_mg_m3  * 0.013,
                                   exp == "SYN22-5" ~ mean_n_mg_m3  * 0.009,
                                   exp == "OL22-3" ~ mean_n_mg_m3  * 0.01),
         sd_pon_optics = case_when(exp == "TP22-14" ~ sd_n_mg_m3  * 0.013,
                                   exp == "SYN22-5" ~ sd_n_mg_m3  * 0.009,
                                   exp == "OL22-3" ~ sd_n_mg_m3  * 0.01)) %>% 
  rename(poc_culture = mean_c_mg_m3,
         pon_culture = mean_n_mg_m3,
         sd_poc_culture = sd_c_mg_m3,
         sd_pon_culture = sd_n_mg_m3) %>% 
  select(exp, tp, poc_culture, sd_poc_culture, pon_culture, sd_pon_culture, poc_optics, sd_poc_optics, pon_optics, sd_pon_optics, mean_cn, sd_cn)

acs_bottle <-
  read_csv("~/Box Sync/Phyto_bbp/DATA/FINAL/diel_optics_phyto_cultures/analyses/acs/processed_acs_bottle_11102023.csv") %>% 
  select(-c(7:11)) %>% 
  left_join(., times) %>% 
  distinct()

fcm_bottle <-
  read_csv(
    "~/Box Sync/Phyto_bbp/DATA/FINAL/diel_optics_phyto_cultures/analyses/fcm/Cultures_and_Optics/processed_fcm_bottle.csv"
  )
bb3_bottle <-
  read_csv("~/Box Sync/Phyto_bbp/DATA/FINAL/diel_optics_phyto_cultures/analyses/bb3/processed_bb3_bottle_11102023.csv") %>% 
  select(-c(ap:bp, date:plot_datetime)) %>%  left_join(., times) %>% 
  distinct()

frr_bottle <-
  read_csv("~/Box Sync/Phyto_bbp/DATA/FINAL/diel_optics_phyto_cultures/analyses/frr/processed_frr_bottle.csv") %>% 
  select(-c(datetime:plot_datetime)) %>% 
  left_join(., times) %>% 
  distinct()

bf_optics <- left_join(acs_bottle , bb3_bottle) %>% 
  left_join(., fcm_bottle %>% filter(source == "Optics")) %>% 
  select(1:6, 15:19, everything()) %>% 
  distinct()

bf_culture <- left_join(fcm_bottle %>% filter(source == "Culture"), frr_bottle) %>% 
  distinct()

bf <- bind_rows(bf_culture, bf_optics) %>% 
  left_join(., poc) %>% 
  arrange(source, datetime) %>% 
  distinct() %>% 
  select(exp:plot_datetime, Fo:Fv_Fo, cells:ssc, poc_culture:sd_cn, everything()) %>% 
  mutate(phyto = case_when(phyto == "italic('O. Lucimarinus')" ~ "italic('O. lucimarinus')",
                           phyto == "italic('Synechococcus ') (WH8102)" ~ "italic('Synechococcus ') (WH8102)",
                           phyto == "italic('T. Pseudonana')" ~ "italic('T. pseudonana')"))
```

    ```

    ```{r,summary file, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
acs_sum <-
  read_csv(
    "~/Box Sync/Phyto_bbp/DATA/FINAL/diel_optics_phyto_cultures/analyses/acs/processed_acs_summary_11102023.csv"
  ) %>% 
  select(-c(6:10)) %>% 
  left_join(., times) %>% 
  distinct()

fcm_sum <-
  read_csv(
    "~/Box Sync/Phyto_bbp/DATA/FINAL/diel_optics_phyto_cultures/analyses/fcm/Cultures_and_Optics/processed_fcm_summary.csv"
  )

bb3_sum <-
  read_csv(
    "~/Box Sync/Phyto_bbp/DATA/FINAL/diel_optics_phyto_cultures/analyses/bb3/processed_bb3_summary_11102023.csv"
  ) %>% 
  select(-c(date:plot_datetime, mean_cp:sd_bb_b)) %>% 
   left_join(., times) %>% 
  distinct()

frr_sum <-
  read_csv(
    "~/Box Sync/Phyto_bbp/DATA/FINAL/diel_optics_phyto_cultures/analyses/frr/processed_frr_summary.csv"
  ) %>% select(-c(datetime:plot_datetime)) %>% 
   left_join(., times) 

sf_optics <- left_join(acs_sum, bb3_sum) %>% 
  left_join(., fcm_sum %>% filter(source == "Optics")) %>% 
  distinct()

sf_culture <- left_join(fcm_sum %>% filter(source == "Culture"), frr_sum) %>% 
  distinct()

sf <- bind_rows(sf_culture, sf_optics) %>% 
  arrange(source, datetime) %>% 
  distinct() %>% 
  left_join(., poc) %>% 
  select(exp:plot_datetime, mean_Fo:sd_Fv_Fo, mean_cells:sd_ssc,  poc_culture:sd_cn,  everything()) %>% 
  mutate(phyto = case_when(phyto == "italic('O. Lucimarinus')" ~ "italic('O. lucimarinus')",
                           phyto == "italic('Synechococcus ') (WH8102)" ~ "italic('Synechococcus ') (WH8102)",
                           phyto == "italic('T. Pseudonana')" ~ "italic('T. pseudonana')"))
```

    ```

########################### 

# SAVE DATA

########################### 

    ```{r,save data}

``` r
write_csv(bf, "FINAL_BOTTLE_11102023.csv")
write_csv(sf, "FINAL_SUMMARY_11102023.csv")
```

    ```

- 

END

- 
