POC
================
Nicholas Baetge
Last compiled on 20 April, 2023

- 

<!-- -->

    ```{r,load packages, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
library(tidyverse)
library(lubridate)
library(hms)
```

    ```

############# 

# IMPORT & WRANGLE DATA

############## 

- 

<!-- -->

    ```{r,import, message=FALSE}

``` r
volumes <- read_csv("poc_sample_list.csv") %>% 
  select(exp, id, vol)

run_list <- "poc_run_list.xlsx" %>% 
  readxl::excel_sheets() %>% 
  set_names() %>% 
  map_dfr(readxl::read_excel, path = "poc_run_list.xlsx") %>% 
  drop_na()

data <-  "CHN_results.xlsx" %>% 
  readxl::excel_sheets() %>% 
  set_names() %>% 
  map_dfr(readxl::read_excel, path = "CHN_results.xlsx") %>% 
  select(-c(1,3, 6)) %>% 
  rename(run_sample = 1,
         c_ug = 2,
         n_ug = 3)

combined <- left_join(run_list, data) %>% 
  left_join(., volumes) %>% 
  arrange(exp, bottle, tp, sample) 
```

    ```

- 

############# 

# BLANK CORRECTIONS & DERIVED VARS

############## 

    ```{r,subtract blank values, message=FALSE}

``` r
blanks <- combined %>%
  filter(sample == "blank") %>% 
  drop_na(c_ug) %>% 
  group_by(exp) %>% 
  summarize(across(
    .cols = c_ug:n_ug,
    .fns = list(mean = mean, sd = sd),
    na.rm = TRUE,
    .names = "{fn}_{'blank'}_{col}"
  ))  %>% 
  ungroup()


blank_corrected <- combined %>% 
  filter(!sample == "blank") %>% 
  drop_na(c_ug) %>% 
  group_by(exp, bottle, tp, sample) %>% 
  filter(!run_sample %in% c("1_13", "1_14", "2_13")) %>% 
  mutate(c_ug = case_when(exp == "TP22-14" ~ c_ug,
                            exp == "SYN22-5" ~ c_ug,
                            exp == "OL22-3" ~ sum(c_ug)),
         n_ug = case_when(exp == "TP22-14" ~ n_ug,
                            exp == "SYN22-5" ~ n_ug,
                            exp == "OL22-3" ~ sum(n_ug)),
         ) %>% 
  ungroup() %>% 
  select(exp, bottle, tp, sample, vol, c_ug, n_ug) %>% 
  distinct() %>% 
  mutate(c_corr = case_when(exp == "TP22-14" ~ c_ug - blanks$mean_blank_c_ug[blanks$exp == "TP22-14"],
                            exp == "SYN22-5" ~ c_ug - blanks$mean_blank_c_ug[blanks$exp == "SYN22-5"],
                            exp == "OL22-3" ~ c_ug - blanks$mean_blank_c_ug[blanks$exp == "OL22-3"]
                            ),
         n_corr = case_when(exp == "TP22-14" ~ n_ug - blanks$mean_blank_n_ug[blanks$exp == "TP22-14"],
                            exp == "SYN22-5" ~ n_ug - blanks$mean_blank_n_ug[blanks$exp == "SYN22-5"],
                            exp == "OL22-3" ~ n_ug - blanks$mean_blank_n_ug[blanks$exp == "OL22-3"]
                            ), 
         c_mg = c_corr / 10^3,
         n_mg = n_corr / 10^3,
         vol_m3 = vol / 10^6,
         c_mg_m3 = c_mg / vol_m3,
         n_mg_m3 = n_mg / vol_m3,
         cn = (c_mg /12000) / (n_mg / 14000)
         ) %>% 
  arrange(exp, tp, bottle, sample) %>% 
  mutate(id = paste(exp, bottle, tp, sample, sep = "_")) %>% 
  filter(!id %in% c( "OL22-3_D_8_a", "OL22-3_B_9_b", "OL22-3_D_9_c", "OL22-3_D_10_a", "OL22-3_D_10_b", "OL22-3_D_11_c", "OL22-3_D_12_c")) %>% 
  select(-id) %>% 
  mutate(id = paste(exp, tp,  sep = "_")) %>% 
  filter(!id == "OL22-3_3") %>% 
  select(-id)
```

    ```

- 

<!-- -->

    ```{r,summarize data , message=FALSE}

``` r
chn <- blank_corrected %>% 
  # filter(exp == "TP22-14") %>% 
  arrange(bottle, tp, sample) %>% 
  select(exp, bottle, tp, c_mg, n_mg, c_mg_m3, n_mg_m3, cn) %>% 
  group_by(exp, tp) %>% 
  summarize(across(
    .cols = c_mg:cn,
    .fns = list(mean = mean, sd = sd),
    na.rm = TRUE,
    .names = "{fn}_{col}"
  ))  %>% 
  ungroup()
```

    ```

- 

########################### 

# SAVE DATA

########################### 

    ```{r,save data}

``` r
write_csv(chn, "processed_chn.csv")
```

    ```

- 

END
