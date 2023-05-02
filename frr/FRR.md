FRR
================
Nicholas Baetge
Last compiled on 13 February, 2023

- 

<!-- -->

    ```{r,load packages, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
library(tidyverse)
library(lubridate)
library(hms)
```

    ```

- 

############# 

# IMPORT & WRANGLE DATA

############## 

- 

<!-- -->

    ```{r,read in frr data, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
#save filenames that are in the folder of interest
tp1_files <- fs::dir_ls(path = "TP22-13/", regexp = "fit\\.csv$")
tp2_files <- fs::dir_ls(path = "TP22-14/", regexp = "fit\\.csv$")
syn1_files <- fs::dir_ls(path = "SYN22-4/", regexp = "fit\\.csv$")
syn2_files <- fs::dir_ls(path = "SYN22-5/", regexp = "fit\\.csv$")
ol1_files <- fs::dir_ls(path = "OL22-2/", regexp = "fit\\.csv$")
ol2_files <- fs::dir_ls(path = "OL22-3/", regexp = "fit\\.csv$")
```

    ```

- 

<!-- -->

    ```{r,wrangle frr data and combine, message=FALSE, warning=FALSE}

``` r
#the first two rows of the bb data output contain header and unit data. 

headers <- names(read_csv(tp2_files[1], n_max = 0))

data <-
  readr::read_csv(
    c(tp2_files, syn1_files, syn2_files, ol1_files, ol2_files),
    col_names = headers,
    skip = 2,
    id = "file_name"
  ) %>%
  select(file_name, DATE:QBP_Size) %>%
  mutate(file_name = gsub("frr/", "", file_name),
    file_name = gsub("_fit.csv", "", file_name)) %>%
  separate(file_name,
           sep = "/",
           into = c("exp", "bottle_tp")) %>%
  separate(bottle_tp,
           sep = "_",
           into = c("bottle", "tp")) %>%
  mutate_at(vars(tp), as.numeric) %>%
  arrange(exp, bottle, tp) %>%
  select(1:5, 19:23) %>%
  rename(
    date = DATE,
    time = TIME,
    Fv_Fm = 9,
    Fv_Fo = 10
  ) %>%
  mutate(exp_no = case_when(
    exp %in% c("TP22-13", "SYN22-4", "OL22-2") ~ "1",
    !exp %in% c("TP22-13", "SYN22-4", "OL22-2") ~ "2",
  ),
  .after = exp) %>%
  mutate(
    phyto = case_when(
      exp %in% c("TP22-13", "TP22-14") ~ "italic('T. Pseudonana')",
      exp %in% c("SYN22-4", "SYN22-5") ~ "italic('Synechococcus ') (WH8102)",
      exp %in% c("OL22-2", "OL22-3") ~ "italic('O. Lucimarinus')"
    ),
    .after = exp_no
  ) %>%
  mutate(datetime = ymd_hms(paste(date, time)), .before = date) %>%
  filter(!datetime == ymd_hms("2022-11-30 09:34:26")) %>%
  mutate(
    datetime = floor_date(datetime, "hour"),
    time = as_hms(datetime),
    plot_date = case_when(
      date %in% c(
        ymd("2022-05-09"),
        ymd("2022-05-18"),
        ymd("2022-08-16"),
        ymd("2022-08-24"),
        ymd("2022-11-30"),
        ymd("2022-12-07")
      ) ~ ymd("2022-11-27"),
      !date %in% c(
        ymd("2022-05-09"),
        ymd("2022-05-18"),
        ymd("2022-08-16"),
        ymd("2022-08-24"),
        ymd("2022-11-30"),
        ymd("2022-12-07")
      ) ~ ymd("2022-11-28")
    ),
    plot_datetime = ymd_hms(paste(plot_date, time)),
    .after = time
  ) %>%
  mutate(source = "Culture", .after = phyto) %>%
  mutate_at(vars(exp_no), as.character) %>% 
  filter(!exp == "TP22-4", !bottle == "C", !tp == 8)
```

    ```

- 

########################### 

# SUMMARIZE

########################### 

    ```{r,summarize and merge with fcm data, message=FALSE}

``` r
frr_sum <- data %>%
  select(exp, tp, 12:16) %>%
  group_by(exp, tp) %>%
  summarize(across(
    .cols = Fo:Fv_Fo, 
    .fns = list(mean = mean, sd = sd),
    na.rm = TRUE,
    .names = "{fn}_{col}"
  )) %>%
  ungroup() %>% 
  left_join(data %>% select(1:4, 6:11) %>% distinct(), .) %>% 
  arrange(datetime) %>% 
  filter(!datetime %in% c("2022-05-19 00:00:00", "2022-08-24 14:00:00", "2022-12-01 01:00:00") )
```

    ```

- 

########################### 

# SAVE DATA

########################### 

    ```{r,save data}

``` r
write_csv(data, "processed_frr_bottle.csv")
write_csv(frr_sum, "processed_frr_summary.csv")
```

    ```

- 

END
