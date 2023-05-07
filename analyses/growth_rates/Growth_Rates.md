Growth Rates
================
Nicholas Baetge
Last compiled on 20 February, 2023

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

    ```{r,import, message=FALSE, warning=FALSE}

``` r
data <-  "cultures.xlsx" %>% 
  readxl::excel_sheets() %>% 
  set_names() %>% 
  map_dfr(readxl::read_excel, path = "cultures.xlsx") %>% 
  mutate(datetime = ymd_hms(paste(date, time)), .after = phyto) %>% 
  select(-date, -time)
```

    ```

    ```{r,import lookup table for division rates, message=FALSE, warning=FALSE}

``` r
lookup <- readxl::read_excel("growth_rate_lookup.xlsx") %>% 
  rename(mu = mew_d)
```

    ```

- 

############# 

# CALCULATE GROWTH RATES

############## 

Here we calculate growth rates as the average daily growth rate from
cell abundances (except for TP22-13 and TP22-14, which were run in
turbidostat mode), from the night time decrease in FSC, and from
dilution volume (for TP22-13 and TP22-14).

    ```{r,calculate growthrates, message=FALSE, warning=FALSE}

``` r
gr <- data %>%
  group_by(phyto, bottle, cycle) %>%
  mutate(
    interv = interval(lag(datetime), datetime),
    duration = as.numeric(interv),
    interv_h = duration / 3600,
    interv_d = as.numeric(interv_h) / 24,
    log_cells = log10(mean_cells_ml),
    mu = ((log_cells - lag(log_cells)) * 2.303) / interv_d) %>% 
  ungroup()
```

    ```

- 

<!-- -->

    ```{r,summarize, message=FALSE, warning=FALSE}

``` r
summary <- gr %>% 
  select(phyto, mu) %>% 
  drop_na(mu) %>% 
  group_by(phyto) %>% 
   summarize(across(
    .cols = mu,
    .fns = list(med = median, mean = mean, sd = sd, count = ~n()),
    # .fns = list(mean = mean, sd = sd),
    na.rm = TRUE,
    .names = "{fn}_{col}"
  ))  

summary
```

    ## # A tibble: 3 × 5
    ##   phyto med_mu mean_mu sd_mu count_mu
    ##   <chr>  <dbl>   <dbl> <dbl>    <int>
    ## 1 ol     0.180   0.362 0.669       22
    ## 2 syn    0.381   0.218 0.658       34
    ## 3 tp     0.902   0.950 0.567       50

    ```

- 

<!-- -->

    ```{r,refer to lookup, message=FALSE, warning=FALSE}

``` r
mean_mu <- summary %>% 
  select(phyto, mean_mu) %>% 
  rename(mu = mean_mu) %>% 
  mutate(mu = round(mu, 2)) 


divisions <- lookup %>% 
  add_row(mu = c(mean_mu$mu)) %>% 
  arrange(mu) %>% 
  mutate(divisions_d = zoo::na.approx(divisions_d, mu),
         doubling_d = zoo::na.approx(doubling_d, mu)) %>% 
  distinct() %>% 
  left_join(mean_mu, .) %>% 
  left_join(., summary %>% select(phyto, count_mu)) %>% 
  mutate(divisions = divisions_d * count_mu)

divisions
```

    ## # A tibble: 3 × 6
    ##   phyto    mu divisions_d doubling_d count_mu divisions
    ##   <chr> <dbl>       <dbl>      <dbl>    <int>     <dbl>
    ## 1 ol     0.36       0.519       1.93       22      11.4
    ## 2 syn    0.22       0.318       3.19       34      10.8
    ## 3 tp     0.95       1.37        0.73       50      68.6

    ```
