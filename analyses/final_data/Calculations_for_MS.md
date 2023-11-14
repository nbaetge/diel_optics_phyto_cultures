Calculations_for_MS
================
Nicholas Baetge
Last compiled on 13 November, 2023

- 

<!-- -->

    ```{r,load packages, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
library(tidyverse)
library(rstatix)
```

    ```

############# 

# IMPORT DATA

############## 

    ```{r,import data,  message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
bf <- read_csv("FINAL_BOTTLE_11102023.csv")
sf <- read_csv("FINAL_SUMMARY_11102023.csv")
gaus <- read_csv("Gaussian_Decompositions.csv")
pigs <- read_csv("Pigment_Ratios.csv")
pf <- read_csv("PAR.csv")
```

    ```

    ```{r,cells data,  message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
cells_data <-  sf %>% 
  select(exp, exp_no, phyto, source, tp, plot_datetime, mean_cells, sd_cells) %>% 
  distinct() %>% 
  filter(!exp == "OL22-3" | !tp == 4) %>% 
  mutate_at(vars(mean_cells, sd_cells), ~ . * 1E6) %>% 
  ungroup()
```

    ```

    ```{r,fcm data, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
fcm_data <- sf %>% 
  select(exp, exp_no, phyto, source, tp, plot_datetime, mean_fsc, mean_ssc, mean_red, mean_yel) %>% 
  distinct() %>% 
  rename_with(~str_remove(., 'mean_')) %>% 
  filter(!exp == "OL22-3" | !tp == 4) %>% 
  filter(!exp == "SYN22-5" | !source == "Culture") %>% 
  group_by(exp, exp_no, phyto, plot_datetime) %>% 
  summarize(across(
    .cols = fsc:yel,
    .fns = list(mean = mean, sd = sd),
    na.rm = TRUE,
    .names = "{fn}_{col}"
  ))  %>% 
  ungroup() %>% 
  group_by(exp) %>% 
  mutate_at(vars(5,7, 9,11), ~ ./(first(.))) %>% 
  ungroup() 
```

    ```

    ```{r,frr data, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
frr_data <- sf %>%
  select(exp, exp_no, phyto, tp, plot_datetime, mean_Fv_Fm) %>% 
  drop_na(mean_Fv_Fm) 
```

    ```

``` r
pigs_data <- pigs %>% 
  filter(spectra %in% c("ppc_chla", "psc_chla", "pub_peb", "pe_chla", "chlb_chla")) %>% 
  mutate(peak = ifelse(peak == Inf, NA, peak)) %>% 
  select(exp:plot_datetime, spectra, peak) %>% 
  pivot_wider(names_from = spectra, values_from = peak)
```

    ```{r,optics data,  message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
optics_data <- bf %>%
  filter(source == "Optics", wl %in% c(470, 532, 660)) %>%
  select(exp, exp_no, phyto, bottle, tp, plot_datetime,  wl, cells, chl_line_height,  poc_optics, pon_optics, cp, ap, bp,  bbp, bb_b) %>% 
  mutate_at(vars(cells), ~ . * 1E6) %>%  # convert cell abundance from cells/ml to m^-3
  mutate_at(vars(poc_optics, pon_optics), ~ . / 1E3) %>% # convert from mg/m3 to g/m3
  mutate_at(vars(wl), as.character) %>%
  distinct() %>% 
  mutate(chl_cell = (chl_line_height/cells) * 10^12, # mg chl/cell to fg chl per cell
         poc_cell = (poc_optics/cells) * 10^12, #g/cell convert to pg c per cell
         pon_cell = (pon_optics/cells) * 10^12, #g/cell convert to pg n per cell
         c_chl = ((poc_optics * 1E3)/chl_line_height), #mg to mg
         a_star = ap / (chl_line_height / 1E3), # m2 per g chl
         sigma_a = (ap / cells)  * 10^12, #m2 per cell * 10^12
         a_poc = ap / poc_optics, #m2 per g C
         c_star = cp /  (chl_line_height / 1E3),
         sigma_c = (cp / cells) * 10^12,
         c_poc = cp / poc_optics,
         bb_star = bbp /  (chl_line_height / 1E3),
         sigma_bb = (bbp / cells) * 10^12,
         bb_poc = bbp / poc_optics
         )

optics_summary <- optics_data %>% 
  group_by(exp, exp_no, phyto, tp, plot_datetime,  wl) %>% 
   summarize(across(
    .cols = cells:bb_poc,
    .fns = list(mean = mean, sd = sd),
    na.rm = TRUE,
    .names = "{fn}_{col}"
  ))  %>% 
     ungroup() %>% 
  arrange(factor(phyto, levels = unique(optics_data$phyto)))

cell_comp_data <-  bf %>%
  filter(source == "Optics") %>%
  select(exp, exp_no, phyto, tp, plot_datetime,  poc_optics, sd_poc_optics,  pon_optics, sd_pon_optics, mean_cn, sd_cn) %>% 
  distinct() %>% 
  left_join(., optics_summary %>% select(exp:plot_datetime,mean_chl_line_height, sd_chl_line_height, mean_chl_cell, sd_chl_cell,mean_c_chl, sd_c_chl, mean_poc_cell, sd_poc_cell, mean_pon_cell, sd_pon_cell)) %>% 
  select(exp:sd_pon_optics, mean_poc_cell, sd_poc_cell, mean_pon_cell, sd_pon_cell, mean_cn:sd_c_chl) %>% 
  distinct()
```

    ```

############# 

# CALCS & STATS

############## 

    ```{r,cell calcs, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
cells_data %>%
  group_by(phyto, exp, source) %>%
  # filter(exp %in% c("OL22-2", "OL22-3")) %>%
  summarize(across(
    .cols = mean_cells,
    .fns = list(min = min, max = max),
    na.rm = TRUE,
    .names = "{fn}_{col}"
  )) %>%
  ungroup() %>%
  view()
```

    ```

    ```{r,frr calcs, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
# frr_data %>% 
#   filter(exp %in% c("OL22-2", "OL22-3")) %>%
#   # view() 
#   group_by(phyto, exp) %>% 
#   summarize(across(
#     .cols = mean_Fv_Fm,
#     .fns = list(min = min, max = max),
#     na.rm = TRUE,
#     .names = "{fn}_{col}"
#   ))  %>% 
#   # ungroup() %>% 
#   group_by(phyto) %>% 
#    summarize(across(
#     .cols = min_mean_Fv_Fm:max_mean_Fv_Fm,
#     .fns = list(mean = mean, sd = sd),
#     na.rm = TRUE,
#     .names = "{fn}_{col}"
#   ))  %>% 
#   ungroup() 
```

    ```

    ```{r,cell composition calcs, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
# cell_comp_data %>%
#   group_by(exp) %>%
#   select(exp:plot_datetime, contains("chl")) %>%
#   distinct() %>%
#   view()
  # filter(tp <= 7) %>%
  # summarize(across(
  #   .cols = mean_chl_cell,
  #   # .fns = list(min = min, max = max),
  #   .fns = list(mean = mean, sd = sd),
  #   na.rm = TRUE,
  #   .names = "{fn}_{col}"
  # ))  %>%
  # ungroup()


# cell_comp_data %>%
#   filter(exp %in% c("SYN22-4", "SYN22-5")) %>%
#   select(exp_no, plot_datetime, mean_chl_cell) %>%
#   distinct() %>%
#   pivot_wider(names_from = exp_no, values_from = mean_chl_cell) %>%
#     drop_na() %>%
#   rename(exp1 = 2,
#          exp2 = 3) %>%
#   mutate(perc_diff = threadr::percentage_change(exp2, exp1)) %>%
#  summarize(across(
#    .cols = perc_diff,
#    # .fns = list(min = min, max = max),
#    .fns = list(mean = mean, sd = sd),
#    na.rm = TRUE,
#    .names = "{fn}_{col}"
#  ))  %>%
#  ungroup()
```

    ```

    ```{r,pigment calcs, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
# pigs_data %>% 
#   group_by(exp) %>% 
#   summarize(across(
#     .cols = ppc_chla:chlb_chla,
#     .fns = list(min = min, max = max),
#     # .fns = list(mean = mean, sd = sd),
#     na.rm = TRUE,
#     .names = "{fn}_{col}"
#   ))  %>%
#   ungroup()
```

    ```

    ```{r, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
optics_summary %>%
  # filter(exp %in% c("OL22-2", "OL22-3")) %>%
  # filter(!between(plot_datetime, as_datetime("2022-11-27 21:00:00 UTC"), as_datetime("2022-11-28 07:00:00 UTC" )) ) %>% 
  select(exp:wl, mean_ap, sd_ap, mean_cp, sd_cp, mean_bbp, sd_bbp, mean_bb_b, sd_bb_b) %>%
   # select(exp:wl,  mean_bbp, sd_bbp, mean_bb_b, sd_bb_b) %>%
  # filter(wl == 470) %>%
  # view()
  group_by(exp, wl) %>%
  # group_by(wl) %>%
  # group_by(tp) %>%
  # filter(between(tp, 5,9)) %>%
  summarize(across(
    .cols = c(mean_ap, mean_cp, mean_bbp, mean_bb_b),
    .fns = list(min = min, max = max),
    # .fns = list(mean = mean, sd = sd),
    na.rm = TRUE,
    .names = "{fn}_{col}"
  ))  %>%
  mutate(min_mean_ap = abs(min_mean_ap)) %>% 
  mutate(ap_fc = max_mean_ap/min_mean_ap,
         cp_fc = max_mean_cp/min_mean_cp,
         bbp_fc = max_mean_bbp/min_mean_bbp,
         ap_pc = threadr::percentage_change(min_mean_ap,max_mean_ap),
         cp_pc = threadr::percentage_change(min_mean_cp,max_mean_cp),
         bbp_pc = threadr::percentage_change(min_mean_bbp,max_mean_bbp))  %>% 
  ungroup() 
```

    ## # A tibble: 15 × 16
    ##    exp     wl    min_mean_ap max_mean_ap min_mean_cp max_mean_cp min_mean_bbp
    ##    <chr>   <chr>       <dbl>       <dbl>       <dbl>       <dbl>        <dbl>
    ##  1 OL22-2  470        0.0979      0.121        0.554       0.731     0.00353 
    ##  2 OL22-2  532        0.0339      0.0427       0.427       0.565     0.00419 
    ##  3 OL22-2  660        0.0268      0.0335       0.291       0.383     0.00292 
    ##  4 OL22-3  470        0.0992      0.175        0.459       0.789     0.00276 
    ##  5 OL22-3  532        0.0330      0.0567       0.345       0.595     0.00322 
    ##  6 OL22-3  660        0.0294      0.0515       0.229       0.393     0.00232 
    ##  7 SYN22-4 470        0.115       0.198        0.649       1.24      0.00164 
    ##  8 SYN22-4 532        0.0505      0.0831       0.501       0.969     0.00237 
    ##  9 SYN22-4 660        0.0271      0.0470       0.330       0.634     0.00261 
    ## 10 SYN22-5 470        0.162       0.278        0.916       1.72      0.00163 
    ## 11 SYN22-5 532        0.0709      0.115        0.696       1.34      0.00264 
    ## 12 SYN22-5 660        0.0409      0.0666       0.441       0.861     0.00299 
    ## 13 TP22-14 470        0.0252      0.0557       0.232       0.570     0.000950
    ## 14 TP22-14 532        0.0141      0.0267       0.215       0.596     0.000850
    ## 15 TP22-14 660        0.0125      0.0266       0.162       0.543     0.000655
    ## # ℹ 9 more variables: max_mean_bbp <dbl>, min_mean_bb_b <dbl>,
    ## #   max_mean_bb_b <dbl>, ap_fc <dbl>, cp_fc <dbl>, bbp_fc <dbl>, ap_pc <dbl>,
    ## #   cp_pc <dbl>, bbp_pc <dbl>

    ```

    ```{r,tidy data for tests, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
stats_data <- optics_summary %>% 
  left_join(., cell_comp_data %>% select(1:5, mean_cn, mean_chl_cell, mean_c_chl)) %>% 
  left_join(., fcm_data %>% select(1:4, mean_fsc, mean_ssc, mean_red)) %>% 
  left_join(., frr_data) %>% 
  left_join(., pigs_data %>% select(1:3, 5, 11:15)) %>% 
  select(-contains("sd")) %>% 
  mutate(phyto = case_when(phyto == "italic('O. lucimarinus')" ~ "o_lucimarinus",
                           phyto ==  "italic('Synechococcus ') (WH8102)" ~ "synechococcus",
                           phyto == "italic('T. pseudonana')" ~ "t_pseudonana"))
```

    ```

    ```{r,summary stats by experiment, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
summary_exp_stats <- stats_data %>% 
  group_by(phyto, exp, wl) %>% 
  get_summary_stats(mean_cells:chlb_chla)
```

    ```

    ```{r,summary stats by phyto and wl, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
summary_phyto_stats <- stats_data %>% 
  group_by(phyto, wl) %>% 
  get_summary_stats(mean_cells:chlb_chla)
```

    ```

``` r
stats_data %>% 
  group_by( wl) %>% 
  get_summary_stats(mean_cells:chlb_chla)
```

    ## # A tibble: 96 × 14
    ##    wl    variable         n     min      max   median       q1       q3      iqr
    ##    <chr> <fct>        <dbl>   <dbl>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl>
    ##  1 470   mean_cells      66 7.25e+9 1.12e+12 1.29e+11 5.04e+10 5.68e+11 5.17e+11
    ##  2 470   mean_chl_li…    67 1.93e+0 5.83e+ 0 3.19e+ 0 2.67e+ 0 4.06e+ 0 1.39e+ 0
    ##  3 470   mean_poc_op…    41 2.8 e-2 4.95e- 1 1.96e- 1 1.08e- 1 3.84e- 1 2.76e- 1
    ##  4 470   mean_pon_op…    41 9   e-3 9   e- 2 3.8 e- 2 1.4 e- 2 7   e- 2 5.5 e- 2
    ##  5 470   mean_cp         67 2.32e-1 1.72e+ 0 7.17e- 1 5.7 e- 1 1.08e+ 0 5.14e- 1
    ##  6 470   mean_ap         67 2.5 e-2 2.78e- 1 1.48e- 1 1.03e- 1 1.75e- 1 7.2 e- 2
    ##  7 470   mean_bp         67 2.07e-1 1.48e+ 0 5.87e- 1 4.86e- 1 9.19e- 1 4.33e- 1
    ##  8 470   mean_bbp        67 1   e-3 5   e- 3 3   e- 3 2   e- 3 4   e- 3 2   e- 3
    ##  9 470   mean_bb_b       67 2   e-3 8   e- 3 4   e- 3 3   e- 3 7   e- 3 4   e- 3
    ## 10 470   mean_chl_ce…    66 3.59e+0 3.03e+ 2 4.13e+ 1 5.43e+ 0 7.14e+ 1 6.59e+ 1
    ## # ℹ 86 more rows
    ## # ℹ 5 more variables: mad <dbl>, mean <dbl>, sd <dbl>, se <dbl>, ci <dbl>

    ```{r,tp stats, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
tp_stats <- stats_data %>% 
   # filter(tp > 2) %>%
   group_by(phyto, wl) %>% 
  filter(phyto == "t_pseudonana") %>%
  cor_test(vars = c("mean_ap", "mean_cp", "mean_bbp", "mean_bb_b"), vars2 = c("mean_ap", "mean_cp", "mean_bbp", "mean_bb_b", "mean_cells", "mean_chl_line_height", "mean_chl_cell", "mean_poc_optics", "mean_poc_cell","mean_cn", "mean_c_chl", "mean_fsc", "mean_ssc", "ppc_chla", "psc_chla", "mean_red", "mean_Fv_Fm"), use = "pairwise.complete.obs", method = "spearman") %>% 
    # filter(p < 0.05) %>%
  arrange(phyto, var1,  var2, wl, cor)
```

    ```

    ```{r,syn stats, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
syn_stats <- stats_data %>% 
  group_by(phyto, wl) %>% 
  filter(phyto == "synechococcus") %>%
  cor_test(vars = c("mean_ap", "mean_cp", "mean_bbp", "mean_bb_b"), vars2 = c("mean_ap", "mean_cp", "mean_bbp", "mean_bb_b", "mean_cells", "mean_chl_line_height", "mean_chl_cell", "mean_poc_optics", "mean_poc_cell","mean_cn", "mean_c_chl", "mean_fsc", "mean_ssc",  "pub_peb", "pe_chla", "mean_red", "mean_Fv_Fm"), use = "pairwise.complete.obs", method = "spearman") %>% 
    # filter(p < 0.05) %>% 
  arrange(phyto, var1,  var2, wl, cor)
```

    ```

    ```{r,ol stats, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
ol_stats <- stats_data %>% 
  group_by(phyto, wl) %>% 
  filter(phyto == "o_lucimarinus") %>%
  cor_test(vars = c("mean_ap", "mean_cp", "mean_bbp", "mean_bb_b"), vars2 = c("mean_ap", "mean_cp", "mean_bbp", "mean_bb_b", "mean_cells", "mean_chl_line_height", "mean_chl_cell", "mean_poc_optics", "mean_poc_cell","mean_cn", "mean_c_chl", "mean_fsc", "mean_ssc", "ppc_chla", "chlb_chla", "mean_red", "mean_Fv_Fm"), use = "pairwise.complete.obs", method = "spearman") %>%
  # filter(p < 0.05) %>% 
  arrange(phyto, var1,  var2, wl, cor)
```

    ```

    ```{r,combine stats, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
levels = c("ap", "cp", "bbp", "bb_b")

levels2 = c("ap", "cp", "bbp", "bb_b",  "cells", "fsc", "ssc", "poc_optics", "poc_cell", "cn", "c_chl", "chl_line_height", "chl_cell", "ppc_chla", "psc_chla", "chlb_chla", "pe_chla", "pub_peb", "Fv_Fm", "red")


stats_table <- bind_rows(tp_stats, syn_stats, ol_stats) %>% 
  # filter(!var1 == var2) %>%
  mutate_if(., 
                is.character, 
                str_remove_all, 
                pattern = "mean_") %>% 
   arrange(phyto, factor(var1, levels = levels), factor(var2, levels = levels2), wl, cor) 
  #   filter(!var1 == "ap" | !var2 %in%c("c_chl", "chl_line_height", "chl_cell", "ppc_chla", "psc_chla", "chlb_chla", "pe_chla", "pub_peb")) %>% 
  # filter(!var1 == "bb_b" | !var2 == "cp") %>% 
  # filter(!var1 == "cp" | !var2 == "bb_b") 
```

    ```

    ```{r,save stats, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
write_csv(stats_table, "Correlations.csv")
```

    ```
