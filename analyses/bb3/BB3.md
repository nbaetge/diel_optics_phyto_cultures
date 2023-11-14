BB3
================
Nicholas Baetge
Last compiled on 14 November, 2023

- 

<!-- -->

    ```{r,load packages, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
library(tidyverse)
library(lubridate)
library(hms)
library(patchwork)
```

    ```

- 

<!-- -->

    ```{r,custom theme for plots, wrapper = TRUE}

``` r
custom.theme <- theme(
  legend.position = "top",
  legend.title = element_text(size = 23),
  legend.key.size = unit(0.7, "cm"),
  legend.text = element_text(size = 23),
  axis.title = element_text(size = 23, face = "bold"),
  panel.spacing.x = unit(2, "cm"),
  strip.text.x = element_text(size = 23, color = "black", face = "bold"),
  strip.text.y = element_text(size = 23, color = "black", face = "bold"),
  strip.background = element_rect(
    color = "black",
    fill = alpha('light grey', 0.4),
    linewidth = 1.5,
    linetype = "solid"
  )
)
```

    ```

- 

############# 

# IMPORT & WRANGLE DATA

############## 

- 

<!-- -->

    ```{r,read in acs data for bb derivations, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
#data for bb wavelengths
acs_bottle <-
  read_csv("~/Box Sync/Phyto_bbp/DATA/FINAL/diel_optics_phyto_cultures/analyses/acs/processed_acs_bottle_11102023.csv") %>% 
  filter(wl %in% c(470, 532, 660)) %>%
  select(-c(c_se, a_se))
```

    ```

    ```{r,read in acs and fcm data, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
#summary files
acs_summary <- read_csv("~/Box Sync/Phyto_bbp/DATA/FINAL/diel_optics_phyto_cultures/analyses/acs/processed_acs_summary_11102023.csv") %>% select(-datetime)
fcm_summary <- read_csv("~/Box Sync/Phyto_bbp/DATA/FINAL/diel_optics_phyto_cultures/analyses/fcm/Cultures_and_Optics/processed_fcm_summary.csv") %>% 
  filter(source == "Optics") %>% 
  select(-datetime)

acs_fcm <- left_join(acs_summary, fcm_summary %>% select(1:5, 10:21))

times <- acs_fcm %>% select(1:9) %>% distinct() %>% 
  group_by(exp, tp) %>% 
  mutate(tf = ifelse(time == first(time), T, F)) %>% 
  filter(tf == T) %>% 
  select(-tf) %>% 
  ungroup()
```

    ```

    ```{r,read in bb3 data, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
#save filenames that are in the folder of interest
tp2_files <- fs::dir_ls(path = "TP22-14/", regexp = "\\.csv$")
syn1_files <- fs::dir_ls(path = "SYN22-4/", regexp = "\\.csv$")
syn2_files <- fs::dir_ls(path = "SYN22-5/", regexp = "\\.csv$")
ol1_files <- fs::dir_ls(path = "OL22-2/", regexp = "\\.csv$")
ol2_files <- fs::dir_ls(path = "OL22-3/", regexp = "\\.csv$")
```

    ```

- 

<!-- -->

    ```{r,wrangle bb data and combine, message=FALSE, warning=FALSE}

``` r
#the first two rows of the bb data output contain header and unit data. 

headers <-
  c("exp",
    "bottle",
    "tp",
    "bb",
    "date",
    "start",
    colnames(read_csv(tp2_files[1])))

tp2 <- tp2_files %>%
  purrr::map_dfr(read.csv,
                 skip = 2,
                 header = F,
                 .id = "source") %>%
  #shorten the filenames
  mutate(source = gsub(".csv", "", source)) %>%
  separate(source,
           sep = "_",
           into = c("exp", "tp", "bb", "date", "start")) %>%
  separate(exp,
           sep = "/",
           into = c("exp", "bottle")) %>%
  arrange(V1) %>%
  rename_all( ~ headers) 

syn1 <- syn1_files %>%
  purrr::map_dfr(read.csv,
                 skip = 2,
                 header = F,
                 .id = "source") %>%
  mutate(source = gsub(".csv", "", source)) %>%
  separate(source,
           sep = "_",
           into = c("exp", "tp", "bb", "date", "start")) %>%
  separate(exp,
           sep = "/",
           into = c("exp", "bottle")) %>%
  arrange(V1) %>%
  rename_all( ~ headers) 

syn2 <- syn2_files %>%
  purrr::map_dfr(read.csv,
                 skip = 2,
                 header = F,
                 .id = "source") %>%
  mutate(source = gsub(".csv", "", source)) %>%
  separate(source,
           sep = "_",
           into = c("exp", "tp", "bb", "date", "start")) %>%
  separate(exp,
           sep = "/",
           into = c("exp", "bottle")) %>%
  arrange(V1) %>%
  rename_all( ~ headers) 

ol1 <- ol1_files %>%
  purrr::map_dfr(read.csv,
                 skip = 2,
                 header = F,
                 .id = "source") %>%
  mutate(source = gsub(".csv", "", source)) %>%
  separate(source,
           sep = "_",
           into = c("exp", "tp", "bb", "date", "start")) %>%
  separate(exp,
           sep = "/",
           into = c("exp", "bottle")) %>%
  arrange(V1) %>%
  rename_all( ~ headers) 

ol2 <- ol2_files %>%
  purrr::map_dfr(read.csv,
                 skip = 2,
                 header = F,
                 .id = "source") %>%
  mutate(source = gsub(".csv", "", source)) %>%
  separate(source,
           sep = "_",
           into = c("exp", "tp", "bb", "date", "start")) %>%
  separate(exp,
           sep = "/",
           into = c("exp", "bottle")) %>%
  arrange(V1) %>%
  rename_all( ~ headers )

data <- bind_rows(tp2, syn1, syn2, ol1, ol2) %>%
  select(-c(bb, date, start, time)) %>%
  rename_all(~stringr::str_replace(.,"beta","")) %>% 
  pivot_longer(cols = c(4:6),
               names_to = "wl",
               values_to = "count") %>% 
  mutate_at(vars(tp, wl), as.numeric) %>% 
  group_by(exp, bottle, tp, wl) %>%
  mutate(med = median(count),
         se = sd(count) / sqrt(length(count))) %>%
  ungroup() %>%
  select(-count) %>% 
  distinct() %>% 
  left_join(times, .) %>% 
  mutate(
    exp2 = case_when(
      exp == "TP22-14" ~ "TP-2",
      exp == "SYN22-4" ~ "SYN-1",
      exp == "SYN22-5" ~ "SYN-2",
      exp == "OL22-2" ~ "OL-1",
      exp == "OL22-3" ~ "OL-2"
    ),
    .after = "exp"
  ) %>% 
  arrange(date, time, wl)
```

    ```

- 

<!-- -->

    ```{r,count plot, fig.height=10, fig.width=16, message=FALSE, warning=FALSE}

``` r
data %>% 
ggplot(aes(
    x = wl,
    y = med,
    color = bottle,
    group = interaction(exp, bottle, tp))) +
 labs(x = expression(bold(Wavelength ~ (nm ^ -1))),
       y = expression(bold(Count)),
       fill = expression(bold(Bottle))) +
  geom_errorbar(aes(ymin = med - se, ymax = med + se),
                 width = 10, alpha = 0.7) +
  geom_line(size = 1, alpha = 0.7) +
  geom_point(shape = 21, size = 4, alpha = 0.7
  ) +
  facet_wrap(~ factor(exp2, levels = unique(data$exp2))) +
  guides(color = "none") +
  theme_linedraw(21) +
  custom.theme
```

![](BB3_files/figure-gfm/count%20plot-1.png)<!-- -->

    ```

- 

### blank subtract

    ```{r,identify sample times, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
sample_times <- data %>%
  filter(!grepl('FSW', bottle)) %>%
  select(exp, plot_datetime) %>%
  distinct() %>%
  mutate(wl = 1) %>%
  group_by(exp, plot_datetime) %>%
  group_modify( ~ add_row(., wl = c(470, 532, 660))) %>%
  ungroup() %>%
  filter(wl != 1)
```

    ```

- 

<!-- -->

    ```{r,interpolate filtrate values for sample times, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
filt <- data %>%
  filter(grepl('FSW', bottle)) %>% 
  select(exp, plot_datetime, wl, med) %>%
  full_join(., sample_times) %>%
  arrange(exp, wl, plot_datetime) %>%
  group_by(exp, wl) %>%
  mutate(med = zoo::na.approx(med, na.rm = T)) %>%
  ungroup() %>%
  rename(filt = med)
```

    ```

- 

<!-- -->

    ```{r,subtract filter values from sample values, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
filt_subtract <- data %>%
  filter(!grepl('FSW', bottle)) %>%
  left_join(., filt) %>%
  mutate(count = med - filt, .after = "wl") %>%
  select(-c(med, filt)) %>%
  rename(count_se = se)
```

    ```

- 

########################### 

# SPECTRA SUBTRACTION

########################### 

- 

<!-- -->

    ```{r,correct spectra to remove imprint of previous sample}

``` r
spectra_subtract <- filt_subtract %>%
  arrange(exp, tp, wl) %>%
  group_by(exp, tp, wl) %>%
  mutate(
    count_corr = ifelse(bottle == "A", count, NA),
    count_corr = ifelse(bottle == "B", count[bottle == "B"] - count[bottle == "A"], count_corr),
    count_corr = ifelse(bottle == "C", count[bottle == "C"] - count[bottle == "B"], count_corr)
  ) %>%
  ungroup() %>%
  select(1:12, count_corr, count_se) %>%
  rename(count = count_corr) %>%
  arrange(date, time,  wl) 
```

    ```

- 

########################### 

# DERIVE BETA FROM COUNTS

########################### 

- 

<!-- -->

    ```{r,import calibration slopes and caculate beta, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
slopes <-
  read_csv("calibrations/processed_slopes_112322.csv") %>%
  mutate(year = year(date)) %>%
  filter(year == 2022) %>% 
  select(wl, calib) %>%
  group_by(wl) %>% 
   summarize(across(
    .cols = calib,
    .fns = list(mean = mean),
    na.rm = TRUE,
    .names = "slope"
  )) 

beta <- spectra_subtract %>%
  left_join(., slopes) %>%
  mutate(beta = count * slope) 
```

    ```

- 

########################### 

# BETA to BBP

########################### 

    ```{r,interpolate chi for 124˚ from Sullivan and Twardowski, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
L = 0.0391
E = 0.4

#values from sullivan and twardowski 2009
angle = c(90, 100, 110, 120, 130, 140, 150, 160, 170)
chi = c(0.684, 0.858, 1.00, 1.097, 1.153, 1.167, 1.156, 1.131, 1.093)
chi_sd = c(0.034, 0.032, 0.026, 0.032, 0.044, 0.049, 0.054, 0.054, 0.057)

chi_table <- tibble(angle, chi, chi_sd) %>%
  add_row(angle = 124) %>%
  arrange(angle) %>%
  mutate(chi = zoo::na.approx(chi, angle))

chi = as.numeric(chi_table %>% filter(angle == 124) %>% select(chi))
```

    ```

- 

<!-- -->

    ```{r,calculate bbp, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
bbp <- left_join(beta, acs_bottle %>% select(exp, bottle, tp, wl, ap, cp, bp)) %>% 
  mutate(
    beta_corr = beta * exp(L * (ap + E * bp)),
    bbp = 2 * pi * chi * beta_corr,
    bb_b = bbp/bp
  )  %>% 
  mutate(id = paste(exp, bottle, tp, sep = "_")) %>%
  filter(!id %in% c("SYN22-4_B_3", "OL22-3_C_5")) %>% 
  select(-exp2, -id)
```

    ```

- 

########################### 

# SUMMARIZE

########################### 

    ```{r,summarize and merge with fcm data, message=FALSE}

``` r
bb3_sum <- bbp %>%
  select(exp, tp, wl, bbp, bb_b) %>%
  mutate_at(vars(tp), as.numeric) %>%
  drop_na(bbp) %>% 
  group_by(exp, tp, wl) %>%
  summarize(across(
    .cols = bbp:bb_b, 
    .fns = list(mean = mean, sd = sd),
    na.rm = TRUE,
    .names = "{fn}_{col}"
  )) %>%
  ungroup() %>% 
  left_join(acs_fcm, .) %>% 
  drop_na(mean_bbp)
```

    ```

########################### 

# SAVE DATA

########################### 

    ```{r,save data}

``` r
write_csv(bbp, "processed_bb3_bottle_11102023.csv")
write_csv(bb3_sum, "processed_bb3_summary_11102023.csv")
```

    ```

- 

END
