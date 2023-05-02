ACS
================
Nicholas Baetge
Last compiled on 21 April, 2023

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

    ```{r,read in acs data, message=FALSE, warning=FALSE, wrapper = TRUE}

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

    ```{r,wrangle acs data and combine, message=FALSE, warning=FALSE}

``` r
#the first two rows of the acs data output contain header and unit data. we'll store those separately...units first and we'll circle back to headers
units <- colnames(read_csv(tp2_files[1], skip = 1))

headers <-
  c("exp",
    "bottle",
    "tp",
    "ac",
    "date",
    "start",
    colnames(read_csv(tp2_files[1])))

#extract the wavelengths from the "a" column and store as a vector (these are also the same for c)
wavelengths <-
  units[4] %>% strsplit(., split = " ") %>% .[[1]] %>%  str_remove(., "1/m  lambda=")
wl <- wavelengths[1:80]

#store all data from all files one data frame (this is what the map_dfr function does). We'll add the filenames to the dataframe as well, in a column called "source.

tp2 <- tp2_files %>%
  purrr::map_dfr(read.csv,
                 skip = 2,
                 header = F,
                 .id = "source") %>%
  #shorten the filenames
  mutate(source = gsub(".csv", "", source)) %>%
  separate(source,
           sep = "_",
           into = c("exp", "tp", "ac", "date", "start")) %>%
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
  #shorten the filenames
  mutate(source = gsub(".csv", "", source)) %>%
  separate(source,
           sep = "_",
           into = c("exp", "tp", "ac", "date", "start")) %>%
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
  #shorten the filenames
  mutate(source = gsub(".csv", "", source)) %>%
  separate(source,
           sep = "_",
           into = c("exp", "tp", "ac", "date", "start")) %>%
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
  #shorten the filenames
  mutate(source = gsub(".csv", "", source)) %>%
  separate(source,
           sep = "_",
           into = c("exp", "tp", "ac", "date", "start")) %>%
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
  #shorten the filenames
  mutate(source = gsub(".csv", "", source)) %>%
  separate(source,
           sep = "_",
           into = c("exp", "tp", "ac", "date", "start")) %>%
  separate(exp,
           sep = "/",
           into = c("exp", "bottle")) %>%
  arrange(V1) %>%
  rename_all( ~ headers)

data <- bind_rows(tp2, syn1, syn2, ol1, ol2)
```

    ```

- 

## attenuation data

    ```{r,extract nested attenuation data, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
c <- data %>%
  select(c) %>%
  mutate(c = str_squish(c)) %>%
  #remove the brackets
  mutate(c = str_replace_all(c, "\\[|\\]", "")) %>%
  mutate(c = str_squish(c)) %>%
  separate(c, sep = " ", into = wavelengths) %>%
  mutate_all(as.numeric)
```

    ```

- 

<!-- -->

    ```{r,wrangle attenuation data, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
c_data <- bind_cols(data, c) %>%
  select(-c, -a) %>%
  group_by(exp, bottle, tp) %>%
  mutate(
    time = ymd_hms(time),
    time = with_tz(time, tzone = "America/Los_Angeles"),
    date = ymd(date),
    start = first((round_date(time, unit = "15 sec"))),
    start_time = as_hms(start), 
    date = as_date(start),
    .after = "ac"
  ) %>%
  ungroup() %>%
  arrange(date, start) %>%
  select(c(1:3, 7, 6, 5, 13:length(.))) %>%
  rename(datetime = 4,
         time = 6)
```

    ```

- 

### filter events

- 

<!-- -->

    ```{r,wrangle attenuation data for filter events, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
c_filt <- c_data %>%
  filter(grepl('FSW', bottle)) %>%
  pivot_longer(cols = c(7:87),
               names_to = "wl",
               values_to = "c") %>%
  mutate_at(vars(wl, c), as.numeric) %>%
  group_by(exp, bottle, tp, datetime, wl) %>%
  mutate(med = median(c),
         se = sd(c) / sqrt(length(c))) %>%
  ungroup() %>%
  select(-c) %>%
  distinct() %>%
  mutate(
    exp2 = case_when(
      exp == "TP22-14" ~ "TP-2",
      exp == "SYN22-4" ~ "SYN-1",
      exp == "SYN22-5" ~ "SYN-2",
      exp == "OL22-2" ~ "OL-1",
      exp == "OL22-3" ~ "OL-2"
    ),
    .after = "exp"
  ) 
```

    ```

- 

<!-- -->

    ```{r,filtrate attenuation plot, fig.height=10, fig.width=16, message=FALSE, warning=FALSE}

``` r
c_filt %>%
  ggplot(aes(
    x = wl,
    y = med,
    color = as.numeric(tp),
    group = interaction(exp, tp)
  )) +
  labs(x = expression(bold(Wavelength ~ (nm ^ -1))),
       y = expression(bold(c[filt] ~ (m ^ -1))),
       color = expression(bold(Time ~ point))) +
  geom_line() +
  geom_errorbar(aes(ymin = med - se, ymax = med + se)) +
  facet_wrap(~ factor(exp2, levels = unique(c_filt$exp2))) +
  scale_color_gradient2(
    low = scales::muted("red"),
    mid = "white",
    high = scales::muted("blue"),
    midpoint = 8,
    na.value = "grey50",
    guide = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black"
    )
  ) +
  
  theme_linedraw(21) +
  custom.theme +
  theme(legend.key.width = unit(2.5, "cm"))
```

![](ACS_files/figure-gfm/filtrate%20attenuation%20plot-1.png)<!-- -->

    ```

- 

### blank subtract

    ```{r,identify sample times, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
sample_times <- c_data %>%
  filter(!grepl('FSW', bottle)) %>%
  select(exp, datetime) %>%
  distinct() %>%
  mutate(wl = 1) %>%
  group_by(exp, datetime) %>%
  group_modify( ~ add_row(., wl = c(as.numeric(wavelengths)))) %>%
  ungroup() %>%
  filter(wl != 1)
```

    ```

- 

<!-- -->

    ```{r,interpolate filtrate values for sample times, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
c_filt_interp <- c_filt %>%
  select(exp, datetime, wl, med) %>%
  full_join(., sample_times) %>%
  arrange(exp, wl, datetime) %>%
  group_by(exp, wl) %>%
  mutate(med = zoo::na.approx(med, na.rm = T)) %>%
  ungroup() %>%
  rename(c_filt = med)
```

    ```

- 

<!-- -->

    ```{r,subtract filter values from sample values, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
c_filt_subtract <- c_data %>%
  filter(!grepl('FSW', bottle)) %>%
  pivot_longer(cols = c(7:87),
               names_to = "wl",
               values_to = "c") %>%
  mutate_at(vars(wl, c), as.numeric) %>%
  group_by(exp, bottle, tp, datetime, wl) %>%
  arrange(exp, bottle, tp, wl) %>%
  mutate(med = median(c),
         se = sd(c) / sqrt(length(c))) %>%
  ungroup() %>%
  select(-c) %>%
  distinct() %>%
  left_join(., c_filt_interp) %>%
  mutate(c = med - c_filt, .after = "wl") %>%
  select(-c(med, c_filt)) %>%
  rename(c_se = se)
```

    ```

- 

## absorption data

    ```{r,extract nested absorption data, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
a <- data %>%
  select(a) %>%
  mutate(a = str_squish(a)) %>%
  #remove the brackets
  mutate(a = str_replace_all(a, "\\[|\\]", "")) %>%
  mutate(a = str_squish(a)) %>%
  separate(a, sep = " ", into = wavelengths) %>%
  mutate_all(as.numeric)
```

    ```

- 

<!-- -->

    ```{r,wrangle absorption data, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
a_data <- bind_cols(data, a) %>%
  select(-c, -a) %>%
  group_by(exp, bottle, tp) %>%
  mutate(
    time = ymd_hms(time),
    time = with_tz(time, tzone = "America/Los_Angeles"),
    date = ymd(date),
    start = first((round_date(time, unit = "15 sec"))),
    start_time = as_hms(start),
    date = as_date(start),
    .after = "ac"
  ) %>%
  ungroup() %>%
  arrange(date, start) %>%
  select(c(1:3, 7, 6, 5, 13:length(.))) %>%
  rename(datetime = 4,
         time = 6)
```

    ```

- 

### filter events

- 

<!-- -->

    ```{r,wrangle absorption data for filter events, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
a_filt <- a_data %>%
  filter(grepl('FSW', bottle)) %>%
  pivot_longer(cols = c(7:87),
               names_to = "wl",
               values_to = "a") %>%
  mutate_at(vars(wl, a), as.numeric) %>%
  group_by(exp, bottle, tp, datetime, wl) %>%
  mutate(med = median(a),
         se = sd(a) / sqrt(length(a))) %>%
  ungroup() %>%
  select(-a) %>%
  distinct() %>%
  mutate(
    exp2 = case_when(
      exp == "TP22-14" ~ "TP-2",
      exp == "SYN22-4" ~ "SYN-1",
      exp == "SYN22-5" ~ "SYN-2",
      exp == "OL22-2" ~ "OL-1",
      exp == "OL22-3" ~ "OL-2"
    ),
    .after = "exp"
  ) 
```

    ```

- 

<!-- -->

    ```{r,filtrate absorption plot, fig.height=10, fig.width=16, message=FALSE, warning=FALSE}

``` r
a_filt %>%
  ggplot(aes(
    x = wl,
    y = med,
    color = as.numeric(tp),
    group = interaction(exp, tp)
  )) +
  labs(x = expression(bold(Wavelength ~ (nm ^ -1))),
       y = expression(bold(a[filt] ~ (m ^ -1))),
       color = expression(bold(Time ~ point))) +
  geom_line() +
  geom_errorbar(aes(ymin = med - se, ymax = med + se)) +
  facet_wrap(~ factor(exp2, levels = unique(c_filt$exp2))) +
  scale_color_gradient2(
    low = scales::muted("red"),
    mid = "white",
    high = scales::muted("blue"),
    midpoint = 8,
    na.value = "grey50",
    guide = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black"
    )
  ) +
  theme_linedraw(21) +
  custom.theme +
  theme(legend.key.width = unit(2.5, "cm"))
```

![](ACS_files/figure-gfm/filtrate%20absorption%20plot-1.png)<!-- -->

    ```

- 

### blank subtract

- 

<!-- -->

    ```{r,interpolate filtrate absorption values for sample times, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
a_filt_interp <- a_filt %>%
  select(exp, datetime, wl, med) %>%
  full_join(., sample_times) %>%
  arrange(exp, wl, datetime) %>%
  group_by(exp, wl) %>%
  mutate(med = zoo::na.approx(med, na.rm = T)) %>%
  ungroup() %>%
  rename(a_filt = med)
```

    ```

- 

<!-- -->

    ```{r,subtract filter absorption values from sample values, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
a_filt_subtract <- a_data %>%
  filter(!grepl('FSW', bottle)) %>%
  pivot_longer(cols = c(7:87),
               names_to = "wl",
               values_to = "a") %>%
  mutate_at(vars(wl, a), as.numeric) %>%
  group_by(exp, bottle, tp, datetime, wl) %>%
  arrange(exp, bottle, tp, wl) %>%
  mutate(med = median(a),
         se = sd(a) / sqrt(length(a))) %>%
  ungroup() %>%
  select(-a) %>%
  distinct() %>%
  left_join(., a_filt_interp) %>%
  mutate(a = med - a_filt, .after = "wl") %>%
  select(-c(med, a_filt)) %>%
  rename(a_se = se)
```

    ```

- 

<!-- -->

    ```{r,combine ap and cp data and calc b, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
acs_filt_subtract <- left_join(c_filt_subtract, a_filt_subtract) %>%
  mutate(b = c - a) %>%
  mutate(id = paste(exp, bottle, sep = "_"), .before = "exp")
```

    ```

- 

########################### 

# TEMP AND SCATTERING CORRECTION

########################### 

This is done using the MATLAB script “rd_data.m” written by Emmanuel
Boss -files needed are ResidualTemperatureAndScatteringCorrection.m and
Sullivan_etal_2006_instrumentspecific.xls

- 

<!-- -->

    ```{r,export data for correction}

``` r
acs_filt_subtract %>% group_by(id) %>%
  group_walk( ~ write_csv(.x, paste0(.y$id, ".csv")))
```

    ```

- 

<!-- -->

    ```{r,import and wrangle corrected data, message=FALSE}

``` r
corrected_files <-
  fs::dir_ls(path = "corrected/", regexp = "\\.csv$")

acs_corrected <- corrected_files %>%
  purrr::map_dfr(read.csv,
                 skip = 1,
                 header = F,
                 .id = "source") %>%
  #shorten the filenames
  mutate(source = gsub(".csv", "", source),
         source = gsub("corrected/", "", source))  %>%
  separate(source,
           sep = "_",
           into = c("exp", "bottle", "coeff"))  %>%
  group_by(exp, bottle, coeff) %>%
  mutate(tp = row_number(), .after = "bottle") %>%
  mutate(tp = ifelse(exp == "TP22-14" & tp == 12, 13, tp),
         tp = ifelse(exp == "TP22-14" & tp == 11, 12, tp),
         tp = ifelse(exp == "TP22-14" & tp == 10, 11, tp)) %>% 
  ungroup() %>%
  rename_all( ~ c("exp", "bottle", "tp", "coeff", wl)) %>%
  pivot_longer(5:84, names_to = "wl", values_to = "value") %>%
  pivot_wider(names_from = "coeff", values_from = "value") %>%
  mutate(bp = cp - ap) %>%
  mutate_at(vars(tp), as.character) %>%
  mutate_at(vars(wl), as.numeric) %>%
  left_join(.,
            acs_filt_subtract %>% select(exp, bottle, tp, datetime, date, time, wl, c_se, a_se)) %>%
  select(1:3, 8:9, 4, 6, 11, 5, 12, 7) %>%
  ungroup() %>% 
  arrange(datetime, wl) %>%
  mutate(
    exp2 = case_when(
      exp == "TP22-14" ~ "TP-2",
      exp == "SYN22-4" ~ "SYN-1",
      exp == "SYN22-5" ~ "SYN-2",
      exp == "OL22-2" ~ "OL-1",
      exp == "OL22-3" ~ "OL-2"
    ),
    .after = "exp"
  ) 
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
spectra_subtract <- acs_corrected %>%
  arrange(exp, tp, wl) %>%
  group_by(exp, tp, wl) %>%
  mutate(
    c_corr = ifelse(bottle == "A", cp, NA),
    c_corr = ifelse(bottle == "B", cp[bottle == "B"] - cp[bottle == "A"], c_corr),
    c_corr = ifelse(bottle == "C", cp[bottle == "C"] - cp[bottle == "B"], c_corr),
    a_corr = ifelse(bottle == "A", ap, NA),
    a_corr = ifelse(bottle == "B", ap[bottle == "B"] - ap[bottle == "A"], a_corr),
    a_corr = ifelse(bottle == "C", ap[bottle == "C"] - ap[bottle == "B"], a_corr),
    b_corr = c_corr - a_corr
  ) %>%
  ungroup() %>%
  select(1:7, c_corr, c_se, a_corr, a_se, b_corr) %>%
  rename(cp = 8,
         ap = 10,
         bp = 12) %>%
  arrange(datetime, wl) %>% 
  mutate(id = paste(exp, bottle, tp, sep = "_")) %>% 
  filter(!id == "SYN22-4_B_3") %>% 
  select(-id)
```

    ```

- 

<!-- -->

    ```{r,corrected ap plot}

``` r
ap.plot <- spectra_subtract %>%
  ggplot(aes(
    x = wl,
    y = ap,
    color = as.numeric(tp),
    group = interaction(exp, bottle, tp)
  )) +
  facet_grid(bottle ~ factor(exp2, levels = unique(spectra_subtract$exp2)),
             labeller = label_parsed,
             scales = "free_y") +
  geom_errorbar(aes(ymin = ap - a_se, ymax = ap + a_se)) +
  geom_line(linewidth = 1) +
  labs(x = expression(bold(Wavelength ~ (nm ^ -1))),
       y = expression(bold(a[p] ~ (m ^ -1))),
       color = "Time point") +
  scale_color_gradient2(
    low = scales::muted("red"),
    mid = "white",
    high = scales::muted("blue"),
    midpoint = 8,
    na.value = "grey50",
    guide = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black"
    )
  ) +
  theme_linedraw(21) +
  custom.theme +
  theme(legend.key.width = unit(2.5, "cm"))
```

    ```

- 

<!-- -->

    ```{r,corrected cp plot}

``` r
cp.plot <- spectra_subtract %>%
  ggplot(aes(
    x = wl,
    y = cp,
    color = as.numeric(tp),
    group = interaction(exp, bottle, tp)
  )) +
  facet_grid(bottle ~ factor(exp2, levels = unique(spectra_subtract$exp2)),
             labeller = label_parsed,
             scales = "free_y") +
  geom_errorbar(aes(ymin = cp - c_se, ymax = cp + c_se)) +
  geom_line(linewidth = 1) +
  labs(x = expression(bold(Wavelength ~ (nm ^ -1))),
       y = expression(bold(c[p] ~ (m ^ -1))),
       color = "Time point") +
  scale_color_gradient2(
    low = scales::muted("red"),
    mid = "white",
    high = scales::muted("blue"),
    midpoint = 8,
    na.value = "grey50",
    guide = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black"
    )
  ) +
  theme_linedraw(21) +
  custom.theme +
  theme(legend.key.width = unit(2.5, "cm"))
```

    ```

- 

<!-- -->

    ```{r,pariculate coefficient plots, fig.asp=0.8, fig.width=24, message=FALSE, warning=FALSE}

``` r
(ap.plot + theme(axis.title.x = element_blank())) /  (cp.plot + guides(color = "none")) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A')  &
  theme(legend.position = "bottom",
        plot.tag = element_text(face = "bold")) 
```

![](ACS_files/figure-gfm/pariculate%20coefficient%20plots-1.png)<!-- -->

    ```

- 

########################### 

# DERIVE CHL FROM ALH

########################### 

- 

calculating chlorophyll concentrations based on Roesler and Barnard
(2015), where chl (mg/m^3) = (alh (m^-1) + 0.012 (m^1)) / 0.0104
(m^2/mg)

    ```{r,calculate chl, message=FALSE}

``` r
acs_alh <- spectra_subtract %>%
  group_by(exp, bottle, tp) %>%
  group_modify(~ add_row(., wl = c(650, 676, 715))) %>%
  arrange(exp, bottle, tp, wl) %>%
  mutate(ap_interp = zoo::na.approx(ap, wl, na.rm = T)) %>%
  fill(c(4:6)) %>%
  filter(wl %in% c(650, 676, 715)) %>%
  select(-c(8:12)) %>%
  rename(ap = ap_interp) %>%
  ungroup() %>%
  group_by(exp, bottle, tp) %>%
  mutate(
    blank = ((ap[wl == 715] - ap[wl == 650]) / (715 - 650)) * (676 - 650) + ap[wl == 650],
    abs_line_height = ap[wl == 676] - blank,
    chl_line_height = (abs_line_height + 0.012) / 0.014
  ) %>%
  ungroup() %>%
  select(exp, exp2, bottle, tp, abs_line_height, chl_line_height) %>%
  distinct() %>%
  left_join(spectra_subtract, .)
```

    ```

    ```{r,interpolate values for bb derivations, message=FALSE, warning=FALSE}

``` r
acs_bottle <- spectra_subtract %>%
  group_by(exp, bottle, tp) %>%
  group_modify(~ add_row(., wl = c(470, 532, 660))) %>% 
  arrange(exp, bottle, tp, wl) %>%
  mutate(ap_interp = zoo::na.approx(ap, wl, na.rm = T),
         cp_interp = zoo::na.approx(cp, wl, na.rm = T),
         bp_interp = zoo::na.approx(bp, wl, na.rm = T)) %>% 
  fill(c(4:6)) %>%
  ungroup() %>% 
  left_join(., acs_alh) %>% 
  group_by(exp, bottle) %>%
  fill(c(16:17)) %>%
  ungroup() %>% 
  mutate(cp = cp_interp,
         ap = ap_interp,
         bp = bp_interp) %>% 
  select(-c(ap_interp:bp_interp)) %>% 
  mutate(datetime = floor_date(datetime, "hour")) %>% 
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
   mutate(
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
    .after = date
  ) %>% 
  mutate(source = "Optics", .after = phyto) %>% 
  mutate_at(vars(exp_no), as.character) %>% 
  arrange(datetime) %>% 
  select(-exp2) 
```

    ```

- 

########################### 

# SUMMARIZE

########################### 

    ```{r,summarize fcm data, message=FALSE}

``` r
acs_sum <- acs_bottle %>%
  select(exp, tp, wl, cp, ap, bp, abs_line_height, chl_line_height) %>%
  mutate_at(vars(tp), as.numeric) %>%
  group_by(exp, tp, wl) %>%
  summarize(across(
    .cols = cp:chl_line_height,
    .fns = list(mean = mean, sd = sd),
    na.rm = TRUE,
    .names = "{fn}_{col}"
  )) %>%
  ungroup() %>% 
  left_join(acs_bottle %>% select(1:4, 6:11) %>% distinct() %>% mutate_at(vars(tp), as.numeric), .)
```

    ```

########################### 

# SAVE DATA

########################### 

    ```{r,save data}

``` r
write_csv(acs_sum, "processed_acs_summary.csv")
write_csv(acs_bottle, "processed_acs_bottle.csv")
```

    ```

- 

END
