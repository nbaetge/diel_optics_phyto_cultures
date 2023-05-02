FCM
================
Nicholas Baetge
Last compiled on 21 April, 2023

- 

<!-- -->

    ```{r,load packages, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
library(flowCore)
library(openCyto)
library(ggcyto)
library(tidyverse)
library(hms)
library(lubridate)
library(purrr)
library(patchwork)
library(janitor)
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

### import fcs files and metadata

    ```{r,read in fcs files, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
tp1_fs <-
  flowCore::read.flowSet(path = "TP22-13_FCS/",
                         pattern = ".fcs",
                         emptyValue = F)
tp2_fs <-
  flowCore::read.flowSet(path = "TP22-14_FCS/",
                         pattern = ".fcs",
                         emptyValue = F)
syn1_fs <-
  flowCore::read.flowSet(path = "SYN22-4_FCS/",
                         pattern = ".fcs",
                         emptyValue = F)
syn2_fs <-
  flowCore::read.flowSet(path = "SYN22-5_FCS/",
                         pattern = ".fcs",
                         emptyValue = F)
ol1_fs <-
  flowCore::read.flowSet(path = "OL22-2_FCS/",
                         pattern = ".fcs",
                         emptyValue = F)
ol2_fs <-
  flowCore::read.flowSet(path = "OL22-3_FCS/",
                         pattern = ".fcs",
                         emptyValue = F)
```

    ```

- 

<!-- -->

    ```{r,extract acquisition times, message=FALSE, wrapper = TRUE}

``` r
tp1_times <-
  flowCore::flowSet_to_list(tp1_fs) %>% map(., pluck, "description", "$BTIM") %>% map(as_tibble) %>% bind_rows(., .id = "column_label") %>% rename(rowname = "column_label", acq_time = value) %>%
  left_join(
    .,
    flowCore::flowSet_to_list(tp1_fs) %>% map(., pluck, "description", "$ETIM") %>% map(as_tibble) %>% bind_rows(., .id = "column_label") %>% rename(rowname = "column_label", end_time = value)
  )

tp2_times <-
  flowCore::flowSet_to_list(tp2_fs) %>% map(., pluck, "description", "$BTIM") %>% map(as_tibble) %>% bind_rows(., .id = "column_label") %>% rename(rowname = "column_label", acq_time = value) %>%
  left_join(
    .,
    flowCore::flowSet_to_list(tp2_fs) %>% map(., pluck, "description", "$ETIM") %>% map(as_tibble) %>% bind_rows(., .id = "column_label") %>% rename(rowname = "column_label", end_time = value)
  )

syn1_times <-
  flowCore::flowSet_to_list(syn1_fs) %>% map(., pluck, "description", "$BTIM") %>% map(as_tibble) %>% bind_rows(., .id = "column_label") %>% rename(rowname = "column_label", acq_time = value) %>%
  left_join(
    .,
    flowCore::flowSet_to_list(syn1_fs) %>% map(., pluck, "description", "$ETIM") %>% map(as_tibble) %>% bind_rows(., .id = "column_label") %>% rename(rowname = "column_label", end_time = value)
  )

syn2_times <-
  flowCore::flowSet_to_list(syn2_fs) %>% map(., pluck, "description", "$BTIM") %>% map(as_tibble) %>% bind_rows(., .id = "column_label") %>% rename(rowname = "column_label", acq_time = value) %>%
  left_join(
    .,
    flowCore::flowSet_to_list(syn2_fs) %>% map(., pluck, "description", "$ETIM") %>% map(as_tibble) %>% bind_rows(., .id = "column_label") %>% rename(rowname = "column_label", end_time = value)
  )

ol1_times <-
  flowCore::flowSet_to_list(ol1_fs) %>% map(., pluck, "description", "$BTIM") %>% map(as_tibble) %>% bind_rows(., .id = "column_label") %>% rename(rowname = "column_label", acq_time = value) %>%
  left_join(
    .,
    flowCore::flowSet_to_list(ol1_fs) %>% map(., pluck, "description", "$ETIM") %>% map(as_tibble) %>% bind_rows(., .id = "column_label") %>% rename(rowname = "column_label", end_time = value)
  )

ol2_times <-
  flowCore::flowSet_to_list(ol2_fs) %>% map(., pluck, "description", "$BTIM") %>% map(as_tibble) %>% bind_rows(., .id = "column_label") %>% rename(rowname = "column_label", acq_time = value) %>%
  left_join(
    .,
    flowCore::flowSet_to_list(ol2_fs) %>% map(., pluck, "description", "$ETIM") %>% map(as_tibble) %>% bind_rows(., .id = "column_label") %>% rename(rowname = "column_label", end_time = value)
  )
```

    ```

- 

<!-- -->

    ```{r,create metadata files, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
tp1_meta <-
  read_csv("TP22-13_IS.csv",
           skip = 7,
           locale = locale(encoding = "latin1")) %>%
  select(1, 2, 11, 12, 6, 9, 53, 58) %>%
  clean_names() %>%
  mutate(
    work_list = gsub(
      "C:/Documents and Settings/GTI/Desktop/Nick Baetge/2022/5/TP2213/",
      "",
      work_list
    ),
    work_list = gsub(".csv", "", work_list)
  ) %>%
  rename(
    well = 1,
    sample = 2,
    acq_date = 3,
    acq_time = 4,
    dil = 5,
    vol = 6,
    flow_rate = 7,
    plate = 8
  ) %>%
  mutate_if(is.character, stringr::str_trim) %>%
  separate(sample,
           into = c("source", "type", "tp"),
           sep = "-") %>%
  group_by(source, type, tp) %>%
  mutate(
    sample_time = round_hms(first(as_hms(acq_time)), digits = -2),
    .after = acq_time,
    acq_date = gsub("MAY", 5, acq_date),
    date = dmy_hms(paste(acq_date, sample_time))
  ) %>%
  ungroup() %>%
  select(plate, everything()) %>%
  mutate(exp = "TP22-13", .before = "plate")


tp2_meta <-
  read_csv("TP22-14_IS.csv",
           skip = 7,
           locale = locale(encoding = "latin1")) %>%
  select(1, 2, 11, 12, 6, 9, 53, 58) %>%
  clean_names() %>%
  mutate(
    work_list = gsub("C:/Documents and Settings/GTI/Desktop/TP2214/",
                     "",
                     work_list),
    work_list = gsub(".csv", "", work_list)
  ) %>%
  rename(
    well = 1,
    sample = 2,
    acq_date = 3,
    acq_time = 4,
    dil = 5,
    vol = 6,
    flow_rate = 7,
    plate = 8
  ) %>%
  mutate_if(is.character, stringr::str_trim) %>%
  separate(sample,
           into = c("source", "type", "tp"),
           sep = "-") %>%
  group_by(source, type, tp) %>%
  mutate(
    sample_time = round_hms(first(as_hms(acq_time)), digits = -2),
    .after = acq_time,
    acq_date = gsub("MAY", 5, acq_date),
    date = dmy_hms(paste(acq_date, sample_time))
  ) %>%
  ungroup() %>%
  select(plate, everything()) %>%
  mutate(exp = "TP22-14", .before = "plate")


syn1_meta <-
  read_csv("SYN22-4_IS.csv",
           skip = 7,
           locale = locale(encoding = "latin1")) %>%
  select(1, 2, 11, 12, 6, 9, 53, 58) %>%
  clean_names() %>%
  mutate(
    work_list = gsub(
      "C:/Documents and Settings/GTI/Desktop/Nick Baetge/Syn/2022/8/SYN22-4/",
      "",
      work_list
    ),
    work_list = gsub(".csv", "", work_list)
  ) %>%
  rename(
    well = 1,
    sample = 2,
    acq_date = 3,
    acq_time = 4,
    dil = 5,
    vol = 6,
    flow_rate = 7,
    plate = 8
  ) %>%
  mutate_if(is.character, stringr::str_trim) %>%
  separate(sample,
           into = c("source", "type", "tp"),
           sep = "-") %>%
  group_by(source, type, tp) %>%
  mutate(
    sample_time = round_hms(first(as_hms(acq_time)), digits = -2),
    .after = acq_time,
    acq_date = gsub("AUG", 8, acq_date),
    date = dmy_hms(paste(acq_date, sample_time))
  ) %>%
  ungroup() %>%
  select(plate, everything()) %>%
  mutate(exp = "SYN22-4", .before = "plate")

syn2_meta <-
  read_csv("SYN22-5_IS.csv",
           skip = 7,
           locale = locale(encoding = "latin1")) %>%
  select(1, 2, 11, 12, 6, 9, 53, 58) %>%
  clean_names() %>%
  mutate(
    work_list = gsub(
      "C:/Documents and Settings/GTI/Desktop/Nick Baetge/Syn/2022/8/SYN22-5/plate_maps/",
      "",
      work_list
    ),
    work_list = gsub(".csv", "", work_list)
  ) %>%
  rename(
    well = 1,
    sample = 2,
    acq_date = 3,
    acq_time = 4,
    dil = 5,
    vol = 6,
    flow_rate = 7,
    plate = 8
  ) %>%
  mutate_if(is.character, stringr::str_trim) %>%
  separate(sample,
           into = c("source", "type", "tp"),
           sep = "-") %>%
  group_by(source, type, tp) %>%
  mutate(
    sample_time = round_hms(first(as_hms(acq_time)), digits = -2),
    .after = acq_time,
    acq_date = gsub("AUG", 8, acq_date),
    date = dmy_hms(paste(acq_date, sample_time))
  ) %>%
  ungroup() %>%
  select(plate, everything()) %>%
  mutate(exp = "SYN22-5", .before = "plate")

ol1_meta <-
  read_csv(
    "OL22-2_IS.csv",
    skip = 7,
    locale = readr::locale(encoding = "latin1")
  ) %>%
  select(1,2,11,12,6,9,53, 58) %>%
  clean_names() %>% 
  mutate(work_list = gsub("C:/Documents and Settings/GTI/Desktop/Nick Baetge/OL/OL22-2/", "", work_list),
         work_list = gsub(".csv", "", work_list)) %>% 
  rename(
    well = 1,
    sample = 2,
    acq_date = 3,
    acq_time = 4,
    dil = 5,
    vol = 6,
    flow_rate = 7,
    plate = 8
  ) %>%
   mutate_if(is.character, stringr::str_trim) %>% 
  separate(sample, into = c("source", "type", "tp"), sep = "-") %>%
  mutate(dil = 1) %>% 
  group_by(source, type, tp) %>%
  mutate(sample_time = round_hms(first(as_hms(acq_time)), digits = -2), .after = acq_time,
         acq_date = gsub("NOV", 11, acq_date),
         acq_date = gsub("DEC", 12, acq_date),
         date = dmy_hms(paste(acq_date, sample_time))) %>%
  ungroup() %>% 
  select(plate, everything()) %>%
  mutate(exp = "OL22-2", .before = "plate")

ol2_meta <-
  read_csv(
    "OL22-3_IS.csv",
    skip = 7,
    # locale = readr::locale(encoding = "latin1")
  ) %>%
  select(1,2,11,12,6,9,53, 58) %>%
  clean_names() %>% 
  mutate(work_list = gsub("C:/Documents and Settings/GTI/Desktop/Nick Baetge/OL/OL22-3/", "", work_list)) %>% 
  mutate(work_list = gsub(".csv", "", work_list)) %>% 
  rename(
    well = 1,
    sample = 2,
    acq_date = 3,
    acq_time = 4,
    dil = 5,
    vol = 6,
    flow_rate = 7,
    plate = 8
  ) %>%
   mutate_if(is.character, stringr::str_trim) %>% 
  separate(sample, into = c("source", "type", "tp"), sep = "-") %>%
  mutate(dil = 1) %>% 
  group_by(source, type, tp) %>%
  mutate(sample_time = round_hms(first(as_hms(acq_time)), digits = -2), .after = acq_time,
         acq_date = gsub("DEC", 12, acq_date),
         date = dmy_hms(paste(acq_date, sample_time))) %>%
  ungroup() %>% 
  select(plate, everything()) %>%
  mutate(exp = "OL22-3", .before = "plate") %>% 
  mutate_at(vars(acq_time, sample_time), as_hms) 

meta <- bind_rows(tp1_meta, tp2_meta, syn1_meta, syn2_meta, ol1_meta, ol2_meta)
```

    ```

- 

<!-- -->

    ```{r,calculate acquisition volume, message=FALSE, wrapper = TRUE}

``` r
tp1_fs_data <- pData(tp1_fs) %>%
  mutate(name = gsub(".fcs", "", name)) %>%
  separate(name, into = c("t1", "name"), sep = "_") %>%
  separate(name, into = c("plate", "well"), sep = "-") %>%
  select(-c(t1)) %>%
  mutate_at(2, as.numeric) %>%
  rownames_to_column() %>%
  left_join(., tp1_times) %>%
  left_join(., tp1_meta %>% select(-acq_time)) %>%
  arrange(plate, well) %>%
  group_by(source, type, tp) %>%
  mutate(rep = seq(1:n()), .after = tp) %>%
  ungroup() %>% 
  column_to_rownames(var = "rowname") %>%
  relocate(c(acq_time, end_time), .after = sample_time) %>%
  mutate_at(vars(acq_time, end_time), as_hms) %>%
  mutate(record_time = as.numeric(end_time - acq_time),
         calc_vol = flow_rate * record_time) %>%
  rename(ins_vol = vol)


tp2_fs_data <- pData(tp2_fs) %>%
  mutate(name = gsub(".fcs", "", name)) %>%
  separate(name, into = c("t1", "name"), sep = "_") %>%
  separate(name, into = c("plate", "well"), sep = "-") %>%
  select(-c(t1)) %>%
  mutate_at(2, as.numeric) %>%
  rownames_to_column() %>%
  left_join(., tp2_times) %>%
  left_join(., tp2_meta %>% select(-acq_time)) %>%
  arrange(plate, well) %>%
  group_by(source, type, tp) %>%
  mutate(rep = seq(1:n()), .after = tp) %>%
  ungroup() %>% 
  column_to_rownames(var = "rowname") %>%
  relocate(c(acq_time, end_time), .after = sample_time) %>%
  mutate_at(vars(acq_time, end_time), as_hms) %>%
  mutate(record_time = as.numeric(end_time - acq_time),
         calc_vol = flow_rate * record_time) %>%
  rename(ins_vol = vol)

syn1_fs_data <- pData(syn1_fs) %>%
  mutate(name = gsub(".fcs", "", name)) %>%
  separate(name, into = c("t1", "name"), sep = "_") %>%
  separate(name, into = c("plate", "well"), sep = "-") %>%
  select(-c(t1)) %>%
  mutate_at(2, as.numeric) %>%
  rownames_to_column() %>%
  left_join(., syn1_times) %>%
  left_join(., syn1_meta %>% select(-acq_time)) %>%
  arrange(plate, well) %>%
  group_by(source, type, tp) %>%
  mutate(rep = seq(1:n()), .after = tp) %>%
  ungroup() %>% 
  column_to_rownames(var = "rowname") %>%
  relocate(c(acq_time, end_time), .after = sample_time) %>%
  mutate_at(vars(acq_time, end_time), as_hms) %>%
  mutate(record_time = as.numeric(end_time - acq_time),
         calc_vol = flow_rate * record_time) %>%
  rename(ins_vol = vol)

syn2_fs_data <- pData(syn2_fs) %>%
  mutate(name = gsub(".fcs", "", name)) %>%
  separate(name, into = c("t1", "name"), sep = "_") %>%
  separate(name, into = c("plate", "well"), sep = "-") %>%
  select(-c(t1)) %>%
  mutate_at(2, as.numeric) %>%
  rownames_to_column() %>%
  left_join(., syn2_times) %>%
  left_join(., syn2_meta %>% select(-acq_time)) %>%
  arrange(plate, well) %>%
  group_by(source, type, tp) %>%
  mutate(rep = seq(1:n()), .after = tp) %>%
  ungroup() %>% 
  column_to_rownames(var = "rowname") %>%
  relocate(c(acq_time, end_time), .after = sample_time) %>%
  mutate_at(vars(acq_time, end_time), as_hms) %>%
  mutate(record_time = as.numeric(end_time - acq_time),
         calc_vol = flow_rate * record_time) %>%
  rename(ins_vol = vol)

ol1_fs_data <- pData(ol1_fs) %>%
  mutate(name = gsub(".fcs", "", name)) %>%
  separate(name, into = c("t1", "name"), sep = "_") %>%
  separate(name, into = c("plate", "well"), sep = "-") %>%
  select(-c(t1)) %>%
  mutate_at(2, as.numeric) %>%
  rownames_to_column() %>%
  left_join(., ol1_times) %>%
  left_join(., ol1_meta %>% select(-acq_time)) %>%
  arrange(plate, well) %>%
  group_by(source, type, tp) %>%
  mutate(rep = seq(1:n()), .after = tp) %>%
  ungroup() %>% 
  column_to_rownames(var = "rowname") %>%
  relocate(c(acq_time, end_time), .after = sample_time) %>%
  mutate_at(vars(acq_time, end_time), as_hms) %>%
  mutate(record_time = as.numeric(end_time - acq_time),
         calc_vol = flow_rate * record_time) %>%
  rename(ins_vol = vol)

ol2_fs_data <- pData(ol2_fs) %>%
  mutate(name = gsub(".fcs", "", name)) %>%
  separate(name, into = c("t1", "name"), sep = "_") %>%
  separate(name, into = c("plate", "well"), sep = "-") %>%
  select(-c(t1)) %>%
  mutate_at(2, as.numeric) %>%
  rownames_to_column() %>%
  left_join(., ol2_times) %>%
  left_join(., ol2_meta %>% select(-acq_time)) %>%
  arrange(plate, well) %>%
  group_by(source, type, tp) %>%
  mutate(rep = seq(1:n()), .after = tp) %>%
  ungroup() %>% 
  column_to_rownames(var = "rowname") %>%
  relocate(c(acq_time, end_time), .after = sample_time) %>%
  mutate_at(vars(acq_time, end_time), as_hms) %>%
  mutate(record_time = as.numeric(end_time - acq_time),
         calc_vol = flow_rate * record_time) %>%
  rename(ins_vol = vol)
```

    ```

- 

<!-- -->

    ```{r,apply wrangled metadata to fs, wrapper = TRUE}

``` r
pData(tp1_fs) <- tp1_fs_data
pData(tp2_fs) <- tp2_fs_data
pData(syn1_fs) <- syn1_fs_data
pData(syn2_fs) <- syn2_fs_data
pData(ol1_fs) <- ol1_fs_data
pData(ol2_fs) <- ol2_fs_data
```

    ```

- 

############# 

# GATE CYTOGRAMS

############## 

- 

<!-- -->

    ```{r,list channels, wrapper = TRUE}

``` r
colnames(tp1_fs)
```

    ##  [1] "FSC-HLin" "SSC-HLin" "GRN-HLog" "YEL-HLog" "RED-HLog" "FSC-HLog"
    ##  [7] "SSC-HLog" "GRN-HLin" "YEL-HLin" "RED-HLin" "TIME"     "GRN-A"   
    ## [13] "GRN-ALog" "GRN-W"

    ```

### gate

#### thalassiosira

    ```{r,tp1 root, wrapper = TRUE}

``` r
tp1_gs <- flowWorkspace::GatingSet(tp1_fs)
tp2_gs <- flowWorkspace::GatingSet(tp2_fs)

tp1_root <- gh_pop_get_data(tp1_gs[[11]], "root", returnType = "flowFrame")

autoplot(tp1_root, "FSC-HLog") ; autoplot(tp1_root, "RED-HLog") ; autoplot(tp1_root, "FSC-HLog", "RED-HLog", bins = 300) ; autoplot(tp1_root, "FSC-HLog", "SSC-HLog", bins = 200) 
```

![](Culture_and_Optics_FCM_files/figure-gfm/tp1%20root-1.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/tp1%20root-2.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/tp1%20root-3.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/tp1%20root-4.png)<!-- -->

    ```

- 

<!-- -->

    ```{r,tp gates, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
gt_tp <- gatingTemplate("TP_Gating.csv")
gt_gating(gt_tp, tp1_gs)
gt_gating(gt_tp, tp2_gs)
```

    ```

- 

<!-- -->

    ```{r,tp1 cells, wrapper = TRUE}

``` r
tp1_cells <- gh_pop_get_data(tp1_gs[[13]], "cells", returnType = "flowFrame")

autoplot(tp1_cells, "FSC-HLog") ; autoplot(tp1_cells, "RED-HLog") ; autoplot(tp1_cells, "FSC-HLog", "RED-HLog", bins = 200) ; autoplot(tp1_cells, "FSC-HLog", "SSC-HLog", bins = 200) 
```

![](Culture_and_Optics_FCM_files/figure-gfm/tp1%20cells-1.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/tp1%20cells-2.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/tp1%20cells-3.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/tp1%20cells-4.png)<!-- -->

    ```

- 

#### synechococcus

    ```{r,syn1 root, wrapper = TRUE}

``` r
syn1_gs <- flowWorkspace::GatingSet(syn1_fs)
syn2_gs <- flowWorkspace::GatingSet(syn2_fs)

syn1_root <- gh_pop_get_data(syn1_gs[[20]], "root", returnType = "flowFrame")

autoplot(syn1_root, "FSC-HLog") ; autoplot(syn1_root, "YEL-HLog") ; autoplot(syn1_root, "FSC-HLog", "YEL-HLog", bins = 300) ; autoplot(syn1_root, "FSC-HLog", "SSC-HLog", bins = 200) 
```

![](Culture_and_Optics_FCM_files/figure-gfm/syn1%20root-1.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/syn1%20root-2.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/syn1%20root-3.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/syn1%20root-4.png)<!-- -->

    ```

- 

<!-- -->

    ```{r,syn gates, message=FALSE, warning=FALSE}

``` r
gt_syn <- gatingTemplate("SYN_Gating.csv")
gt_gating(gt_syn, syn1_gs)
gt_gating(gt_syn, syn2_gs)
```

    ```

- 

<!-- -->

    ```{r,syn1 cells, wrapper = TRUE}

``` r
syn1_cells <- gh_pop_get_data(syn1_gs[[20]], "cells", returnType = "flowFrame")

autoplot(syn1_cells, "FSC-HLog") ; autoplot(syn1_cells, "YEL-HLog") ; autoplot(syn1_cells, "RED-HLog") ; autoplot(syn1_cells, "FSC-HLog", "YEL-HLog", bins = 200) ; autoplot(syn1_cells, "FSC-HLog", "SSC-HLog", bins = 200) 
```

![](Culture_and_Optics_FCM_files/figure-gfm/syn1%20cells-1.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/syn1%20cells-2.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/syn1%20cells-3.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/syn1%20cells-4.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/syn1%20cells-5.png)<!-- -->

    ```

- 

<!-- -->

    ```{r,syn1 low yel, wrapper = TRUE}

``` r
syn1_l <- gh_pop_get_data(syn1_gs[[20]], "YEL-HLog-", returnType = "flowFrame")

autoplot(syn1_l, "FSC-HLog") ; autoplot(syn1_l, "YEL-HLog") ; autoplot(syn1_l, "FSC-HLog", "YEL-HLog", bins = 200) ; autoplot(syn1_l, "FSC-HLog", "SSC-HLog", bins = 200) 
```

![](Culture_and_Optics_FCM_files/figure-gfm/syn1%20low%20yel-1.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/syn1%20low%20yel-2.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/syn1%20low%20yel-3.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/syn1%20low%20yel-4.png)<!-- -->

    ```

- 

<!-- -->

    ```{r,syn1 high yel, wrapper = TRUE}

``` r
syn1_h <- gh_pop_get_data(syn1_gs[[20]], "YEL-HLog+", returnType = "flowFrame")

autoplot(syn1_h, "FSC-HLog") ; autoplot(syn1_h, "YEL-HLog") ; autoplot(syn1_h, "FSC-HLog", "YEL-HLog", bins = 200) ; autoplot(syn1_h, "FSC-HLog", "SSC-HLog", bins = 200) 
```

![](Culture_and_Optics_FCM_files/figure-gfm/syn1%20high%20yel-1.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/syn1%20high%20yel-2.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/syn1%20high%20yel-3.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/syn1%20high%20yel-4.png)<!-- -->

    ```

- 

#### ostreococcus

    ```{r,ol1 root, wrapper = TRUE}

``` r
ol1_gs <- flowWorkspace::GatingSet(ol1_fs)
ol2_gs <- flowWorkspace::GatingSet(ol2_fs)

ol1_root <- gh_pop_get_data(ol1_gs[[127]], "root", returnType = "flowFrame")

autoplot(ol1_root, "FSC-HLog") ; autoplot(ol1_root, "RED-HLog") ;  autoplot(ol1_root, "GRN-HLog") ;  autoplot(ol1_root, "YEL-HLog") ; autoplot(ol1_root, "GRN-HLog", "YEL-HLog", bins = 300) ; autoplot(ol1_root, "FSC-HLog", "SSC-HLog", bins = 200)
```

![](Culture_and_Optics_FCM_files/figure-gfm/ol1%20root-1.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/ol1%20root-2.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/ol1%20root-3.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/ol1%20root-4.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/ol1%20root-5.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/ol1%20root-6.png)<!-- -->

    ```

``` r
g <- openCyto:::gate_flowclust_2d(ol1_root, xChannel = "GRN-HLog", yChannel = "YEL-HLog", K=1, quantile=0.9)
```

    ## The prior specification has no effect when usePrior=no

    ## Using the serial version of flowClust

``` r
p <- autoplot(ol1_root, "GRN-HLog", "YEL-HLog", bins = 200)
p+geom_gate(g)
```

![](Culture_and_Optics_FCM_files/figure-gfm/ol%20gate%20plot-1.png)<!-- -->

- 

<!-- -->

    ```{r,ol gates, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
gt_ol <- gatingTemplate("OL_Gating.csv")
gt_gating(gt_ol, ol1_gs)
gt_gating(gt_ol, ol2_gs)
```

    ```

- 

<!-- -->

    ```{r,ol1 nondebris, wrapper = TRUE}

``` r
# ol1_nd <- gh_pop_get_data(ol1_gs[[127]], "nondebris", returnType = "flowFrame")
# 
# autoplot(ol1_nd, "FSC-HLog") ; autoplot(ol1_nd, "SSC-HLog") ;  autoplot(ol1_nd, "RED-HLog") ; autoplot(ol1_nd, "GRN-HLog") ;  autoplot(ol1_nd, "YEL-HLog") ; autoplot(ol1_nd, "GRN-HLog", "YEL-HLog", bins = 300) ; autoplot(ol1_nd, "FSC-HLog", "SSC-HLog", bins = 200)
```

    ```

``` r
# g <- openCyto:::.boundary(ol1_nb, channels = c("FSC-HLog", "SSC-HLog"), min = c(0.25, 0.25), max = c(2.5, 2.5))
# p <- autoplot(ol1_nb, "FSC-HLog", "SSC-HLog", bins = 200)
# p+geom_gate(g)
```

    ```{r,ol1 cells, wrapper = TRUE}

``` r
ol1_cells <- gh_pop_get_data(ol1_gs[[12]], "cells", returnType = "flowFrame")

autoplot(ol1_cells, "FSC-HLog") ; autoplot(ol1_cells, "SSC-HLog") ;  autoplot(ol1_cells, "RED-HLog") ; autoplot(ol1_cells, "GRN-HLog") ;  autoplot(ol1_cells, "YEL-HLog") ; autoplot(ol1_cells, "GRN-HLog", "YEL-HLog", bins = 300) ; autoplot(ol1_cells, "FSC-HLog", "SSC-HLog", bins = 200)
```

![](Culture_and_Optics_FCM_files/figure-gfm/ol1%20cells-1.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/ol1%20cells-2.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/ol1%20cells-3.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/ol1%20cells-4.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/ol1%20cells-5.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/ol1%20cells-6.png)<!-- -->![](Culture_and_Optics_FCM_files/figure-gfm/ol1%20cells-7.png)<!-- -->

    ```

- 

############# 

# CALCULATE ABUNDANCES

############## 

- 

<!-- -->

    ```{r,extract gate data, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
tp1_pop <-
  CytoExploreR::cyto_extract(tp1_gs, parent = "/cells", raw = T) %>%
  map(as_tibble) %>%
  map(., ~ select(., contains("HLog"))) %>%
  bind_rows(., .id = "column_label") %>%
  rename(name = "column_label") %>%
  rename_all( ~ str_replace_all(., "-HLog", "")) %>% 
  group_by(name) %>%
  summarize(across(
    .cols = GRN:SSC,
    .fns = list(mean = mean),
    na.rm = TRUE,
    .names = "{col}"
  )) %>%
  ungroup()

tp2_pop <-
  CytoExploreR::cyto_extract(tp2_gs, parent = "/cells", raw = T) %>%
  map(as_tibble) %>%
  map(., ~ select(., contains("HLog"))) %>%
  bind_rows(., .id = "column_label") %>%
  rename(name = "column_label") %>%
  rename_all( ~ str_replace_all(., "-HLog", "")) %>% 
  group_by(name) %>%
  summarize(across(
    .cols = GRN:SSC,
    .fns = list(mean = mean),
    na.rm = TRUE,
    .names = "{col}"
  )) %>%
  ungroup()

syn1_pop <-
  CytoExploreR::cyto_extract(syn1_gs, parent = "/cells", raw = T) %>%
  map(as_tibble) %>%
  map(., ~ select(., contains("HLog"))) %>%
  bind_rows(., .id = "column_label") %>%
  rename(name = "column_label") %>%
  rename_all( ~ str_replace_all(., "-HLog", "")) %>% 
  group_by(name) %>%
  summarize(across(
    .cols = GRN:SSC,
    .fns = list(mean = mean),
    na.rm = TRUE,
    .names = "{col}"
  )) %>%
  ungroup()

syn2_pop <-
  CytoExploreR::cyto_extract(syn2_gs, parent = "/cells", raw = T) %>%
  map(as_tibble) %>%
  map(., ~ select(., contains("HLog"))) %>%
  bind_rows(., .id = "column_label") %>%
  rename(name = "column_label") %>%
  rename_all( ~ str_replace_all(., "-HLog", "")) %>% 
  group_by(name) %>%
  summarize(across(
    .cols = GRN:SSC,
    .fns = list(mean = mean),
    na.rm = TRUE,
    .names = "{col}"
  )) %>%
  ungroup()

ol1_pop <-
  CytoExploreR::cyto_extract(ol1_gs, parent = "nondebris/cells", raw = T) %>%
  map(as_tibble) %>%
  map(., ~ select(., contains("HLog"))) %>%
  bind_rows(., .id = "column_label") %>%
  rename(name = "column_label") %>%
  rename_all( ~ str_replace_all(., "-HLog", "")) %>% 
  group_by(name) %>%
  summarize(across(
    .cols = GRN:SSC,
    .fns = list(mean = mean),
    na.rm = TRUE,
    .names = "{col}"
  )) %>%
  ungroup()

ol2_pop <-
  CytoExploreR::cyto_extract(ol2_gs, parent = "nondebris/cells", raw = T) %>%
  map(as_tibble) %>%
  map(., ~ select(., contains("HLog"))) %>%
  bind_rows(., .id = "column_label") %>%
  rename(name = "column_label") %>%
  rename_all( ~ str_replace_all(., "-HLog", "")) %>% 
  group_by(name) %>%
  summarize(across(
    .cols = GRN:SSC,
    .fns = list(mean = mean),
    na.rm = TRUE,
    .names = "{col}"
  )) %>%
  ungroup()
```

    ```

- 

<!-- -->

    ```{r,merge data and calculate abundances for TP and OL, message=FALSE, wrapper = TRUE}

``` r
#should be able to use gs_pop_get_count_with_metadata, but doesn't work sometimes

tp1_abund <- flowWorkspace::gs_pop_get_count_fast(tp1_gs)  %>%
  left_join(., flowWorkspace::pData(tp1_gs)) %>%
  rename(bottle = type) %>%
  mutate(bottle = "A") %>%
  select(
    name,
    plate,
    well,
    source,
    bottle,
    tp,
    rep,
    date,
    acq_date,
    acq_time,
    end_time,
    sample_time,
    record_time,
    flow_rate,
    calc_vol,
    ParentCount,
    Count
  ) %>%
  mutate_at(vars(3, 6, 7, 14:17), as.numeric) %>%
  filter(!ParentCount < 100) %>%
  mutate_at(vars(sample_time, acq_time), as_hms) %>%
  mutate_at(vars(acq_date), dmy) %>%
  mutate_at(vars(date), ymd_hms) %>%
  arrange(source, bottle, tp, rep) %>%
  mutate(dilution = 1) %>% 
  mutate(cells = round((Count / calc_vol) * 1000) / (dilution)) %>%
  mutate(
    event_s = ParentCount / as.numeric(record_time),
    gate_event_s = Count / as.numeric(record_time),
    gate_percent = (Count / ParentCount) * 100
  ) %>%
  group_by(source, bottle, tp) %>%
  mutate(mean_cells = round(mean(cells))) %>%
  ungroup() %>%
  left_join(., tp1_pop)

tp2_abund <- flowWorkspace::gs_pop_get_count_fast(tp2_gs)  %>%
  left_join(., flowWorkspace::pData(tp2_gs)) %>%
  rename(bottle = type) %>%
  mutate(bottle = case_when(bottle == "X" ~ "A",
                            bottle == "Y" ~ "B",
                            bottle == "Z" ~ "C")) %>%
  select(
    name,
    plate,
    well,
    source,
    bottle,
    tp,
    rep,
    date,
    acq_date,
    acq_time,
    end_time,
    sample_time,
    record_time,
    flow_rate,
    calc_vol,
    ParentCount,
    Count
  ) %>%
  mutate_at(vars(3, 6, 7, 14:17), as.numeric) %>%
  filter(!ParentCount < 100) %>%
  mutate_at(vars(sample_time, acq_time), as_hms) %>%
  mutate_at(vars(acq_date), dmy) %>%
  mutate_at(vars(date), ymd_hms) %>%
  arrange(source, bottle, tp, rep) %>%
  mutate(dilution = 1) %>% 
  mutate(cells = round((Count / calc_vol) * 1000) * dilution) %>%
  mutate(
    event_s = ParentCount / as.numeric(record_time),
    gate_event_s = Count / as.numeric(record_time),
    gate_percent = (Count / ParentCount) * 100
  ) %>%
  group_by(source, bottle, tp) %>%
  mutate(mean_cells = round(mean(cells))) %>%
  mutate(
    a_cells = ifelse(source == "Optics" & bottle == "A", mean_cells, NA),
    b_cells = ifelse(source == "Optics" &
                       bottle == "B", mean_cells, NA),
    c_cells = ifelse(source == "Optics" &
                       bottle == "C", mean_cells, NA),
  ) %>%
  group_by(source, tp) %>%
  fill(a_cells:c_cells, .direction = "downup") %>%
  mutate(
    cells_added = case_when(
      source == "Optics" & bottle == "A" ~ a_cells,
      source == "Optics" &
        bottle == "B" ~ b_cells - a_cells,
      source == "Optics" &
        bottle == "C" ~ c_cells - b_cells,
    )
  ) %>%
  ungroup() %>%
  select(-c(a_cells:c_cells)) %>%
  left_join(., tp2_pop)


ol1_abund <- gs_pop_get_count_fast(ol1_gs)  %>%
  left_join(., pData(ol1_gs)) %>%
  filter(Population == "/nondebris/cells") %>% 
  rename(bottle = type) %>%
  select(
    name,
    plate,
    well,
    source,
    bottle,
    tp,
    rep,
    date,
    acq_date,
    acq_time,
    end_time,
    sample_time,
    record_time,
    flow_rate,
    calc_vol,
    ParentCount,
    Count
  ) %>%
  mutate_at(vars(3, 6, 7, 14:17), as.numeric) %>%
  filter(!ParentCount < 100) %>%
  mutate_at(vars(sample_time, acq_time), as_hms) %>%
  mutate_at(vars(acq_date), dmy) %>%
  mutate_at(vars(date), ymd_hms) %>%
  arrange(source, bottle, tp, rep) %>%
  mutate(
    dilution = case_when(
      tp == 1 & source == "Optics" ~ 0.2,
      tp > 1 & source == "Optics" ~ 0.1,
      source == "Culture" ~ 0.002
    )
  ) %>%
  mutate(cells = round((Count / calc_vol) * 1000) / (dilution)) %>%
  mutate(
    event_s = ParentCount / as.numeric(record_time),
    gate_event_s = Count / as.numeric(record_time),
    gate_percent = (Count / ParentCount) * 100
  ) %>%
  group_by(source, bottle, tp) %>%
  mutate(mean_cells = round(mean(cells))) %>%
  mutate(
    a_cells = ifelse(source == "Optics" & bottle == "A", mean_cells, NA),
    b_cells = ifelse(source == "Optics" &
                       bottle == "B", mean_cells, NA),
    c_cells = ifelse(source == "Optics" &
                       bottle == "C", mean_cells, NA),
  ) %>%
  group_by(source, tp) %>%
  fill(a_cells:c_cells, .direction = "downup") %>%
  mutate(
    cells_added = case_when(
      source == "Optics" & bottle == "A" ~ a_cells,
      source == "Optics" &
        bottle == "B" ~ b_cells - a_cells,
      source == "Optics" &
        bottle == "C" ~ c_cells - b_cells,
    )
  ) %>%
  ungroup() %>%
  select(-c(a_cells:c_cells)) %>%
  left_join(., ol1_pop) %>% 
  filter(!name %in% c("Wells_Plate2-69.fcs", "Wells_Plate2-70.fcs"))


ol2_abund <- gs_pop_get_count_fast(ol2_gs)  %>%
  left_join(., pData(ol2_gs)) %>%
  filter(Population == "/nondebris/cells") %>% 
  rename(bottle = type) %>%
  select(
    name,
    plate,
    well,
    source,
    bottle,
    tp,
    rep,
    date,
    acq_date,
    acq_time,
    end_time,
    sample_time,
    record_time,
    flow_rate,
    calc_vol,
    ParentCount,
    Count
  ) %>%
  mutate_at(vars(3, 6, 7, 14:17), as.numeric) %>%
  filter(!ParentCount < 100) %>%
  mutate_at(vars(sample_time, acq_time), as_hms) %>%
  mutate_at(vars(acq_date), dmy) %>%
  mutate_at(vars(date), ymd_hms) %>%
  arrange(source, bottle, tp, rep) %>%
  mutate(dilution = case_when(source == "Optics" ~ 0.1,
                              source == "Culture" ~ 0.002)) %>%
  mutate(cells = round((Count / calc_vol) * 1000) / (dilution)) %>%
  mutate(
    event_s = ParentCount / as.numeric(record_time),
    gate_event_s = Count / as.numeric(record_time),
    gate_percent = (Count / ParentCount) * 100
  ) %>%
  group_by(source, bottle, tp) %>%
  mutate(mean_cells = round(mean(cells))) %>%
  mutate(
    a_cells = ifelse(source == "Optics" & bottle == "A", mean_cells, NA),
    b_cells = ifelse(source == "Optics" &
                       bottle == "B", mean_cells, NA),
    c_cells = ifelse(source == "Optics" &
                       bottle == "C", mean_cells, NA),
  ) %>%
  group_by(source, tp) %>%
  fill(a_cells:c_cells, .direction = "downup") %>%
  mutate(
    cells_added = case_when(
      source == "Optics" & bottle == "A" ~ a_cells,
      source == "Optics" &
        bottle == "B" ~ b_cells - a_cells,
      source == "Optics" &
        bottle == "C" ~ c_cells - b_cells,
    )
  ) %>%
  ungroup() %>%
  select(-c(a_cells:c_cells)) %>%
  left_join(., ol2_pop) 
```

    ```

- 

Cell abundances for SYN22-5 are low quality due to instrumentation
saturation: - concentrations in the optics system will be corrected
using the relationship between cell abundance and chlorophyll (from ACS)
later - the average ratio of cells in the optics system relative to the
cultures from SYN22-4 will be used to estimate the SYN22-5 culture
concentrations

    ```{r,import acs chl line height data, message=FALSE, wrapper = TRUE}

``` r
# acs <- read_csv("~/Box Sync/Phyto_bbp/DATA/FINAL/2022_Experiments/acs/processed_acs_bottle.csv") %>% 
#   select(exp, bottle, tp, chl_line_height) %>% 
#   distinct() %>% 
#   drop_na() %>% 
#   arrange(exp, bottle, tp) %>% 
#   group_by(exp, tp) %>% 
#   mutate(chl = mean(chl_line_height)) %>% 
#   select(-c(bottle, chl_line_height)) %>% 
#   distinct() %>% 
#   ungroup()
```

    ```

- 

<!-- -->

    ```{r,merge data and calculate abundances for SYN22-4, message=FALSE, wrapper = TRUE}

``` r
syn1_abund <- flowWorkspace::gs_pop_get_count_fast(syn1_gs)  %>%
  left_join(., flowWorkspace::pData(syn1_gs)) %>%
  filter(Population == "/cells") %>% 
  rename(bottle = type) %>%
  mutate(bottle = case_when(bottle == "X" ~ "A",
                            bottle == "Y" ~ "B")) %>%
  select(
    name,
    plate,
    well,
    source,
    bottle,
    tp,
    rep,
    date,
    acq_date,
    acq_time,
    end_time,
    sample_time,
    record_time,
    flow_rate,
    calc_vol,
    ParentCount,
    Count
  ) %>%
  mutate_at(vars(3, 6, 7, 14:17), as.numeric) %>%
  filter(!ParentCount < 100) %>%
  mutate_at(vars(sample_time, acq_time), as_hms) %>%
  mutate_at(vars(acq_date), dmy) %>%
  mutate_at(vars(date), ymd_hms) %>%
  arrange(source, bottle, tp, rep) %>%
  mutate(dilution = ifelse(source == "Culture", 0.1, 1)) %>%
  mutate(cells = round((Count / calc_vol) * 1000) / (dilution)) %>%
  mutate(
    event_s = ParentCount / as.numeric(record_time),
    gate_event_s = Count / as.numeric(record_time),
    gate_percent = (Count / ParentCount) * 100
  ) %>%
  group_by(source, bottle, tp) %>%
  mutate(mean_cells = round(mean(cells))) %>%
  mutate(
    a_cells = ifelse(source == "Optics" & bottle == "A", mean_cells, NA),
    b_cells = ifelse(source == "Optics" &
                       bottle == "B", mean_cells, NA),
  ) %>%
  group_by(source, tp) %>%
  fill(a_cells:b_cells, .direction = "downup") %>%
  mutate(
    cells_added = case_when(
      source == "Optics" & bottle == "A" ~ a_cells,
      source == "Optics" &
        bottle == "B" ~ b_cells - a_cells
    )
  ) %>%
  ungroup() %>%
  select(-c(a_cells:b_cells)) %>%
  left_join(., syn1_pop) %>% 
  filter(!name == "Wells_Plate1-9.fcs")
```

    ```

- 

<!-- -->

    ```{r,lmodel2 of cells and chl for syn22-4 optics, message=FALSE}

``` r
# syn1_chl_cell <- syn1_abund %>% 
#   filter(source == "Optics") %>% 
#   select(tp, cells_added) %>% 
#   distinct() %>% 
#   group_by(tp) %>% 
#   mutate(cells_added = round(mean(cells_added))) %>%
#   distinct() %>% 
#   ungroup() %>% 
#   left_join(., acs %>% filter(exp == "SYN22-4") %>% select(-exp)) %>% 
#   rename(cells = 2, 
#          chl = 3) %>% 
#   summarize(lm = list(lmodel2::lmodel2(formula = cells ~ chl)))
```

    ```

- 

<!-- -->

    ```{r,calculate ratio of optics cells v culture cells for syn22-4, message=FALSE}

``` r
syn224_fraction <- syn1_abund  %>% 
  select(source, bottle, tp, mean_cells, cells_added) %>% 
  distinct() %>% 
  mutate(cells = case_when(source == "Optics" ~ cells_added,
                           source != "Optics" ~ mean_cells),
         cells = ifelse(tp > 13, NA, cells)) %>% 
  select(1,3, 6) %>% 
  distinct() %>% 
  arrange(tp) %>%
  group_by(source, tp) %>% 
  mutate(cells = mean(cells)) %>% 
  distinct() %>% 
  ungroup() %>% 
  group_by(tp) %>% 
  mutate(fraction = last(cells)/first(cells)) %>% 
  ungroup() %>% 
  select(-cells) %>% 
  summarize(across(
    .cols = fraction,
    .fns = list(mean = mean),
    na.rm = TRUE,
    .names = "{col}"
  )) %>%
  ungroup()
```

    ```

- 

<!-- -->

    ```{r,estimate syn22-5 optics cell counts using syn22-4 lmodel2, message=FALSE}

``` r
# sma_slope <- syn1_chl_cell[[1]][[1]]$regression.results$Slope[[3]]
# sma_intercept <- syn1_chl_cell[[1]][[1]]$regression.results$Intercept[[3]]
# 
# syn2_optics_cell_cor <- read_csv("~/Box Sync/Phyto_bbp/DATA/FINAL/2022_Experiments/acs/processed_acs_bottle.csv") %>% 
#   select(exp, bottle, tp, chl_line_height) %>% 
#   distinct() %>% 
#   drop_na() %>% 
#   arrange(exp, bottle, tp) %>% 
#   filter(exp ==  "SYN22-5") %>% 
#   rename(chl = 4) %>% 
#   mutate(corrected_cells =  (sma_slope * chl) - sma_intercept) %>% 
#   select(bottle, tp, corrected_cells) %>% 
#   mutate(source = "Optics")
```

    ```

- 

<!-- -->

    ```{r,merge data and calculate estimated abundances for SYN22-5, message=FALSE, wrapper = TRUE}

``` r
syn2_abund <- gs_pop_get_count_fast(syn2_gs)  %>%
  left_join(., pData(syn2_gs)) %>%
  filter(Population == "/cells") %>% 
  rename(bottle = type) %>%
  select(
    name,
    plate,
    well,
    source,
    bottle,
    tp,
    rep,
    date,
    acq_date,
    acq_time,
    end_time,
    sample_time,
    record_time,
    flow_rate,
    calc_vol,
    ParentCount,
    Count,
    dil
  ) %>%
  mutate_at(vars(3, 6, 7, 14:18), as.numeric) %>%
  filter(!ParentCount < 100) %>%
  filter(bottle != "D" | source != "Optics") %>%
  mutate_at(vars(sample_time, acq_time), as_hms) %>%
  mutate_at(vars(acq_date), dmy) %>%
  mutate_at(vars(date), ymd_hms) %>%
  arrange(source, bottle, tp, rep) %>%
  rename(dilution = dil) %>%
  mutate(cells = round((Count / calc_vol) * 1000) / (dilution)) %>%
  mutate(
    event_s = ParentCount / as.numeric(record_time),
    gate_event_s = Count / as.numeric(record_time),
    gate_percent = (Count / ParentCount) * 100
  ) %>%
  group_by(source, bottle, tp) %>%
  mutate(mean_cells = round(mean(cells))) %>%
  mutate(
    a_cells = ifelse(source == "Optics" & bottle == "A", mean_cells, NA),
    b_cells = ifelse(source == "Optics" &
                       bottle == "B", mean_cells, NA),
    c_cells = ifelse(source == "Optics" &
                       bottle == "C", mean_cells, NA),
  ) %>%
  group_by(source, tp) %>%
  fill(a_cells:c_cells, .direction = "downup") %>%
  mutate(
    cells_added = case_when(
      source == "Optics" & bottle == "A" ~ a_cells,
      source == "Optics" &
        bottle == "B" ~ b_cells - a_cells,
      source == "Optics" &
        bottle == "C" ~ c_cells - b_cells,
    )
  ) %>%
  ungroup() %>%
  select(-c(a_cells:c_cells)) %>%
  left_join(., syn2_pop) %>% 
  # left_join(., syn2_optics_cell_cor) %>%
  # mutate(across(c(cells, mean_cells:cells_added), ~ ifelse(source %in% c("Optics", "Culture"), NA, .))) %>%
  mutate(across(c(cells, mean_cells:cells_added), ~ ifelse(source %in% c("Culture"), NA, .))) %>%
  # mutate(culture_cells = corrected_cells / syn224_fraction$fraction[[1]]) %>% 
  mutate(culture_cells = cells_added / syn224_fraction$fraction[[1]]) %>% 
  arrange(bottle, tp) %>% 
  group_by(bottle, tp) %>% 
  # fill(corrected_cells:culture_cells, .direction = "updown") %>% 
  fill(cells_added:culture_cells, .direction = "updown") %>% 
  ungroup() %>% 
  # mutate(mean_cells = case_when(source == "Culture" ~ culture_cells),
  #        cells_added = case_when(source == "Optics" ~ corrected_cells)) %>% 
  mutate(mean_cells = case_when(source == "Culture" ~ culture_cells)) %>%
  select(1:29) %>% 
  arrange(source, tp)
```

    ```

- 

<!-- -->

    ```{r,combine bottle data, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
dates <- bind_rows(
    tp1_abund %>% mutate(exp = "TP22-13", .before = source),
    tp2_abund %>% mutate(exp = "TP22-14", .before = source),
    syn1_abund %>% mutate(exp = "SYN22-4", .before = source),
    syn2_abund %>% mutate(exp = "SYN22-5", .before = source),
    ol1_abund %>% mutate(exp = "OL22-2", .before = source),
    ol2_abund %>% mutate(exp = "OL22-3", .before = source),
  ) %>%
  select(exp, source, tp, date) %>%
  mutate(date = floor_date(date, "hour")) %>% 
  filter(source == "Culture") %>% 
  distinct() %>% 
  filter(!date %in% c(
    ymd_hms("2022-05-19 08:00:00"),
    ymd_hms("2022-08-25 00:00:00"),
    ymd_hms("2022-12-01 08:00:00"),
    ymd_hms("2022-12-07 12:00:00")
  )) %>% 
  select(-source)

bottle <-
  bind_rows(
    tp1_abund %>% mutate(exp = "TP22-13", .before = source),
    tp2_abund %>% mutate(exp = "TP22-14", .before = source),
    syn1_abund %>% mutate(exp = "SYN22-4", .before = source),
    syn2_abund %>% mutate(exp = "SYN22-5", .before = source),
    ol1_abund %>% mutate(exp = "OL22-2", .before = source),
    ol2_abund %>% mutate(exp = "OL22-3", .before = source),
  ) %>%
  select(exp, source, bottle, tp, date, mean_cells, cells_added, GRN:SSC) %>%
  mutate(date = floor_date(date, "hour")) %>%
  distinct() %>%
  mutate(mean_cells = ifelse(source == "Optics", cells_added, mean_cells)) %>%
  select(-7) %>%
  arrange(exp, tp) %>%
  rename(
    grn = 7,
    yel = 8,
    red = 9,
    fsc = 10,
    ssc = 11
  ) %>%
  group_by(exp, source, bottle, tp, date, mean_cells) %>%
  summarize(across(
    .cols = grn:ssc,
    .fns = list(mean = mean),
    na.rm = TRUE,
    .names = "{col}"
  )) %>%
  ungroup() %>%
  filter(!date %in% c(
    ymd_hms("2022-05-19 08:00:00"),
    ymd_hms("2022-08-25 00:00:00"),
    ymd_hms("2022-12-01 08:00:00"),
    ymd_hms("2022-12-07 12:00:00")
  )) %>%
  select(-date) %>% 
  left_join(., dates) %>% 
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
  rename(datetime = date,
         cells = mean_cells) %>%
  mutate(
    date = as_date(datetime),
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
    .after = datetime
  ) %>%
  mutate_at(vars(exp_no), as.character) %>% 
  arrange(source, datetime) %>% 
  select(exp:tp, datetime:plot_datetime, cells:ssc) %>% 
  mutate(id = paste(exp, source, bottle, tp, sep = "_"),
         cells = ifelse(id %in% c("TP22-14_Optics_B_11", "OL22-3_Optics_C_13"), NA, cells)) %>% 
  select(-id)
```

    ```

- 

<!-- -->

    ```{r,summarize data for tp, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
summary <- bottle %>% 
  select(-bottle) %>% 
  group_by(exp, exp_no, phyto, source, tp, datetime, date, time, plot_date, plot_datetime) %>% 
  summarize(across(
    .cols = cells:ssc,
    .fns = list(mean = mean, sd = sd),
    na.rm = TRUE,
    .names = "{fn}_{col}"
  )) %>%
  ungroup() %>% 
  arrange(source, datetime)
```

    ```

- 

############# 

# CALCULATE POC FILTER ABUNDANCES

############## 

Here we subtract the cell concentration in the POC filtrates from the
culture cell abundance to estimate the abundance on the POC filters

    ```{r,calculate filter abundance, message=FALSE, warning=FALSE, wrapper = TRUE}

``` r
poc <- read_csv("~/Box Sync/Phyto_bbp/DATA/FINAL/2022_Experiments/fcm/POC Filtrates/processed_poc_filtrate_data.csv") %>% 
  mutate(source = "Culture") %>% 
  drop_na() %>% 
  left_join(., bottle %>% select(exp, source, bottle, tp, cells)) %>% 
  mutate(filter_cells = cells - filtrate_cells,
         retention = filter_cells/cells,
         slip = filtrate_cells/cells) %>% 
   summarize(across(
    .cols = retention:slip,
    .fns = list(mean = mean, sd = sd),
    na.rm = TRUE,
    .names = "{fn}_{col}"
  ))  
```

    ```

Even when cells slip through and are detected by the flow cytometer,
98.9 +/- 0.9% cells are captured on the filter

- 

########################### 

# SAVE DATA

########################### 

- 

<!-- -->

    ```{r,write csv, wrapper = TRUE}

``` r
write_csv(bottle, "processed_fcm_bottle.csv")
write_csv(summary, "processed_fcm_summary.csv")
```

    ```

- 

END
