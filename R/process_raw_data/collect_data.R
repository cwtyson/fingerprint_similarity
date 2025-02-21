## Collect fingerprint data from birds

## Housekeeping
library(tidyverse)
library(janitor)
library(lubridate)
library(sf)
library(duckplyr)
library(duckdb)
library(rstatix)
library(text2vec)

## New tag file path
newest_tag_path = "~/Library/CloudStorage/GoogleDrive-cwtyson@gmail.com/My Drive/Zebby_tracking_field_data/tags/zebby_tag_log_20250102.xlsx"

## Compare fingerprints between intervals ##########

## Tags
tag_log  <- readxl::read_excel(newest_tag_path) %>% 
  janitor::clean_names() %>% 
  transmute(species,
            tag = gsub("NA",NA, tag),
            sex,
            bird_band,
            section,
            # date,
            # time,
            tag_start_time = parse_date_time(paste(date, time),
                                             tz = "Australia/Broken_Hill",
                                             orders = c("%d.%m.%Y %H:%M:%S",
                                                        "%d.%m.%Y %H:%M")),
            tag_end_time = parse_date_time(paste(end_date, end_time),
                                           tz = "Australia/Broken_Hill",
                                           orders = c("%d.%m.%Y %H:%M:%S",
                                                      "%d.%m.%Y %H:%M")),
            year = format(tag_start_time, "%Y"))  %>% 
  filter(!is.na(tag)) %>% 
  filter(year == 2024)

## Date filter
date_filter <- as.Date("2024-09-01")

## Connection to database
conn <- DBI::dbConnect(duckdb::duckdb((dir = "/Users/tyson/Library/CloudStorage/GoogleDrive-cwtyson@gmail.com/Other computers/My Mac (1)/databases/australia.duckdb")))

## Checkout connection
# DBI::dbListTables(conn)

## Read in detections from database
dets <- dplyr::tbl(conn, "raw") %>%
  
  ## Keep only station ids matching the specified filter
  duckplyr::filter(station_id %in% c("D82AA0A12259", 
                                     "4BA80216EAEB")) %>%
  
  ## Date filter
  duckplyr::filter(time > date_filter) %>%
  
  # duckplyr::filter(tag_id %in% tag_f) %>%
  
  ## Distinct
  duckplyr::distinct(tag_id,
                     node_id,
                     time,
                     .keep_all = T) %>%
  
  duckplyr::collect() %>%
  
  ## Select and rename
  dplyr::transmute(node = toupper(node_id),
                   date_time = lubridate::with_tz(time, tz = "Australia/Broken_Hill"),
                   tag = tag_id,
                   rssi = tag_rssi) %>%
  duckplyr::arrange(date_time)

## Filter
dets_f <- dets %>%
  
  ## Remove incorrect RSSI
  filter(rssi < 0)

## Node codes
node_codes <- readxl::read_excel("/Users/tyson/Google Drive/My Drive/Zebby_tracking_field_data/nodes/node_codes_20230906.xlsx") %>% 
  clean_names() %>% 
  transmute(node = toupper(node),
            node_number = as.character(node_number))

## Node deployment log
node_log <- readxl::read_excel("/Users/tyson/Google Drive/My Drive/Zebby_tracking_field_data/nodes/node_deployment_log_20241204.xlsx") %>% 
  dplyr::transmute(grid_point,
                   node_number = as.character(round(as.numeric(node_number),0)),
                   deployment_time = lubridate::dmy_hm(paste(start_date, start_time), 
                                                       tz = "Australia/Broken_Hill"),
                   removal_time = lubridate::dmy_hm(paste(end_date, end_time),
                                                    tz = "Australia/Broken_Hill")) %>% 
  ## Join node node
  dplyr::left_join(node_codes,
                   by  = "node_number") %>% 
  
  dplyr::select(node,
                grid_point,
                node_number,
                date_time = deployment_time,
                removal_time)  

## Convert to data.table
nodes <- data.table::data.table(node_log, key = c("node", "date_time"))
dets_f <- data.table::data.table(dets_f, key = c("node", "date_time"))

## Rolling join node log to node records
dets_f <- nodes[dets_f, roll = Inf]

## Remove locations without a grid point and where date time is after removal time
dets_f1 <- dets_f %>%
  dplyr::filter(!is.na(grid_point),
                (date_time < removal_time | is.na(removal_time))) %>%
  
  ## Remove impossible RSSI values
  dplyr::filter(rssi < 0) %>%
  dplyr::arrange(tag,
                 date_time) %>%
  dplyr::select(grid_point,
                tag,
                date_time,
                rssi) %>%
  data.frame() 

## Associate tag with correct band. Convert to data.table and do a roiling join.
tag_log <- tag_log %>% 
  rename(date_time = tag_start_time) %>% 
  mutate(tag_start_time = date_time)
tag_log <- data.table::data.table(tag_log, key = c("tag", "date_time"))
dets_f2 <- data.table::data.table(dets_f1, key = c("tag", "date_time"))

## Rolling join tag log to detections
dets_f2 <- tag_log[dets_f2, roll = Inf]

## Remove tags without a band and where detection time is after removal time
dets_f3 <- dets_f2 %>%
  dplyr::filter(date_time > tag_start_time) %>% 
  dplyr::filter(!is.na(bird_band)) %>% 
  dplyr::filter(date_time < tag_end_time | is.na(tag_end_time)) %>%
  dplyr::arrange(tag,
                 date_time) %>%
  dplyr::select(grid_point,
                bird_band,
                tag,
                species,
                sex,
                section,
                date_time,
                rssi) %>%
  data.frame()

## Get fingerprint each 15 seoncds
dets_int <- dets_f3 %>%
  
  ## Keep tagged birds
  filter(tag %in% tag_log$tag) %>% 
  arrange(date_time) %>% 
  # slice(1:100000) %>%
  
  ## 30 second fingerprint
  mutate(date_time_r = floor_date(date_time,unit="30 seconds")) %>% 
  group_by(bird_band,date_time_r,grid_point) %>% 
  summarise(mean_rssi = mean(rssi)) %>% 
  group_by(bird_band,
           date_time_r) %>%
  arrange(bird_band,
          date_time_r,
          desc(mean_rssi)) %>% 
  ungroup()

## Process to keep top 5 nodes in each interval
dets_int_sum <- dets_int %>% 
  
  ## Assign interval values
  group_by(date_time_r) %>% 
  mutate(int = cur_group_id(),
         inds = n()
  ) %>% 
  filter(inds > 1) %>% 
  select(-inds) %>%
  
  ## Keep sample data
  # filter(int < 10000) %>% 
  arrange(int)  %>% 
  ungroup(date_time_r) %>% 
  select(int,date_time_r,bird_band,grid_point,mean_rssi) %>% 
  group_by(int,
           date_time_r,
           bird_band) %>% 
  pivot_wider(names_from=grid_point,
              values_from=mean_rssi,
              names_prefix = "gp_") %>% 
  group_by(bird_band) %>% 
  mutate(ind_id = cur_group_id()) %>% 
  select(int,
         bird_band,
         ind_id,
         everything()) %>% 
  ungroup()

dets_int_sum <- dets_int_sum %>% 
  arrange(int)

saveRDS(dets_int_sum, "./data/fingerprint_dets_wide.RDS")


