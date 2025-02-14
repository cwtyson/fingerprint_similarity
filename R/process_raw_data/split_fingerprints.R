
## Split fingerprint files

## Tag checks 

## Housekeeping
library(tidyverse)

## Process detections
dets <- readRDS("./data/fingerprint_dets_wide.RDS")

## Summarize detections
dets_int_sum <- dets %>% 
  group_by(int) %>% 
  mutate(inds = n_distinct(bird_band))

dets_int_sum_f <- dets_int_sum %>% 
  filter(inds > 1) %>% 
  select(-inds)

## Split into days 
dets_int_sum_f_daily <- dets_int_sum_f %>% 
  mutate(day = floor_date(date_time_r, unit= "day"))  %>% 
  select(day,int,date_time_r,bird_band,ind_id,everything())

dl <- split(dets_int_sum_f_daily, dets_int_sum_f_daily$day)
names(dl) = lapply(dl, function(x) as.character(unique(x$day)))

## Save as separte file
purrr::iwalk(dl, function(dat, name) saveRDS(dat, file = paste0("./data/fingerprint_similarity/days/", name, ".RDS")))
