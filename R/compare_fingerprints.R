## Combine and assess similarity
library(tidyverse)

newest_tag_path = "~/Library/CloudStorage/GoogleDrive-cwtyson@gmail.com/My Drive/Zebby_tracking_field_data/tags/zebby_tag_log_20250102.xlsx"

files <- list.files("./data/fingerprint_similarity/similarity_estimates/",full.names = T)


sim_df <- lapply(files,readRDS) %>% 
  do.call(rbind,.)

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


## Add group information
sim_df_j <- sim_df %>%
  left_join(tag_log %>%
              select(ind = bird_band,
                     ind_sex = sex,
                     ind_section = section) %>%
              distinct(.keep_all = T)) %>%
  left_join(tag_log %>%
              select(partner = bird_band,
                     partner_sex  = sex,
                     partner_section = section) %>%
              distinct(.keep_all = T)
  ) %>%
  
  mutate(same_group = ifelse(ind_section == partner_section, "yes","no"),
         section = substring(ind_section,1,1),
         hour = as.numeric(format(lubridate::round_date(date_time_r, unit = "hour"),"%H")),
         group_type = case_when((ind_sex %in% c("Male","Female") | partner_sex %in% c("Male","Female")) ~ "Pair",
                                TRUE ~ "Adult + Juv"))


ggplot(sim_df_j) +
  geom_violin(aes(x=same_group,
                   y=score)) 


## Summarise scores
sim_sum <- sim_df_j %>% 
  group_by(same_group,
           section,
           group_type,
           ind_section) %>% 
  summarise(mean_sim = mean(score))


ggplot(sim_sum) +
  
  geom_boxplot(aes(x=same_group,
                   y=mean_sim)) +
  geom_jitter(aes(x=same_group,
                  y=mean_sim)) +
  facet_grid(section~group_type) +
  theme_minimal()


