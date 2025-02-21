## Combine and assess similarity
library(tidyverse)

## Tag path
newest_tag_path = "~/Library/CloudStorage/GoogleDrive-cwtyson@gmail.com/My Drive/Zebby_tracking_field_data/tags/zebby_tag_log_20250102.xlsx"
files <- list.files("./data/fingerprint_similarity/similarity_estimates/", full.names = T)

## Combine
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
  
  ## Adjust order of individual and partner
  rowwise() %>%
  mutate(ind = min(c(ind[1],partner[1])),
         partner = max(c(ind[1],partner[1])))  %>% 
  ungroup() %>% 
  
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
  mutate(family = ifelse(ind_section == partner_section, "yes","no"),
         section = substring(ind_section,1,1),
         hour = as.numeric(format(lubridate::round_date(date_time_r, unit = "hour"),"%H")),
         group_type = case_when((ind_sex %in% c("Male","Female") & partner_sex %in% c("Male","Female")) ~ "Adults",
                                ind_sex == "Juvenile" & partner_sex == "Juvenile" ~ "Juveniles",
                                TRUE ~ "Adult + Juv"),
         groupID = paste(ind,partner,sep="-")) %>% 
  select(groupID,everything())

## Summarise scores
sim_sum <- sim_df_j %>% 
  dplyr::filter(nodes > 5) %>% 
  dplyr::group_by(groupID,
                  family,
                  section,
                  group_type,
                  ind_section,
                  partner_section) %>% 
  dplyr::summarise(mean_sim = 2-mean(score,na.rm = T),
                   counts = n()) %>% 
  dplyr::mutate(group_type = factor(group_type,
                                    levels=c("Juveniles",
                                             "Adult + Juv",
                                             "Adults")))

## Plot similarity between groups
ggplot(sim_sum %>% 
         filter(section %in% c("D","E"))) +
  geom_violin(aes(x=family,
                  y=mean_sim)) +
  ggbeeswarm::geom_quasirandom(
    aes(x=family,
        y=mean_sim,
        alpha=counts,
        color=counts)) +
  facet_grid(section~group_type) +
  scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1", 100, type = "continuous")) + 
  theme_minimal() +
  labs(x="family",y="similarity index") +
  theme(legend.position = "none")

## Plot similarity over time
d2p <- sim_df_j %>% 
  filter(section %in% c("D","E")) %>% 
  filter(family == "yes") %>% 
  filter(group_type == "Adult + Juv")

ggplot(d2p) +
  
  stat_smooth(aes(x=date_time_r,
                  y=score,
                  color=groupID,
                  group=groupID),
              method="lm",
              se = FALSE) +
  facet_grid(group_type~.) +
  scale_color_manual(values = as.character(wesanderson::wes_palette("Zissou1", n_distinct(d2p$groupID), type = "continuous"))) +
  theme_minimal() +
  labs(x="family",y="similarity index") +
  theme(legend.position = "none")

