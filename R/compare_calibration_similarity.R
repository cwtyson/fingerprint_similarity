## Calibration data comparison
library(tidyverse)
library(janitor)


## Read in calcualte similarity values and process ####### 
sim_df <- readRDS("./data/fingerprint_similarity/calibration/lf/cp_sim_parallel_all_metrics_node_filter.RDS")

## Pivot longer
sim_l <- sim_df %>% 
  pivot_longer(angular:Yule2)

## Stop point coordinates
pt_coords <- readRDS("/Users/tyson/Documents/git/zebby_tracking_field/data/2024/walking_validation/stop_points.RDS") %>% 
  sf::st_transform(3308) %>% 
  mutate(x = sf::st_coordinates(.)[,1],
         y = sf::st_coordinates(.)[,2]) %>% 
  select(pt = pt_name,
         x,
         y) %>% 
  sf::st_drop_geometry()

## Calibration point data
cp_data <- readRDS("/Users/tyson/Documents/git/zebby_tracking_field/data/2024/walking_validation/dets_rounds1-2.RDS")

## Pivot wider
cpw <- cp_data %>% 
  
  filter(!is.na(point)) %>% 
  
  ## 30 second fingerprint
  mutate(date_time_r = floor_date(date_time,unit="30 seconds")) %>% 
  group_by(point,place,type,tag,date_time_r,gp) %>% 
  summarise(mean_rssi = mean(rssi)) %>% 
  pivot_wider(names_from=gp,
              values_from=mean_rssi) %>% 
  
  ## Create interval
  group_by(date_time_r) %>% 
  mutate(int = cur_group_id()) %>%
  
  ## Create interval
  group_by(tag) %>% 
  mutate(ind_id = cur_group_id()) %>% 
  ungroup() %>%
  mutate(row = 1:n()) %>% 
  select(row,point,int,ind_id,everything())  

## Join cp information
sim_df_j <- sim_l %>% 
  left_join(cpw %>% 
              select(row_ind = row,
                     point_ind = point,
                     type_ind = type,
                     place_ind = place),
            by = join_by(row_ind)) %>% 
  left_join(cpw %>% 
              select(row_partner = row,
                     point_partner = point,
                     type_partner = type,
                     place_partner = place),
            by = join_by(row_partner)) %>% 
  
  mutate(point_type = ifelse(point_ind==point_partner,
                             "same_point",
                             "diff_point"),
         
         place_type = case_when((place_ind==place_partner & place_ind == "high") ~ "high",
                                (place_ind==place_partner & place_ind == "low") ~ "low",
                                TRUE ~  "diff_place"),
         tag_type = case_when((type_ind==type_partner & type_ind == "hybrid") ~ "hybrid",
                              (type_ind==type_partner & type_ind == "life") ~ "life",
                              TRUE ~ "diff_tag")) %>% 
  left_join(pt_coords,
            by = join_by(point_ind == pt))  %>%
  rename(x_ind = x,
         y_ind = y) %>% 
  left_join(pt_coords,
            by = join_by(point_partner == pt)) %>% 
  rename(x_partner = x,
         y_partner = y) %>% 
  mutate(dist = sqrt((x_ind - x_partner)^2 + (y_ind - y_partner)^2)) 

## Reformat
sim_ref <- sim_df_j %>% 
  transmute(vec_type,
            point_type,
            num_gp,
            name,
            value,
            dist = as.numeric(dist))  %>% 
  filter(value != "NaN") %>% 
  filter(value != "-Inf") %>% 
  
  ## Remove measures without 
  group_by(vec_type,name) %>% 
  mutate(all_same = ifelse(var(value) == 0,TRUE,FALSE)) %>% 
  filter(all_same == FALSE)  %>% 
  
  ## Scale
  group_by(vec_type,name) %>% 
  mutate(value_s = scale(value)[,1]) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(node_num = as.numeric(unlist(strsplit(vec_type, "_"))[5])) 


## Save
saveRDS(sim_ref, "./data/fingerprint_similarity/calibration/lf/similarity_scores.RDS")


## Read in similarity values and analyze #######
sim_ref <- readRDS("./data/fingerprint_similarity/calibration/lf/similarity_scores.RDS")

## Process
sim_p <-sim_ref %>% 
  
  ## Keep only fingerprints with 20 or fewer
  filter(is.na(node_num) | node_num <= 20) %>% 
  group_by(vec_type,
           name) %>% 
  mutate(group = cur_group_id())

## Fit models for each subset
sim_nested <- sim_p %>% 
  arrange(group) %>% 
  ungroup() %>% 
  ## Keep only records with fewer than 20 nodes
  
  nest_by(group,
          vec_type,
          name) 

## Models

## Get correlation for each group
s_cor <- sim_nested %>% 
  mutate(cor_sp =  list(cor(data$dist,
                            data$value_s,
                            method = "spearman")))

## Select relevant columns
s_cor_f <- s_cor %>%
  transmute(group,vec_type,name,cor=unlist(cor_sp)) %>% 
  arrange(desc(abs(cor))) %>% 
  ungroup()

## Keep top metrics

## Filter data to only keep top models
sim_p_f <- sim_p %>% 
  filter(group %in% (s_cor_f %>% 
                       slice(1:9) %>% 
                       pull(group))) %>% 
  mutate(group_name = paste(name, vec_type)) %>% 
  left_join(s_cor_f) %>% 
  ungroup() %>% 
  arrange(desc(cor)) %>% 
  
  ## Reorder 
  mutate(group_name = factor(group_name, levels=rev(unique(group_name[order(cor)])), ordered=TRUE))


str(sim_p_f$group_name)

## Plot top metrics
ggplot(aes(y=dist,
           x=value_s),
       data = sim_p_f) +
  
  geom_jitter(alpha=0.2) +
  stat_smooth(se = FALSE,
              method="lm") +
  
  ggplot2::scale_colour_manual(values = wesanderson::wes_palette("Zissou1", 
                                                                 9, 
                                                                 type = "continuous"), 
                               name = "# Receivers") +
  facet_wrap(group_name~.)  +
  theme_minimal() +
  theme()

# 
# ## Fit LM for each 
# sim_lm <- sim_nested %>% 
#   ungroup() %>%
#   slice(1:3) %>%
#   rowwise() %>% 
#   mutate(model = list(lm(dist ~ value_s, 
#                          data = data)),
#          cor_sp =  list(cor(data$dist,data$value_s,method = "spearman"))) %>% 
#   summarise(rsq = summary(model)$r.squared) %>% 
#   arrange(desc(rsq))


## Show top 9 metrics
top_metrics <- sim_mods %>% 
  ungroup() %>% 
  slice(1:9)

## Filter data to only keep top models
sim_p_f <- sim_p %>% 
  filter(group %in% top_metrics$group) %>% 
  mutate(group_name = paste(name, vec_type))

## Plot top metrics
ggplot(aes(y=dist,
           x=value_s),
       data = sim_p_f) +
  
  geom_jitter(alpha=0.2) +
  stat_smooth(se = FALSE,
              method="lm") +
  
  ggplot2::scale_colour_manual(values = wesanderson::wes_palette("Zissou1", 
                                                                 9, 
                                                                 type = "continuous"), 
                               name = "# Receivers") +
  facet_wrap(group_name~.)  +
  theme_minimal() +
  theme()


## Compare similarity based on node proportion ########

## Node similarity (already long format)
sim_df_l <- readRDS("./data/fingerprint_similarity/calibration/lf/cp_node_similarity.RDS")

## Stop point coordinates
pt_coords <- readRDS("/Users/tyson/Documents/git/zebby_tracking_field/data/2024/walking_validation/stop_points.RDS") %>% 
  sf::st_transform(3308) %>% 
  mutate(x = sf::st_coordinates(.)[,1],
         y = sf::st_coordinates(.)[,2]) %>% 
  select(pt = pt_name,
         x,
         y) %>% 
  sf::st_drop_geometry()

## Calibration point data
cp_data <- readRDS("/Users/tyson/Documents/git/zebby_tracking_field/data/2024/walking_validation/dets_rounds1-2.RDS")

## Pivot wider
cpw <- cp_data %>% 
  
  filter(!is.na(point)) %>% 
  
  ## 30 second fingerprint
  mutate(date_time_r = floor_date(date_time,unit="30 seconds")) %>% 
  group_by(point,place,type,tag,date_time_r,gp) %>% 
  summarise(mean_rssi = mean(rssi)) %>% 
  pivot_wider(names_from=gp,
              values_from=mean_rssi) %>% 
  
  ## Create interval
  group_by(date_time_r) %>% 
  mutate(int = cur_group_id()) %>%
  
  ## Create interval
  group_by(tag) %>% 
  mutate(ind_id = cur_group_id()) %>% 
  ungroup() %>%
  mutate(row = 1:n()) %>% 
  select(row,point,int,ind_id,everything())  

## Join cp information
sim_df_j <- sim_df_l %>% 
  left_join(cpw %>% 
              select(row_ind = row,
                     point_ind = point,
                     type_ind = type,
                     place_ind = place),
            by = join_by(row_ind)) %>% 
  left_join(cpw %>% 
              select(row_partner = row,
                     point_partner = point,
                     type_partner = type,
                     place_partner = place),
            by = join_by(row_partner)) %>% 
  
  mutate(point_type = ifelse(point_ind==point_partner,
                             "same_point",
                             "diff_point"),
         
         place_type = case_when((place_ind==place_partner & place_ind == "high") ~ "high",
                                (place_ind==place_partner & place_ind == "low") ~ "low",
                                TRUE ~  "diff_place"),
         tag_type = case_when((type_ind==type_partner & type_ind == "hybrid") ~ "hybrid",
                              (type_ind==type_partner & type_ind == "life") ~ "life",
                              TRUE ~ "diff_tag")) %>% 
  left_join(pt_coords,
            by = join_by(point_ind == pt))  %>%
  rename(x_ind = x,
         y_ind = y) %>% 
  left_join(pt_coords,
            by = join_by(point_partner == pt)) %>% 
  rename(x_partner = x,
         y_partner = y) %>% 
  mutate(dist = sqrt((x_ind - x_partner)^2 + (y_ind - y_partner)^2)) 

## Lines
ggplot(aes(x=dist,
           y=n_vec_score,
           color=n_filter,
           group=n_filter),
       data = sim_df_j) +
  geom_point() +
  stat_smooth(method="lm") +
  facet_wrap(n_filter~.)

## Boxplot
ggplot(aes(x=n_filter,
           y=n_vec_score,
           group=n_filter,
           color=n_filter),
       data = sim_df_j) +
  geom_boxplot() +
  facet_grid(.~point_type)

## Alpha values ######

## Read in calcualte similarity values and process 
alpha_df <- readRDS("./data/fingerprint_similarity/calibration/lf/alpha.RDS")

ggplot(alpha_df) +
  geom_point(aes(x=num_gp,y=Podani))

## Pivot longer
sim_l <- alpha_df %>% 
  pivot_longer(Podani)

## Stop point coordinates
pt_coords <- readRDS("/Users/tyson/Documents/git/zebby_tracking_field/data/2024/walking_validation/stop_points.RDS") %>% 
  sf::st_transform(3308) %>% 
  mutate(x = sf::st_coordinates(.)[,1],
         y = sf::st_coordinates(.)[,2]) %>% 
  select(pt = pt_name,
         x,
         y) %>% 
  sf::st_drop_geometry()

## Calibration point data
cp_data <- readRDS("/Users/tyson/Documents/git/zebby_tracking_field/data/2024/walking_validation/dets_rounds1-2.RDS")

## Pivot wider
cpw <- cp_data %>% 
  
  filter(!is.na(point)) %>% 
  
  ## 30 second fingerprint
  mutate(date_time_r = floor_date(date_time,unit="30 seconds")) %>% 
  group_by(point,place,type,tag,date_time_r,gp) %>% 
  summarise(mean_rssi = mean(rssi)) %>% 
  pivot_wider(names_from=gp,
              values_from=mean_rssi) %>% 
  
  ## Create interval
  group_by(date_time_r) %>% 
  mutate(int = cur_group_id()) %>%
  
  ## Create interval
  group_by(tag) %>% 
  mutate(ind_id = cur_group_id()) %>% 
  ungroup() %>%
  mutate(row = 1:n()) %>% 
  select(row,point,int,ind_id,everything())  

## Join cp information
sim_df_j <- sim_l %>% 
  left_join(cpw %>% 
              select(row_ind = row,,
                     point_ind = point,
                     type_ind = type,
                     place_ind = place),
            by = join_by(row_ind)) %>% 
  left_join(cpw %>% 
              select(row_partner = row,
                     point_partner = point,
                     type_partner = type,
                     place_partner = place),
            by = join_by(row_partner)) %>% 
  
  mutate(point_type = ifelse(point_ind==point_partner,
                             "same_point",
                             "diff_point"),
         
         place_type = case_when((place_ind==place_partner & place_ind == "high") ~ "high",
                                (place_ind==place_partner & place_ind == "low") ~ "low",
                                TRUE ~  "diff_place"),
         tag_type = case_when((type_ind==type_partner & type_ind == "hybrid") ~ "hybrid",
                              (type_ind==type_partner & type_ind == "life") ~ "life",
                              TRUE ~ "diff_tag")) %>% 
  left_join(pt_coords,
            by = join_by(point_ind == pt))  %>%
  rename(x_ind = x,
         y_ind = y) %>% 
  left_join(pt_coords,
            by = join_by(point_partner == pt)) %>% 
  rename(x_partner = x,
         y_partner = y) %>% 
  mutate(dist = sqrt((x_ind - x_partner)^2 + (y_ind - y_partner)^2)) 

## Reformat
sim_ref <- sim_df_j %>% 
  transmute(vec_type,
            alpha,
            point_ind,
            point_partner,
            point_type,
            num_gp,
            name,
            value,
            dist = as.numeric(dist))  %>% 
  filter(value != "NaN") %>% 
  filter(value != "-Inf") %>% 
  
  ## Remove measures without 
  group_by(vec_type,name) %>% 
  mutate(all_same = ifelse(var(value) == 0,TRUE,FALSE)) %>% 
  filter(all_same == FALSE)  %>% 
  
  ## Scale
  group_by(vec_type,name) %>% 
  mutate(value_s = scale(value)[,1]) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(node_num = as.numeric(unlist(strsplit(vec_type, "_"))[5])) 


## Save
saveRDS(sim_ref, "./data/fingerprint_similarity/calibration/lf/alpha_ref.RDS")


## Read in similarity values and analyze 
sim_ref <- readRDS("./data/fingerprint_similarity/calibration/lf/alpha_ref.RDS")

## Plot
sim_ref %>% 
  filter(vec_type == "exp_NA_remove") %>% 
  ggplot() +
  geom_boxplot(aes(x=point_type,
                   y=value,
                   # group=as.factor(alpha),
                   fill=as.factor(alpha)))


## Process
sim_p <-sim_ref %>% 
  
  ## Keep only fingerprints with 20 or fewer
  filter(is.na(node_num) | node_num <= 20) %>% 
  group_by(vec_type,
           name) %>% 
  mutate(group = cur_group_id())

## Fit models for each subset
sim_nested <- sim_p %>% 
  arrange(group) %>% 
  ungroup() %>% 
  ## Keep only records with fewer than 20 nodes
  
  nest_by(group,
          vec_type,
          name) 

## Models

## Get correlation for each group
s_cor <- sim_nested %>% 
  mutate(cor_sp =  list(cor(data$dist,
                            data$value_s,
                            method = "spearman")))

## Select relevant columns
s_cor_f <- s_cor %>%
  transmute(group,vec_type,name,cor=unlist(cor_sp)) %>% 
  arrange(desc(abs(cor))) %>% 
  ungroup()

## Keep top metrics

## Filter data to only keep top models
sim_p_f <- sim_p %>% 
  filter(group %in% (s_cor_f %>% 
                       slice(1:9) %>% 
                       pull(group))) %>% 
  mutate(group_name = paste(name, vec_type)) %>% 
  left_join(s_cor_f) %>% 
  ungroup() %>% 
  arrange(desc(cor)) %>% 
  
  ## Reorder 
  mutate(group_name = factor(group_name, levels=rev(unique(group_name[order(cor)])), ordered=TRUE))


str(sim_p_f$group_name)

## Plot top metrics
ggplot(aes(y=dist,
           x=value_s),
       data = sim_p_f) +
  
  geom_jitter(alpha=0.2) +
  stat_smooth(se = FALSE,
              method="lm") +
  
  ggplot2::scale_colour_manual(values = wesanderson::wes_palette("Zissou1", 
                                                                 9, 
                                                                 type = "continuous"), 
                               name = "# Receivers") +
  facet_wrap(group_name~.)  +
  theme_minimal() +
  theme()

# 
# ## Fit LM for each 
# sim_lm <- sim_nested %>% 
#   ungroup() %>%
#   slice(1:3) %>%
#   rowwise() %>% 
#   mutate(model = list(lm(dist ~ value_s, 
#                          data = data)),
#          cor_sp =  list(cor(data$dist,data$value_s,method = "spearman"))) %>% 
#   summarise(rsq = summary(model)$r.squared) %>% 
#   arrange(desc(rsq))


## Show top 9 metrics
top_metrics <- sim_mods %>% 
  ungroup() %>% 
  slice(1:9)

## Filter data to only keep top models
sim_p_f <- sim_p %>% 
  filter(group %in% top_metrics$group) %>% 
  mutate(group_name = paste(name, vec_type))

## Plot top metrics
ggplot(aes(y=dist,
           x=value_s),
       data = sim_p_f) +
  
  geom_jitter(alpha=0.2) +
  stat_smooth(se = FALSE,
              method="lm") +
  
  ggplot2::scale_colour_manual(values = wesanderson::wes_palette("Zissou1", 
                                                                 9, 
                                                                 type = "continuous"), 
                               name = "# Receivers") +
  facet_wrap(group_name~.)  +
  theme_minimal() +
  theme()


## Compare similarity based on node proportion ########

## Node similarity (already long format)
sim_df_l <- readRDS("./data/fingerprint_similarity/calibration/lf/cp_node_similarity.RDS")

## Stop point coordinates
pt_coords <- readRDS("/Users/tyson/Documents/git/zebby_tracking_field/data/2024/walking_validation/stop_points.RDS") %>% 
  sf::st_transform(3308) %>% 
  mutate(x = sf::st_coordinates(.)[,1],
         y = sf::st_coordinates(.)[,2]) %>% 
  select(pt = pt_name,
         x,
         y) %>% 
  sf::st_drop_geometry()

## Calibration point data
cp_data <- readRDS("/Users/tyson/Documents/git/zebby_tracking_field/data/2024/walking_validation/dets_rounds1-2.RDS")

## Pivot wider
cpw <- cp_data %>% 
  
  filter(!is.na(point)) %>% 
  
  ## 30 second fingerprint
  mutate(date_time_r = floor_date(date_time,unit="30 seconds")) %>% 
  group_by(point,place,type,tag,date_time_r,gp) %>% 
  summarise(mean_rssi = mean(rssi)) %>% 
  pivot_wider(names_from=gp,
              values_from=mean_rssi) %>% 
  
  ## Create interval
  group_by(date_time_r) %>% 
  mutate(int = cur_group_id()) %>%
  
  ## Create interval
  group_by(tag) %>% 
  mutate(ind_id = cur_group_id()) %>% 
  ungroup() %>%
  mutate(row = 1:n()) %>% 
  select(row,point,int,ind_id,everything())  

## Join cp information
sim_df_j <- sim_df_l %>% 
  left_join(cpw %>% 
              select(row_ind = row,
                     point_ind = point,
                     type_ind = type,
                     place_ind = place),
            by = join_by(row_ind)) %>% 
  left_join(cpw %>% 
              select(row_partner = row,
                     point_partner = point,
                     type_partner = type,
                     place_partner = place),
            by = join_by(row_partner)) %>% 
  
  mutate(point_type = ifelse(point_ind==point_partner,
                             "same_point",
                             "diff_point"),
         
         place_type = case_when((place_ind==place_partner & place_ind == "high") ~ "high",
                                (place_ind==place_partner & place_ind == "low") ~ "low",
                                TRUE ~  "diff_place"),
         tag_type = case_when((type_ind==type_partner & type_ind == "hybrid") ~ "hybrid",
                              (type_ind==type_partner & type_ind == "life") ~ "life",
                              TRUE ~ "diff_tag")) %>% 
  left_join(pt_coords,
            by = join_by(point_ind == pt))  %>%
  rename(x_ind = x,
         y_ind = y) %>% 
  left_join(pt_coords,
            by = join_by(point_partner == pt)) %>% 
  rename(x_partner = x,
         y_partner = y) %>% 
  mutate(dist = sqrt((x_ind - x_partner)^2 + (y_ind - y_partner)^2)) 

## Lines
ggplot(aes(x=dist,
           y=n_vec_score,
           color=n_filter,
           group=n_filter),
       data = sim_df_j) +
  geom_point() +
  stat_smooth(method="lm") +
  facet_wrap(n_filter~.)

## Boxplot
ggplot(aes(x=n_filter,
           y=n_vec_score,
           group=n_filter,
           color=n_filter),
       data = sim_df_j) +
  geom_boxplot() +
  facet_grid(.~point_type)

