## Calibration data comparison
library(tidyverse)
library(janitor)

## Read in similarity values and process ####### 
sim_df <- readRDS("./data/fingerprint_similarity/calibration/cp_sim_parallel_all_metrics_node_filter.RDS")

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
saveRDS(sim_ref, "./data/fingerprint_similarity/calibration/similarity_scores.RDS")


## Read in similarity values and visualize #######
sim_ref <- readRDS("./data/fingerprint_similarity/calibration/similarity_scores.RDS")

## Compare node X filters

## Prepare data
diff_df <- sim_ref %>% 
  
  ## Keep node filters
  filter(grepl("vec_list_exp", vec_type)) %>% 
  # filter(num_gp < 30) %>%
  
  ## Keep distances less than 300 m
  filter(dist < 300) 

## Comparison of distance and scores
ggplot(aes(x=dist,
           y=value_s),
       data = diff_df %>% 
         filter(node_num < 20)) +
  
  stat_smooth(aes(group=node_num,
                  color=node_num),
              se = TRUE,
              method="lm") +
  ggplot2::scale_colour_gradientn(colours = wesanderson::wes_palette("Zissou1", 
                                                                     100, 
                                                                     type = "continuous"), 
                                  name = "# Receivers") +
  facet_wrap(name~.)  +
  theme_minimal() +
  
  labs(x = "Distance (m)",
       y = "Scaled value")


## Boxplot comparing same and different points
ggplot(aes(x=as.factor(point_type),
           fill=vec_type,
           y=value_s), 
       data = sim_ref %>% 
         slice_sample(prop=1)) +
  geom_boxplot(outliers = FALSE)  +
  facet_wrap(name~.,
             scales="free")  +
  ggplot2::scale_fill_manual(values = wesanderson::wes_palette("Zissou1", n_distinct(sim_ref$vec_type), 
                                                               type = "continuous")) +
  theme_minimal()# ggsave("./plots/sim_score_boxplots.JPG",scale = 2)

## Prepare data
diff_df <- sim_ref %>% 
  
  ## Keep node filters
  filter(grepl("vec_list_exp", vec_type)) %>% 
  # filter(num_gp < 30) %>%
  
  ## Keep distances less than 300 m
  filter(dist < 300)

## Comparison of distance and scores
ggplot(aes(x=dist,
           y=value_s),
       data = diff_df) +
  
  stat_smooth(aes(group=node_num,
                  color=node_num),
              se = FALSE,
              method="lm") +
  ggplot2::scale_colour_gradientn(colours = wesanderson::wes_palette("Zissou1", 
                                                                     100, 
                                                                     type = "continuous"), 
                                  name = "# Receivers") +
  facet_wrap(name~.)  +
  theme_minimal() +
  theme()

# ggsave("./plots/sim_score_dist.JPG",scale = 2)

## Read in similarity values and analyze #######
sim_ref <- readRDS("./data/fingerprint_similarity/calibration/similarity_scores.RDS") %>% 
  
  ## Keep only fingerprints with 20 or fewer
  filter(is.na(node_num) | node_num <= 20) %>% 
  group_by(vec_type,
           name) %>% 
  mutate(group = cur_group_id())

## Fit models for each subset
sim_nested <- sim_ref %>% 
  arrange(group) %>% 
  ungroup() %>% 
  ## Keep only records with fewer than 20 nodes
  
  nest_by(group,
          vec_type,
          name) 

## Models

## Fit LM for each 
sim_mods <- sim_nested %>% 
  # slice(1:5) %>% 
  mutate(model = list(lm(dist ~ value_s, 
                         data = data))) %>% 
  summarise(rsq = summary(model)$r.squared) %>% 
  arrange(desc(rsq))


## Show top 9 metrics
top_metrics <- sim_mods %>% 
  ungroup() %>% 
  slice(1:9)

sim_ref_f <- sim_ref %>% 
  filter(group %in% top_metrics$group) %>% 
  mutate(group_name = paste(name, vec_type))

## Plot top metrics
ggplot(aes(y=dist,
           x=value_s),
       data = sim_ref_f) +
  
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















