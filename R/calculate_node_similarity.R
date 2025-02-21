## Fingerprint similarity for calibration data

## Housekeeping
library(tidyverse)
remotes::install_github("vh-d/RPortfolioSimilarity")
library(RPortfolioSimilarity)
library(text2vec)
library(foreach)
library(proxy)

cores = parallel::detectCores()
cl <- parallel::makeForkCluster(cores-1, outfile = "")
doParallel::registerDoParallel(cl)

## Hardcoded values

# ## Calibration point data #######
# cp_data <- readRDS("/Users/tyson/Documents/git/zebby_tracking_field/data/2024/walking_validation/dets_rounds1-2.RDS")
# 
# ## Pivot wider
# cpw <- cp_data %>% 
#   
#   filter(!is.na(point)) %>% 
#   
#   ## 30 second fingerprint
#   mutate(date_time_r = floor_date(date_time,unit="30 seconds")) %>% 
#   group_by(point,place,type,tag,date_time_r,gp) %>% 
#   summarise(mean_rssi = mean(rssi)) %>% 
#   pivot_wider(names_from=gp,
#               values_from=mean_rssi,
#               names_prefix = "gp_") %>% 
#   
#   ## Create interval
#   group_by(date_time_r) %>% 
#   mutate(int = cur_group_id()) %>%
#   
#   ## Create interval
#   group_by(tag) %>% 
#   mutate(ind_id = cur_group_id()) %>% 
#   ungroup() %>%
#   mutate(row = 1:n()) %>% 
#   select(row,point,int,ind_id,everything()) %>% 
#   slice_sample(prop=1)
# 
# saveRDS(cpw, "./data/fingerprint_similarity/calibration/calibration_data_wide.RDS")

## Read in fingerprint data
cpw <- readRDS("./data/fingerprint_similarity/calibration/calibration_data_wide.RDS") %>% 
  arrange(row)

## Get pairwise combinations for each fingerprint
combos <- expand.grid(x = unique(cpw$row),
                      y =unique(cpw$row)) %>% 
  filter(x != y) %>% 
  rowwise() %>% 
  mutate(order = paste(min(x,y),max(x,y))) %>% 
  distinct(order,.keep_all = T) %>% 
  arrange(x,
          y) %>% 
  select(x,
         y) %>% 
  ungroup() %>% 
  mutate(row = 1:n())

## Nodes to keep
ns <- 1:10

## Fitting function #######
fitting_function <- function(row_i)
{
  
  # row_i = 616605
  
  cat(row_i,"\n")
  
  rows <- c(combos[row_i,]$x,combos[row_i,]$y)
  
  ## Keep only matching
  int_sub_inds_node_f <- cpw %>% 
    filter(row %in% rows)
  
  ## Process
  int_sub_inds_node <- int_sub_inds_node_f %>% 
    select(where(~!all(is.na(.)))) %>% 
    pivot_longer(contains("gp_"),
                 names_to = "grid_point") %>% 
    arrange(desc(value))
  
  ## Calculate score based on n nodes
  n_vec <- c()
  for(n in ns){
    
    ## Score
    score_n <- int_sub_inds_node %>% 
      group_by(int) %>% 
      slice(1:n) %>% 
      group_by(grid_point) %>% 
      summarise(count = n()-1) %>% 
      pull(count) %>% 
      sum() / n
    
    n_vec <- c(n_vec,score_n)
    
  }
  
  
  ## Data frame for proportional scores
  sim_n_df <- data.frame(n_vec_score = n_vec,
                         n_filter = ns)
  
  ## Combine 
  sim_n_df_full = sim_n_df %>% 
    mutate(row = row_i,
           ind = int_sub_inds_node_f$tag[1],
           partner = int_sub_inds_node_f$tag[2],
           row_ind = int_sub_inds_node_f$row[1],
           row_partner = int_sub_inds_node_f$row[2])
  
} 


## For each
cp_sim <- foreach(row_i=1:nrow(combos),
                  .packages='tidyverse',
                  .verbose = FALSE) %dopar%
  { fitting_function(row_i) }


## Combine
alld <- do.call(rbind, cp_sim)

## Save
saveRDS(alld, 
        "./data/fingerprint_similarity/calibration/lf/cp_node_similarity.RDS")
