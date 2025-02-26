## Fingerprint similarity for calibration data

## Housekeeping
library(tidyverse)
remotes::install_github("vh-d/RPortfolioSimilarity")
library(RPortfolioSimilarity)
library(text2vec)
library(foreach)
library(proxy)

cores = parallel::detectCores()
cl <- parallel::makeForkCluster(cores-3, outfile = "")
doParallel::registerDoParallel(cl)

## Hardcoded values

## RSSI cutoff to turn into NA
rssi_cutoff = -90

## Alpha for exponential representation
alpha = 24

## Minimum rssi
min_rssi = -117

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
cpw <- readRDS("./data/fingerprint_similarity/calibration/calibration_data_wide.RDS")

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


## Similarity and distance measures - 
sim_measures = c("angular","Braun-Blanquet","Chi-squared","correlation","cosine","Cramer",
                 "Dice","eDice","eJaccard","Fager","Hamman","Jaccard",
                 "Kulczynski1","Kulczynski2","Michael","Mountford","Mozley","Ochiai",
                 "Pearson","Phi","Phi-squared","Russel","simple matching","Simpson","Stiles","Tanimoto",
                 "Tschuprow","Yule","Yule2")

dist_measures = c("Bhjattacharyya","Bray","Canberra","Chord","divergence","Euclidean",
                  "Geodesic","Hellinger","Kullback","Manhattan",
                  "Podani","Soergel","supremum","Wave","Whittaker")

## Combine
measures <- sort(c(sim_measures,dist_measures))

## Fitting function #######
fitting_function <- function(row_i)
{
  
  sim_df_output <- list()
  
  # row_i = 52312
  
  cat(row_i,'\n')
  
  rows <- c(combos[row_i,]$x,combos[row_i,]$y)
  
  ## Keep fingerprints where at least one value is not NA
  int_sub_inds <- cpw %>% 
    filter(row %in% rows) %>% 
    select(where(~!any(is.na(.)))) %>% 
    ungroup()
  
  ## Number of grid points
  num_gp <- ncol(int_sub_inds)-8
  
  ## Calculate similarity if at least 1 node
  if(num_gp >= 1){
    
    ## Remove NAs
    linear_NA_remove<- int_sub_inds %>% 
      select(where(~!any(is.na(.))))
    
    ## Replace NA with minimum value
    linear_NA_replace <- int_sub_inds %>% 
      mutate(across(starts_with("gp_"), ~case_when(is.na(.x) ~ -120,
                                                   TRUE ~ .x)))
    
    ## Remove NAS, then create exponential representation 
    exp_NA_remove <- int_sub_inds %>% 
      
      select(where(~!any(is.na(.)))) %>% 
      
      ## Change to positive representation
      mutate(across(starts_with("gp_"), ~case_when(!is.na(.x) ~ (.x - -117),
                                                   TRUE ~ 0))) %>% 
      
      ## Change to exponeitial representation
      mutate(across(starts_with("gp_"), ~(exp(.x/24))/(exp(--117/24))))
    
    ## Replace NAs, then create exponential representation
    exp_NA_replace <- int_sub_inds %>% 
      
      ## Replace NAs  
      mutate(across(starts_with("gp_"), ~case_when(is.na(.x) ~ -120,
                                                   TRUE ~ .x))) %>% 
      
      ## Change to positive representation
      mutate(across(starts_with("gp_"), ~case_when(!is.na(.x) ~ (.x - -117),
                                                   TRUE ~ 0))) %>% 
      
      ## Change to exponeitial representation
      mutate(across(starts_with("gp_"), ~(exp(.x/24))/(exp(--117/24))))
    
    
    ## Get max number of nodes to use given number of columns in vec_NA_repl
    nodes <- ncol(exp_NA_remove)-8
    
    ## List for vectors
    vec_exp_n_list <- list()
    for(nn in 1:nodes){
      
      ## Create exponential representation of full vector with only top X nodes
      vec_exp_n_list[[nn]] <- exp_NA_remove %>% 
        group_by(row) %>% 
        pivot_longer(starts_with("gp_")) %>% 
        arrange(desc(value)) %>% 
        ungroup() %>% 
        mutate(gp_rank = 1:n()) %>% 
        group_by(name) %>% 
        mutate(gp_rank = min(gp_rank)) %>% 
        filter(gp_rank <=  nn) %>% 
        select(-gp_rank) %>% 
        ungroup() %>% 
        pivot_wider(names_from = name) 
      
    }
    
    vec_exp_n_list_names <- paste0("vec_list_exp_n_", 1:nodes)
    
    ## Minimum vector length will be vec_NA_remove, if at least 1 column:
    if((ncol(linear_NA_remove)-8)>=1){
      
      ## List of vectors
      vec_list <- c(list(linear_NA_remove,
                         linear_NA_replace,
                         exp_NA_replace,
                         exp_NA_remove),
                    vec_exp_n_list)
      
      ## Matrices
      vec_list_mat <- lapply(vec_list,
                             function(x) x %>% 
                               select(-int,-tag,-ind_id,-date_time_r,-point,-place,-type,-row) %>% 
                               as.matrix() 
      )
      
      num_gp = unlist(lapply(vec_list_mat,
                             ncol))
      
      # cat("num_gp:", num_gp[4],"\n")
      
      
      ## Replicate matrices and methods equal number of times
      vec_list_mat_rep <- rep(vec_list_mat, 
                              each = length(measures))
      
      measures_rep <- rep(measures, 
                          length(vec_list_mat))
      
      ## Calculate scores from similarity
      scores_df <- map2(.x = vec_list_mat_rep,
                        .y = measures_rep,
                        .f = function(x,y) suppressWarnings(proxy::dist(x,
                                                                        method = y))) %>% 
        unlist() %>% 
        matrix(nrow = length(vec_list),
               byrow = T) %>% 
        data.frame() %>% 
        mutate(vec_type = c("linear_NA_remove",
                            "linear_NA_replace",
                            "exp_NA_replace",
                            "exp_NA_remove",
                            vec_exp_n_list_names),
               num_gp = num_gp)
      names(scores_df) <- c(measures, "vec_type","num_gp")
      
      sim_df_output[[1]] <- scores_df %>% 
        mutate(row = row_i,
               ind = int_sub_inds$tag[1],
               partner = int_sub_inds$tag[2],
               row_ind = int_sub_inds$row[1],
               row_partner = int_sub_inds$row[2]) %>% 
        select(
          row,
          ind,
          partner,
          row_ind,
          row_partner,
          num_gp,
          vec_type,
          everything()
        )
      
    } 
    
  } 

}

## For each
cp_sim <- foreach(row_i=1:nrow(combos),
                  .packages='tidyverse',
                  .verbose = FALSE) %dopar%
  { fitting_function(row_i) }


alld <- do.call(rbind, cp_sim) %>% 
  dplyr::select(row,ind,partner,row_ind,row_partner,num_gp,everything())

saveRDS(alld, "./data/fingerprint_similarity/calibration/lf/cp_sim_parallel_all_metrics_node_filter.RDS")
