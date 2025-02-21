## Explore alpha values for exponential repesentation

## Housekeeping #######

library(tidyverse)
library(foreach)
library(proxy)

cores = parallel::detectCores()
cl <- parallel::makeForkCluster(cores-3, outfile = "")
doParallel::registerDoParallel(cl)

# ## Compare two fingerprints with different alpha values ######
# 
# na_replace = -120
# 
# ## Increasing alpha decreases disimilarity between values - larger values increases disimilarity for strong values relative to weak values
# alpha_i = 24
# 
# ## Read in fingerprint data
# cpw <- readRDS("./data/fingerprint_similarity/calibration/calibration_data_wide.RDS")
# 
# # ## Get count of non-NA values per fingerprint
# # ints_f <- cpw %>% 
# #   select(row, contains("gp_")) %>% 
# #   group_by(row) %>%
# #   summarise_all(list(~sum(!is.na(.))))  %>% 
# #   rowwise() %>% 
# #   mutate(sum = sum(c_across(contains("gp_")), na.rm = T)) %>% 
# #   select(row,sum) %>% 
# #   left_join(cpw,
# #             by = join_by(row)) %>% 
# #   arrange(desc(sum)) %>% 
# #   select(sum,row,point,place,type,tag, everything()) %>% 
# #   ungroup()
# 
# ## Reformat for plotting
# 
# ftest <- cpw %>% 
#   arrange(row) %>% 
#   slice(c(1,169)) 
# 
# fref_raw <- ftest %>% 
#   arrange(row) %>% 
#   slice(1:2) %>% 
#   
#   select(where(~!any(is.na(.)))) %>% 
#   
#   # ## Change to positive representation
#   # mutate(across(starts_with("gp_"), ~case_when(!is.na(.x) ~ (.x - -117),
#   #                                              TRUE ~ 0))) %>% 
#   # 
#   # ## Change to exponential representation
#   # mutate(across(starts_with("gp_"), ~(exp(.x/alpha_i))/(exp(--117/alpha_i)))) %>%  
#   
#   pivot_longer(contains("gp_")) %>% 
#   ungroup() %>% 
#   select(row,name,value) %>%
#   mutate(value = scale(value)[,1]) %>% 
#   # replace(is.na(.),na_replace) %>% 
#   pivot_wider(names_from = row, 
#               values_from = value) %>% 
#   rename_at(2, ~"row_1") %>% 
#   rename_at(3, ~"row_2") 
# 
# alpha_i = 20
# fref_trans <- ftest %>% 
#   arrange(row) %>% 
#   slice(1:2) %>% 
#   
#   select(where(~!any(is.na(.)))) %>% 
#   
#   ## Change to positive representation
#   mutate(across(starts_with("gp_"), ~case_when(!is.na(.x) ~ (.x - -117),
#                                                TRUE ~ 0))) %>% 
#   
#   ## Change to exponeitial representation
#   mutate(across(starts_with("gp_"), ~(exp(.x/alpha_i))/(exp(--117/alpha_i)))) %>% 
#   pivot_longer(contains("gp_")) %>% 
#   select(row,name,value) %>%
#   mutate(value = scale(value)[,1]) %>% 
#   # replace(is.na(.),na_replace) %>% 
#   pivot_wider(names_from = row, 
#               values_from = value) %>% 
#   rename_at(2, ~"row_1") %>% 
#   rename_at(3, ~"row_2") 
# 
# ## Visualize top fingerprints
# ggplot() +
#   geom_abline(slope = 1,intercept=0) +
#   geom_point(aes(x=row_1,
#                  y=row_2),
#              data= fref_trans) +
#   geom_point(aes(x=row_1,
#                  y=row_2),
#              data= fref_raw,
#              color="gold")


## Hardcoded values ######

## RSSI value to replace NA values
na_rssi = -120

## Alpha for exponential representation
alpha = 1:30

## Read in fingerprint data
cpw <- readRDS("./data/fingerprint_similarity/calibration/calibration_data_wide.RDS") 

## Get pairwise combinations for each fingerprint
combos <- expand.grid(x = unique(cpw$row),
                      y = unique(cpw$row)) %>% 
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
# sim_measures = c("angular","Braun-Blanquet","Chi-squared","correlation","cosine","Cramer",
#                  "Dice","eDice","eJaccard","Fager","Hamman","Jaccard",
#                  "Kulczynski1","Kulczynski2","Michael","Mountford","Mozley","Ochiai",
#                  "Pearson","Phi","Phi-squared","Russel","simple matching","Simpson","Stiles","Tanimoto",
#                  "Tschuprow","Yule","Yule2")
# 
# dist_measures = c("Bhjattacharyya","Bray","Canberra","Chord","divergence","Euclidean",
#                   "Geodesic","Hellinger","Kullback","Manhattan",
#                   "Podani","Soergel","supremum","Wave","Whittaker")
# ## Combine
# measures <- sort(c(sim_measures,dist_measures))


measures = "Podani"

## Fitting function #######
fitting_function <- function(row_i)
{
  
  sim_df_output <- list()
  
  # row_i = 2
  
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
    
    ## List for vectors
    exp_NA_remove_alpha <- list()
    
    ## Alphas
    for(alpha_i in alpha){
      
      ## Remove NAS, then create exponential representation 
      exp_NA_remove_alpha[[alpha_i]] <- int_sub_inds %>% 
        
        select(where(~!any(is.na(.)))) %>% 
        
        ## Change to positive representation
        mutate(across(starts_with("gp_"), ~case_when(!is.na(.x) ~ (.x - -117),
                                                     TRUE ~ 0))) %>% 
        
        ## Change to exponeitial representation
        mutate(across(starts_with("gp_"), ~(exp(.x/alpha_i))/(exp(--117/alpha_i))))
    }
    
    ## Names
    exp_NA_remove_alpha_names <- paste0("vec_list_exp_alpha_1", alpha)
    
    
    ## Create exponential representation of full vector with only top X nodes
    vec_exp_n_list[[nn]] <- exp_NA_remove
    
    
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
    
    ## Names
    vec_exp_n_list_names <- paste0("vec_list_exp_n_", 1:nodes)
    
    ## All names
    vec_names <- list(exp_NA_remove_alpha_names,vec_exp_n_list_names)
    
    
    ## Minimum vector length will be vec_NA_remove, if at least 1 column:
    if((ncol(exp_NA_remove)-8)>=1){
      
      ## List of vectors
      vec_list <- c(exp_NA_remove_alpha,vec_exp_n_list)
      
      ## Matrices
      vec_list_mat <- lapply(vec_list,
                             function(x) x %>% 
                               select(-int,-tag,-ind_id,-date_time_r,-point,-place,-type,-row) %>% 
                               as.matrix() 
      )
      
      num_gp = unlist(lapply(vec_list_mat,
                             ncol))
      
      
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
        mutate(vec_type = c(unlist(vec_names)),
               num_gp = num_gp)
      
      
      names(scores_df) <- c(measures, "vec_type","num_gp")
      
      sim_df_output[[row_i]] <- scores_df %>% 
        mutate(row = row_i,
               alpha = alpha_i,
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

saveRDS(alld, "./data/fingerprint_similarity/calibration/lf/alpha.RDS")


