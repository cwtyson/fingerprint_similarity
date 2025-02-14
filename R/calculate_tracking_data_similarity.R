## Compare fingerprints in parallel

## Housekeeping
library(tidyverse)
library(rstatix)
library(text2vec)
library(foreach)
library(proxy)


cores = parallel::detectCores()
cl <- parallel::makeForkCluster(cores-1, outfile = "")
doParallel::registerDoParallel(cl)

## Read in prepared fingerprinting files
finger_files = list.files("./data/fingerprint_similarity/days/")

## Done files
done <- list.files("./data/fingerprint_similarity/similarity_estimates/")

## Get files to do
finger_files_to_do <- list.files("./data/fingerprint_similarity/days/",full.names = T)[!(finger_files %in% done)]

## Similarity and distance measures - 
sim_measures = c("Cramer")

## Combine
measures <- sim_measures 

## Fitting function
fitting_function <- function(f)
{
  
  
  ## Read in file
  dets_int_sum <- readRDS(finger_files_to_do[f])%>% 
    group_by(int) %>% 
    mutate(inds = n()) %>% 
    select(inds,everything()) %>% 
    filter(inds>1)
  
  
  ## Get ints that are not all NA
  ints_f <- dets_int_sum %>% 
    group_by(int) %>%
    summarise_all(list(~sum(!is.na(.))))  %>% 
    group_by(int) %>% 
    mutate(max = max(c_across(contains("gp_")), na.rm = T)) %>% 
    filter(max > 1) %>% 
    pull(int)  
  
  ## Remove groups with all NA
  dets_int_sum_f <- dets_int_sum %>% 
    filter(int %in% ints_f )

  ## Day to process
  day = unlist(strsplit(finger_files_to_do[f],"/"))[6]
  
  cat("Day\n", day)
  
  ## Set progress bar
  pb <- txtProgressBar(min = 0, max = length(unique(dets_int_sum_f$int)), style = 3)
  cat("\n")
  
  ## For each interval
  sim_df_day <- data.frame()
  for(i in unique(dets_int_sum_f$int)){
    
    ## Progress bar
    Sys.sleep(0.1)
    setTxtProgressBar(pb, which(i == unique(dets_int_sum_f$int)))
    
    # i = unique(dets_int_sum_f$int)[1]
    
    int_sub <- dets_int_sum_f[dets_int_sum_f$int==i,] 
    
    ## Get overlap for each 
    combos <- expand.grid(x = unique(int_sub$ind_id),
                          y = unique(int_sub$ind_id)) %>% 
      filter(x != y) %>% 
      rowwise() %>% 
      mutate(order = paste(min(x,y),max(x,y))) %>% 
      distinct(order,.keep_all = T) %>% 
      arrange(x,y) %>% 
      select(x,y)
    
    ## For each row
    for(row_i in 1:nrow(combos)){
      
      # row_i = 3
      
      ind_ids <- c(combos[row_i,]$x,combos[row_i,]$y)
      
      ## Keep only matching
      int_sub_inds <- int_sub %>% 
        filter(ind_id %in% ind_ids) %>% 
        select(where(~!any(is.na(.)))) %>% 
        ungroup()
      
      ## Calculate similarity if at least 3 grid points detecting
      num_gp <- ncol(int_sub_inds)-6
      
      if(num_gp >= 1){
        
        ## Create exponential representaiton of full vector
        vec_exp <- int_sub_inds %>% 
          
          ## Change to positive representation
          mutate(across(starts_with("gp_"), ~case_when(!is.na(.x) ~ (.x - -117),
                                                       TRUE ~ 0))) %>% 
          
          ## Change to exponeitial representation
          mutate(across(starts_with("gp_"), ~(exp(.x/24))/(exp(--117/24))))
        
        
        ## Change to matrix and calculate similarity score
        mat <- vec_exp%>% 
          select(-inds,-day,-int,-date_time_r,-bird_band,-ind_id) %>% 
          as.matrix()
        
        sim <- proxy::dist(mat,
                           method = "Cramer")
        
        sim_df <- data.frame(int = int_sub_inds$int[1],
                             row = row_i,
                             score = sim[1],
                             ind = int_sub_inds$bird_band[1],
                             partner = int_sub_inds$bird_band[2],
                             row_ind = int_sub_inds$ind_id[1],
                             row_partner = int_sub_inds$ind_id[2]) 
        
        sim_df_day <- rbind(sim_df_day,sim_df)
        
      } 
    }
    cat("\nFinished interval\n")
    
  }
  close(pb)
  
  if(nrow(sim_df_day)>0){
    
    sim_df_day_j <- sim_df_day %>% 
      left_join(dets_int_sum_f %>% 
                  distinct(int,date_time_r),
                by = join_by(int)) %>% 
      select(int,date_time_r,everything())
    
    saveRDS(sim_df_day_j,
            paste0("./data/fingerprint_similarity/similarity_estimates/",
                   day))
    
    cat("Saved file")
    
  }
}

foreach(f=1:length(finger_files_to_do),.packages='tidyverse',
        .verbose = TRUE) %dopar%
  { fitting_function(f) }
