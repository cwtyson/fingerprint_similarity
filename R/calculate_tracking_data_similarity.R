## Compare fingerprints in parallel

## Housekeeping
library(tidyverse)
library(rstatix)
library(text2vec)
library(foreach)


cores = parallel::detectCores()
cl <- parallel::makeForkCluster(cores-3, outfile = "")
doParallel::registerDoParallel(cl)

## Read in track list without outliers
finger_files_to_do = list.files("./data/fingerprint_similarity/days/")

## Done files
done <- list.files("./data/fingerprint_similarity/similarity_estimates/")

finger_files <- list.files("./data/fingerprint_similarity/days/",full.names = T)[!(finger_files_to_do %in% done)]

## Fitting function
fitting_function <- function(f)
{
  
  ## Read in file
  dets_int_sum_f <- readRDS(finger_files[f]) 
  
  day = unlist(strsplit(finger_files[f],"/"))[6]
  
  ## Set progress bar
  pb <- txtProgressBar(min = 0, max = length(unique(dets_int_sum_f$int)), style = 3)
  
  ## For each interval
  sim_df <- data.frame()
  for(i in unique(dets_int_sum_f$int)){
    
    ## Progress bar
    Sys.sleep(0.1)
    setTxtProgressBar(pb, which(i == unique(dets_int_sum_f$int)))
    
    i# i = unique(dets_int_sum_f$int)[1]
    
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
    
    for(row in 1:nrow(combos)){
      
      # row = 1
      
      int_sub_inds <- int_sub %>% 
        filter(ind_id %in% collapse(combos[row,])) %>% 
        select(where(~!any(is.na(.)))) %>% 
        ungroup()
      
      ## Calculate similarity if at least 4 grid points detecting
      num_gp <- ncol(int_sub_inds)-5
      
      if(num_gp > 3){
        
        sim_score = int_sub_inds %>% 
          select(-int,-day,-bird_band,-ind_id,-date_time_r) %>% 
          as.matrix() %>% 
          sim2(.,method = "cosine") %>% 
          data.frame()
        
        sim_df <- rbind(sim_df,
                        data.frame(int = i,
                                   ind = int_sub_inds$bird_band[1],
                                   partner = int_sub_inds$bird_band[2],
                                   score = sim_score[1,2],
                                   num_gp = num_gp))
        
      }
      
    }
    
  }
  close(pb)
  
  if(nrow(sim_df)>0){
    
    sim_df_dt <- sim_df %>% 
      left_join(dets_int_sum_f %>% 
                  distinct(int,date_time_r)) %>% 
      select(int,date_time_r,everything())
    
    saveRDS(sim_df_dt,
            paste0("./data/fingerprint_similarity/similarity_estimates/",
                   day))
    
    cat("Saved file")
    
  }
}

foreach(f=1:length(finger_files),.packages='tidyverse',
        .verbose = TRUE) %dopar%
  { fitting_function(f) }


## Combine estimates  
ests <- list.files("./data/2024/fingerprint_sim/daily_similiarity_ests/",full.names = T)


# ## Combine date time info
# 
# 
# ## 
# 
# ## Add group information
# sim_df_j <- sim_df_dt %>%
#   left_join(tag_log %>%
#               select(ind = bird_band,
#                      ind_sex = sex,
#                      ind_section = section) %>%
#               distinct(.keep_all = T)) %>%
#   left_join(tag_log %>%
#               select(partner = bird_band,
#                      partner_sex   = sex,
#                      partner_section = section) %>%
#               distinct(.keep_all = T)
#   ) %>%
#   
#   mutate(same_group = ifelse(ind_section == partner_section, "yes","no"))
# 
# ggplot(sim_df_j %>% 
#          mutate(hour = format(date_time_r, unit = "%H"))) +
#   geom_jitter(aes(x=hour,
#                   y=score)) +
#   facet_grid(~same_group)
# 
# 
# 



# 
# dets_int_sum_ref_na <- dets_int_sum%>% 
#   
#   ## Replace NA values with random value
#   mutate(across(everything(), 
#                 .fns = ~ifelse(!is.na(.),
#                                ., 
#                                sample(-1000:1000)))) %>% 
#   arrange(int,bird_band)
# 
# sim_df_f <- dets_int_sum %>% 
#   select(-date_time_r) %>% 
#   ## Nest by time interval
#   nest_by(int) %>% 
#   mutate(mat = list(data %>% 
#                       select(-bird_band) %>% 
#                       as.matrix()),
#          bird_band = list(data %>% pull(bird_band)),
#          sim = list(data.frame(sim2(mat,method = "cosine")) %>% 
#                       rename_with(~bird_band)%>% 
#                       mutate(bird_band = bird_band) %>% 
#                       column_to_rownames("bird_band") %>% 
#                       rstatix::cor_gather() %>% 
#                       filter(var1 != var2) %>% 
#                       rename(ind = var1,
#                              partner = var2))) %>% 
#   select(int,
#          sim) %>% 
#   reframe(sim)
# 

# ggplot(sim_df_j) +
#   geom_(aes(x=same_group,
#             y=score)
#   )
# 
# 
# 
# dets <- readRDS("./data/2024/fingerprint_sim/raw_dets.RDS")
# 
# ## Get fingerprint each 15 seoncds
# dets_int <- dets_f3 %>%
#   
#   ## Keep tagged birds
#   filter(tag %in% tag_log$tag) %>% 
#   arrange(date_time) %>% 
#   # slice(1:100000) %>%
#   
#   ## 30 second fingerprint
#   mutate(date_time_r = floor_date(date_time,unit="30 seconds")) %>% 
#   group_by(bird_band,date_time_r,grid_point) %>% 
#   summarise(mean_rssi = mean(rssi)) %>% 
#   group_by(bird_band,
#            date_time_r) %>%
#   arrange(bird_band,
#           date_time_r,
#           desc(mean_rssi)) %>% 
#   ungroup()
# 
# ## Process to keep top 5 nodes in each interval
# dets_int_sum <- dets_int %>% 
#   
#   ## Assign interval values
#   group_by(date_time_r) %>% 
#   mutate(int = cur_group_id(),
#          inds = n()
#   ) %>% 
#   filter(inds > 1) %>% 
#   select(-inds) %>%
#   
#   ## Keep sample data
#   # filter(int < 10000) %>% 
#   arrange(int)  %>% 
#   ungroup(date_time_r) %>% 
#   select(int,date_time_r,bird_band,grid_point,mean_rssi) %>% 
#   group_by(int,
#            date_time_r,
#            bird_band) %>% 
#   pivot_wider(names_from=grid_point,
#               values_from=mean_rssi) %>% 
#   group_by(bird_band) %>% 
#   mutate(ind_id = cur_group_id()) %>% 
#   select(int,
#          bird_band,
#          ind_id,
#          everything()) %>% 
#   ungroup()
# 
# 
# ## Set progress bar
# pb <- txtProgressBar(min = 0, max = length(unique(dets_int_sum$int)), style = 3)
# 
# ## For each interval
# sim_df <- data.frame()
# for(i in unique(dets_int_sum$int)){
#   
#   ## Progress bar
#   Sys.sleep(0.1)
#   setTxtProgressBar(pb, which(i == 1:nrow(combos)))
#   
#   # i = 1
#   
#   int_sub <- dets_int_sum[dets_int_sum$int==i,]
#   
#   ## Get overlap for each 
#   combos <- expand.grid(x = unique(int_sub$ind_id),
#                         y = unique(int_sub$ind_id)) %>% 
#     filter(x != y) %>% 
#     rowwise() %>% 
#     mutate(order = paste(min(x,y),max(x,y))) %>% 
#     distinct(order,.keep_all = T) %>% 
#     arrange(x,y) %>% 
#     select(x,y)
#   
#   for(row in 1:nrow(combos)){
#     
#     # row = 1122
#     
#     
#     int_sub_inds <- int_sub %>% 
#       filter(ind_id %in% collapse(combos[row,])) %>% 
#       select(where(~!any(is.na(.)))) %>% 
#       ungroup()
#     
#     ## Calculate similarity if at least 4 grid points detecting
#     num_gp <- ncol(int_sub_inds)-4
#     
#     if(num_gp > 3){
#       
#       sim_score = int_sub_inds %>% 
#         select(-int,-bird_band,-ind_id,-date_time_r) %>% 
#         as.matrix() %>% 
#         sim2(.,method = "cosine") %>% 
#         data.frame()
#       
#       sim_df <- rbind(sim_df,
#                       data.frame(int = i,
#                                  ind = int_sub_inds$bird_band[1],
#                                  partner = int_sub_inds$bird_band[2],
#                                  score = sim_score[1,2],
#                                  num_gp = num_gp))
#       
#     }
#     
#   }
#   
# }
# close(pb)
# 
# 
# 
# 
# 
# 
# 
# dets_int_sum_ref_na <- dets_int_sum%>% 
#   
#   ## Replace NA values with random value
#   mutate(across(everything(), 
#                 .fns = ~ifelse(!is.na(.),
#                                ., 
#                                sample(-1000:1000)))) %>% 
#   arrange(int,bird_band)
# 
# sim_df_f <- dets_int_sum %>% 
#   select(-date_time_r) %>% 
#   ## Nest by time interval
#   nest_by(int) %>% 
#   mutate(mat = list(data %>% 
#                       select(-bird_band) %>% 
#                       as.matrix()),
#          bird_band = list(data %>% pull(bird_band)),
#          sim = list(data.frame(sim2(mat,method = "cosine")) %>% 
#                       rename_with(~bird_band)%>% 
#                       mutate(bird_band = bird_band) %>% 
#                       column_to_rownames("bird_band") %>% 
#                       rstatix::cor_gather() %>% 
#                       filter(var1 != var2) %>% 
#                       rename(ind = var1,
#                              partner = var2))) %>% 
#   select(int,
#          sim) %>% 
#   reframe(sim)
# 
# ## Add group information
# sim_df_j <- sim_df %>% 
#   left_join(tag_log %>% 
#               select(ind = bird_band,
#                      ind_sex = sex,
#                      ind_section = section) %>% 
#               distinct(.keep_all = T)) %>% 
#   left_join(tag_log %>% 
#               select(partner = bird_band,
#                      partner_sex   = sex,
#                      partner_section = section) %>% 
#               distinct(.keep_all = T)
#   ) %>% 
#   
#   mutate(same_group = ifelse(ind_section == partner_section, "yes","no"))
# 
# ggplot(sim_df_j) +
#   geom_(aes(x=same_group,
#             y=score)
#   )
