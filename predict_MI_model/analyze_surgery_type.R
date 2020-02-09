
pr_class_fun <- function(x,code_df){
  if(x %in% code_df$code){
    pr_type <- code_df %>% filter(code == x) %>% 
      select(type) %>% 
      unlist()
  }
  if(!(x %in% code_df$code)){
    pr_type = NA
  }
  return(as.character(pr_type))
}

pr_list <- list(c(1:9),c(10:12),
                c(13:21), 
                c(22:35), 
                c(36:42),
                c(43:50, 62:63),
                c(51:61),
                c(66:99), 
                c(100:118),
                c(119:125,129:132),
                c(126:128, 133:141),
                c(142:164),
                c(165:167),
                c(168:175),
                176)
names(pr_list) = c("Neuro", "Endo","Optho","ENT","Thoracic","Cardiac","Vascular","General","GU", "Gynecologic", "Obstetric", "Orthopedic","Breast","Skin/burn", "Transplant")

pr_df <- stack(pr_list) %>% 
  as.data.frame() %>% 
  stats::setNames(c("code","type"))

data$surgery_type = NA

t1 <- Sys.time()
data$surgery_type <- sapply(as.numeric(as.character(data$prccs1)), pr_class_fun, code_df = pr_df)
save(data, file = "data_raw_with_surgery_type.rda")
t2 <- Sys.time()
print(t2 - t1)

prop_df <- data.frame(type = names(sum), frac = round(as.numeric(sum)/sum(as.numeric(sum)),3))

#prop_df_by_MI <- data %>% 
  #group_by(MI) %>% 
  #summarize()
