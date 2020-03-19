## Exploratory plots ## ------------------------------------------------------------------------------------------------
dfm <- melt(data_all, id.vars = c("ind","year"))

p_RCRI <- ggplot(dfm %>% filter(variable == "RCRI_pt"), aes(x = value )) + 
  geom_bar(aes(fill = year), stat = "count", position = "dodge", alpha = .5) +
  xlab("RCRI") + ylab("") +
  ggtitle("RCRI distribution by year") + 
  plot_themes# + 
# facet_grid(year~., scales = "free_y")

if(save_plots){
  save_plot("RCRI_distribution.pdf", p_RCRI, base_width = 8, base_height = 4)
}

M <- cor(cor_data)
diag(M) <- 0
p_cor <- corrplot(M,
                  method = "shade",
                  type = "lower", 
                  order = "hclust",
                  tl.col = "black",
                  tl.cex = .75)

cor_data_acute <- data %>% select(c(age_factor, 
                                    NSTEMI,
                                    severe_MI,
                                    Bleed, 
                                    PE,
                                    DVT,
                                    Afib,
                                    PNA,
                                    sepsis)) %>% 
  mutate(STEMI = 1- as.numeric(NSTEMI)) %>% 
  apply(.,2,as.numeric) 

M_acute <- cor(cor_data_acute)
diag(M_acute) <- 0
p_cor_acute <- corrplot(M_acute,
                        method = "shade",
                        type = "lower",
                        diag=T,
                        order = "hclust",
                        tl.col = "black",
                        tl.cex = .75)

