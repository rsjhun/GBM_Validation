# 1. Install relevant packages -------------------------------------------------
list_of_packages <- c(
   "tidyverse", "survival", "Hmisc", "caret", "ggsci"
)

new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if (length(new_packages)) install.packages(new_packages)
lapply(list_of_packages, require, character.only = TRUE)

# 2. Import AFT model ----------------------------------------------------------
model_aft <- readRDS("C:/Users/Ray Jhun/Documents/CNOC_Projects/GBM Validation/model_aft.rds")

# 3. Import example test set ---------------------------------------------------
library(readxl)
df.test <- read_excel("C:/Users/Ray Jhun/Documents/CNOC_Projects/GBM Validation/TrialTEST4.xlsx")
df_testorig<-readRDS("C:/Users/Ray Jhun/Documents/CNOC_Projects/GBM Validation/GBMDATA.rds")


df.test$Sex<-factor(df.test$Sex
#                       ,levels=c(1,2),labels="Female","Male"
                    )
df.test$Race<-factor(df.test$Race
#                       ,levels=c(1,2,3,4),labels="White","Black","Asian","Other"
                     )
df.test$Hispanic<-as.factor(df.test$Hispanic
#                            ,levels=c(1,2),labels="No","Yes"
                            )
df.test$Married<-as.factor(df.test$Married
#                           ,levels=c(1,2),labels="No","Yes"
                           )
df.test$Insurance<-as.factor(df.test$Insurance
#                             ,levels=c(1,2),labels="Insured","Uninsured/Medicaid"
                             )
df.test$Laterality<-as.factor(df.test$Laterality
#                              ,levels=c(1,2,3),labels="Left","Right", "Midline"
                              )
df.test$Location<-as.factor(df.test$Location
#                            ,levels=c(1:8),labels="Frontal Lobe", "Temporal Lobe", "Parietal Lobe", "Occipital Lobe","Ventricle, NOS", "Cerebellum, NOS", "Brain Stem", "Overlapping Region of Brain"
                            )
df.test$TumorExtension<-as.factor(df.test$TumorExtension
#                                  ,levels=c(1,2,3),labels="Confined to Primary Location","Ventricles", "Midline Crossing"
                                  )

df.test$EOR<-as.factor(df.test$EOR
#                       ,levels=c(1:3),labels="Biopsy","Sub-Total Resection", "Gross-Total Resection"
                       )
df.test$Radiotherapy<-as.factor(df.test$Radiotherapy
#                                ,levels=c(1,2),labels="Yes","No"
                                )
df.test$Chemotherapy<-as.factor(df.test$Chemotherapy
#                                ,levels=c(1,2),labels="Yes","No"
                                )

str(df.test)
  
# 4. Compute binary survival probability and individualized kaplan meier curves-
pred_obj <- predict(model_aft, newdata = df.test, se.fit = T, type = "quantile", p = seq(0.001, 0.999, by = .001)) 
df.pred <- pred_obj$fit %>% t() %>% as_data_frame()

convert_to_prob_per_time <- function(column) {
   find_closest_to_n <- function(y,n) {
   above <- max(which(y < n))
   above <- ifelse(above == -Inf, 1,above)
   below <- min(which(y > n))
   below <- ifelse(below == Inf, 999,below)
   diff_above <- n - (y[above])
   diff_below <- (y[below]) - n
   diff_total <- abs(diff_above) + abs(diff_below)
   rel_diff_above <- abs(diff_above)/diff_total
   rel_diff_below <- abs(diff_below)/diff_total
   weight_above <- (1 - rel_diff_above)
   weight_below <- (1 - rel_diff_below)
   inv.prob <- (weight_above*above + weight_below*below)/1000
   1 - inv.prob
   }
   map_dbl(c(1:120), ~find_closest_to_n(column, .))}

df.curves <- map_df(df.pred, ~convert_to_prob_per_time(.)) %>% rbind(rep(1, ncol(.)),.)

one_year_survival_prob <- df.curves[12,] %>% as.vector() %>% as.numeric()

df.test <- df.test %>% mutate(one_year_survival_prob = one_year_survival_prob)

# 5. Compute performance metrics -----------------------------------------------

#5A: C-index
rcorr.object <- rcorr.cens(one_year_survival_prob, Surv(df.test$Observation, df.test$Dead))
rcorr.object

df.cal <- df.test %>% dplyr::select(one_year_survival_prob, one_year_survival) %>%
   mutate_at(vars(one_year_survival), funs(factor(., levels = c("Yes", "No")))) %>% na.omit()

cal_obj <- calibration(one_year_survival ~ one_year_survival_prob,
   data = df.cal, cuts = 10); ggplot(cal_obj)
plot(cal_obj, type = "l", auto.key = list(columns = 3, lines = TRUE, points = FALSE))

#5B: Calibration plot
col4 = pal_nejm(palette = c("default"), alpha = 0.8)(4)

(calibration_plot <- ggplot(cal_obj) + 
      geom_line(aes(linetype = calibModelVar), alpha = 1) + theme_bw() +
   theme(legend.position = c(0.75,0.1),
         axis.title = element_text(size = 20),
         axis.text = element_text(size = 18),
         legend.title = element_text(size = 22),
         legend.text = element_text(size = 18),
      legend.box.background = element_rect(colour = "black"),
      plot.margin = margin(1, 1, 0.7, 0.7, "cm")) +
   scale_linetype_manual(name = "Model",
      labels = c("AFT model"),
      values = c(1)) +
   scale_color_manual(name = "Model",
      labels = c("AFT model"),
      values = col4[1]) +
   labs(x = "Predicted Probability (%)",y = "Observed Event Rate (%)") +
   scale_y_continuous(expand = c(0.05,0.05)) +
   scale_x_continuous(expand = c(0.05,0.05)))

