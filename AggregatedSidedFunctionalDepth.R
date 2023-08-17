library(reticulate)
library(fdaoutlier)

pickle = import("pickle")

X_ = py_load_object("/Users/Guille/Dropbox/ProcessedDataTexas/vatic_output/Texas-7k/FunctionalDepthsTexas/ProcessedAggregatedData.pkl", pickle = "pickle")
L_ = X_[[1]]
S_ = X_[[2]]
W_ = X_[[3]]
G_ = X_[[4]]
N_ = X_[[5]]
X_ = list(L_, S_, W_, G_, N_)

N_days      = dim(L_)[3]
N_features  = length(X_)
N_metrics   = 3
N_scenarios = 1000
N_hours     = 24

FD_ = array(rep(1, N_metrics*N_days*N_features*N_scenarios), 
            dim = c(N_metrics, N_features, N_days, N_scenarios))
# Loop over asset features
for (i_feature in 1:N_features){
  print(i_feature)
  # Loop over available days
  for (i_day in 1:N_days){
    # Get day data asset features
    x_ = X_[[i_feature]][,,i_day]
    x_ = t(diff(t(x_)))
    # Compute functional depth metrics
    FD_[1, i_feature, i_day,]  = extreme_rank_length(x_, type = "two_sided")
    FD_[2, i_feature, i_day,]  = extreme_rank_length(x_, type = "one_sided_left")
    FD_[3, i_feature, i_day,]  = extreme_rank_length(x_, type = "one_sided_right")
  }}

#py_save_object(FD_, "/Users/Guille/Dropbox/ProcessedDataTexas/vatic_output/Texas-7k/FunctionalDepthsTexas/AreaAggSidedFuntionalDepths.pkl", pickle = "pickle")
py_save_object(FD_, "/Users/Guille/Dropbox/ProcessedDataTexas/vatic_output/Texas-7k/FunctionalDepthsTexas/DiffAggSidedFuntionalDepths.pkl", pickle = "pickle")