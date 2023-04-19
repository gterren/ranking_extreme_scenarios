library(reticulate)
library(fdaoutlier)
#library(ddalpha)
library(fda.usc)

pickle = import("pickle")

X_ = py_load_object("/Users/Guille/Dropbox/ProcessedDataTexas/vatic_output/Texas-7k/FunctionalDepthsTexas/ProcessedZonalTexas.pkl", pickle = "pickle")
L_ = X_[[1]]
S_ = X_[[2]]
W_ = X_[[3]]
G_ = X_[[4]]
N_ = X_[[5]]
X_ = list(L_, S_, W_, G_, N_)

N_days      = dim(L_)[4]
N_zones     = dim(L_)[3]
N_features  = length(X_)
N_metrics   = 10
N_scenarios = 1000

FD_ = array(rep(1, N_metrics*N_days*N_features*N_zones*N_scenarios), 
            dim = c(N_metrics, N_features, N_zones, N_days, N_scenarios))
# Loop over asset features
for (i_feature in 1:N_features){
  print(i_feature)
  # Loop over grid zones
  for (i_zone in 1:N_zones){
    print(i_zone)
    # Loop over available days
    for (i_day in 1:N_days){
      # Day data
      x_ = t(X_[[i_feature]][,, i_zone, i_day])
      #x_ = t(diff(t(x_)))
      # Compute functional depth metrics
      FD_[1, i_feature, i_zone, i_day,]  = band_depth(x_)
      FD_[2, i_feature, i_zone, i_day,]  = modified_band_depth(x_)
      FD_[3, i_feature, i_zone, i_day,]  = directional_quantile(x_, quantiles = c(0.025, 0.975))
      FD_[3, i_feature, i_zone, i_day,]  = max(FD_[3, i_feature, i_zone, i_day,]) - FD_[3, i_feature, i_zone, i_day,] + min(FD_[3, i_feature, i_vatic, i_day,]) 
      FD_[4, i_feature, i_zone, i_day,]  = extremal_depth(x_)
      FD_[5, i_feature, i_zone, i_day,]  = linfinity_depth(x_)
      FD_[6, i_feature, i_zone, i_day,]  = extreme_rank_length(x_, type = "two_sided")
      #FD_[7, i_feature, i_day,]  = projection_depth(x_)
      FD_[7, i_feature, i_zone, i_day,]  = depth.mode(fdata(x_))$dep
      FD_[8, i_feature, i_zone, i_day,]  = depth.RT(fdata(x_))$dep
      FD_[9, i_feature, i_zone, i_day,]  = depth.FM(fdata(x_))$dep
      FD_[10, i_feature, i_zone, i_day,] = depth.RP(fdata(x_))$dep
      #FD_[12, i_feature, i_zone, i_day,] = depth.RPD(fdata(x_))$dep
      #FD_[13, i_feature, i_zone, i_day,] = depth.FSD(fdata(x_))$dep
      #FD_[14, i_feature, i_zone, i_day,] = depth.KFSD(fdata(x_))$dep
    }}}

py_save_object(FD_, "/Users/Guille/Dropbox/ProcessedDataTexas/vatic_output/Texas-7k/FunctionalDepthsTexas/AreaZonalFuntionalDepths.pkl", pickle = "pickle")
#py_save_object(FD_, "/Users/Guille/Dropbox/ProcessedDataTexas/vatic_output/Texas-7k/FunctionalDepthsTexas/DiffZonalFuntionalDepths.pkl", pickle = "pickle")