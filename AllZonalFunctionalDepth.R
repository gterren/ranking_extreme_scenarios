library(reticulate)
library(fdaoutlier)
#library(ddalpha)
library(fda.usc)
library(abind)

pickle = import("pickle")

X_ = py_load_object("/Users/Guille/Dropbox/ProcessedDataTexas/vatic_output/Texas-7k/FunctionalDepthsTexas/ProcessedZonalTexas.pkl", pickle = "pickle")
L_ = X_[[1]]
S_ = X_[[2]]
W_ = X_[[3]]
G_ = X_[[4]]
N_ = X_[[5]]
X_ = list(L_, S_, W_, G_, N_)

i_feature = 3
i_zone    = 3
i_day     = 1

x_ = t(X_[[i_feature]][,, i_zone, i_day])
dim(x_)

matplot(t(x_), type="l", col = 1)

days_ = c(25, 1, 3, 6)

#N_days      = dim(L_)[4]
N_days      = length(days_)
N_zones     = dim(L_)[3]
N_features  = length(X_)
N_metrics   = 10
N_scenarios = 1000


FD_ = array(rep(1, N_metrics*N_days*N_zones*N_scenarios), 
            dim = c(N_metrics, N_zones, N_days*N_scenarios))

  # Loop over grid zones
for (i_zone in 1:N_zones){
  print(i_zone)
  x_ = t(X_[[i_feature]][,, i_zone, days_[1]])
  # Loop over available days
  for (i_day in 2:N_days){
    y_ = t(X_[[i_feature]][,, i_zone, days_[i_day]])
    dim(y_)
    # Day data
    x_ = abind(x_, y_, along = 1)
  }
  print(dim(x_))
  # Compute functional depth metrics
  FD_[1, i_zone,]  = band_depth(x_)
  FD_[2, i_zone,]  = modified_band_depth(x_)
  FD_[3, i_zone,]  = directional_quantile(x_, quantiles = c(0.025, 0.975))
  FD_[3, i_zone,]  = max(FD_[3, i_zone,]) - FD_[3, i_zone,] + min(FD_[3, i_zone,]) 
  FD_[4, i_zone,]  = extremal_depth(x_)
  FD_[5, i_zone,]  = linfinity_depth(x_)
  FD_[6, i_zone,]  = extreme_rank_length(x_, type = "two_sided")
  FD_[7, i_zone,]  = depth.mode(fdata(x_))$dep
  FD_[8, i_zone,]  = depth.RT(fdata(x_))$dep
  FD_[9, i_zone,]  = depth.FM(fdata(x_))$dep
  FD_[10, i_zone,] = depth.RP(fdata(x_))$dep
}

py_save_object(FD_, "/Users/Guille/Dropbox/ProcessedDataTexas/vatic_output/Texas-7k/FunctionalDepthsTexas/AllAreaZonalFuntionalDepths.pkl", pickle = "pickle")
