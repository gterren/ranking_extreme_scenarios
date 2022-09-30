library(fdapace)
library(refund)
library(latex2exp)
library(QuantPsyc)
library(ks)
library(refund)
library(mvtnorm)
library(ddalpha)
library(fdaoutlier)

load_files = function(path, date, source) {
  main_dir   = paste(path, date, source, sep = '/')
  dir_list   = list.dirs(main_dir, recursive = FALSE) 
  list_files = list.files(path = main_dir,
                          pattern = "*.csv", 
                          full.names = T) 
  N = length(list_files)
  y_ac_ = list("list", N)
  y_fc_ = list("list", N)
  Y_sc_ = list("list", N)
  name_ = list("list", N)
  date = sub(date, pattern = "SimDat_", replacement = "")
  for (i in 1:N) {
    data_ = read.csv(list_files[i], sep = ',', header = FALSE)
    name_[[i]] = sub(sub(".*/", "", list_files[i]), pattern = ".csv$", replacement = "")
    #name_[[i]] = gsub(name, pattern = "_", replacement = " ")
    y_ac_[[i]] = matrix(as.numeric(unlist(data_[2, -c(1:2)])), ncol = 1, nrow = 24)
    y_fc_[[i]] = matrix(as.numeric(unlist(data_[3, -c(1:2)])), ncol = 1, nrow = 24)
    Y_sc_[[i]] = matrix(as.numeric(unlist(data_[-c(1:3),-c(1:2)])), ncol = 24, nrow = 1000)
  }
  return(list(y_ac_, y_fc_, Y_sc_, name_, date))
}

# SimDat_20180613
# SimDat_20180911
# SimDat_20181126

# SimDat_20180702
# SimDat_20180722
# SimDat_20180718

simulation = 'SimDat_20180722'
source     = 'Wind'

path_to_data   = '/Users/Guille/Desktop/extreme_scenarios/data'
path_to_images = '/Users/Guille/Desktop/extreme_scenarios/images'

X_    = load_files(path_to_data, simulation, source)
y_ac_ = X_[[1]]
y_fc_ = X_[[2]]
Y_sc_ = X_[[3]]
name_ = X_[[4]]
date  = X_[[5]]
print(name_)
print(date)

N_assets   = length(Y_sc_) 
N_features = dim(Y_sc_[[1]])[2]
N_samples  = dim(Y_sc_[[1]])[1]
print(N_assets)
print(N_features)
print(N_samples)

i_asset = 27
Y_sc_asset_ = Y_sc_[[i_asset]]
y_ac_asset_ = y_ac_[[i_asset]]
y_fc_asset_ = y_fc_[[i_asset]]
name_[i_asset]

#depth_ = modified_band_depth(Y_sc_asset_)
depth_ = directional_quantile(Y_sc_asset_)
#depth_ = extreme_rank_length(Y_sc_asset_)

dim(Y_sc_asset_)
Col_ = viridis::viridis(1000)

ranking_ = order(depth_)
print(ranking_)

matplot(t(Y_sc_asset_[ranking_,]), type = 'l', col = Col_, lwd = 1, lty = 1)

# Define the continuum
t_   = seq(1, N_features, length.out = N_features)
T_   = rep(t_, N_samples)
IDs_ = rep(1:N_samples, each = N_features)
length(t_)
length(T_)
length(IDs_)

# Compute Function PCA
fun  = MakeFPCAInputs(IDs = IDs_, tVec = T_, t(Y_sc_asset_))
fPCA = FPCA(fun$Ly, fun$Lt, optns = list(maxK = 5))

mu_  = fPCA$mu
Xi_  = fPCA$xiEst
Phi_ = fPCA$phi

matplot(t(Y_sc_asset_), xlab = 'Time [hh]', 
        ylab = 'Energy [kWh]',
        ylim = range(min(Y_sc_asset_), max(Y_sc_asset_)),
        type = 'l', col = 'gray75', lwd = 0.125, lty = 1)
par(new=TRUE)
matplot(y_fc_asset_, xlab = 'Time [hh]', 
        ylab = 'Energy [kWh]',
        ylim = range(min(Y_sc_asset_), max(Y_sc_asset_)),
        type = 'b', col = 'blue', lwd = 1., lty = 1, pch = 1)
par(new=TRUE)
matplot(y_ac_asset_, xlab = 'Time [hh]', 
        ylab = 'Energy [kWh]',
        ylim = range(min(Y_sc_asset_), max(Y_sc_asset_)),
        type = 'b', col = 'red', lwd = 1., lty = 1, pch = 2)
par(new=TRUE)
matplot(mu_, xlab = 'Time [hh]', 
        ylab = 'Energy [kWh]',
        main = paste(date, paste(source, name_[i_asset], sep = ' | '), sep = ' | '),
        ylim = range(min(Y_sc_asset_), max(Y_sc_asset_)),
        type = 'b', col = 'green', lwd = 1., lty = 1, pch = 3)
grid(lty = 2,      
     col = 'black', 
     lwd = 0.5) 
legend(1, max(Y_sc_asset_), legend = c(TeX('$f_i^{(sc)}$'), TeX('$f^{(ac)}$'), TeX('$f^{(fc)}$'), TeX('$mu(t)$')), 
       col = c('gray75', 'red', 'blue', 'green'), lty = 1, cex = 0.75)

col_ = c('red', 'blue', 'green', 'yellow', 'magenta')

matplot(Phi_, xlab = 'Time [hh]', 
        ylab = 'Score', 
        ylim = range(c(min(Phi_), max(Phi_))),
        main = paste(date, paste(source, name_[i_asset], sep = ' | '), sep = ' | '),
        type = 'l', col = col_, lwd = 1, lty = 1)
grid(lty = 2,      
     col = 'black', 
     lwd = 0.5)   
legend(1, max(Phi_), legend = c(TeX('$phi_1$'), TeX('$phi_2$'), TeX('$phi_3$'), TeX('$phi_4$'), TeX('$phi_5$')), 
       col = col_, lty = 1, cex = 0.75)

matplot(t(Xi_), xlab = 'PC Factor', 
        ylab = 'Score',
        ylim = range(min(Xi_), max(Xi_)),
        type = 'l', col = 'gray75', lwd = 0.125, lty = 1)
grid(lty = 2,      
     col = 'black', 
     lwd = 0.5)   

N_scenarios = 100
Z_ = array(NA, dim = c(N_scenarios, N_features))
dim(Z_)
length(Z_[[i,]])
length( f(Y_sc_asset_[i,]))
for (i in 1:N_scenarios) {
  f = ecdf(Y_sc_asset_[i,])
  Z_[i,] = f(Y_sc_asset_[i,])
}

d_ = colSum(abs(.5 - Z_))
d_ = 1. - colSum(.5 - F_), axis = 0)/24.




plot(z_)

datafA$vals = t(Y_sc_asset_)
datafA$args = c(1:24)
dim(datafA$vals)
length(datafA$args)

depthf.BD(datafA, datafA)


z_1_ = dmvnorm(Xi_, colMeans(Xi_), cov(Xi_), FALSE)

Z_1_ = array(NA, dim = c(1000, 5))
for (i in 1:5) {
  Z_1_[,i] = dnorm(Xi_[,i], mean(Xi_[,i]), sd(Xi_[,i]), FALSE)
}

Z_2_ = array(NA, dim = c(1000, 24))
for (i in 1:24) {
  Z_2_[,i] = dnorm(Y_[,i], mean(Y_[,i]), sd(Y_[,i]), FALSE)
}

z_1_ = colprods(t(Z_1_))
length(z_1_)
z_2_ = colprods(t(Z_2_))
length(z_2_)

xiCol_ = viridis::viridis(1000)

idx_1_ = order(z_1_)
idx_2_ = order(z_2_)

matplot(t(Y_[idx_1_,]), type = 'l',col = xiCol_, lwd = 1, lty = 1)
matplot(t(Y_[idx_2_,]), type = 'l',col = xiCol_, lwd = 1, lty = 1)

matplot(t(Xi_[idx_1_,]), type = 'l',col = xiCol_, lwd = 1, lty = 1)
matplot(t(Xi_[idx_2_,]), type = 'l',col = xiCol_, lwd = 1, lty = 1)

idx_1_ = order(Xi_[,1])
idx_2_ = order(Xi_[,2])
idx_3_ = order(Xi_[,3])
idx_4_ = order(Xi_[,4])
idx_5_ = order(Xi_[,5])

plot(Xi_[idx_1_,3], Xi_[idx_1_,2], col = xiCol_, lwd = 1, lty = 1)

matplot(t(Y_[idx_1_,]), type = 'l',col = xiCol_, lwd = 1, lty = 1)
matplot(t(Y_[idx_2_,]), type = 'l',col = xiCol_, lwd = 1, lty = 1)
matplot(t(Y_[idx_3_,]), type = 'l',col = xiCol_, lwd = 1, lty = 1)
matplot(t(Y_[idx_4_,]), type = 'l',col = xiCol_, lwd = 1, lty = 1)
matplot(t(Y_[idx_5_,]), type = 'l',col = xiCol_, lwd = 1, lty = 1)

matplot(t(Xi_[idx_1_,]), type = 'l',col = xiCol_, lwd = 1, lty = 1)
matplot(t(Xi_[idx_2_,]), type = 'l',col = xiCol_, lwd = 1, lty = 1)
matplot(t(Xi_[idx_3_,]), type = 'l',col = xiCol_, lwd = 1, lty = 1)
matplot(t(Xi_[idx_4_,]), type = 'l',col = xiCol_, lwd = 1, lty = 1)
matplot(t(Xi_[idx_5_,]), type = 'l',col = xiCol_, lwd = 1, lty = 1)


matplot(phis_[,2,], xlab = 'Time [hh]',
        ylab = TeX('$phi_{2i}(t)$'), 
        main = paste(date, paste(source, 'PC2', sep = ': '), sep = ' - '),
        type = 'l', col = 'gray75', lwd = 0.5, lty = 1)



col_ = c('red', 'blue', 'green', 'yellow', 'magenta')
par(mfrow = c(1,1), mar = c(4,5,2,1), family = 'serif')
matplot(fPCA$phi, xlab = 'Time [hh]', 
        ylab = 'Score', 
        ylim = range(c(min(fPCA$phi), max(fPCA$phi))),
        main = paste(date, paste(source, name_[i], sep = ' | '), sep = ' | '),
        type = 'l', col = col_, lwd = 1, lty = 1)


matplot(t(Y_sc_prime_), xlab = 'Time [hh]', 
        ylab = 'Energy [kWh]',
        ylim = range(min(Y_sc_prime_), max(Y_sc_prime_)),
        type = 'l', col = 'gray75', lwd = 0.125, lty = 1)
par(new=TRUE)
matplot(y_fc_prime_, xlab = 'Time [hh]', 
        ylab = 'Energy [kWh]',
        ylim = range(min(Y_sc_prime_), max(Y_sc_prime_)),
        type = 'b', col = 'blue', lwd = 1., lty = 1, pch = 1)
par(new=TRUE)
matplot(y_ac_prime_, xlab = 'Time [hh]', 
        ylab = 'Energy [kWh]',
        ylim = range(min(Y_sc_prime_), max(Y_sc_prime_)),
        type = 'b', col = 'red', lwd = 1., lty = 1, pch = 2)