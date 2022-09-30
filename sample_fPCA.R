library(fdapace)
library(refund)
library(latex2exp)
library(QuantPsyc)
library(ks)
library(refund)
library(mvtnorm)

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

simulation = 'SimDat_20180718'
source     = 'load'

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

i = 1
Y_sc_prime_ = Y_sc_[[i]]
y_ac_prime_ = y_ac_[[i]]
y_fc_prime_ = y_fc_[[i]]
length(Y_sc_)
dim(Y_sc_prime_)

# Define the continuum
t_   = seq(1, N_features, length.out = N_features)
T_   = rep(t_, N_samples)
IDs_ = rep(1:N_samples, each = N_features)
length(t_)
length(T_)
length(IDs_)

fun  = MakeFPCAInputs(IDs = IDs_, tVec = T_, t(Y_sc_[[i]]))
fPCA = FPCA(fun$Ly, fun$Lt)

Xi_ = fPCA$xiEst
Y_  = Y_sc_[[i]]

dim(Xi_)
dim(Y_)
dim(cov(Y_))

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

matplot(t(Y_sc_[[i]]), xlab = 'Time [hh]', 
        ylab = 'Energy [kWh]',
        ylim = range(min(Y_sc_[[i]]), max(Y_sc_[[i]])),
        type = 'l', col = 'gray75', lwd = 0.125, lty = 1)
par(new=TRUE)
matplot(y_fc_[[i]], xlab = 'Time [hh]', 
        ylab = 'Energy [kWh]',
        ylim = range(min(Y_sc_[[i]]), max(Y_sc_[[i]])),
        type = 'b', col = 'blue', lwd = 1., lty = 1, pch = 1)
par(new=TRUE)
matplot(y_ac_[[i]], xlab = 'Time [hh]', 
        ylab = 'Energy [kWh]',
        ylim = range(min(Y_sc_[[i]]), max(Y_sc_[[i]])),
        type = 'b', col = 'red', lwd = 1., lty = 1, pch = 2)
par(new=TRUE)
matplot(fPCA$mu, xlab = 'Time [hh]', 
        ylab = 'Energy [kWh]',
        main = paste(date, paste(source, name_[i], sep = ' | '), sep = ' | '),
        ylim = range(min(Y_sc_[[i]]), max(Y_sc_[[i]])),
        type = 'b', col = 'green', lwd = 1., lty = 1, pch = 3)

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



N_components_max = 5

mus_   = array(NA, dim = c(N_features, N_assets))
cums_  = array(NA, dim = c(N_components_max, N_assets))
phis_  = array(NA, dim = c(N_features, N_components_max, N_assets))
fPCA_  = list("list", N_assets)

optns = list(maxK = N_components_max)

# Create Paths for Saving Results
dir.create(paste(path_to_images, simulation, sep = '/'), showWarnings = TRUE)
path_to_save = paste(path_to_images, simulation, source, sep = '/')
print(path_to_save)
dir.create(path_to_save, showWarnings = TRUE)

for (i in 1:N_assets) {
  #
  fun  = MakeFPCAInputs(IDs = IDs_, tVec = T_, t(Y_sc_[[i]]))
  fPCA = FPCA(fun$Ly, fun$Lt, optns)
  
  # Get Results
  idx_           = c(1:length(fPCA$cumFVE))
  fPCA_[[i]]     = fPCA
  mus_[,i]       = fPCA$mu
  cums_[idx_,i]  = fPCA$cumFVE
  phis_[,idx_,i] = fPCA$phi
  
  # Saving Image 

  y_ac_prime_ = matrix(unlist(y_ac_), ncol = length(y_ac_), nrow = 24)
  y_fc_prime_ = matrix(unlist(y_fc_), ncol = length(y_fc_), nrow = 24)
  
  image_name = paste(path_to_save,  paste(name_[i], '_sc.png', sep = ''), sep = '/')
  png(file = image_name, units = 'in', width = 6, height = 5, res = 150)
  par(mfrow = c(1,1), mar = c(4,5,2,1), family = 'serif')
  matplot(t(Y_sc_[[i]]), xlab = 'Time [hh]', 
                         ylab = 'Energy [kWh]',
                         ylim = range(min(Y_sc_[[i]]), max(Y_sc_[[i]])),
                         type = 'l', col = 'gray75', lwd = 0.125, lty = 1)
  par(new=TRUE)
  matplot(y_fc_[[i]], xlab = 'Time [hh]', 
                      ylab = 'Energy [kWh]',
                      ylim = range(min(Y_sc_[[i]]), max(Y_sc_[[i]])),
                      type = 'b', col = 'blue', lwd = 1., lty = 1, pch = 1)
  par(new=TRUE)
  matplot(y_ac_[[i]], xlab = 'Time [hh]', 
                      ylab = 'Energy [kWh]',
                      ylim = range(min(Y_sc_[[i]]), max(Y_sc_[[i]])),
                      type = 'b', col = 'red', lwd = 1., lty = 1, pch = 2)
  par(new=TRUE)
  matplot(fPCA$mu, xlab = 'Time [hh]', 
                   ylab = 'Energy [kWh]',
                   main = paste(date, paste(source, name_[i], sep = ' | '), sep = ' | '),
                   ylim = range(min(Y_sc_[[i]]), max(Y_sc_[[i]])),
                   type = 'b', col = 'green', lwd = 1., lty = 1, pch = 3)
  grid(lty = 2,      
       col = 'black', 
       lwd = 0.5) 
  legend(1, max(Y_sc_[[i]]), legend = c(TeX('$f_i^{(sc)}$'), TeX('$f^{(ac)}$'), TeX('$f^{(fc)}$'), TeX('$mu(t)$')), 
                             col = c('gray75', 'red', 'blue', 'green'), lty = 1, cex = 0.75)
  dev.off()
  
  image_name = paste(path_to_save,  paste(name_[i], '_eigen.png', sep = ''), sep = '/')
  png(file = image_name, units = 'in', width = 6, height = 5, res = 150)
  col_ = c('red', 'blue', 'green', 'yellow', 'magenta')
  par(mfrow = c(1,1), mar = c(4,5,2,1), family = 'serif')
  matplot(fPCA$phi, xlab = 'Time [hh]', 
                    ylab = 'Score', 
                    ylim = range(c(min(fPCA$phi), max(fPCA$phi))),
                    main = paste(date, paste(source, name_[i], sep = ' | '), sep = ' | '),
                    type = 'l', col = col_, lwd = 1, lty = 1)
  grid(lty = 2,      
       col = 'black', 
       lwd = 0.5)   
  legend(1, max(fPCA$phi), legend = c(TeX('$phi_1$'), TeX('$phi_2$'), TeX('$phi_3$'), TeX('$phi_4$'), TeX('$phi_5$')), 
         col    = col_, lty = 1, cex = 0.75)
  dev.off()
  
}

image_name = paste(path_to_save,  paste('mus.png', sep = ''), sep = '/')
png(file = image_name, units = 'in', width = 6, height = 5, res = 150)
par(mfrow = c(1,1), mar = c(4,5,2,1), family = 'serif')
matplot(mus_, xlab = 'Time [hh]', 
              ylab = TeX('$mu_i(t)$'), 
              main = paste(date, paste(source, TeX('$mu(t)$'), sep = ': '), sep = ' - '),
              type = 'l', col = 'gray75', lwd = 0.5, lty = 1)
grid(lty = 2,      
     col = 'black', 
     lwd = 0.5) 
dev.off()

#image_name = paste(path_to_save,  paste('FVE.png', sep = ''), sep = '/')
#png(file = image_name, units = 'in', width = 6, height = 5, res = 150)
#par(mfrow = c(1,1), mar = c(4,5,2,1), family = 'serif')
#matplot(cums_, xlab = 'Principal Component #', 
#               ylab = 'Score', 
#               main = paste(date, paste(source, 'Variance Explained', sep = ': '), sep = ' - '),
#               ylim = range(c(0, 1)),
#               type = 'l', col = 'gray75', lwd = 0.5, lty = 1)
#grid(lty = 2,      
#     col = 'black', 
#     lwd = 0.5) 
#dev.off()

image_name = paste(path_to_save,  paste('PC1.png', sep = ''), sep = '/')
png(file = image_name, units = 'in', width = 6, height = 5, res = 150)
par(mfrow = c(1,1), mar = c(4,5,2,1), family = 'serif')
matplot(phis_[,1,], xlab = 'Time [hh]',
                    ylab = TeX('$phi_{1i}(t)$'), 
                    main = paste(date, paste(source, 'PC1', sep = ': '), sep = ' - '),
                    type = 'l', col = 'gray75', lwd = 0.5, lty = 1)
grid(lty = 2,      
     col = 'black', 
     lwd = 0.5) 
dev.off()

image_name = paste(path_to_save,  paste('PC2.png', sep = ''), sep = '/')
par(mfrow = c(1,1), mar = c(4,5,2,1), family = 'serif')
png(file = image_name, units = 'in', width = 6, height = 5, res = 150)
matplot(phis_[,2,], xlab = 'Time [hh]',
                    ylab = TeX('$phi_{2i}(t)$'), 
                    main = paste(date, paste(source, 'PC2', sep = ': '), sep = ' - '),
                    type = 'l', col = 'gray75', lwd = 0.5, lty = 1)
grid(lty = 2,      
     col = 'black', 
     lwd = 0.5) 
dev.off()

image_name = paste(path_to_save,  paste('PC3.png', sep = ''), sep = '/')
png(file = image_name, units = 'in', width = 6, height = 5, res = 150)
par(mfrow = c(1,1), mar = c(4,5,2,1), family = 'serif')
matplot(phis_[,3,], xlab = 'Time [hh]',
        ylab = TeX('$phi_{3i}(t)$'), 
        main = paste(date, paste(source, 'PC3', sep = ': '), sep = ' - '),
        type = 'l', col = 'gray75', lwd = 0.5, lty = 1)
grid(lty = 2,      
     col = 'black', 
     lwd = 0.5) 
dev.off()

image_name = paste(path_to_save,  paste('PC4.png', sep = ''), sep = '/')
png(file = image_name, units = 'in', width = 6, height = 5, res = 150)
par(mfrow = c(1,1), mar = c(4,5,2,1), family = 'serif')
matplot(phis_[,4,], xlab = 'Time [hh]',
        ylab = TeX('$phi_{4i}(t)$'), 
        main = paste(date, paste(source, 'PC4', sep = ': '), sep = ' - '),
        type = 'l', col = 'gray75', lwd = 0.5, lty = 1)
grid(lty = 2,      
     col = 'black', 
     lwd = 0.5) 
dev.off()

image_name = paste(path_to_save,  paste('PC5.png', sep = ''), sep = '/')
png(file = image_name, units = 'in', width = 6, height = 5, res = 150)
par(mfrow = c(1,1), mar = c(4,5,2,1), family = 'serif')
matplot(phis_[,5,], xlab = 'Time [hh]',
        ylab = TeX('$phi_{5i}(t)$'), 
        main = paste(date, paste(source, 'PC5', sep = ': '), sep = ' - '),
        type = 'l', col = 'gray75', lwd = 0.5, lty = 1)
grid(lty = 2,      
     col = 'black', 
     lwd = 0.5) 
dev.off()

i = 3
N_cols = 20
# Kernel density estimation
fPCA  = fPCA_[[i]]
xiKDS = kde(fPCA$xiEst, eval.points = fPCA$xiEst)
length(xiKDS$estimate)
xiQ_ = quantile(xiKDS$estimate, probs = seq(0, 1, l = N_cols + 1))
print(xiQ_)
xiCol_ = viridis::viridis(N_cols)[cut(xiKDS$estimate, breaks = xiQ_)]

par(mfrow = c(1,1), mar = c(4,5,2,1), family = 'serif')
matplot(t(Y_sc_[[i]]), xlab = 'Time [hh]',
                       ylab = 'Energy [kWh]',
                       type = 'l', col = xiCol_, lwd = 0.5, lty = 1)
grid(lty = 2,      
     col = 'black', 
     lwd = 0.5) 
