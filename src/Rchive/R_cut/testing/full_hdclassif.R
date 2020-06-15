library(SNPRelate)
library(plyr)
# install.packages("HDclassif")
library(HDclassif)
library(scrime)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

## loading the gds of the data and pullling some attributes out
gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds"
source("Analysis\\R\\functions\\data_loading.R")

genotypes <- replace(genotypes, genotypes == 0, 1)
genotypes <- replace(genotypes, genotypes == 3, NA)

genotypesImputed <- knncatimpute(t(genotypes))

# res <- hddc(genotypesImputed, model = "ALL")
# res <- hddc(genotypesImputed, model = "ABkQkDk")

# install.packages("dbscan")
library(dbscan)
res <- optics(genotypesImputed, minPts = 10)
for (i in 1:30) {
  print(extractDBSCAN(res, eps_cl = i))
}

res2 <- extractDBSCAN(res, eps_cl = 28)
str(res)
res$cluster
# ### MAPDP
# library(ClustMAPDP)
# library(matrixStats)
# 
# ## PCA
# wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds")
# pca <- snpgdsPCA(wheat, num.thread = 4, autosome.only = F, eigen.cnt = -1)
# snpgdsClose(wheat)
# 
# 
# N0 <- 1;                           # Prior count (concentration parameter)
# m0 <- rowMeans(t(pca$eigenvect[,1:10]));                     # Normal-Wishart prior mean
# a0 <- 10;                            # Normal-Wishart prior scale
# c0 <- 10/365;                          # Normal-Wishart prior degrees of freedom
# B0 <- diag(1./(0.05*rowVars(t(pca$eigenvect[,1:10]))));    # Normal-Wishart prior precision
# 
# start.time <- Sys.time()
# Rprof(tmp <- tempfile())
# r <- clustMapDP(t(pca$eigenvect[,1:10]),N0,m0,a0,c0,B0);
# Rprof()
# summaryRprof(tmp)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# print(time.taken)
# 
# r$z
# 
# # D <- 2; # Data dimensionality
# # K <- 5; # Number of clusters
# # 
# # #Generate random cluster locations and covariance matrices
# # manyZeros <- rep(0,D*D*K);
# # Sigma <- array(manyZeros,c(D,D,K));
# # Mu <- matrix(0,D,K);
# # 
# # for(k in 1:K){
# #   s <- matrix(rnorm(D*D),D);
# #   s <- t(s)%*%s;
# #   Sigma[,,k]<-0.15*s;
# #   Mu[,k] <- 2.8*matrix(rnorm(D),D)
# # }
# # # Generate random cluster mixture weights
# # Pi <- matrix(runif(K),1);
# # Pi <- Pi/sum(Pi);
# # 
# # N <- 4000;
# # 
# # # Generate categorical data for cluster assignments z
# # Z <- sample(K,N,replace = TRUE,prob=Pi);
# # 
# # #Generate multivariate Gaussian data for X
# # X <- matrix(0,D,N);
# # for(k in 1:K){
# #   i <- (Z==k);
# #   M <- sum(i);
# #   X[,i] <- sqrtm(Sigma[,,k])%*%matrix(rnorm(D*M),D)+repmat(Mu[,k],1,M);
# # }
# # i <- sample(1:N,N,replace=FALSE);
# # Z <- Z[i];
# # X <- X[,i];
# 
# # Set up Normal-Wishart MAP-DP prior parameters
# N0 <- 1;                           # Prior count (concentration parameter)
# m0 <- rowMeans(X);                     # Normal-Wishart prior mean
# a0 <- 10;                            # Normal-Wishart prior scale
# c0 <- 10/N;                          # Normal-Wishart prior degrees of freedom
# B0 <- diag(1./(0.05*rowVars(X)));    # Normal-Wishart prior precision
# start.time <- Sys.time()
# Rprof(tmp <- tempfile())
# r <- clustMapDP(X,N0,m0,a0,c0,B0);
# Rprof()
# summaryRprof(tmp)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# print(time.taken)

