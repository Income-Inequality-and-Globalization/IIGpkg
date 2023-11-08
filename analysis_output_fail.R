out <- readRDS("data/output/checkme.rds")
dim(out$Gibbs2_SM_SA$B)

out$Gibbs2_SM_SA$B[,,1]
min(out$Gibbs2_SM_SA$B[1, 1,])
max(out$Gibbs2_SM_SA$B[1, 1,])
checks <- sort(out$Gibbs2_SM_SA$B[1, 1,], decreasing = TRUE)[1:7]
id_fail_outlier <- which(out$Gibbs2_SM_SA$B[1, 1,] %in% checks)

test_mcmc_seq <- out$Gibbs2_SM_SA$B[1,1, ]
test_mcmc_seq <- test_mcmc_seq[-c(id_fail_outlier)]

check <- IIGpkg:::get_fails(out, iter_fail = id_fail_outlier[1])
dim(out$Gibbs2_SM_SA$f)
dim(out$Gibbs2_SM_SA$B)
dim(out$Gibbs2_SM_SA$D)
dim(out$Gibbs2_SM_SA$A)
dim(out$Gibbs2_SM_SA$V)




upper <- readRDS("./inst/legacy-output/upper_fail.rds")
lower <- readRDS("./inst/legacy-output/lower_fail.rds")
Sigma <- readRDS("./inst/legacy-output/sigma_fail.rds")
mean <- readRDS("./inst/legacy-output/mean_fail.rds")

num_samples <- 25000

set.seed(123)
out_test <- tmvnsim::tmvnsim(num_samples, length(upper), lower = lower, upper = upper, means = mean, sigma = Sigma)$samp
# set.seed(123)
# out_test2 <- tmvnsim::tmvnsim(25000, length(upper), lower = lower, upper = upper, means = mean, sigma = Sigma)$samp
# 
# identical(out_test, out_test2)
# plot(out_test[, 3], type = "l")

out_test_new_01 <- tmvtnorm::rtmvnorm(num_samples, mean = as.vector(mean), sigma = Sigma, lower = lower, upper = upper)
plot(out_test_new_01[, 3], type = "l")

out_test_new_02 <- tmvtnsim::rtmvnorm(mean = as.vector(mean), sigma = Sigma, lower = lower, upper = upper, n = num_samples)
plot(out_test_new_02[, 3], type = "l")


bench::mark(old = tmvnsim::tmvnsim(num_samples, length(upper), lower = lower, upper = upper, means = mean, sigma = Sigma),
            new_01 = tmvtnorm::rtmvnorm(num_samples, mean = as.vector(mean), sigma = Sigma, lower = lower, upper = upper),
            new_02 = tmvtnsim::rtmvnorm(mean = as.vector(mean), sigma = Sigma, lower = lower, upper = upper, n = num_samples),
            check = FALSE)

set.seed(123)
tmvtnsim::rtmvnorm(mean = as.vector(mean), sigma = Sigma, lower = lower, upper = upper, n = 1)
set.seed(123)
tmvtnsim:::rtmvnormcpp(mean = matrix(mean, ncol = ncol(Sigma)),
                       sigma = Sigma,
                       blc = diag(ncol(Sigma)),
                       lower = matrix(lower, ncol = ncol(Sigma)),
                       upper = matrix(upper, ncol = ncol(Sigma)),
                       init = matrix(0, ncol = ncol(Sigma)),
                       burn = 10)



colMeans(out_test2)
apply(out_test2, 2, sd)

colMeans(out_test_new_01)
apply(out_test_new_01, 2, sd)

colMeans(out_test_new_02)
apply(out_test_new_01, 2, sd)

# 
# 
# 
# d = 3;
# rho = 0.9;
# sigma = matrix(0, d, d);
# sigma = rho^abs(row(sigma) - col(sigma));
# blc = diag(1,d);
# n = 1000;
# mean = matrix(rep(1:d,n), nrow=n, ncol=d, byrow=TRUE);
# set.seed(1203)
# result = tmvtnsim::rtmvnorm(mean, sigma, blc, -1, 1, burn=50)
# 
# set.seed(42)
# test1 <- tmvtnsim::rtmvnorm(mean=1:d, sigma, blc, c(-1, -1, -1), c(1, 1, 1), burn=50, n=1000)
# set.seed(42)
# test2 <- tmvtnsim::rtmvnorm(mean=1:d, sigma, blc, -1, 1, burn=50, n=1000)



