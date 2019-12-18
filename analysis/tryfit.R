library(rstan)
library(updog)
data("snpdat")

## prepare data
subsnp <- snpdat[snpdat$snp == "SNP3", ]
ploidy <- 6
datalist <- list(K = ploidy,
                 N = nrow(subsnp),
                 x = subsnp$counts,
                 n = subsnp$size)

## run STAN
vbmodel <- rstan::stan_model(file = "../code/vargeno.stan", )
vbout <- rstan::vb(object = vbmodel, data = datalist)

## Get posterior genotypes
get_postprob <- function(x) {
  postprobmat <- matrix(summary(x, pars = "logpostprobmat")$summary[, "mean"], ncol = ploidy + 1, byrow = TRUE)
  rowmax <- apply(postprobmat, 1, max)
  postprobmat <- postprobmat - rowmax
  postprobmat <- exp(postprobmat - log(rowSums(exp(postprobmat))))
  return(postprobmat)
}

get_maxgeno <- function(postprobmat) {
  apply(postprobmat, 1, which.max) - 1
}

postprobmat <- get_postprob(vbout)

genoest <- get_maxgeno(postprobmat)

## expit function
expit <- function(x) {
  1 / (exp(-x) + 1)
}

## get parameter estimates
get_pars <- function(x) {
  tparvec <- summary(x, pars = c("logit_eps", "logit_tau", "log_h"))$summary[, "mean"]
  parvec <- vector(mode = "numeric", length = 3)
  names(parvec) <- c("eps", "tau", "h")
  parvec[["eps"]] <- expit(tparvec[["logit_eps"]])
  parvec[["tau"]] <- expit(tparvec[["logit_tau"]])
  parvec[["h"]] <- exp(tparvec[["log_h"]])
  return(parvec)
}

parvec <- get_pars(vbout)

plot_geno(refvec = datalist$x,
          sizevec = datalist$n, ploidy = 6,
          geno = genoest,
          seq = parvec[[1]],
          bias = parvec[[3]])

## Compare to updog
uout <- flexdog(refvec = datalist$x, sizevec = datalist$n, ploidy = datalist$K, model = "hw")

plot(uout)

uout$seq
parvec[["eps"]]

uout$par$alpha
summary(vbout, pars = "alpha")$summary[, "mean"]

uout$od
parvec[["tau"]]

uout$bias
parvec[["h"]]

## Run good old MCMC
sout <- rstan::stan(file = "../code/vargeno.stan", data = datalist)
sout

mcmc_postprob <- get_postprob(sout)
mcmc_maxpost  <- get_maxgeno(mcmc_postprob)

plot(mcmc_postprob, postprobmat)
abline(0, 1)

sum(mcmc_maxpost != genoest)






