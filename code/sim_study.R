#####################
## Example single iteration
#####################

library(updog)

## Set possible parameter values here that we iterate over.

## Set the simulation parameters for this iteration
## vary these each iteration.
## Do same parameter settings for 100 iteration each setting.

## For loop begins here
nsamp          <- 30
seq_error      <- 0.01
ploidy         <- 6
bias_val       <- 1
overdispersion <- 0
allele_freq    <- 0.1

## Generate genotypes
genovec <- rgeno(n           = nsamp,
                 ploidy      = ploidy,
                 model       = "hw",
                 allele_freq = allele_freq)

## Generate read-count data data.
sizevec <- rep(50, length.out = nsamp)
refvec <- rflexdog(sizevec = sizevec,
                   geno    = genovec,
                   ploidy  = ploidy,
                   seq     = seq_error,
                   bias    = bias_val,
                   od      = overdispersion)

## Fit stan on sizevec and refvec

## Evaluate how close our genotype estimates are to genovec

## Use flexdog() to estimate genotypes

## Evaluate how close flexdog() genotype estimates are to genovec.

## For loop ends here

## Result: misclassification error rates for updog and stan.

