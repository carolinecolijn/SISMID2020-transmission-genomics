## Exercise 2 - outbreaker

library(ape)
library(seqinr)
library(outbreaker2)
library(adegenet)

# for this to run as-is, make sure your fasta and csv data files are in your working directory e.g.
# setwd("location of your files")


## Read in the files
dna_fasta_sub_B.1.13<-fasta2DNAbin("dna_fasta_B.1.13.FASTA")
date_non_na_B.1.13 <- read.csv("metadata_B.1.13.csv",header=T)
metadata_b113 <- read.csv("metadata_B.1.13.csv",header=T)
# test the labels match:
identical(labels(dna_fasta_sub_B.1.13), as.vector(date_non_na_B.1.13$sequence_name)) #True!

dates_B113 <- as.Date(metadata_b113 [,8])



## wdens and fdens
#gen time ~ gamma(mean = 5.2 days, sd = 1.72 days) 
#samp time ~ Weibull(shape = 1.73,  scale = 9.85) 
# we truncate to a max of 15 days - this is a bit of a hack but will suffice here
x<- seq(1, 15)
# scale to 1 function
scaleto1 <- function(x){x/sum(x)}
# define the distributions
w <- scaleto1(sapply(x, dgamma, shape=9.1, rate=1/0.56))
f <- scaleto1(sapply(x, dweibull, shape=1.73, scale=9.85))

#plot the generation/sampling time
col <- "#6666cc"
plot(w, type = "h", xlim = c(0, 15), 
     lwd = 30, col = col, lend = 2, 
     xlab = "Days after infection", 
     ylab = "p(new case)", 
     main = "Generation time distribution")
plot(f, type = "h", xlim = c(0, 15), 
     lwd = 30, col = col, lend = 2, 
     xlab = "Days after infection", 
     ylab = "p(new case)", 
     main = "Sampling time distribution")




## Run outbreaker

dates_B113 <- as.Date(dates_B113)
# take the name from the metadata file or from the dna labels, either works just the same:
names(dates_B113) <- as.vector(metadata_b113$sequence_name)
#names(dates) = labels(dna_fasta_sub_B.1.13)

seq_outbreaker_B.1.13 <- outbreaker_data(dates = dates_B113 , dna = dna_fasta_sub_B.1.13, w_dens = w, f_dens = f)
# we don't have any contact tracing data

# for testing, a configuration where we run very few iterations
config <- create_config(n_iter = 100,
                        sample_every = 1)

# set a seed to make the results reproducible
set.seed(999)
res <- outbreaker(data = seq_outbreaker_B.1.13, config)




## Outputs
summary(res)
class(res)
dim(res)
plot(res)
plot(res,"prior")
# adjust burn in accordingly to your number of iterations
plot(res,"prior", burnin=20)
plot(res,"mu", burnin=20)
plot(res,burnin=20)
plot(res,"mu","hist",burnin=20)
plot(res,"mu","density",burnin=20)
plot(res,type="alpha",burnin=20)
plot(res,type="t_inf",burnin=20) 
plot(res,type="kappa",burnin=20)
plot(res,type="network",burnin=20,min_support=0.01)
summary(res)
# by the look of these results, we definitely need more than 100 iterations!!


