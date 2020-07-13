
## Exercise 2 - outbreaker

library(ape)
library(outbreaker2)
#install.packages("Hmisc")
library(Hmisc)
library(lubridate)

# for this to run as-is, make sure your fasta and txt data files are in your working directory e.g.
# setwd("location of your files")



#read data and set them as dna and dates
dna <- read.dna(file = "roetzer2013.fasta",format = "fasta")
dates <- read.table(file="roetzer_dates.txt")
dates <- as.Date(dates[,1])
dates<-(year(dates) - 1997)*12 + month(dates) + day(dates)/monthDays(dates)



#sample w and f by gamma distribution
#TB: gen time ~ gamma(shape = 1.3, rate = 0.3), samp time ~ gamma(shape = 1.1, rate = 0.4) 
#  = gen time ~ gamma(shape = 1.3, rate = 0.025) months, samp time ~ gamma(shape = 1.1, rate = 0.033) months

x<- seq(1, 60) # we truncate to a max time of 60 months for both distributions
scaleto1 <- function(x){x/sum(x)}
w <- scaleto1(sapply(x, dgamma, shape=1.3, rate=0.025))
f <- scaleto1(sapply(x, dgamma, shape=1.1, rate=0.033))



#plot the generation/sampling time
col <- "#6666cc"
plot(w, type = "h", xlim = c(0, 60), 
     lwd = 30, col = col, lend = 2, 
     xlab = "Months after infection", 
     ylab = "p(new case)", 
     main = "Generation time distribution")
plot(f, type = "h", xlim = c(0, 60), 
     lwd = 30, col = col, lend = 2, 
     xlab = "Months after infection", 
     ylab = "p(new case)", 
     main = "Sampling time distribution")


# set up the names correctly
names(dates) = labels(dna)
data <- outbreaker_data(dna = dna, dates = dates, w_dens = w, f_dens = f)

# for testing, a configuration where we run very few iterations
config2 <- create_config(n_iter = 100,
                         sample_every = 1,
                         move_kappa = FALSE)

# set a seed to make the results reproducible
set.seed(47463)

res <- outbreaker(data, config2)


# View the results
class(res)
dim(res)
res
names(res)

plot(res)
plot(res, "prior")
plot(res, "mu")
plot(res, "t_inf_15")
# adjust burn in accordingly for your number of iterations (say about 20% as a rough)
plot(res, burnin = 20)
plot(res, "mu", "density", burnin = 20)
plot(res, type = "alpha", burnin = 20)
plot(res, type = "t_inf", burnin = 20)
plot(res, type = "kappa", burnin = 20)
plot(res, type = "network", burnin = 20, min_support = 0.01)
summary(res)
