##################################################
### Code for generating correct generation and ###
### sampling time distributions for outbreaker ###
##################################################

# outbreaker requires the generation time w_dens and sampling time f_dens
# distributions in a specific format: they should be numeric vectors of length M 
# indicating the generation time/sampling time distribution t= 1, 2, ... time steps
# after infection.

# For example, if the collection dates are provided at a daily resolution, then the
# time scale is 1 day and a vector w_dens should be provided where the first element 
# corresponds to the probability that the generation time is 1 day long, the second
# element corresponds to the probability that the generation time is 2 days long etc. 

# This vector should therefore sum to 1, or if it does not then outbreaker will 
# scale it as such. 

# If instead we desire the generation/sampling time accordinging to some known 
# distribution - say a gamma distribution with known shape and scale, we can set
# them up as follows, truncating the distribution to a reasonable number of days:

# we truncate to a max of 15 days - this is a bit of a hack but will suffice here
x<- seq(1, 15)

# scale to 1 function
scaleto1 <- function(x){x/sum(x)}

# Let's say we want gamma with shape  = 5, rate = 1 => mean = 2, for both w and f

w <- scaleto1(sapply(x, dgamma, shape=5, rate=1))


f <- scaleto1(sapply(x, dgamma, shape=5, rate=1))

