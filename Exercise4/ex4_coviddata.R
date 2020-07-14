## Exercise 4 - TransPhylo

library(TransPhylo)
library(ape)
library(coda)
library(lubridate)
library(phytools)

set.seed(0)
ph <- read.nexus("covid_datedtree.nex")# a timed phylogeny for the COVID data (made in iqtree2)
metadata_b113 <- read.csv("metadata_B.1.13.csv",header=T)
# this phylogeny has the samples in a different order to before, so lets reorder the metadata to match:
metadata_b113 <- metadata_b113[order(match(metadata_b113$sequence_name,ph$tip.label)),]
# test the labels match:
identical(ph$tip.label, as.vector(metadata_b113$sequence_name)) #True - nice!

# get dates into the right format - we only need a date to set the date of the phylogeny in time,
# so let's label the first sample as 'day 0' and work from there. 
dates <- as.Date(metadata_b113[,8])
length_time <- as.numeric(max(dates) -  min(dates))

# we see that this phylogeny doesn't contain much detail, since all the covid genomes sampled are very similar
plot(ph) 

# remove any multifurcations:
ph <- multi2di(ph) 
# do we have any negative branch lengths?
sum(ph$edge.length<0) # no! good news
#ph$edge.length <- pmax(ph$edge.length,1/365) # so we don't need this

plot(ph)
axisPhylo(backward = F)
ptree <- ptreeFromPhylo(ph,dateLastSample=length_time)
plot(ptree)

# COVID-19: gen time ~ gamma(mean = 5.2 days, sd = 1.72 days), samp time ~ Weibull(shape = 1.73,  scale = 9.85)
# But, TransPhylo needs gamma distributions only. 
# Let's approximate with: 
# gen time ~ gamma(mean = 5.2 days, sd = 1.72 days), samp time ~ gamma(shape = 1.73,  scale = 5.5)
w.shape=9.1
w.scale=1/0.56
ws.shape=1.73
ws.scale=5.5


# set dateT - the time observation stopped. 
# As with the TB data, we set this to a a little while after the last sampling date to avoid
# lots of unsampled cases towards the tips
dateT=length_time + 10


#-----------------------------------------------------------
##MCMC

# I'm using a small number of iterations, just for a preliminary run. You will want to increase
# this in order to get better results - we should use our MCMC diagnostics to confirm if we have enough
# iterations (trace plots, ESS)
res<-inferTTree(ptree,mcmcIterations=5000,w.shape=w.shape,w.scale=w.scale,
                ws.shape=ws.shape,ws.scale=ws.scale, dateT=dateT, startPi=0.9, updatePi=F)
# We can turn on/off the estimation of different parameters depending on how good the
# estimation is, and what we are interested to learn/how much prior knowledge we have
# I have turned off updating Pi here, but you could update it if you wish

plot(res)
mcmc=convertToCoda(res)
effectiveSize(mcmc)
med=medTTree(res)
plot(med)

##transmission
ttree=extractTTree(med)
plot(ttree,type='detailed',w.shape,w.scale)

mat=computeMatTDist(res)
lattice::levelplot(mat,xlab='',ylab='')

a=getIncidentCases(res,show.plot = T)
a=getGenerationTimeDist(res,show.plot = T)
a=getInfectionTimeDist(res,k=c('1','2'),show.plot = T)

a=getOffspringDist(res,k=c('1','2'),show.plot = T)#might get error due to zero length replacement

# I have not included code for it here, but it would be interesting to look at if transmission 
# seems to be clustered by location (BRIS - Bristol, CAMB = Cambridge)

# VisNetwork

library(visNetwork)
library(RColorBrewer)
source("transphylo_extras.R") 
# this .R file needs to be in your working directory, else you will need to change the path

mywiw=computeMatWIW(res)

mynp = networkPlot(res[[1]]$ctree,showTimes = TRUE,shapefactor = 3)
modnp=mynp
modnp$edges$width=1

modnp$nodes$label=as.character(modnp$nodes$label)
modnp$nodes$label[which(modnp$nodes$groups=="unsampled")]="unsamp"
modnp$nodes$font.size=ifelse(modnp$nodes$groups == "sampled", 20, 10)
visNetwork(modnp$nodes,modnp$edges,width = "900px",height="600px") %>% 
  visLegend(width=0.2,addNodes=mynp$lnodes,useGroups=F) %>% visSave(file="demo.html")
# uncomment the last section to save the network to html
# the tip names are really long which causes the bubbles to be giant!
# You can edit the labels by changing modnp$nodes$label to something shorter







