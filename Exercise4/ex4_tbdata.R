
## Exercise 4 - TransPhylo

library(TransPhylo)
library(ape)
library(coda)
library(lubridate)

set.seed(0)
ph <- read.nexus("roetzer2013.nex")# the TB .nex tree provided in Ex4
dates <- as.Date(read.table("roetzer_dates.txt")[,1])
dates <- decimal_date(dates)

# normally we might run these, but actually we don't have any in this tree:
#ph <- multi2di(ph) # remove multifurcations
#ph$edge.length <- pmax(ph$edge.length,1/365) # make sure all branch lengths are at least 1 day

plot(ph)
axisPhylo(backward = F)
ptree <- ptreeFromPhylo(ph,dateLastSample=(max(dates)))
plot(ptree)

# TB: gen time ~ gamma(shape = 1.3, rate = 0.3) years, samp time ~ gamma(shape = 1.1, rate = 0.4) years
w.shape=1.3
w.scale=1/0.3 # scale is 1/rate 
ws.shape=1.1
ws.scale=1/0.4

dateT=max(dates)+ 30/365 # let's use about 1 month later

#-----------------------------------------------------------
##MCMC

# I'm using a small number of iterations, just for a preliminary run. You will want to increase
# this in order to get better results - we should use our MCMC diagnostics to confirm if we have enough
# iterations (trace plots, ESS)
res<-inferTTree(ptree,mcmcIterations=5000,w.shape=w.shape,w.scale=w.scale,
                ws.shape=ws.shape,ws.scale=ws.scale, dateT=dateT, startPi=0.9, 
                updatePi=F, updateNeg = FALSE)
# We can turn on/off the estimation of different parameters depending on how good the
# estimation is, and what we are interested to learn/how much prior knowledge we have
# I have turned off updating Neg here

plot(res)
mcmc=convertToCoda(res)
effectiveSize(mcmc)
med=medTTree(res)
plot(med)
# We are getting a lot of unsampled cases near the tips, so we might want to increase dateT, 
# but we should also just do more iterations

##transmission
ttree=extractTTree(med)
plot(ttree,type='detailed',w.shape,w.scale)
# it's a little hard to read this with so many cases, it's better for smaller outbreaks

mat=computeMatTDist(res)
lattice::levelplot(mat,xlab='',ylab='')

a=getIncidentCases(res,show.plot = T)
a=getGenerationTimeDist(res,show.plot = T)
a=getInfectionTimeDist(res,k=c('1','2'),show.plot = T)

a=getOffspringDist(res,k=c('1','2'),show.plot = T)#might get error due to zero length replacement



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
  visLegend(width=0.2,addNodes=mynp$lnodes,useGroups=F)   #%>% visSave(file="demo.html")
# uncomment the last section to save the network to html




