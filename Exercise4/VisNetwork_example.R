#############################
## VisNetwork example code ##
#############################

# This example is set up to take output from TransPhylo, and create a VisNetwork
# For this code to run as-is, you need a TransPhylo output object named 'res' in 
# your environment

install.packages("visNetwork")

library(TransPhylo)
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



