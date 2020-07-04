#' @param record MCMC output produced by inferTTree
#' @param k Case whose posterior infection times are to be extracted. Either an integer or a string matching one of the case names in the data
#' @return A vector of posterior infection times for case k 
#' @examples
#' gettimes(record[501:length(record)],"Case1") # posterior infection times for case "Case1", disregarding the first 500 record entries
#' @author Caroline Colijn
#' gettimes(record,2) # all posterior infection times for case 2
getInfectionTimes <- function(record,k) {
  if (is.numeric(k)) {
    mytimes= vapply(1:length(record), function(x) { tt=extractTTree(record[[x]]$ctree); return(tt$ttree[k,1])},FUN.VALUE=1);
    return(mytimes)}
  else {
    mytimes= vapply(1:length(record), function(x) { tt=extractTTree(record[[x]]$ctree); ii=which(tt$nam==k); return(tt$ttree[ii,1])},FUN.VALUE=1);
    return(mytimes)}
}

#' @param ctree Combined tree object 'ctree' from the record produced by inferTTree, for example record[[1]]$ctree where record is the output of inferTTree
#' @return Vector of times between becoming infected and infecting others (generation times) in the posterior
#' @author Caroline Colijn
#' @examples 
#' getGenerationTimes(record[[1]]$ctree)
getGenerationTimes <-  function(ctree) { tt=extractTTree(ctree)$ttree;
# 3rd column of ttree lists the infectors; exclude source
infectors=tt[,3]; infectors=infectors[infectors!=0];
# times at which each infector infected others:
infothers=tt[tt[,3]!=0,1]; 
# times at which each infector was herself infected:
gotinfd=tt[infectors,1];
return(infothers-gotinfd);}

#' @param ctree Combined tree object 'ctree' from the record produced by inferTTree, for example record[[1]]$ctree where record is the output of inferTTree
#' @return Vector of times between becoming infected and getting sampled
#' @author Caroline Colijn
#' @examples 
#' getGenerationTimes(record[[1]]$ctree)
getTimesToSampling <-  function(ctree) { tt=extractTTree(ctree)$ttree;
ns=sum(!is.na(tt[,2]))
return(tt[1:ns,2]-tt[1:ns,1])
}

# number unsampled cases 
#' @param ctree Combined tree object 'ctree' from the record produced by inferTTree, for example record[[1]]$ctree where record is the output of inferTTree
#' @return The number of unsampled cases in ctree
#' @examples 
#' getNumberUnsampled(ctree)
#' @author Caroline Colijn
getNumberUnsampled <- function(ctree) {tt=extractTTree(ctree)$ttree; return(sum(is.na(tt[,2])))}

# simple visualisation of transmission tree with visNetwork
library(visNetwork)
#' @param ctree Combined tree object 'ctree' from the record produced by inferTTree, for example record[[1]]$ctree where record is the output of inferTTree
#' @author Caroline Colijn
#' @example 
#' vninfo=networkPlot(record[[10]]$ctree)
#' visNetwork(vninfo$nodes,vninfo$edges) %>% visLegend(width=0.2,addNodes=vninfo$lnodes,useGroups=F)
networkPlot <- function(ctree,showTimes=T,shapefactor=3) {
    tt=extractTTree(ctree)
    info1=tt$ttree  
    numCases=nrow(info1)
    numSamp=length(tt$nam)  # number sampled; first group
    numUnsamp=nrow(info1)-numSamp; # number unsampled 
    SimpleWiw <- cbind(info1[,3],1:numCases) # infector, infectee
    if (!showTimes)  nodes <- data.frame(id = 1:nrow(info1), label=c(tt$nam, (numSamp+1):numCases) )
    if (showTimes) {
      infTimes=tt$ttree[,1]
      labs=paste(c(tt$nam, (numSamp+1):numCases)," (",0.1*round(10*infTimes),") ",sep="") # should really convert to dd-mm-yy
      nodes <- data.frame(id = 1:numCases,label=labs)
    }
   nodes$groups=c(rep("sampled",numSamp), rep("unsampled",numUnsamp))
  nodes$value=1; nodes$value[nodes$groups=="sampled"]=shapefactor
  colors=brewer.pal(4,"Spectral")
  pal<-colorRampPalette(colors) 
  # early cases have higher values in this ordering: 
  nodes$color=pal(numCases)[numCases-rank(infTimes)+1]
  nodes$shape="circle"

  
  # 3 demonstrative colours for the colour key: 
  lnodes <- data.frame(label = c(round(min(infTimes)), round(median(infTimes)),round(max(infTimes))),
                       shape = c( "cirle"), color = pal(3)[c(3,2,1)],size=1)
  edges <- data.frame(from = SimpleWiw[,1], to = SimpleWiw[,2],
                      arrows="to")
  thesource=edges$to[which(edges$from==0)]
#   nodes$shape[thesource]="star" # thought this made the source look too different
  nodes$color[thesource]="darkgrey"
  # more demonst cols
  lnodes<- data.frame(label=round(quantile(infTimes, seq(0,by=0.1,1))),shape=c("circle"), color=pal(11)[seq(11,1)],size=1);
#   visNetwork(nodes, edges) %>% visLegend(width=0.3,addNodes=lnodes,useGroups = F)
  return(list(nodes=nodes,edges=edges,lnodes=lnodes))
  }

# requires treespace. 
# create list of wiw information in order to compute transmission tree distances
#' @param record  MCMC output produced by inferTTree
#' @return list of MRCI information required by wiwTreeDist, one entry for each transmission tree that is included
 getTTreeDistInfo <- function(record) {
  matList <- lapply(1:length(record), function(x) {
  info <- extractTTree(record[[x]]$ctree)$ttree
  wiw <- cbind(info[,3],1:length(info[,1]))
  findMRCIs(wiw)$mrciDepths
})
return(matList)
}

 
 
 
 
 
