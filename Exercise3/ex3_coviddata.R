## Exercise 3 - phylogeny building

library(stats)
library(ade4)
library(ape)
library(adegenet)
library(phangorn)


#-----------------------------------------------------------
## Read in the files
dna<-fasta2DNAbin("dna_fasta_B.1.13.FASTA")
metadata_b113 <- read.csv("metadata_B.1.13.csv",header=T)
# test the labels match:
identical(labels(dna), as.vector(metadata_b113$sequence_name)) #True!

dates_B113 <- metadata_b113 [,8]


dna
class(dna)
object.size(as.character(dna))/object.size(dna)
as.character(dna)[1:5, 1:10]
unclass(dna)[1:5, 1:10]
typeof(unclass(dna)[1:5, 1:10])
annot <- dates_B113
head(annot)

#-----------------------------------------------------------
##Distance-based phlogenies
D <- dist.dna(dna, model = "TN93")
class(D)
length(D)

#-----------------------------------------------------------
##Plot the pairwise distance directly
temp <- as.data.frame(as.matrix(D))
table.paint(temp, cleg = 0, clabel.row = 0.5, clabel.col = 0.5)
View(temp)


temp <- t(as.matrix(D))
temp <- temp[, ncol(temp):1]
par(mar = c(1, 5, 5, 1))
image(x = 1:40, y = 1:40, temp, col = rev(heat.colors(100)),
      xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(side = 2, at = 1:40, lab = rownames(dna), las = 2, cex.axis = 0.5)
axis(side = 3, at = 1:40, lab = rownames(dna), las = 3, cex.axis = 0.5)

#-----------------------------------------------------------
## Building trees

# TREE 1 - NJ
tre1 <- nj(D)
class(tre1)
tre1
plot(tre1, cex = 0.6)
title("A simple NJ tree for the COVID data")
write.tree(tre1, file = "COVIDnjTree")
# Notice that all the distances are really short, since these were sampled over a small time frame
# We do see quite a lot of clustering by location (CAMB, BRIS)

plot(tre1, show.tip = FALSE)
title("Unrooted NJ tree for the COVID data")
myPal <- colorRampPalette(c("red", "yellow", "green", "blue"))
tiplabels(annot, bg = fac2col(annot, col.pal = myPal),
          cex = 0.5) # we use fac2col instead because the dates default to a factor
temp <- pretty(seq.Date(from = as.Date("2020-03-19"), to = as.Date("2020-04-07"), by = 1), 5)
legend("topright", fill = any2col(temp, col.pal = myPal)$col,
       leg = temp, ncol = 2)

# is it rooted?
is.rooted(tre1)

plot(tre1, type = "unrooted", show.tip = FALSE)
title("Unrooted NJ tree for the COVID data")
tiplabels(annot, bg = fac2col(annot, col.pal = myPal),
          cex = 0.5)

# Let's root it

tre1.2 <- root(tre1, out = 1)
plot(tre1.2, show.tip = FALSE, edge.width = 2)
title("Rooted NJ tree for the COVID data")
tiplabels(annot, bg = transp(fac2col(annot, col.pal = myPal),                                   + 0.7), cex = 0.5, fg = "transparent")
axisPhylo()
temp <- pretty(seq.Date(from = as.Date("2020-03-19"), to = as.Date("2020-04-07"), by = 1), 5)
legend("topright", fill = transp(any2col(temp, col.pal = myPal)$col),
       leg = temp, ncol = 2)

# Assess the quality
x <- as.vector(D)
y <- as.vector(as.dist(cophenetic(tre1.2)))
plot(x, y, xlab = "original distance", ylab = "distance in the tree",
     main = "Is NJ appropriate?", pch = 20, col = transp("black",
                                                         0.1), cex = 3)
abline(lm(y ~ x), col = "red")
cor(x, y)^2

tre1.3 <- as.phylo(hclust(D, method = "average"))
y <- as.vector(as.dist(cophenetic(tre1.3)))
plot(x, y, xlab = "original distance", ylab = "distance in the tree",
     main = "Is UPGMA appropriate?", pch = 20, col = transp("black",
                                                            0.1), cex = 3)
abline(lm(y ~ x), col = "red")
cor(x, y)^2
plot(tre1.3, cex = 0.5)
title("UPGMA tree") # is UPGMA more appropriate here than in the example/ TB data?

# Bootstrapping
myBoots <- boot.phylo(tre1.2, dna, function(e) root(njs(dist.dna(e,model = "TN93")), 1))
myBoots

plot(tre1.2, show.tip = FALSE, edge.width = 2)
title("NJ tree + bootstrap values for the COVID data")
tiplabels(frame = "none", pch = 20, col = transp(fac2col(annot,
                                                         col.pal = myPal), 0.7), cex = 3, fg = "transparent")
axisPhylo()
temp <- pretty(seq.Date(from = as.Date("2020-03-19"), to = as.Date("2020-04-07"), by = 1), 5)
legend("topright", fill = transp(any2col(temp, col.pal = myPal)$col),
       leg = temp, ncol = 2)
nodelabels(myBoots, cex = 0.6)

# Collapse poorly supported nodes
temp <- tre1.2
N <- length(tre1.2$tip.label)
toCollapse <- match(which(myBoots < 70) + N, temp$edge[, 2])
temp$edge.length[toCollapse] <- 0
tre1.3 <- di2multi(temp, tol = 1e-05)

plot(tre1.3, show.tip = FALSE, edge.width = 2)
title("NJ tree for COVID data after collapsing weak nodes")
tiplabels(annot, bg = transp(fac2col(annot, col.pal = myPal),
                             0.7), cex = 0.5, fg = "transparent")
axisPhylo()
temp <- pretty(seq.Date(from = as.Date("2020-03-19"), to = as.Date("2020-04-07"), by = 1), 5)
legend("topright", fill = transp(any2col(temp, col.pal = myPal)$col),
       leg = temp, ncol = 2)



# TREE 2 - BIONJ
tre2 <- bionj(D)
class(tre2)
tre2
plot(tre2, cex = 0.6)
title("An improved NJ tree for the COVID data")
write.tree(tre2, file = "TBnjTree2")

# TREE 3 - FAST ME
tre3 <- fastme.bal(D)
class(tre3)
tre3
plot(tre3, cex = 0.6)
title("Minimum evolution tree  for the COVID data")
write.tree(tre3, file = "MinEvolTree")

# Notice we still keep seeing this one out-cluster

# TREE 3 - HCLUST
tre4 <- hclust(D)
class(tre4)
tre4
plot(tre4, cex = 0.6)
#title("Hierarchical clustering for the COVID data")
write.phyDat (tre4, file = "HclustTree")

# We can create very different looking trees from different methods!

#-----------------------------------------------------------
## Maximum parsimony

set.seed(5737292)
dna2 <- as.phyDat(dna)
class(dna2)
dna2
tre.ini <- nj(dist.dna(dna, model = "raw"))
tre.ini
parsimony(tre.ini, dna2)
tre.pars <- optim.parsimony(tre.ini, dna2) # apparently it's already optimal
tre.pars

plot(tre.pars, type = "unr", show.tip = FALSE, edge.width = 2)
title("Maximum-parsimony tree for the COVID data")
tiplabels(annot, bg = transp(fac2col(annot, col.pal = myPal),
                             + 0.7), cex = 0.5, fg = "transparent")
temp <- pretty(seq.Date(from = as.Date("2020-03-19"), to = as.Date("2020-04-07"), by = 1), 5)
legend("bottomleft", fill = transp(any2col(temp, col.pal = myPal)$col),
       leg = temp, ncol = 2)

#-----------------------------------------------------------
## Maximum likelihood phylogeny

dna2 <- as.phyDat(dna)
class(dna2)
dna2

tre.ini <- nj(dist.dna(dna, model = "TN93"))
tre.ini
pml(tre.ini, dna2, k = 4)
na.posi <- which(apply(as.character(dna), 2, function(e) any(!e %in%
                                                               c("a", "t", "g", "c")))) 
# We have some missing data - let's fix this
temp <- apply(as.character(dna), 2, function(e) sum(!e %in% c("a",
                                                              "t", "g", "c")))
plot(temp, type = "l", col = "blue", xlab = "Position in HA segment",
     ylab = "Number of NAs")
# remove the locations of NAs:
dna3 <- dna[, -na.posi]
dna3
table(as.character(dna3))
dna4 <- as.phyDat(dna3)

tre.ini <- nj(dist.dna(dna3, model = "TN93"))
fit.ini <- pml(tre.ini, dna4, k = 4)
fit.ini

fit <- optim.pml(fit.ini, optNni = TRUE, optBf = TRUE, optQ = TRUE,
                 optGamma = TRUE)


fit
class(fit)
names(fit)
anova(fit.ini, fit)
AIC(fit.ini)
AIC(fit) # new tree is better!

tre4 <- root(fit$tree, 1)
plot(tre4, show.tip = FALSE, edge.width = 2)
title("Maximum-likelihood tree for the COVID data")
tiplabels(annot, bg = transp(fac2col(annot, col.pal = myPal),
                             0.7), cex = 0.5, fg = "transparent")
axisPhylo()
temp <- pretty(seq.Date(from = as.Date("2020-03-19"), to = as.Date("2020-04-07"), by = 1), 5)
legend("topright", fill = transp(any2col(temp, col.pal = myPal)$col),
       leg = temp, ncol = 2)

write.nexus(tre4, file = "ML_covidtree.NEX")

