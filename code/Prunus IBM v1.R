library(raster)
library(sp)

#################
# GROWTH PARAMS #
#################

(grow.T.parm <- t(read.csv("data/tree growth parameters Apr21.csv", header=T, row.names=1)))
(grow.T.parm.unscale <- read.csv("data/tree growth model unscale params.csv", header=T, row.names=1))
(grow.S.parm <- read.csv("data/seedling lmer growth params CAA Apr21.csv", header=T, row.names=1))
(surv.parm <- t(read.csv("data/surv params Apr21.csv", header=T, row.names=1)))
(surv.parm.unscale <- read.csv("data/survival model unscale params.csv", header=T, row.names=1))
(fruit.parm <- read.csv("data/fruiting parameters Apr21.csv", header=T, row.names=1))
(tran.parm <- read.csv("data/transition params CAA Jun21.csv", header=T, row.names=1))
(rec.parm <- t(read.csv("data/dispersal kernel parameters Apr21.csv", header=T, row.names=1)))
(HM.parm <- t(read.csv("data/HM values.csv", header=T, row.names=1)))

# standardize some spp names
colnames(surv.parm) <- gsub(" ", ".", colnames(surv.parm))
colnames(grow.T.parm) <- gsub(" ", ".", colnames(grow.T.parm))


#################################
# DEFINING VITAL RATE FUNCTIONS #
#################################

# Note on notations:
# Seedling phase size units are log-transformed height, h
# Adult phase size units are log-transformed DBH, z
# At the end of each year interval, sizes are denoted by h1 and z1

# intraS, intraA denote symmetric and asymmetric intraspecific competition, respectively
# interS, interA, interspecific
# CAE (conspecific adult effect) and HAE (heterospecific) denote similar, but are used for seedling growth only
# also, for adult tree growth, inter/intra are calculated with distance from focal tree incorporated


## TREE Survival function, logistic regression
# note: need to add interS effect eventually

sT_z <- function(terrain, trees, sp, intraA, interAHM)
{
  m.par <- surv.parm
  z <- trees[,"logdbh"]
  plot_type <- ifelse(extract(terrain, trees[,c("x","y")])==2, "wet", "dry")
  plot_type[is.na(plot_type)] <- "dry"
  
  # scale intraS and interS, then calculate p1 and p2 from them
  intraA.scaled <- (log(intraS + surv.parm.unscale["intraA.unscale.logoffset",]) -
                      surv.parm.unscale["intraA.unscale.mu",]) / surv.parm.unscale["intraA.unscale.sigma",]
  interAHM.scaled <- (log(interAHM + surv.parm.unscale["interAHM.unscale.logoffset",]) 
                      - surv.parm.unscale["interAHM.unscale.mu",]) /  surv.parm.unscale["interAHM.unscale.sigma",]
  
  p1 <- m.par[paste0("p1.", plot_type), sp] + (intraS.scaled * m.par["p1.intraS", sp])
  p2 <- m.par["p2", sp] + (intraS.scaled * m.par["p2.intraS", sp])
  
  surv.prob <- ifelse( z <= 1,
                       m.par["K",sp] / ( 1 + exp(-m.par["r1",sp] * (z - p1)) ),
                       m.par["K",sp] / ( 1 + exp(-m.par["r2",sp] * (z - p2)) )
  )
  
  # return binary outcome
  surv <- rbinom(n = length(z), prob = surv.prob, size = 1)
  surv[is.na(surv)] <- 0
  return(surv)
}
#sT_z(nssf.m, ppoT, "Prunus.polystachya", 100)
#sT_z(nssf.m, ppoT, "Prunus.polystachya", intra.calc(ppoT, ppoS, "Prunus.polystachya")[[1]])

## TREE Growth function
# this is the only function that uses dbh instead of logdbh
# note: need to add interS effect eventually

GT_z1z <- function(terrain, trees, sp, intraA)
{
  m.par <- grow.T.parm
  z <- trees[,"logdbh"]
  dbh <- exp(z)
  plot_type <- ifelse(extract(terrain, trees[,c("x","y")])==2, "wet", "dry")
  plot_type[is.na(plot_type)] <- "dry"
  
  # scale intraA and interS
  intraA.scaled <- (log(intraA + grow.T.parm.unscale["grow.intraA.unscale.logoffset",]) -
                      grow.T.parm.unscale["grow.intraA.unscale.mu",]) / grow.T.parm.unscale["grow.intraA.unscale.sigma",]
  #interS.scaled <- (log(interS) - grow.T.parm.unscale["grow.interS.unscale.mu",]) / 
  #  grow.T.parm.unscale["grow.interS.unscale.sigma",]
  loggrowth <- (
    m.par["a",sp] * 
      (dbh ^ (m.par["b", sp] + intraA.scaled * m.par["b.intraA", sp])) *
      exp(-dbh * (m.par["c", sp] + intraA.scaled * m.par["c.intraA", sp]))
  ) + rnorm(n = length(z), mean = 0, sd = m.par["sigma",sp]) # add the noise in (sigma)
  dbh1 <- dbh + (exp(loggrowth) - 1)
  dbh1[which(dbh1 < 1)] <- 1
  z1 <- log(dbh1)
  names(z1) <- NULL
  return(z1)
}
#plot(GT_z1z(nssf.m, ppoT, "Prunus.polystachya", intra.dist.calc(ppoT)) ~ ppoT[,"logdbh"], cex=3)
#abline(0,1, lty=2)

## FLOWERING function, logistic regression

p_bz <- function(terrain, trees, sp)
{
  m.par <- fruit.parm
  z <- trees[,"logdbh"]
  plot_type <- ifelse(extract(terrain, trees[,c("x","y")])==2, "wet", "dry")
  plot_type[is.na(plot_type)] <- "dry"
  
  linear.p <- m.par[paste0("fruit.int.",plot_type),sp] + m.par["fruit.z",sp] * z      # linear predictor
  
  # return binary outcome
  surv <- rbinom(n = length(linear.p), prob = 1/(1+exp(-linear.p)), size = 1)	# logistic transformation to probability
  surv[is.na(surv)] <- 0
  return(surv)
}
#p_bz(nssf.m, ppoT, "Prunus.polystachya")

## SEEDLING production function (this is seedling recruitment, not seed production)

b_z <- function(terrain, trees, sp)
{
  # first few lines can be modified for later expansion
  m.par <- fruit.parm
  z <- trees[,"logdbh"]
  plot_type <- ifelse(extract(terrain, trees[,c("x","y")])==2, "wet", "dry")
  plot_type[is.na(plot_type)] <- "dry"
  
  N <- m.par["rec.int",sp] + m.par["fruited.z",sp] * z 	# seed production of a size z plant
  N[is.na(N)] <- 0
  return(round(exp(N),0))						# convert log(recruitment) to actual recruitment no.
}
#b_z(nssf.m, ppoT, "Prunus.polystachya")
#p_bz(nssf.m, ppoT, "Prunus.polystachya")*b_z(nssf.m, ppoT, "Prunus.polystachya")

## SEEDLING recruit size pdf

c_0h1 <- function(n, sp)
{
  # first few lines can be modified for later expansion
  m.par <- rec.parm
  h1 <- rnorm(n, mean = m.par["rcsz",sp], sd = m.par["rcsd",sp])
  return(h1)
}
#c_0h1(100, "Prunus.polystachya")

## SEEDLING Growth function
# note: need to add HAE effect eventually

GS_h1h <- function(terrain, seedlings, sp, CAE)
{
  # first few lines can be modified for later expansion
  m.par <- grow.S.parm
  h <- seedlings[,"logheight"]
  plot_type <- ifelse(extract(terrain, seedlings[,c("x","y")])==2, "wet", "dry")
  plot_type[is.na(plot_type)] <- "dry"
  CAE.scaled <- ( log(CAE+10) - m.par["grow.S.CAE.unscale.mu",sp] ) / m.par["grow.S.CAE.unscale.sigma",sp]
  
  mu <- m.par["grow.S.int",sp] + 					# intercept
    m.par["grow.S.CAE", sp] * CAE.scaled +			# CAE
    m.par["grow.S.h", sp] * h +        				# height
    m.par["grow.S.CAExh", sp] * h * CAE.scaled			# interaction
  
  h1 <- mu + rnorm(n = length(mu), mean = 0, sd = m.par["grow.S.sigma",sp])
  
  names(h1) <- NULL
  return(h1)
}
#plot(GS_h1h(nssf.m, ppoS, "Prunus.polystachya", 100) ~ ppoS[,"logheight"])
#points(GS_h1h(nssf.m, ppoS, "Prunus.polystachya", 1000) ~ ppoS[,"logheight"], col="red")
#abline(0,1,lty=2,lwd=2)

## SEEDLING Survival function, logistic regression

sS_h <- function(terrain, seedlings, sp, intraS)
{
  m.par <- surv.parm
  h <- seedlings[,"logheight"]
  plot_type <- ifelse(extract(terrain, seedlings[,c("x","y")])==2, "wet", "dry")
  plot_type[is.na(plot_type)] <- "dry"
  
  # transform logheight to logdbh
  z <- tran.parm["tran.int",sp] + h*tran.parm["tran.h",sp]
  
  # scale intra, then calculate p1 from it
  intraS.scaled <- (log(intraS + surv.parm.unscale["surv.intraS.unscale.logoffset",]) -
                      surv.parm.unscale["surv.intraS.unscale.mu",]) / surv.parm.unscale["surv.intraS.unscale.sigma",]
  
  p1 <- m.par[paste0("p1.", plot_type), sp] + (intraS.scaled * m.par["p1.intraS", sp])
  
  surv.prob <- m.par["K",sp] / ( 1 + exp(-m.par["r1",sp] * (z - p1)) )
  
  # return binary outcome
  surv <- rbinom(n = length(z), prob = surv.prob, size = 1)
  surv[is.na(surv)] <- 0
  return(surv)
}
#sS_h(nssf.m, ppoS, "Prunus.polystachya", 100) 
#sS_h(nssf.m, ppoS, "Prunus.polystachya", intra.calc(ppoT, ppoS, "Prunus.polystachya")[[2]]) * GS_h1h(nssf.m, ppoS, "Prunus.polystachya", 100)
#plot(jitter(sS_h(nssf.m, ppoS, "Prunus.polystachya")) ~ sqrt(ppoS[,"logheight"]))

## SEEDLING-TO-SAPLING transition function

# takes h1 from seedlings, thus seedlings df must have already completed growth/mortality
T_z1h1 <- function(trees, seedlings, sp)
{
  m.par <- tran.parm
  h1 <- seedlings[,"logheight"]
  
  # compute logdbh from logheight
  mu <- (m.par["tran.int",sp] + m.par["tran.h",sp] * h1)
  
  # only initiate transitions if there is at least one seedling larger than 1cm DBH
  if(sum(mu>0) > 0){
    tran.index <- which(mu > 0) # mu > 0 because log(1) = 0
    tran.loc <- ppoS[tran.index,c("x","y")]
    # add sapling to tree df
    newsaplings <- cbind(tran.loc, mu[tran.index])
    names(newsaplings) <- names(trees)
    trees <-  rbind(trees, newsaplings)
    # remove seedling from seedling stack
    seedlings <- seedlings[-tran.index,]
  }
  return(list(trees, seedlings))
}
# ppoTS <- T_z1h1(ppoT, ppoS, "Prunus.polystachya")
# nrow(ppoTS[[1]]); nrow(ppoTS[[2]])

## DISPERSAL KERNEL function

# define a 2DT pdf
twoDT.pdf <- function(dist, S, L) { S / 
    (pi * L * (1 + (dist^2/L))^(S+1))
}

# this function takes input of germinant number, generates vector of distances from mother tree 
# up to 100m from mother tree, and to a spatial resolution of 0.1m (because of discretization of pdf in "support")
twoDT.sample <- function(n, sp) {
  support <- seq(0.1, 100, 0.1)
  m.par <- rec.parm
  probs <- twoDT.pdf(support, m.par["S",sp], m.par["L",sp])
  probs <- probs/sum(probs) # normalize probability 
  sample(
    support, 
    n, # n draws
    TRUE,  # with replacement
    probs # using these probabilities
  )
}
#hist(twoDT.sample(200, "Gironniera.nervosa"), breaks=100)

######################################
# FUNCTIONS FOR EXTRACTING VARIABLES #
######################################

## Conspecific adult density around seedlings (CAE)
# CAE is sum BA of all adult stems in 20x20 plot
# we want to convert this 400m^2 into a circle about each focal seedling
# 400 = pi*r^2
CAE.calc <- function(trees, seedlings, r = sqrt(400/pi)){
  # spDists returns distance matrix with rows as trees and cols as seedlings
  Ts.in.r <- ifelse(spDists(as.matrix(trees), as.matrix(seedlings)) < r, 1, 0)
  raw.CAE <- as.vector(pi*(exp(trees$logdbh)/2)^2) %*% Ts.in.r
  return(as.vector(raw.CAE))
}
#CAE.calc(ppoT, ppoS)

## Symmetric intraspecific competition (intraS)
# almost same as above, except includes seedlings
# returns intraS for adults [[1]] and seedlings [[2]]
intra.calc <- function(trees, seedlings, sp, r = sqrt(400/pi)){
  # convert seedling heights to dbh
  h <- seedlings$logheight
  z <- tran.parm["tran.int",sp] + h*tran.parm["tran.h",sp]
  # concatenate adult and seedling zs
  z <- c(trees$logdbh, z)
  
  # rbind coordinates of trees and seedlings
  ALL <- as.matrix(rbind(trees[,1:2], seedlings[,1:2]))
  # calc distance and find all within 400m2 radius of focal
  dists <- spDists(ALL)
  Ts.in.r <- ifelse(dists < sqrt(400/pi), 1, 0)
  diag(Ts.in.r) <- 0
  # for asymmetric competition, we only want pairs in which focal is the larger of the two, so need to compute another matrix
  zz <- matrix(rep(z, each=length(z)), ncol=length(z))
  big.mat <- ifelse(t(zz)>zz, 1, 0)
  pairs.to.count <- Ts.in.r * big.mat
  # use matrix multiplication to obtain sum intraA for each individual
  raw.intraA <- pi*(exp(z)/2)^2 %*% pairs.to.count
  return(list(as.vector(raw.intraA)[1:nrow(trees)], as.vector(raw.intraA)[(nrow(trees)+1):(nrow(trees)+nrow(seedlings))]))
}
#intra.calc(ppoT, ppoS, sp="Prunus.polystachya")
#intraA <- intra.calc(ppoT, ppoS, sp="Prunus.polystachya")[[2]]
#hist((log(intraS + surv.parm.unscale["surv.intraS.unscale.logoffset",]) -
#      surv.parm.unscale["surv.intraS.unscale.mu",]) / surv.parm.unscale["surv.intraS.unscale.sigma",])

## Distance-explicit asymmetric intraspecific competition (intraA.dist)
# almost same as above, except only for adult trees, and more spatially sensitive
intra.dist.calc <- function(trees, r=sqrt(1600/pi)){
  # calc distance
  dists <- spDists(trees)
  dbh <- exp(trees$logdbh)
  # diagonal elements and trees which are > sqrt(400/pi) m from focal should not be computed (note: 1600m2 is the area of a 40x40m plot)
  # to get rid simply set to a very large value
  dists[which(dists > r)] <- 1e10
  diag(dists) <- 1e10
  dists <- dists^-1
  raw.intraS.dist <- (dists %*% dbh)
  # intraA is simply intraS divided by dbh
  return(raw.intraS.dist/dbh)
}
#intra.dist.calc(ppoT)
#hist(intraA <- intra.dist.calc(ppoT))
#hist((log(intraA + grow.T.parm.unscale["grow.intraA.unscale.logoffset",]) -
#      grow.T.parm.unscale["grow.intraA.unscale.mu",]) / grow.T.parm.unscale["grow.intraA.unscale.sigma",])

## Heterospecific adult density around seedlings (HAE)
# HAE is sum BA of all non focal species adult stems in 20x20 plot
# we want to convert this 400m^2 into a circle about each focal seedling
# 400 = pi*r^2
HAE.calc <- function(trees.hetero, seedlings.con, r = sqrt(400/pi)){
  # spDists returns distance matrix with rows as trees and cols as seedlings
  Ts.in.r <- ifelse(spDists(as.matrix(trees.hetero), as.matrix(seedlings.con)) < r, 1, 0)
  raw.CAE <- as.vector(pi*(exp(trees$logdbh)/2)^2) %*% Ts.in.r
  return(as.vector(raw.CAE))
}
#CAE.calc(sceT, ppoS)

## Symmetric interspecific competition (interS)
# almost same as above, except includes seedlings
# returns interS for adults [[1]] and seedlings [[2]]
inter.calc <- function(trees.hetero, trees.con, seedlings.con, sp, r = sqrt(400/pi)){
  # convert seedling heights to dbh
  h <- seedlings.con$logheight
  z <- tran.parm["tran.int",sp] + h*tran.parm["tran.h",sp]
  
  # rbind coordinates of trees and seedlings
  con.ALL <- as.matrix(rbind(trees.con[,1:2], seedlings.con[,1:2]))
  # calc distance and find all within 400m2 radius of focal
  dists <- spDists(trees.hetero[,1:2], con.ALL)
  Ts.in.r <- ifelse(dists < sqrt(400/pi), 1, 0)
  trees.hetero.ba <- matrix(pi*(exp(trees.hetero$logdbh)/2)^2, nrow=1)
  raw.interS <- trees.hetero.ba %*% as.matrix(Ts.in.r)
  return(list(as.vector(raw.interS)[1:nrow(trees.con)], as.vector(raw.interS)[(nrow(trees.con)+1):(nrow(trees.con)+nrow(seedlings.con))]))
}
#inter.calc(sceT, ppoT, ppoS, sp="Prunus.polystachya")
#interS <- inter.calc(sceT, ppoT, ppoS, sp="Prunus.polystachya")[[2]]
#hist((log(interS) -
#      surv.parm.unscale["surv.interS.unscale.mu",]) / surv.parm.unscale["surv.interS.unscale.sigma",])

## Distance-explicit symmetric interspecific competition (intraS.dist)
# almost same as above, except only for adult trees, and more spatially sensitive
inter.dist.calc <- function(trees.hetero, trees.con, r=sqrt(1600/pi)){
  # calc distance
  dists <- spDists(trees.hetero, trees.con)
  dbh <- exp(trees.hetero$logdbh)
  # diagonal elements and trees which are > sqrt(400/pi) m from focal should not be computed (note: 1600m2 is the area of a 40x40m plot)
  # to get rid simply set to a very large value
  dists[which(dists > r)] <- 1e10
  diag(dists) <- 1e10
  dists <- dists^-1
  raw.interS.dist <- matrix(dbh, nrow=1) %*% dists
  return(raw.interS.dist)
}
#inter.dist.calc(sceT, ppoT)
#hist(interS <- inter.dist.calc(sceT, ppoT))
#hist((log(interS) -
#      grow.T.parm.unscale["grow.interS.unscale.mu",]) / grow.T.parm.unscale["grow.interS.unscale.sigma",])


#################
# SETUP BASEMAP #
#################

nssf.m <- raster("data/basemap 10m-res.tif")
nssf.m <- crop(nssf.m, extent(nssf.m, 220, 260, 120, 160))
plot(nssf.m)

library(poweRlaw)
## initialize PPO
init.ppo.n <- 400

# 20cm seedling height is the size of newly germinated seedling. calculate the theoretical DBH of this
ppo.minlogdbh <- tran.parm["tran.int","Prunus.polystachya"] + tran.parm["tran.h","Prunus.polystachya"] * log(20)
ppo.minDBH <- exp(ppo.minlogdbh)
# use this minimum and the power law to create the size distribution (in untransformed DBH)
init.ppo.dbh <- rplcon(init.ppo.n, ppo.minDBH, 2)

# subset the data into seedlings and adults
ppoS.logdbh <- log(init.ppo.dbh[which(init.ppo.dbh < 1)])
ppoT.logdbh <- log(init.ppo.dbh[which(init.ppo.dbh >= 1)])

# translate theoretical log DBH of seedlings into logheight (units used for all operations)
ppoS.logheight <- (ppoS.logdbh - tran.parm["tran.int","Prunus.polystachya"]) / tran.parm["tran.h","Prunus.polystachya"]

# create random points for ppoT, assign them logdbh values
ppoT.ras <- sampleRandom(nssf.m, size=length(ppoT.logdbh), na.rm=T, asRaster=T)
ppoT.ras[!is.na(ppoT.ras)] <- ppoT.logdbh
ppoT <- data.frame(rasterToPoints(ppoT.ras, spatial=F))
names(ppoT) <- c("x", "y", "logdbh")

# create random points for ppoS, assign them logheight values
ppoS.ras <- sampleRandom(nssf.m, size=length(ppoS.logheight), na.rm=T, asRaster=T)
ppoS.ras[!is.na(ppoS.ras)] <- ppoS.logheight
ppoS <- data.frame(rasterToPoints(ppoS.ras, spatial=F))
names(ppoS) <- c("x", "y", "logheight")

## initialize SCE
init.sce.n <- 500

# 10cm seedling height is the size of newly germinated seedling. calculate the theoretical DBH of this
sce.minlogdbh <- tran.parm["tran.int","Strombosia.ceylanica"] + tran.parm["tran.h","Strombosia.ceylanica"] * log(10)
sce.minDBH <- exp(sce.minlogdbh)
# use this minimum and the power law to create the size distribution (in untransformed DBH)
init.sce.dbh <- rplcon(init.sce.n, sce.minDBH, 2)

# subset the data into seedlings and adults
sceS.logdbh <- log(init.sce.dbh[which(init.sce.dbh < 1)])
sceT.logdbh <- log(init.sce.dbh[which(init.sce.dbh >= 1)])

# translate theoretical log DBH of seedlings into logheight (units used for all operations)
sceS.logheight <- (sceS.logdbh - tran.parm["tran.int","Strombosia.ceylanica"]) / tran.parm["tran.h","Strombosia.ceylanica"]

# create random points for sceT, assign them logdbh values
sceT.ras <- sampleRandom(nssf.m, size=length(sceT.logdbh), na.rm=T, asRaster=T)
sceT.ras[!is.na(sceT.ras)] <- sceT.logdbh
sceT <- data.frame(rasterToPoints(sceT.ras, spatial=F))
names(sceT) <- c("x", "y", "logdbh")

# create random points for sceS, assign them logheight values
sceS.ras <- sampleRandom(nssf.m, size=length(sceS.logheight), na.rm=T, asRaster=T)
sceS.ras[!is.na(sceS.ras)] <- sceS.logheight
sceS <- data.frame(rasterToPoints(sceS.ras, spatial=F))
names(sceS) <- c("x", "y", "logheight")

## initialize "All other species" for background interspecific competition
init.aos.n <- 2000

# 10cm seedling height is the size of newly germinated seedling. calculate the theoretical DBH of this
aos.minlogdbh <- tran.parm["tran.int","All.other.spp"] + tran.parm["tran.h","All.other.spp"] * log(10)
aos.minDBH <- exp(aos.minlogdbh)
# use this minimum and the power law to create the size distribution (in untransformed DBH)
init.aos.dbh <- rplcon(init.aos.n, aos.minDBH, 2)

# subset the data into seedlings and adults
aosS.logdbh <- log(init.aos.dbh[which(init.aos.dbh < 1)])
aosT.logdbh <- log(init.aos.dbh[which(init.aos.dbh >= 1)])

# translate theoretical log DBH of seedlings into logheight (units used for all operations)
aosS.logheight <- (aosS.logdbh - tran.parm["tran.int","All.other.spp"]) / tran.parm["tran.h","All.other.spp"]

# create random points for aosT, assign them logdbh values
aosT.ras <- sampleRandom(nssf.m, size=length(aosT.logdbh), na.rm=T, asRaster=T)
aosT.ras[!is.na(aosT.ras)] <- aosT.logdbh
aosT <- data.frame(rasterToPoints(aosT.ras, spatial=F))
names(aosT) <- c("x", "y", "logdbh")

# create random points for aosS, assign them logheight values
aosS.ras <- sampleRandom(nssf.m, size=length(aosS.logheight), na.rm=T, asRaster=T)
aosS.ras[!is.na(aosS.ras)] <- aosS.logheight
aosS <- data.frame(rasterToPoints(aosS.ras, spatial=F))
names(aosS) <- c("x", "y", "logheight")



# Define a colour palette
col.pal <- c("#E3C16F", "#946846", "#FAFF70", "#6D213C", "#65DEF1", "#F3F5F6", "#BAAB68", "#CBC5EA")
col.t <- adjustcolor(col.pal, alpha.f=0.6)

plot(nssf.m, legend=F, col=col.pal[c(6,5)])
# PPO
points(ppoT, cex=ppoT$logdbh, col=col.t[1], pch=16)
points(data.frame(ppoS), col=col.t[1], pch=4, cex=0.1)
# SCE
points(sceT, cex=sceT$logdbh, col=col.t[2], pch=16)
points(data.frame(sceS), col=col.t[2], pch=4, cex=0.1)
# All other species (static--just for a background interspecific competition pressure)
points(aosT, cex=aosT$logdbh,  pch=1)
scalebar(100, xy=c(367200, 152900), type="bar", lonlat=F, below="metres", divs=4)

#######################
# COMMENCE SIMULATION #
#######################

ppoT.init <- ppoT
ppoS.init <- ppoS

sceT.init <- sceT
sceS.init <- sceS

time <- 0
maxTime <- 100

(n.ppoT <- nrow(ppoT))
(n.ppoS <- nrow(ppoS))

(n.sceT <- nrow(sceT))
(n.sceS <- nrow(sceS))

(z.ppoT <- mean(ppoT$logdbh))
(h.ppoS <- mean(ppoS$logheight))

(z.sceT <- mean(sceT$logdbh))
(h.sceS <- mean(sceS$logheight))

while (time < maxTime) {
  
  # Recruitment: initiate fruiting, but don't add new recruits to ppoS yet (let them grow/die first)
  fruiting.index <- which(p_bz(nssf.m, ppoT, "Prunus.polystachya")==1)
  prod.vec <- b_z(nssf.m, ppoT[fruiting.index,], "Prunus.polystachya")
  parent.loc <- ppoT[fruiting.index, 1:2]
  
  # Seedling growth and survival
  seedling.CAE <- CAE.calc(ppoT, ppoS)
  intraS <- intra.calc(ppoT, ppoS, "Prunus.polystachya")
  ppoS$logheight <- sS_h(nssf.m, ppoS, "Prunus.polystachya", intraS[[2]]) * GS_h1h(nssf.m, ppoS, "Prunus.polystachya", seedling.CAE)
  if(sum(ppoS$logheight==0) > 0)  ppoS <- ppoS[-which(ppoS$logheight==0),]
  rm(seedling.CAE)
  
  # Tree growth and survival
  intraA <- intra.dist.calc(ppoT)
  ppoT$logdbh <- sT_z(nssf.m, ppoT, "Prunus.polystachya", intraS[[1]]) * GT_z1z(nssf.m, ppoT, "Prunus.polystachya", intraA)
  if(sum(ppoT$logdbh==0) > 0)  ppoT <- ppoT[-which(ppoT$logdbh==0),]
  rm(intraA)
  
  # Seedling-sapling transition
  ppoTS <- T_z1h1(ppoT, ppoS, "Prunus.polystachya")
  ppoT <- ppoTS[[1]]
  ppoS <- ppoTS[[2]]
  rm(intraS)
  
  # This is the current number of seedlings:
  n.old.ppoS <- nrow(ppoS)
  
  # Recruitment: now add new recruits to ppoS
  if(length(prod.vec) > 0){
    for(i in 1:length(prod.vec)){
      distances <- twoDT.sample(prod.vec[i], "Prunus.polystachya")
      x <- parent.loc[i,1] + (distances * cos(runif(length(distances), min=1, max=360)))
      y <- parent.loc[i,2] + (distances * sin(runif(length(distances), min=1, max=360)))
      logheight <- c_0h1(length(distances), "Prunus.polystachya")
      ppoS <- rbind(ppoS, cbind(x, y, logheight))
    }
  }
  # Kill off recruits that are located too close (< 10cm) to each other
  inter.rec.dists <- spDists(ppoS[(n.old.ppoS+1):nrow(ppoS),])
  diag(inter.rec.dists) <- NA
  clustered.recs <- match(names(which(apply(inter.rec.dists, 1, min, na.rm=T) < 0.2)), rownames(ppoS))
  rm(inter.rec.dists)
  if(length(clustered.recs >0)){
    dying.recs <- sample(clustered.recs, size=round(length(clustered.recs)*0.5,0), replace=F)
    ppoS <- ppoS[-dying.recs,]
  }
  rm(clustered.recs)
  
  # Take stock of all individuals
  n1.ppoT <- nrow(ppoT)
  n.ppoT <- c(n.ppoT, n1.ppoT) 
  n1.ppoS <- nrow(ppoS)
  n.ppoS <- c(n.ppoS, n1.ppoS) 
  
  # calculate mean sizes of adults and saplings and add to string  
  z.ppoT <- c(z.ppoT, mean(ppoT$logdbh))
  h.ppoS <- c(h.ppoS, mean(ppoS$logheight))
  
  # move to next time point
  time <- time + 1
  print(time) # slow the model
}

par(mfrow=c(1,2), mar=c(4,4,2,2))
plot(nssf.m, main="Initial condition", legend=F)
points(ppoT.init, cex=ppoT.init$logdbh)
points(data.frame(ppoS.init), col="grey", pch=4, cex=0.1)
legend('bottomleft', bg="white", legend=c("Adult tree", "Seedling"), 
       pch=c(1,4), pt.cex=c(3,0.1), col=c("black","grey"))
plot(nssf.m, main=paste0("After ", time, " years"), legend=F)
points(data.frame(ppoS), col="grey", pch=4, cex=0.1)
points(ppoT, cex=ppoT$logdbh)
scalebar(100, xy=c(367200, 152900), type="bar", lonlat=F, below="metres", divs=4)


par(mfrow=c(1,2), mar=c(5.5,5.5,2,2))
plot(n.ppoT ~ c(1:length(n.ppoT)), lwd=5, col="forestgreen", type="l",
     ylab="Adult population size", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5)
plot(n.ppoS ~ c(1:length(n.ppoS)), lwd=5, col="brown", type="l",
     ylab="Seedling population size", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5)

par(mfrow=c(1,2))
plot(z.ppoT ~ c(1:length(z.ppoT)), lwd=5, col="forestgreen", type="l",
     ylab="Mean adult DBH (log-transformed)", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5)
plot(h.ppoS ~ c(1:length(h.ppoS)), lwd=5, col="brown", type="l",
     ylab="Mean seedling height (log-transformed)", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5)

par(mfrow=c(1,2))
hist(ppoS$logheight)
hist(ppoT$logdbh)

