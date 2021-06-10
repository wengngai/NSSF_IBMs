library(raster)
library(sp)

#################
# GROWTH PARAMS #
#################

(grow.T.parm <- t(read.csv("data/tree growth parameters Apr21.csv", header=T, row.names=1)))
(grow.T.parm.unscale <- read.csv("data/tree growth model unscale params.csv", header=T, row.names=1))
(grow.S.parm <- read.csv("data/seedling lmer growth params CAA Apr21.csv", header=T, row.names=1))
(surv.parm <- t(read.csv("data/surv params Jun21.csv", header=T, row.names=1)))
(surv.parm.unscale <- read.csv("data/survival model unscale params.csv", header=T, row.names=1))
(fruit.parm <- read.csv("data/fruiting parameters Apr21.csv", header=T, row.names=1))
(tran.parm <- read.csv("data/transition params CAA Jun21.csv", header=T, row.names=1))
(rec.parm <- t(read.csv("data/dispersal kernel parameters Apr21.csv", header=T, row.names=1)))
(HM.parm <- t(read.csv("data/HM values.csv", header=T, row.names=1)))

# standardize some spp names
colnames(surv.parm) <- gsub(" ", ".", colnames(surv.parm))
colnames(grow.T.parm) <- gsub(" ", ".", colnames(grow.T.parm))
colnames(HM.parm) <- gsub(" ", ".", colnames(HM.parm))



##################
# LOAD FUNCTIONS #
##################
source("code/functions.R")



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
init.sce.n <- 600

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

# Use 60cm min seedling height for all other spp. (just need competition effect, and young seedlings contribute little)
aos.minlogdbh <- tran.parm["tran.int","All.other.spp"] + tran.parm["tran.h","All.other.spp"] * log(60)
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
  # PPO ("Prunus.polystachya")
  fruiting.index.ppo <- which(p_bz(nssf.m, ppoT, "Prunus.polystachya")==1)
  prod.vec.ppo <- b_z(nssf.m, ppoT[fruiting.index.ppo,], "Prunus.polystachya")
  parent.loc.ppo <- ppoT[fruiting.index.ppo, 1:2]
  # SCE ("Strombosia.ceylanica")
  fruiting.index.sce <- which(p_bz(nssf.m, sceT, "Strombosia.ceylannica")==1)
  prod.vec.sce <- b_z(nssf.m, sceT[fruiting.index.sce,], "Strombosia.ceylannica")
  parent.loc.sce <- sceT[fruiting.index.sce, 1:2]
  
  # Extract competition measures
  # PPO ("Prunus.polystachya")
  inter.on.PPO <- inter.calc(rbind(sceT, aosT), rbind(sceS, aosS), ppoT, ppoS, sp="Prunus.polystachya")
  intra.on.PPO <- intra.calc(ppoT, ppoS, sp="Prunus.polystachya")
  CAE.on.ppoS <- CAE.calc(ppoT, ppoS)
  HAE.on.ppoS <- HAE.calc(rbind(sceT, aosT), ppoS)
  intra.dist.on.PPO <- intra.dist.calc(ppoT)
  inter.dist.on.PPO <- inter.dist.calc(rbind(sceT, aosT), ppoT)
  # SCE ("Strombosia.ceylanica")
  inter.on.SCE <- inter.calc(rbind(ppoT, aosT), rbind(ppoS, aosS), sceT, sceS, sp="Strombosia.ceylanica")
  intra.on.SCE <- intra.calc(sceT, sceS, sp="Strombosia.ceylanica")
  CAE.on.sceS <- CAE.calc(sceT, sceS)
  HAE.on.sceS <- HAE.calc(rbind(ppoT, aosT), sceS)
  intra.dist.on.SCE <- intra.dist.calc(sceT)
  inter.dist.on.SCE <- inter.dist.calc(rbind(ppoT, aosT), sceT)
  
  # Seedling growth and survival
  # PPO ("Prunus.polystachya")
  ppoS$logheight <- sS_h(nssf.m, ppoS, "Prunus.polystachya", intra.on.PPO[[4]], inter.on.PPO[[4]]) * 
    GS_h1h(nssf.m, ppoS, "Prunus.polystachya", CAE.on.ppoS, HAE.on.ppoS)
  if(sum(ppoS$logheight==0) > 0)  ppoS <- ppoS[-which(ppoS$logheight==0),]
  # SCE ("Strombosia.ceylanica")
  sceS$logheight <- sS_h(nssf.m, sceS, "Strombosia.ceylanica", intra.on.SCE[[4]], inter.on.SCE[[4]]) * 
    GS_h1h(nssf.m, sceS, "Strombosia.ceylanica", CAE.on.sceS, HAE.on.sceS)
  if(sum(sceS$logheight==0) > 0)  sceS <- sceS[-which(sceS$logheight==0),]
  
  # Tree growth and survival
  # PPO ("Prunus.polystachya")
  ppoT$logdbh <- sT_z(nssf.m, ppoT, "Prunus.polystachya", intra.on.PPO[[1]], intra.on.PPO[[3]], inter.on.PPO[[1]], inter.on.PPO[[3]]) * 
    GT_z1z(nssf.m, ppoT, "Prunus.polystachya", intra.dist.on.PPO, inter.dist.on.PPO)
  if(sum(ppoT$logdbh==0) > 0)  ppoT <- ppoT[-which(ppoT$logdbh==0),]
  # SCE ("Strombosia.ceylanica")
  sceT$logdbh <- sT_z(nssf.m, sceT, "Strombosia.ceylanica", intra.on.SCE[[1]], intra.on.SCE[[3]], inter.on.SCE[[1]], inter.on.SCE[[3]]) * 
    GT_z1z(nssf.m, sceT, "Strombosia.ceylanica", intra.dist.on.SCE, inter.dist.on.SCE)
  if(sum(sceT$logdbh==0) > 0)  sceT <- sceT[-which(sceT$logdbh==0),]

  # Conserve memory by deleting competition indices (don't need anymore)
  #rm()
  
  # Seedling-sapling transition
  # PPO ("Prunus.polystachya")
  ppoTS <- T_z1h1(ppoT, ppoS, "Prunus.polystachya")
  ppoT <- ppoTS[[1]]
  ppoS <- ppoTS[[2]]
  # SCE ("Strombosia.ceylanica")
  sceTS <- T_z1h1(sceT, sceS, "Strombosia.ceylanica")
  sceT <- sceTS[[1]]
  sceS <- sceTS[[2]]
  
  # This is the current number of seedlings:
  # PPO ("Prunus.polystachya")
  n.old.ppoS <- nrow(ppoS)
  # SCE ("Strombosia.ceylanica")
  n.old.sceS <- nrow(sceS)
  
  # Recruitment: now add new recruits to ppoS
  # PPO ("Prunus.polystachya")
  if(length(prod.vec.ppo) > 0){
    for(i in 1:length(prod.vec.ppo)){
      distances <- twoDT.sample(prod.vec.ppo[i], "Prunus.polystachya")
      angles <- runif(length(distances), min=0, max=360)
      x <- parent.loc.ppo[i,1] + (distances * cos(angles))
      y <- parent.loc.ppo[i,2] + (distances * sin(angles))
      logheight <- c_0h1(length(distances), "Prunus.polystachya")
      ppoS <- rbind(ppoS, cbind(x, y, logheight))
    }
  }
  # SCE ("Strombosia.ceylanica")
  if(length(prod.vec.sce) > 0){
    for(i in 1:length(prod.vec.sce)){
      distances <- twoDT.sample(prod.vec.sce[i], "Strombosia.ceylanica")
      angles <- runif(length(distances), min=0, max=360)
      x <- parent.loc.sce[i,1] + (distances * cos(angles))
      y <- parent.loc.sce[i,2] + (distances * sin(angles))
      logheight <- c_0h1(length(distances), "Strombosia.ceylanica")
      sceS <- rbind(sceS, cbind(x, y, logheight))
    }
  }
  
  # Kill off 50% of all recruits that are located < 20cm from each other
  # PPO ("Prunus.polystachya")
  inter.rec.dists <- spDists(ppoS[(n.old.ppoS+1):nrow(ppoS),])
  diag(inter.rec.dists) <- NA
  clustered.recs <- match(names(which(apply(inter.rec.dists, 1, min, na.rm=T) < 0.2)), rownames(ppoS))
  if(length(clustered.recs >0)){
    dying.recs <- sample(clustered.recs, size=round(length(clustered.recs)*0.5,0), replace=F)
    ppoS <- ppoS[-dying.recs,]
  }
  # SCE ("Strombosia.ceylanica")
  inter.rec.dists <- spDists(sceS[(n.old.sceS+1):nrow(sceS),])
  diag(inter.rec.dists) <- NA
  clustered.recs <- match(names(which(apply(inter.rec.dists, 1, min, na.rm=T) < 0.2)), rownames(ppoS))
  if(length(clustered.recs >0)){
    dying.recs <- sample(clustered.recs, size=round(length(clustered.recs)*0.5,0), replace=F)
    sceS <- sceS[-dying.recs,]
  }
  
  # Take stock of all individuals
  # PPO ("Prunus.polystachya")
  n1.ppoT <- nrow(ppoT)
  n.ppoT <- c(n.ppoT, n1.ppoT) 
  n1.ppoS <- nrow(ppoS)
  n.ppoS <- c(n.ppoS, n1.ppoS) 
  # SCE ("Strombosia.ceylanica")
  n1.sceT <- nrow(sceT)
  n.sceT <- c(n.sceT, n1.sceT) 
  n1.sceS <- nrow(sceS)
  n.sceS <- c(n.sceS, n1.sceS) 
  
  # calculate mean sizes of adults and saplings and add to string  
  # PPO ("Prunus.polystachya")
  z.ppoT <- c(z.ppoT, mean(ppoT$logdbh))
  h.ppoS <- c(h.ppoS, mean(ppoS$logheight))
  # SCE ("Strombosia.ceylanica")
  z.sceT <- c(z.sceT, mean(sceT$logdbh))
  h.sceS <- c(h.sceS, mean(sceS$logheight))
  
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

