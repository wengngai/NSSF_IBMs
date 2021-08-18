library(raster)
library(sp)
library(poweRlaw)
library(doParallel)
library(rgdal)

# Settings for parallerisation
registerDoParallel(cores = 8)


#################
# GROWTH PARAMS #
#################

grow.T.parm <- t(read.csv("data/tree growth parameters Aug21.csv", header=T, row.names=1))
grow.T.parm.unscale <- read.csv("data/tree growth model unscale params.csv", header=T, row.names=1)
grow.S.parm <- read.csv("data/seedling lmer growth params CAA Apr21.csv", header=T, row.names=1)
surv.parm <- t(read.csv("data/surv params Aug21.csv", header=T, row.names=1))
surv.parm.unscale <- read.csv("data/survival model unscale params Aug21.csv", header=T, row.names=1)
fruit.parm <- read.csv("data/fruiting parameters Apr21.csv", header=T, row.names=1)
tran.parm <- read.csv("data/transition params CAA Jun21.csv", header=T, row.names=1)
rec.parm <- t(read.csv("data/dispersal kernel parameters Apr21.csv", header=T, row.names=1))
HM.parm <- t(read.csv("data/HM values.csv", header=T, row.names=1))

# standardize some spp names
colnames(surv.parm) <- gsub(" ", ".", colnames(surv.parm))
colnames(grow.T.parm) <- gsub(" ", ".", colnames(grow.T.parm))
colnames(HM.parm) <- gsub(" ", ".", colnames(HM.parm))

# add a "min size" row to rec.parm
rec.parm <- rbind(rec.parm, rec.parm["rcsz",] - 2*rec.parm["rcsd",])
rownames(rec.parm)[5] <- "minsize"

# add "population" parameters to HM.parm
HM.parm <- data.frame(HM.parm, population=c(0.5,0.5,0.5,0.5,0.5,0.5))

##################
# LOAD FUNCTIONS #
##################
source("code/functions.R")


#################
# SETUP BASEMAP #
#################

# set scenario (usual or extreme)
scenario <- "usual"
#scenario <- "extreme"

source("code/setup_map_three spp.R")


#######################
# COMMENCE SIMULATION #
#######################

time <- 1
maxTime <- 10

# take a snapshot of initial conditions
ppoT.init <- ppoT
ppoS.init <- ppoS
sceT.init <- sceT
sceS.init <- sceS
ppiT.init <- ppiT
ppiS.init <- ppiS
aosT.init <- aosT

# record population sizes
n.ppoT <- nrow(ppoT)
n.ppoS <- nrow(ppoS)
n.sceT <- nrow(sceT)
n.sceS <- nrow(sceS)
n.ppiT <- nrow(ppiT)
n.ppiS <- nrow(ppiS)
n.aosT <- nrow(aosT)

# record mean sizes
z.ppoT <- mean(ppoT$logdbh)
h.ppoS <- mean(ppoS$logheight)
z.sceT <- mean(sceT$logdbh)
h.sceS <- mean(sceS$logheight)
z.ppiT <- mean(ppiT$logdbh)
h.ppiS <- mean(ppiS$logheight)
z.aosT <- mean(aosT$logdbh)

# create table for total BA in swamp vs non-swamp
ba.ppoT <- data.frame(time=seq(time, maxTime, 1), nonswamp=NA, swamp=NA, landscape.swamp.prop=NA)
ba.sceT <- data.frame(time=seq(time, maxTime, 1), nonswamp=NA, swamp=NA, landscape.swamp.prop=NA)
ba.ppiT <- data.frame(time=seq(time, maxTime, 1), nonswamp=NA, swamp=NA, landscape.swamp.prop=NA)
ba.aosT <- data.frame(time=seq(time, maxTime, 1), nonswamp=NA, swamp=NA, landscape.swamp.prop=NA)

# create vectors to record deaths, recruitment
deaths.ppoT <- c()
deaths.ppoS <- c()
deaths.sceT <- c()
deaths.sceS <- c()
deaths.ppiT <- c()
deaths.ppiS <- c()
deaths.aosT <- c()

recs.ppoS <- c()
recs.sceS <- c()
recs.ppiS <- c()
recs.aosT <- c()

# progress bar
pb <- txtProgressBar(min = 1, max = maxTime, style = 3)

while (time <= maxTime) {
  
  # Choose landscape based on scenario and time (beyond 22 years, just use the last time point (=2042))
  if(scenario=="usual")  nssf.m <- nssf.usual[[ifelse(time>22, 22, time)]]
  if(scenario=="extreme")  nssf.m <- nssf.extreme[[ifelse(time>22, 22, time)]]
  
  # Recruitment: initiate fruiting, but don't add new recruits to seedling dfs yet (let them grow/die first)
  # PPO ("Prunus.polystachya")
  fruiting.index.ppo <- which(p_bz(nssf.m, ppoT, "Prunus.polystachya")==1)
  prod.vec.ppo <- b_z(nssf.m, ppoT[fruiting.index.ppo,], "Prunus.polystachya")
  parent.loc.ppo <- ppoT[fruiting.index.ppo, 1:2]
  # SCE ("Strombosia.ceylanica")
  fruiting.index.sce <- which(p_bz(nssf.m, sceT, "Strombosia.ceylanica")==1)
  prod.vec.sce <- b_z(nssf.m, sceT[fruiting.index.sce,], "Strombosia.ceylanica")
  parent.loc.sce <- sceT[fruiting.index.sce, 1:2]
  # PPI ("Pometia pinnata")
  fruiting.index.ppi <- which(p_bz(nssf.m, ppiT, "Pometia.pinnata")==1)
  prod.vec.ppi <- b_z(nssf.m, ppiT[fruiting.index.ppi,], "Pometia.pinnata")
  parent.loc.ppi <- ppiT[fruiting.index.ppi, 1:2]
  
  # Extract competition measures
  # PPO ("Prunus.polystachya")
  inter.on.PPO <- inter.calc(rbind(sceT, ppiT, aosT), rbind(sceS, ppiS), ppoT, ppoS, 
                             sp="Prunus.polystachya", grid.neighbours = grid.neighbours)
  intra.on.PPO <- intra.calc(ppoT, ppoS, sp="Prunus.polystachya", grid.neighbours = grid.neighbours)
  CAE.on.ppoS <- CAE.calc(ppoT, ppoS, grid.neighbours = grid.neighbours)
  HAE.on.ppoS <- HAE.calc(rbind(sceT, ppiT, aosT), ppoS, grid.neighbours = grid.neighbours)
  intra.dist.on.PPO <- intra.dist.calc(ppoT, grid.neighbours = grid.neighbours)
  inter.dist.on.PPO <- inter.dist.calc(rbind(sceT, ppiT, aosT), ppoT, grid.neighbours = grid.neighbours)
  # SCE ("Strombosia.ceylanica")
  inter.on.SCE <- inter.calc(rbind(ppoT, ppiT, aosT), rbind(ppoS, ppiS), sceT, sceS, 
                             sp="Strombosia.ceylanica", grid.neighbours = grid.neighbours)
  intra.on.SCE <- intra.calc(sceT, sceS, sp="Strombosia.ceylanica", grid.neighbours = grid.neighbours)
  CAE.on.sceS <- CAE.calc(sceT, sceS, grid.neighbours = grid.neighbours)
  HAE.on.sceS <- HAE.calc(rbind(ppoT, ppiT, aosT), sceS, grid.neighbours = grid.neighbours)
  intra.dist.on.SCE <- intra.dist.calc(sceT, grid.neighbours = grid.neighbours)
  inter.dist.on.SCE <- inter.dist.calc(rbind(ppoT, ppiT, aosT), sceT, grid.neighbours = grid.neighbours)
  # PPI ("Pometia.pinnata")
  inter.on.PPI <- inter.calc(rbind(sceT, ppoT, aosT), rbind(sceS, ppoS), ppiT, ppiS, 
                             sp="Pometia.pinnata", grid.neighbours = grid.neighbours)
  intra.on.PPI <- intra.calc(ppiT, ppiS, sp="Pometia.pinnata", grid.neighbours = grid.neighbours)
  CAE.on.ppiS <- CAE.calc(ppiT, ppiS, grid.neighbours = grid.neighbours)
  HAE.on.ppiS <- HAE.calc(rbind(sceT, ppoT, aosT), ppiS, grid.neighbours = grid.neighbours)
  intra.dist.on.PPI <- intra.dist.calc(ppiT, grid.neighbours = grid.neighbours)
  inter.dist.on.PPI <- inter.dist.calc(rbind(sceT, ppoT, aosT), ppiT, grid.neighbours = grid.neighbours)
  # AOS ("population")
  inter.on.AOS <- inter.calc(rbind(ppoT, sceT, ppiT), rbind(ppoS, sceS, ppiS), aosT, aosS, 
                             sp="All.other.spp", grid.neighbours = grid.neighbours)
  #intra.on.AOS <- intra.calc(aosT, aosS, sp="population", grid.neighbours = grid.neighbours)
  #intra.dist.on.AOS <- intra.dist.calc(aosT, grid.neighbours = grid.neighbours)
  inter.dist.on.AOS <- inter.dist.calc(rbind(ppoT, sceT, ppiT), aosT, grid.neighbours = grid.neighbours)
  
  # Seedling growth and survival
  # PPO ("Prunus.polystachya")
  #ppoS$logheight <- sS_h(nssf.m, ppoS, "Prunus.polystachya", intra.on.PPO[[4]], inter.on.PPO[[4]]) * 
  #test using mean inter and intra values to see if these are causing the problem
  ppoS$logheight <- sS_h(nssf.m, ppoS, "Prunus.polystachya", exp(surv.parm.unscale["intraA.unscale.mu",]), exp(surv.parm.unscale["interAHM.unscale.mu",])) * 
    GS_h1h(nssf.m, ppoS, "Prunus.polystachya", CAE.on.ppoS, HAE.on.ppoS)
  deaths.ppoS <- c(deaths.ppoS, sum(ppoS$logheight < rec.parm["minsize", "Prunus.polystachya"]) )
  if(sum(ppoS$logheight < rec.parm["minsize", "Prunus.polystachya"], na.rm = TRUE) > 0){
    ppoS <- ppoS[-which(ppoS$logheight < rec.parm["minsize", "Prunus.polystachya"]),] }
  # SCE ("Strombosia.ceylanica")
  #sceS$logheight <- sS_h(nssf.m, sceS, "Strombosia.ceylanica", intra.on.SCE[[4]], inter.on.SCE[[4]]) *
  #test using mean inter and intra values to see if these are causing the problem
  sceS$logheight <- sS_h(nssf.m, sceS, "Strombosia.ceylanica", exp(surv.parm.unscale["intraA.unscale.mu",]), exp(surv.parm.unscale["interAHM.unscale.mu",])) * 
    GS_h1h(nssf.m, sceS, "Strombosia.ceylanica", CAE.on.sceS, HAE.on.sceS)
  deaths.sceS <- c(deaths.sceS, sum(sceS$logheight < rec.parm["minsize", "Strombosia.ceylanica"]) )
  if(sum(sceS$logheight < rec.parm["minsize", "Strombosia.ceylanica"], na.rm = TRUE) > 0){
    sceS <- sceS[-which(sceS$logheight < rec.parm["minsize", "Strombosia.ceylanica"]),] }
  # PPI ("Pometia.pinnata")
  #ppiS$logheight <- sS_h(nssf.m, ppiS, "Pometia.pinnata", intra.on.PPI[[4]], inter.on.PPI[[4]]) * 
  #test using mean inter and intra values to see if these are causing the problem
  ppiS$logheight <- sS_h(nssf.m, ppiS, "Pometia.pinnata", exp(surv.parm.unscale["intraA.unscale.mu",]), exp(surv.parm.unscale["interAHM.unscale.mu",])) *
    GS_h1h(nssf.m, ppiS, "Pometia.pinnata", CAE.on.ppiS, HAE.on.ppiS)
  deaths.ppiS <- c(deaths.ppiS, sum(ppiS$logheight < rec.parm["minsize", "Pometia.pinnata"]) )
  if(sum(ppiS$logheight < rec.parm["minsize", "Pometia.pinnata"], na.rm = TRUE) > 0){
    ppiS <- ppiS[-which(ppiS$logheight < rec.parm["minsize", "Pometia.pinnata"]),] }
  
  # Tree growth and survival
  # PPO ("Prunus.polystachya")
  ppoT$logdbh <- sT_z(nssf.m, ppoT, "Prunus.polystachya", intra.on.PPO[[1]], intra.on.PPO[[3]], inter.on.PPO[[1]], inter.on.PPO[[3]]) * 
    GT_z1z(nssf.m, ppoT, "Prunus.polystachya", intra.dist.on.PPO, inter.dist.on.PPO)
  deaths.ppoT <- c(deaths.ppoT, sum(ppoT$logdbh==0, na.rm = TRUE))
  if(sum(ppoT$logdbh==0, na.rm = TRUE) > 0)  ppoT <- ppoT[-which(ppoT$logdbh==0),]
  # SCE ("Strombosia.ceylanica")
  sceT$logdbh <- sT_z(nssf.m, sceT, "Strombosia.ceylanica", intra.on.SCE[[1]], intra.on.SCE[[3]], inter.on.SCE[[1]], inter.on.SCE[[3]]) * 
    GT_z1z(nssf.m, sceT, "Strombosia.ceylanica", intra.dist.on.SCE, inter.dist.on.SCE)
  deaths.sceT <- c(deaths.sceT, sum(sceT$logdbh==0, na.rm = TRUE))
  if(sum(sceT$logdbh==0, na.rm = TRUE) > 0)  sceT <- sceT[-which(sceT$logdbh==0),]
  # PPI ("Pometia.pinnata")
  ppiT$logdbh <- sT_z(nssf.m, ppiT, "Pometia.pinnata", intra.on.PPI[[1]], intra.on.PPI[[3]], inter.on.PPI[[1]], inter.on.PPI[[3]]) * 
    GT_z1z(nssf.m, ppiT, "Pometia.pinnata", intra.dist.on.PPI, inter.dist.on.PPI)
  deaths.ppiT <- c(deaths.ppiT, sum(ppiT$logdbh==0, na.rm = TRUE))
  if(sum(ppiT$logdbh==0, na.rm = TRUE) > 0)  ppiT <- ppiT[-which(ppiT$logdbh==0),]
  
  # AOS growth, survival and recruitment
  # 0.003544643 AOS stems recruit per m2 per year (rate excludes the rate of the 3 focal spp.)
  # use a fixed, small value for intra
  aosT$logdbh <- sT_z(nssf.m, aosT, "population", 100, 100, inter.on.AOS[[1]], inter.on.AOS[[3]]) * 
    GT_z1z(nssf.m, aosT, "population", 100, inter.dist.on.AOS)
  deaths.aosT <- c(deaths.aosT, sum(aosT$logdbh==0, na.rm = TRUE))
  if(sum(aosT$logdbh==0, na.rm = TRUE) > 0)  aosT <- aosT[-which(aosT$logdbh==0),]
  n.aos.rec <- round(crop.area * rnorm(n=1, mean=0.003544643, sd=0.0005), 0) / dilution
  recs.aosT <- c(recs.aosT, n.aos.rec)
  aos.loc.rec <- spsample(crop.poly, n.aos.rec, type = 'random')
  aos.grid.rec <- extract(nssf.100, aos.loc.rec)
  aosT <- rbind(aosT, data.frame(coordinates(aos.loc.rec), logdbh = log(5), grid = aos.grid.rec))
  
  # Conserve memory by deleting competition indices (don't need anymore)
  rm(inter.on.PPO, intra.on.PPO, CAE.on.ppoS, HAE.on.ppoS, intra.dist.on.PPO, inter.dist.on.PPO,
     inter.on.SCE, intra.on.SCE, CAE.on.sceS, HAE.on.sceS, intra.dist.on.SCE, inter.dist.on.SCE,
     inter.on.PPI, intra.on.PPI, CAE.on.ppiS, HAE.on.ppiS, intra.dist.on.PPI, inter.dist.on.PPI,
     inter.on.AOS, inter.dist.on.AOS)
  
  # Seedling-sapling transition
  # PPO ("Prunus.polystachya")
  ppoTS <- T_z1h1(ppoT, ppoS, "Prunus.polystachya")
  ppoT <- ppoTS[[1]]
  ppoS <- ppoTS[[2]]
  rm(ppoTS)
  # SCE ("Strombosia.ceylanica")
  sceTS <- T_z1h1(sceT, sceS, "Strombosia.ceylanica")
  sceT <- sceTS[[1]]
  sceS <- sceTS[[2]]
  rm(sceTS)
  # PPI ("Pometia.pinnata")
  ppiTS <- T_z1h1(ppiT, ppiS, "Pometia.pinnata")
  ppiT <- ppiTS[[1]]
  ppiS <- ppiTS[[2]]
  rm(ppiTS)
  
  # This is the current number of seedlings:
  n.old.ppoS <- nrow(ppoS)
  n.old.sceS <- nrow(sceS)
  n.old.ppiS <- nrow(ppiS)
  
  # Recruitment: now add new recruits to ppoS
  # PPO ("Prunus.polystachya")
  if(length(prod.vec.ppo) > 0){
    for(i in 1:length(prod.vec.ppo)){
      distances <- twoDT.sample(prod.vec.ppo[i], "Prunus.polystachya", max_dist = 410)
      angles <- runif(length(distances), min=0, max=360)
      x <- parent.loc.ppo[i,1] + (distances * cos(angles))
      y <- parent.loc.ppo[i,2] + (distances * sin(angles))
      grid <- extract(nssf.100, cbind(x,y))
      logheight <- c_0h1(length(distances), "Prunus.polystachya")
      recruits <- cbind(x, y, logheight, grid)
      # use overlap function to filter out recruits which fail because they land on spots already occupied by trees/seedlings
      if (nrow(recruits) > 0) {
        recruits.success <- rm.overlap(recruits=recruits, 
                                       trees=rbind(ppoT, sceT, ppiT, aosT), 
                                       seedlings=rbind(ppoS, sceS, ppiS),
                                       grid.neighbours=grid.neighbours)
        # if there are still successfully recruiting seedlings, rbind them to existing seedlings
        if(nrow(recruits.success) > 0) {
          ppoS <- rbind(ppoS, recruits.success)
          recs.ppoS <- c(recs.ppoS, nrow(recruits.success))
        }
      }
      if(nrow(recruits)==0 | nrow(recruits.success)==0) recs.ppoS <- c(recs.ppoS, 0)
    }
  }
  # SCE ("Strombosia.ceylanica")
  if(length(prod.vec.sce) > 0){
    for(i in 1:length(prod.vec.sce)){
      distances <- twoDT.sample(prod.vec.sce[i], "Strombosia.ceylanica", max_dist = 410)
      angles <- runif(length(distances), min=0, max=360)
      x <- parent.loc.sce[i,1] + (distances * cos(angles))
      y <- parent.loc.sce[i,2] + (distances * sin(angles))
      grid <- extract(nssf.100, cbind(x,y))
      logheight <- c_0h1(length(distances), "Strombosia.ceylanica")
      recruits <- cbind(x, y, logheight, grid)
      if (nrow(recruits) > 0) {
        recruits.success <- rm.overlap(recruits, 
                                       trees=rbind(ppoT, sceT, ppiT, aosT), 
                                       seedlings=rbind(ppoS, sceS, ppiS),
                                       grid.neighbours=grid.neighbours)
        if(nrow(recruits.success) > 0) {
          sceS <- rbind(sceS, recruits.success)
          recs.sceS <- c(recs.sceS, nrow(recruits.success))
        } 
      }
      if(nrow(recruits)==0 | nrow(recruits.success)==0) recs.sceS <- c(recs.sceS, 0)
    }
  }
  # PPI ("Pometia.pinnata")
  if(length(prod.vec.ppi) > 0){
    for(i in 1:length(prod.vec.ppi)){
      distances <- twoDT.sample(prod.vec.ppi[i], "Pometia.pinnata", max_dist = 410)
      angles <- runif(length(distances), min=0, max=360)
      x <- parent.loc.ppi[i,1] + (distances * cos(angles))
      y <- parent.loc.ppi[i,2] + (distances * sin(angles))
      grid <- extract(nssf.100, cbind(x,y))
      logheight <- c_0h1(length(distances), "Pometia.pinnata")
      recruits <- cbind(x, y, logheight, grid)
      # use overlap function to filter out recruits which fail because they land on spots already occupied by trees/seedlings
      if (nrow(recruits) > 0) {
        recruits.success <- rm.overlap(recruits=recruits, 
                                       trees=rbind(ppoT, sceT, ppiT, aosT), 
                                       seedlings=rbind(ppoS, sceS, ppiS),
                                       grid.neighbours=grid.neighbours)
        # if there are still successfully recruiting seedlings, rbind them to existing seedlings
        if(nrow(recruits.success) > 0) {
          ppiS <- rbind(ppiS, recruits.success)
          recs.ppiS <- c(recs.ppiS, nrow(recruits.success))
        }
      }
      if(nrow(recruits)==0 | nrow(recruits.success)==0) recs.ppiS <- c(recs.ppiS, 0)
    }
  }
  
  # remove seedlings which recruit out of landscape boundary (temporary solution--need to torus-ize landscape eventually)
  if(sum(is.na(ppoS$grid)) > 0)  ppoS <- ppoS[-which(is.na(ppoS$grid)),]
  if(sum(is.na(sceS$grid)) > 0)  sceS <- sceS[-which(is.na(sceS$grid)),]
  if(sum(is.na(ppiS$grid)) > 0)  ppiS <- ppiS[-which(is.na(ppiS$grid)),]
  
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
  # PPI ("Pometia.pinnata")
  n1.ppiT <- nrow(ppiT)
  n.ppiT <- c(n.ppiT, n1.ppiT) 
  n1.ppiS <- nrow(ppiS)
  n.ppiS <- c(n.ppiS, n1.ppiS) 
  # AOS ("All.other.spp")
  n1.aosT <- nrow(aosT)
  n.aosT <- c(n.aosT, n1.aosT) 
  

  # calculate mean sizes of adults and saplings and add to string  
  # PPO ("Prunus.polystachya")
  z.ppoT <- c(z.ppoT, mean(ppoT$logdbh))
  h.ppoS <- c(h.ppoS, mean(ppoS$logheight))
  # SCE ("Strombosia.ceylanica")
  z.sceT <- c(z.sceT, mean(sceT$logdbh))
  h.sceS <- c(h.sceS, mean(sceS$logheight))
  # PPI ("Pometia.pinnata")
  z.ppiT <- c(z.ppiT, mean(ppiT$logdbh))
  h.ppiS <- c(h.ppiS, mean(ppiS$logheight))
  # AOS ("All.other.spp")
  z.aosT <- c(z.aosT, mean(aosT$logdbh))
  
  # record total basal areas of adults in swamp versus non-swamp areas
  ba.ppoT[which(ba.ppoT$time==time), 2:3] <- tapply(pi*(exp(ppoT[,"logdbh"])/2)^2, extract(nssf.m, ppoT[,1:2]), sum)
  ba.sceT[which(ba.sceT$time==time), 2:3] <- tapply(pi*(exp(sceT[,"logdbh"])/2)^2, extract(nssf.m, sceT[,1:2]), sum)
  ba.ppiT[which(ba.ppiT$time==time), 2:3] <- tapply(pi*(exp(ppiT[,"logdbh"])/2)^2, extract(nssf.m, ppiT[,1:2]), sum)
  ba.aosT[which(ba.aosT$time==time), 2:3] <- tapply(pi*(exp(aosT[,"logdbh"])/2)^2, extract(nssf.m, aosT[,1:2]), sum)
  # to calculate SSI need to know area of swamp versus non-swamp across landscape at that time point
  ba.ppoT[which(ba.ppoT$time==time), "landscape.swamp.prop"] <- sum(values(nssf.m), na.rm = T) / length(na.omit(values(nssf.m)))
  ba.sceT[which(ba.sceT$time==time), "landscape.swamp.prop"] <- sum(values(nssf.m), na.rm = T) / length(na.omit(values(nssf.m)))
  ba.ppiT[which(ba.ppiT$time==time), "landscape.swamp.prop"] <- sum(values(nssf.m), na.rm = T) / length(na.omit(values(nssf.m)))
  ba.aosT[which(ba.aosT$time==time), "landscape.swamp.prop"] <- sum(values(nssf.m), na.rm = T) / length(na.omit(values(nssf.m)))
  
  # print message to check status during runs
  message(paste(n1.ppoT + n1.ppoS + n1.sceT + n1.sceS + n1.ppiT + n1.ppiS, "stems at year", time, "  "),
          appendLF = TRUE)
    
  # move to next time point
  time <- time + 1
  # update progress bar
  setTxtProgressBar(pb, time)
}

close(pb)

# output to be saved
out3.usual <- 
  list(nssf.m = nssf.m, 
       ppoT.init = ppoT.init, sceT.init = sceT.init, ppiT.init = ppiT.init, 
       ppoS.init = ppoS.init, sceS.init = sceS.init, ppiS.init = ppiS.init,
       ppoT = ppoT, sceT = sceT, ppiT = ppiT,
       ppoS = ppoS, sceS = sceS, ppiS = ppiS,
       n.ppoT = n.ppoT, n.sceT = n.sceT, n.ppiT = n.ppiT,
       n.ppoS = n.ppoS, n.sceS = n.sceS, n.ppiS = n.ppiS,
       z.ppoT = z.ppoT, z.sceT = z.sceT, z.ppiT = z.ppiT,
       h.ppoS = h.ppoS, h.sceS = h.sceS, h.ppiS = h.ppiS,
       ba.ppoT = ba.ppoT, ba.sceT = ba.sceT, ba.ppiT = ba.ppiT,
       deaths.ppoT = deaths.ppoT, deaths.sceT = deaths.sceT, deaths.ppiT = deaths.ppiT, deaths.aosT = deaths.aosT,
       deaths.ppoS = deaths.ppoS, deaths.sceS = deaths.sceS, deaths.ppiS = deaths.ppiS, 
       recs.ppoS = recs.ppoS, recs.sceS = recs.sceS, recs.ppiS = recs.ppiS, recs.aosT = recs.aosT,
       aosT = aosT, n.aosT = n.aosT, z.aosT = z.aosT, ba.aosT = ba.aosT)
saveRDS(out3.usual, file = "out/sim_out3_usual.rds")


