library(raster)
library(sp)
library(poweRlaw)
library(doParallel)
library(rgdal)

# Settings for parallerisation
registerDoParallel(cores = 64)


#################
# GROWTH PARAMS #
#################

grow.T.parm <- t(read.csv("data/tree growth parameters Apr21.csv", header=T, row.names=1))
grow.T.parm.unscale <- read.csv("data/tree growth model unscale params.csv", header=T, row.names=1)
grow.S.parm <- read.csv("data/seedling lmer growth params CAA Apr21.csv", header=T, row.names=1)
surv.parm <- t(read.csv("data/surv params Jun21.csv", header=T, row.names=1))
surv.parm.unscale <- read.csv("data/survival model unscale params.csv", header=T, row.names=1)
fruit.parm <- read.csv("data/fruiting parameters Apr21.csv", header=T, row.names=1)
tran.parm <- read.csv("data/transition params CAA Jun21.csv", header=T, row.names=1)
rec.parm <- t(read.csv("data/dispersal kernel parameters Apr21.csv", header=T, row.names=1))
HM.parm <- t(read.csv("data/HM values.csv", header=T, row.names=1))

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
source("code/setup_map_four spp.R")


#######################
# COMMENCE SIMULATION #
#######################

time <- 1
maxTime <- 23

# take a snapshot of initial conditions
ppoT.init <- ppoT
ppoS.init <- ppoS
sceT.init <- sceT
sceS.init <- sceS
ppiT.init <- ppiT
ppiS.init <- ppiS
gneT.init <- gneT
gneS.init <- gneS

# record population sizes
n.ppoT <- nrow(ppoT)
n.ppoS <- nrow(ppoS)
n.sceT <- nrow(sceT)
n.sceS <- nrow(sceS)
n.ppiT <- nrow(ppiT)
n.ppiS <- nrow(ppiS)
n.gneT <- nrow(gneT)
n.gneS <- nrow(gneS)

# record mean sizes
z.ppoT <- mean(ppoT$logdbh)
h.ppoS <- mean(ppoS$logheight)
z.sceT <- mean(sceT$logdbh)
h.sceS <- mean(sceS$logheight)
z.ppiT <- mean(ppiT$logdbh)
h.ppiS <- mean(ppiS$logheight)
z.gneT <- mean(gneT$logdbh)
h.gneS <- mean(gneS$logheight)

# create table for total BA in swamp vs non-swamp
ba.ppoT <- data.frame(time=seq(time, maxTime, 1), nonswamp=NA, swamp=NA, landscape.swamp.prop=NA)
ba.sceT <- data.frame(time=seq(time, maxTime, 1), nonswamp=NA, swamp=NA, landscape.swamp.prop=NA)
ba.ppiT <- data.frame(time=seq(time, maxTime, 1), nonswamp=NA, swamp=NA, landscape.swamp.prop=NA)
ba.gneT <- data.frame(time=seq(time, maxTime, 1), nonswamp=NA, swamp=NA, landscape.swamp.prop=NA)

# progress bar
pb <- txtProgressBar(min = 1, max = maxTime, style = 3)
while (time <= maxTime) {
  
  # Choose landscape based on scenario and time (beyond 22 years, just use the last time point (=2042))
  nssf.m <- nssf.extreme[[ifelse(time>22, 22, time)]]
  
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
  # GNE ("Gironniera nervosa")
  fruiting.index.gne <- which(p_bz(nssf.m, gneT, "Gironniera.nervosa")==1)
  prod.vec.gne <- b_z(nssf.m, gneT[fruiting.index.gne,], "Gironniera.nervosa")
  parent.loc.gne <- gneT[fruiting.index.gne, 1:2]
  
  # Extract competition measures
  # PPO ("Prunus.polystachya")
  inter.on.PPO <- inter.calc(rbind(sceT, ppiT, gneT, aosT), rbind(sceS, ppiS, gneS), ppoT, ppoS, 
                             sp="Prunus.polystachya", grid.neighbours = grid.neighbours)
  intra.on.PPO <- intra.calc(ppoT, ppoS, sp="Prunus.polystachya", grid.neighbours = grid.neighbours)
  CAE.on.ppoS <- CAE.calc(ppoT, ppoS, grid.neighbours = grid.neighbours)
  HAE.on.ppoS <- HAE.calc(rbind(sceT, ppiT, gneT, aosT), ppoS, grid.neighbours = grid.neighbours)
  intra.dist.on.PPO <- intra.dist.calc(ppoT, grid.neighbours = grid.neighbours)
  inter.dist.on.PPO <- inter.dist.calc(rbind(sceT, ppiT, gneT, aosT), ppoT, grid.neighbours = grid.neighbours)
  # SCE ("Strombosia.ceylanica")
  inter.on.SCE <- inter.calc(rbind(ppoT, ppiT, gneT, aosT), rbind(ppoS, ppiS, gneS), sceT, sceS, 
                             sp="Strombosia.ceylanica", grid.neighbours = grid.neighbours)
  intra.on.SCE <- intra.calc(sceT, sceS, sp="Strombosia.ceylanica", grid.neighbours = grid.neighbours)
  CAE.on.sceS <- CAE.calc(sceT, sceS, grid.neighbours = grid.neighbours)
  HAE.on.sceS <- HAE.calc(rbind(ppoT, ppiT, gneT, aosT), sceS, grid.neighbours = grid.neighbours)
  intra.dist.on.SCE <- intra.dist.calc(sceT, grid.neighbours = grid.neighbours)
  inter.dist.on.SCE <- inter.dist.calc(rbind(ppoT, ppiT, gneT, aosT), sceT, grid.neighbours = grid.neighbours)
  # PPI ("Pometia.pinnata")
  inter.on.PPI <- inter.calc(rbind(sceT, ppoT, gneT, aosT), rbind(sceS, ppoS, gneS), ppiT, ppiS, 
                             sp="Pometia.pinnata", grid.neighbours = grid.neighbours)
  intra.on.PPI <- intra.calc(ppiT, ppiS, sp="Pometia.pinnata", grid.neighbours = grid.neighbours)
  CAE.on.ppiS <- CAE.calc(ppiT, ppiS, grid.neighbours = grid.neighbours)
  HAE.on.ppiS <- HAE.calc(rbind(sceT, ppoT, gneT, aosT), ppiS, grid.neighbours = grid.neighbours)
  intra.dist.on.PPI <- intra.dist.calc(ppiT, grid.neighbours = grid.neighbours)
  inter.dist.on.PPI <- inter.dist.calc(rbind(sceT, ppoT, gneT, aosT), ppiT, grid.neighbours = grid.neighbours)
  # GNE ("Gironniera.nervosa")
  inter.on.GNE <- inter.calc(rbind(ppoT, ppiT, sceT, aosT), rbind(ppoS, ppiS, sceS), gneT, gneS, 
                             sp="Gironniera.nervosa", grid.neighbours = grid.neighbours)
  intra.on.GNE <- intra.calc(gneT, gneS, sp="Gironniera.nervosa", grid.neighbours = grid.neighbours)
  CAE.on.gneS <- CAE.calc(gneT, gneS, grid.neighbours = grid.neighbours)
  HAE.on.gneS <- HAE.calc(rbind(ppoT, ppiT, sceT, aosT), gneS, grid.neighbours = grid.neighbours)
  intra.dist.on.GNE <- intra.dist.calc(gneT, grid.neighbours = grid.neighbours)
  inter.dist.on.GNE <- inter.dist.calc(rbind(ppoT, ppiT, sceT, aosT), gneT, grid.neighbours = grid.neighbours)
  
  # Seedling growth and survival
  # PPO ("Prunus.polystachya")
  ppoS$logheight <- sS_h(nssf.m, ppoS, "Prunus.polystachya", intra.on.PPO[[4]], inter.on.PPO[[4]]) * 
    GS_h1h(nssf.m, ppoS, "Prunus.polystachya", CAE.on.ppoS, HAE.on.ppoS)
  if(sum(ppoS$logheight==0, na.rm = TRUE) > 0)  ppoS <- ppoS[-which(ppoS$logheight==0),]
  # SCE ("Strombosia.ceylanica")
  sceS$logheight <- sS_h(nssf.m, sceS, "Strombosia.ceylanica", intra.on.SCE[[4]], inter.on.SCE[[4]]) * 
    GS_h1h(nssf.m, sceS, "Strombosia.ceylanica", CAE.on.sceS, HAE.on.sceS)
  if(sum(sceS$logheight==0, na.rm = TRUE) > 0)  sceS <- sceS[-which(sceS$logheight==0),]
  # PPI ("Pometia.pinnata")
  ppiS$logheight <- sS_h(nssf.m, ppiS, "Pometia.pinnata", intra.on.PPI[[4]], inter.on.PPI[[4]]) * 
    GS_h1h(nssf.m, ppiS, "Pometia.pinnata", CAE.on.ppiS, HAE.on.ppiS)
  if(sum(ppiS$logheight==0, na.rm = TRUE) > 0)  ppiS <- ppiS[-which(ppiS$logheight==0),]
  # GNE ("Gironniera.nervosa")
  gneS$logheight <- sS_h(nssf.m, gneS, "Gironniera.nervosa", intra.on.GNE[[4]], inter.on.GNE[[4]]) * 
    GS_h1h(nssf.m, gneS, "Gironniera.nervosa", CAE.on.gneS, HAE.on.gneS)
  if(sum(gneS$logheight==0, na.rm = TRUE) > 0)  gneS <- gneS[-which(gneS$logheight==0),]
  
  # Tree growth and survival
  # PPO ("Prunus.polystachya")
  ppoT$logdbh <- sT_z(nssf.m, ppoT, "Prunus.polystachya", intra.on.PPO[[1]], intra.on.PPO[[3]], inter.on.PPO[[1]], inter.on.PPO[[3]]) * 
    GT_z1z(nssf.m, ppoT, "Prunus.polystachya", intra.dist.on.PPO, inter.dist.on.PPO)
  if(sum(ppoT$logdbh==0, na.rm = TRUE) > 0)  ppoT <- ppoT[-which(ppoT$logdbh==0),]
  # SCE ("Strombosia.ceylanica")
  sceT$logdbh <- sT_z(nssf.m, sceT, "Strombosia.ceylanica", intra.on.SCE[[1]], intra.on.SCE[[3]], inter.on.SCE[[1]], inter.on.SCE[[3]]) * 
    GT_z1z(nssf.m, sceT, "Strombosia.ceylanica", intra.dist.on.SCE, inter.dist.on.SCE)
  if(sum(sceT$logdbh==0, na.rm = TRUE) > 0)  sceT <- sceT[-which(sceT$logdbh==0),]
  # PPI ("Pometia.pinnata")
  ppiT$logdbh <- sT_z(nssf.m, ppiT, "Pometia.pinnata", intra.on.PPI[[1]], intra.on.PPI[[3]], inter.on.PPI[[1]], inter.on.PPI[[3]]) * 
    GT_z1z(nssf.m, ppiT, "Pometia.pinnata", intra.dist.on.PPI, inter.dist.on.PPI)
  if(sum(ppiT$logdbh==0, na.rm = TRUE) > 0)  ppiT <- ppiT[-which(ppiT$logdbh==0),]
  # GNE ("Gironniera.nervosa")
  gneT$logdbh <- sT_z(nssf.m, gneT, "Gironniera.nervosa", intra.on.GNE[[1]], intra.on.GNE[[3]], inter.on.GNE[[1]], inter.on.GNE[[3]]) * 
    GT_z1z(nssf.m, gneT, "Gironniera.nervosa", intra.dist.on.GNE, inter.dist.on.GNE)
  if(sum(gneT$logdbh==0, na.rm = TRUE) > 0)  gneT <- gneT[-which(gneT$logdbh==0),]
  
  # Conserve memory by deleting competition indices (don't need anymore)
  rm(inter.on.PPO, intra.on.PPO, CAE.on.ppoS, HAE.on.ppoS, intra.dist.on.PPO, inter.dist.on.PPO,
     inter.on.SCE, intra.on.SCE, CAE.on.sceS, HAE.on.sceS, intra.dist.on.SCE, inter.dist.on.SCE,
     inter.on.PPI, intra.on.PPI, CAE.on.ppiS, HAE.on.ppiS, intra.dist.on.PPI, inter.dist.on.PPI,
     inter.on.GNE, intra.on.GNE, CAE.on.gneS, HAE.on.gneS, intra.dist.on.GNE, inter.dist.on.GNE)
  
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
  # GNE ("Gironniera.nervosa")
  gneTS <- T_z1h1(gneT, gneS, "Gironniera.nervosa")
  gneT <- gneTS[[1]]
  gneS <- gneTS[[2]]
  rm(gneTS)
  
  # This is the current number of seedlings:
  n.old.ppoS <- nrow(ppoS)
  n.old.sceS <- nrow(sceS)
  n.old.ppiS <- nrow(ppiS)
  n.old.gneS <- nrow(gneS)
  
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
                                       trees=rbind(ppoT, sceT, ppiT, gneT, aosT), 
                                       seedlings=rbind(ppoS, sceS, ppiS, gneS),
                                       grid.neighbours=grid.neighbours)
        # if there are still successfully recruiting seedlings, rbind them to existing seedlings
        if(nrow(recruits.success) > 0) ppoS <- rbind(ppoS, recruits.success)
      }
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
                                       trees=rbind(ppoT, sceT, ppiT, gneT, aosT), 
                                       seedlings=rbind(ppoS, sceS, ppiS, gneS),
                                       grid.neighbours=grid.neighbours)
        if(nrow(recruits.success) > 0) sceS <- rbind(sceS, recruits.success)
      }
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
                                       trees=rbind(ppoT, sceT, ppiT, gneT, aosT), 
                                       seedlings=rbind(ppoS, sceS, ppiS, gneS),
                                       grid.neighbours=grid.neighbours)
        # if there are still successfully recruiting seedlings, rbind them to existing seedlings
        if(nrow(recruits.success) > 0) ppiS <- rbind(ppiS, recruits.success)
      }
    }
  }
  # GNE ("Gironniera.nervosa")
  if(length(prod.vec.gne) > 0){
    for(i in 1:length(prod.vec.gne)){
      distances <- twoDT.sample(prod.vec.gne[i], "Gironniera.nervosa", max_dist = 410)
      angles <- runif(length(distances), min=0, max=360)
      x <- parent.loc.gne[i,1] + (distances * cos(angles))
      y <- parent.loc.gne[i,2] + (distances * sin(angles))
      grid <- extract(nssf.100, cbind(x,y))
      logheight <- c_0h1(length(distances), "Gironniera.nervosa")
      recruits <- cbind(x, y, logheight, grid)
      if (nrow(recruits) > 0) {
        recruits.success <- rm.overlap(recruits, 
                                       trees=rbind(ppoT, sceT, ppiT, gneT, aosT), 
                                       seedlings=rbind(ppoS, sceS, ppiS, gneS),
                                       grid.neighbours=grid.neighbours)
        if(nrow(recruits.success) > 0) gneS <- rbind(gneS, recruits.success)
      }
    }
  }
  # Kill off 50% of all recruits that are located < 10cm from each other
  # temporarily omitted
  #ppoS <- kill.rec(ppoS, n.old.ppoS)
  #sceS <- kill.rec(sceS, n.old.sceS)
  #ppiS <- kill.rec(ppiS, n.old.ppiS)
  #gneS <- kill.rec(gneS, n.old.gneS)

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
  # GNE ("Gironniera.nervosa")
  n1.gneT <- nrow(gneT)
  n.gneT <- c(n.gneT, n1.gneT) 
  n1.gneS <- nrow(gneS)
  n.gneS <- c(n.gneS, n1.gneS) 
  
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
  # GNE ("Gironniera.nervosa")
  z.gneT <- c(z.gneT, mean(gneT$logdbh))
  h.gneS <- c(h.gneS, mean(gneS$logheight))
  
  # record total basal areas of adults in swamp versus non-swamp areas
  ba.ppoT[which(ba.ppoT$time==time), 2:3] <- tapply(pi*(exp(ppoT[,"logdbh"])/2)^2, extract(nssf.m, ppoT[,1:2]), sum)
  ba.sceT[which(ba.sceT$time==time), 2:3] <- tapply(pi*(exp(sceT[,"logdbh"])/2)^2, extract(nssf.m, sceT[,1:2]), sum)
  ba.ppiT[which(ba.ppiT$time==time), 2:3] <- tapply(pi*(exp(ppiT[,"logdbh"])/2)^2, extract(nssf.m, ppiT[,1:2]), sum)
  ba.gneT[which(ba.gneT$time==time), 2:3] <- tapply(pi*(exp(gneT[,"logdbh"])/2)^2, extract(nssf.m, gneT[,1:2]), sum)
  # to calculate SSI need to know area of swamp versus non-swamp across landscape at that time point
  ba.ppoT[which(ba.ppoT$time==time), "landscape.swamp.prop"] <- sum(values(nssf.m), na.rm = T) / length(na.omit(values(nssf.m)))
  ba.sceT[which(ba.sceT$time==time), "landscape.swamp.prop"] <- sum(values(nssf.m), na.rm = T) / length(na.omit(values(nssf.m)))
  ba.ppiT[which(ba.ppiT$time==time), "landscape.swamp.prop"] <- sum(values(nssf.m), na.rm = T) / length(na.omit(values(nssf.m)))
  ba.gneT[which(ba.gneT$time==time), "landscape.swamp.prop"] <- sum(values(nssf.m), na.rm = T) / length(na.omit(values(nssf.m)))
          
          
  # print message to check status during runs
  message(paste(n1.ppoT + n1.ppoS + n1.sceT + n1.sceS + n1.ppiT + n1.ppiS + n1.gneT + n1.gneS , "stems at year", time, "  "),
          appendLF = TRUE)
    
  # move to next time point
  time <- time + 1
  # update progress bar
  setTxtProgressBar(pb, time)
}

close(pb)

# output to be saved
out.extreme <- 
  list(nssf.m = nssf.m, 
       ppoT.init = ppoT.init, sceT.init = sceT.init, ppiT.init = ppiT.init, gneT.init = gneT.init,
       ppoS.init = ppoS.init, sceS.init = sceS.init, ppiS.init = ppiS.init, gneS.init = gneS.init,
       ppoT = ppoT, sceT = sceT, ppiT = ppiT, gneT = gneT,
       ppoS = ppoS, sceS = sceS, ppiS = ppiS, gneS = gneS,
       n.ppoT = n.ppoT, n.sceT = n.sceT, n.ppiT = n.ppiT, n.gneT = n.gneT,
       n.ppoS = n.ppoS, n.sceS = n.sceS, n.ppiS = n.ppiS, n.gneS = n.gneS,
       z.ppoT = z.ppoT, z.sceT = z.sceT, z.ppiT = z.ppiT, z.gneT = z.gneT,
       h.ppoS = h.ppoS, h.sceS = h.sceS, h.ppiS = h.ppiS, h.gneS = h.gneS,
       ba.ppoT = ba.ppoT, ba.sceT = ba.sceT, ba.ppiT = ba.ppiT, ba.gneT = ba.gneT)
saveRDS(out.extreme, file = "out/sim_out_extreme.rds")
#saveRDS(out, file = "D:\\National University of Singapore\\Chong Kwek Yan - CRSF\\Data\\IBM\\out\\Extreme3 scenario 50 years.rds")

