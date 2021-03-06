library(raster)
library(sp)
library(poweRlaw)
library(doParallel)


# Settings for parallerisation
registerDoParallel(cores = 24)


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
source("code/setup_map.R")


#######################
# COMMENCE SIMULATION #
#######################

time <- 1
maxTime <- 100

ppoT.init <- ppoT
ppoS.init <- ppoS

sceT.init <- sceT
sceS.init <- sceS

n.ppoT <- nrow(ppoT)
n.ppoS <- nrow(ppoS)

n.sceT <- nrow(sceT)
n.sceS <- nrow(sceS)

z.ppoT <- mean(ppoT$logdbh)
h.ppoS <- mean(ppoS$logheight)

z.sceT <- mean(sceT$logdbh)
h.sceS <- mean(sceS$logheight)

ba.ppoT <- sum(pi * (exp(ppoT$logdbh) / 2)^2)
ba.sceT <- sum(pi * (exp(sceT$logdbh) / 2)^2)

# progress bar
pb <- txtProgressBar(min = 1, max = maxTime, style = 3)
while (time < maxTime) {
  
  # Choose landscape based on scenario and time (beyond 22 years, just use the last time point (=2042))
  nssf.m <- nssf.extreme[[ifelse(time>22, 22, time)]]
  
  # Recruitment: initiate fruiting, but don't add new recruits to ppoS yet (let them grow/die first)
  # PPO ("Prunus.polystachya")
  fruiting.index.ppo <- which(p_bz(nssf.m, ppoT, "Prunus.polystachya")==1)
  prod.vec.ppo <- b_z(nssf.m, ppoT[fruiting.index.ppo,], "Prunus.polystachya")
  parent.loc.ppo <- ppoT[fruiting.index.ppo, 1:2]
  # SCE ("Strombosia.ceylanica")
  fruiting.index.sce <- which(p_bz(nssf.m, sceT, "Strombosia.ceylanica")==1)
  prod.vec.sce <- b_z(nssf.m, sceT[fruiting.index.sce,], "Strombosia.ceylanica")
  parent.loc.sce <- sceT[fruiting.index.sce, 1:2]
  
  # Extract competition measures
  # PPO ("Prunus.polystachya")
  inter.on.PPO <- inter.calc(rbind(sceT, aosT), rbind(sceS, aosS), ppoT, ppoS, 
                             sp="Prunus.polystachya", grid.neighbours = grid.neighbours)
  intra.on.PPO <- intra.calc(ppoT, ppoS, sp="Prunus.polystachya", grid.neighbours = grid.neighbours)
  CAE.on.ppoS <- CAE.calc(ppoT, ppoS, grid.neighbours = grid.neighbours)
  HAE.on.ppoS <- HAE.calc(rbind(sceT, aosT), ppoS, grid.neighbours = grid.neighbours)
  intra.dist.on.PPO <- intra.dist.calc(ppoT, grid.neighbours = grid.neighbours)
  inter.dist.on.PPO <- inter.dist.calc(rbind(sceT, aosT), ppoT, grid.neighbours = grid.neighbours)
  # SCE ("Strombosia.ceylanica")
  inter.on.SCE <- inter.calc(rbind(ppoT, aosT), rbind(ppoS, aosS), sceT, sceS, 
                             sp="Strombosia.ceylanica", grid.neighbours = grid.neighbours)
  intra.on.SCE <- intra.calc(sceT, sceS, sp="Strombosia.ceylanica", grid.neighbours = grid.neighbours)
  CAE.on.sceS <- CAE.calc(sceT, sceS, grid.neighbours = grid.neighbours)
  HAE.on.sceS <- HAE.calc(rbind(ppoT, aosT), sceS, grid.neighbours = grid.neighbours)
  intra.dist.on.SCE <- intra.dist.calc(sceT, grid.neighbours = grid.neighbours)
  inter.dist.on.SCE <- inter.dist.calc(rbind(ppoT, aosT), sceT, grid.neighbours = grid.neighbours)
  
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
  rm(inter.on.PPO, intra.on.PPO, CAE.on.ppoS, HAE.on.ppoS, intra.dist.on.PPO, inter.dist.on.PPO,
     inter.on.SCE, intra.on.SCE, CAE.on.sceS, HAE.on.sceS, intra.dist.on.SCE, inter.dist.on.SCE)
  
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
  
  # This is the current number of seedlings:
  # PPO ("Prunus.polystachya")
  n.old.ppoS <- nrow(ppoS)
  # SCE ("Strombosia.ceylanica")
  n.old.sceS <- nrow(sceS)
  
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
                                       trees=rbind(ppoT, sceT, aosT), 
                                       seedlings=rbind(ppoS, sceS, aosS),
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
                                       trees=rbind(ppoT, sceT, aosT), 
                                       seedlings=rbind(ppoS, sceS, aosS),
                                       grid.neighbours=grid.neighbours)
        if(nrow(recruits.success) > 0) sceS <- rbind(sceS, recruits.success)
      }
    }
  }

  # Kill off 50% of all recruits that are located < 10cm from each other
  # PPO ("Prunus.polystachya")
  ppoS <- kill.rec(ppoS, n.old.ppoS)

  # SCE ("Strombosia.ceylanica")
  sceS <- kill.rec(sceS, n.old.sceS)

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

  # calculate total basal areas of adults 
  # not doing saplings for now, because they only have height (can convert later?)
  ba.ppoT <- c(ba.ppoT, sum(pi * (exp(ppoT$logdbh) / 2)^2))
  ba.sceT <- c(ba.sceT, sum(pi * (exp(sceT$logdbh) / 2)^2))

  # print message to check status during runs
  message(paste(n1.ppoT + n1.ppoS + n1.sceT + n1.sceS, "stems at year", time, "  "),
          appendLF = TRUE)
    
  # move to next time point
  time <- time + 1
  # update progress bar
  setTxtProgressBar(pb, time)
}

close(pb)

# output to be saved
out <- 
  list(nssf.m = nssf.m, 
       ppoT.init = ppoT.init, sceT.init = sceT.init,
       ppoS.init = ppoS.init, sceS.init = sceS.init,
       ppoT = ppoT, sceT = sceT,
       ppoS = ppoS, sceS = sceS,
       n.ppoT = n.ppoT, n.sceT = n.sceT,
       n.ppoS = n.ppoS, n.sceS = n.sceS,
       z.ppoT = z.ppoT, z.sceT = z.sceT,
       h.ppoS = h.ppoS, h.sceS = h.sceS,
       ba.ppoT = ba.ppoT, ba.sceT = ba.sceT)
saveRDS(out, file = "out/sim_out.rds")
#saveRDS(out, file = "D:\\National University of Singapore\\Chong Kwek Yan - CRSF\\Data\\IBM\\out\\Extreme3 scenario 50 years.rds")

