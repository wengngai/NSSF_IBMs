#################
# PLOT OUTPUTS  #
#################

usual <- readRDS(file = "out/sim_out3_usual_Aug21_11years.rds")
summary(usual)

# Define a colour palette
col.pal <- c("#E3C16F", "#946846", "#FAFF70", "#6D213C", "#65DEF1", "#F3F5F6", "#BAAB68", "#CBC5EA")
col.t <- adjustcolor(col.pal, alpha.f=0.6)

# read in hydrological model predictions (each with 22 years: 2021-2042)
# Scenarios
nssf.usual <- list()
nssf.extreme <- list()
nssf.low <- list()
nssf.high <- list()

# read in and reset values to 0/1 (resampling issues)
for(i in 1:22){
        nssf.usual[[i]] <- raster("data/As usual scenario.tif", band=i)
        values(nssf.usual[[i]]) <- ifelse(values(nssf.usual[[i]]) < 0.5, 0, 1)
        nssf.extreme[[i]] <- raster("data/Extreme scenario.tif", band=i)
        values(nssf.extreme[[i]]) <- ifelse(values(nssf.extreme[[i]]) < 0.5, 0, 1)
        nssf.low[[i]] <- raster("data/Low RF scenario.tif", band=i)
        values(nssf.low[[i]]) <- ifelse(values(nssf.low[[i]]) < 0.5, 0, 1)
        nssf.high[[i]] <- raster("data/High RF scenario.tif", band=i)
        values(nssf.high[[i]]) <- ifelse(values(nssf.high[[i]]) < 0.5, 0, 1)
}

# need to crop landscape to a central 500*500-m simulation range
crop.extent <- extent(367000, 367500, 153000, 153500)
crop.area <- 500*500
#plot(nssf.usual[[1]])
#rect(367000, 153000, 367500, 153500)

for(i in 1:22){
        nssf.usual[[i]] <- crop(nssf.usual[[i]], crop.extent)
        nssf.extreme[[i]] <- crop(nssf.extreme[[i]], crop.extent)
        nssf.low[[i]] <- crop(nssf.low[[i]], crop.extent)
        nssf.high[[i]] <- crop(nssf.high[[i]], crop.extent)
}

#############
# LOCATIONS #
#############

#pdf(file = "D:\\Dropbox\\twn idiwn\\Post doc\\IBM temp\\distributions 3spp Aug21.pdf", width=16, height=9)
prop = 0.25
# Initial condition
par(mfrow=c(1,2), mar=c(5,5,3,2))
plot(nssf.usual[[1]], legend=F, col=col.pal[c(6,5)], main="2021")
# PPO
i.ppoT <- sample(1:nrow(usual$ppoT), prop*nrow(usual$ppoT))
points(usual$ppoT.init[i.ppoT,], cex=usual$ppoT.init$logdbh[i.ppoT]/2, col=col.t[1], pch=16)
# SCE
i.sceT <- sample(1:nrow(usual$sceT), prop*nrow(usual$sceT))
points(usual$sceT.init[i.sceT,], cex=usual$sceT.init$logdbh[i.sceT]/2, col=col.t[2], pch=16)
# PPI
i.ppiT <- sample(1:nrow(usual$ppiT), prop*nrow(usual$ppiT))
points(usual$ppiT.init[i.ppiT,], cex=usual$ppiT.init$logdbh[i.ppiT]/2, col=col.t[3], pch=16)
# All other species (static--just for a background iproperspecific competition pressure)
i.aosT <- sample(1:nrow(usual$aosT), prop*nrow(usual$aosT))
#points(aosT[i.aosT,], cex=aosT$logdbh[i.aosT]/2,  pch=1, col="grey70")
legend('topleft', col=col.t[c(1,2,3)], pch=16, bty="n", pt.cex=2,
       legend=c(
               "P. polystachya (non-swamp generalist)",
               "S. ceylanica (swamp generalist)",
               "P. pinnata (swamp specialist)"
               ))

# Final condition
plot(nssf.usual[[22]], legend=F, col=col.pal[c(6,5)], main="2042")
# PPO
points(usual$ppoT[i.ppoT,], cex=usual$ppoT$logdbh[i.ppoT]/2, col=col.t[1], pch=16)
# SCE
points(usual$sceT[i.sceT,], cex=usual$sceT$logdbh[i.sceT]/2, col=col.t[2], pch=16)
# PPI
points(usual$ppiT[i.ppiT,], cex=usual$ppiT$logdbh[i.ppiT]/2, col=col.t[3], pch=16)
# All other species (static--just for a background interspecific competition pressure)
#points(aosT[i.aosT,], cex=aosT$logdbh[i.aosT]/2,  pch=1, col="grey70")
scalebar(400, xy=c(367800, 151550), type="bar", lonlat=F, below="metres", divs=2)
dev.off()

####################
# POPULATION SIZES #
####################

#pdf(file = "D:\\Dropbox\\twn idiwn\\Post doc\\IBM temp\\populations 3spp Aug21.pdf", width=16, height=9)
par(mfrow=c(1,2), mar=c(5.5,5.5,2,2))
plot(log(usual$n.ppoT) ~ c(1:length(usual$n.ppoT)), lwd=5, col=col.pal[2], type="l",
     ylab="Adult population size (log-transformed)", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=log(c(min(c(usual$n.ppoT, usual$n.sceT, usual$n.ppiT, usual$n.aosT))+1, 
            max(c(usual$n.ppoT, usual$n.sceT, usual$n.ppiT, usual$n.aosT))))
     )
lines(log(usual$n.sceT) ~ c(1:length(usual$n.sceT)), lwd=5, col=col.pal[1])
lines(log(usual$n.ppiT) ~ c(1:length(usual$n.ppiT)), lwd=5, col=col.pal[3])
lines(log(usual$n.aosT) ~ c(1:length(usual$n.aosT)), lwd=5, col="grey")
legend('topleft', col=c(col.pal[c(2,1,3)], "grey"), bty="n", lwd=3, lty=1, title="Species",
       legend=c(
               "P. polystachya (non-swamp generalist)",
               "S. ceylanica (swamp generalist)",
               "P. pinnata (swamp specialist)",
               "All other species"
       ))
plot(log(usual$n.ppoS+1) ~ c(1:length(usual$n.ppoS)), lwd=5, col=col.pal[2], type="l",
     ylab="Seedling population size (log-transformed)", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=log(c(min(c(usual$n.ppoS, usual$n.sceS, usual$n.ppiS))+1, 
            max(c(usual$n.ppoS, usual$n.sceS, usual$n.ppiS))))
     )
lines(log(usual$n.sceS+1) ~ c(1:length(usual$n.sceS)), lwd=5, col=col.pal[1])
lines(log(usual$n.ppiS+1) ~ c(1:length(usual$n.ppiS)), lwd=5, col=col.pal[3])
dev.off()

#####################
# DEMOGRAPHIC RATES #
#####################

ppoT.mort <- usual$deaths.ppoT / usual$n.ppoT[-length(usual$n.ppoT)]
ppoS.mort <- usual$deaths.ppoS / usual$n.ppoS[-length(usual$n.ppoS)]
sceT.mort <- usual$deaths.sceT / usual$n.sceT[-length(usual$n.sceT)]
sceS.mort <- usual$deaths.sceS / usual$n.sceS[-length(usual$n.sceS)]
ppiT.mort <- usual$deaths.ppiT / usual$n.ppiT[-length(usual$n.ppiT)]
ppiS.mort <- usual$deaths.ppiS / usual$n.ppiS[-length(usual$n.ppiS)]
aosT.mort <- usual$deaths.aosT / usual$n.aosT[-length(usual$n.aosT)]

par(mfrow=c(1,2), mar=c(5.5,5.5,2,2))
plot(ppoT.mort ~ c(1:length(ppoT.mort)), lwd=5, col=col.pal[2], type="l",
     ylab="Annual mortality rate", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=c(0, max(c(ppoT.mort, ppoS.mort, sceT.mort, sceS.mort, ppiT.mort, ppiS.mort), na.rm=T)))
lines(sceT.mort ~ c(1:length(sceT.mort)), lwd=5, col=col.pal[1])
lines(ppiT.mort ~ c(1:length(ppiT.mort)), lwd=5, col=col.pal[3])
lines(aosT.mort ~ c(1:length(aosT.mort)), lwd=5, col="grey")

lines(ppoS.mort ~ c(1:length(ppoS.mort)), lwd=5, col=col.pal[2], lty=2)
lines(sceS.mort ~ c(1:length(sceS.mort)), lwd=5, col=col.pal[1], lty=2)
lines(ppiS.mort ~ c(1:length(ppiS.mort)), lwd=5, col=col.pal[3], lty=2)

legend('topleft', col=col.pal[c(2,1,3)], bty="n", lwd=3, lty=1, title="Species",
       legend=c(
               "P. polystachya (non-swamp generalist)",
               "S. ceylanica (swamp generalist)",
               "P. pinnata (swamp specialist)"
       ))

plot(usual$recs.ppoS ~ c(1:length(usual$recs.ppoS)), lwd=5, col=col.pal[2], type="l",
     ylab="Annual recruitment volume", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=c(0, 
            max(c(usual$recs.ppoS, usual$recs.sceS, usual$recs.ppiS))))
lines(usual$recs.sceS ~ c(1:length(usual$recs.sceS)), lwd=5, col=col.pal[1])
lines(usual$recs.ppiS ~ c(1:length(usual$recs.ppiS)), lwd=5, col=col.pal[3])

####################
# INDIVIDUAL SIZES #
####################

#pdf(file = "D:\\Dropbox\\twn idiwn\\Post doc\\IBM temp\\sizes 3spp Aug21.pdf", width=16, height=9)
par(mfrow=c(1,2), mar=c(5.5,5.5,2,2))
plot(usual$z.ppoT ~ c(1:length(usual$z.ppoT)), lwd=5, col=col.pal[2], type="l",
     ylab="Mean adult DBH (log-transformed)", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=c(min(c(usual$z.ppoT, usual$z.sceT, usual$z.ppiT), na.rm=T), 
            max(c(usual$z.ppoT, usual$z.sceT, usual$z.ppiT), na.rm=T)))
lines(usual$z.sceT ~ c(1:length(usual$z.sceT)), lwd=5, col=col.pal[1])
lines(usual$z.ppiT ~ c(1:length(usual$z.ppiT)), lwd=5, col=col.pal[3])
legend('topleft', col=col.pal[c(2,1,3)], bty="n", lwd=3, lty=1, title="Species",
       legend=c(
               "P. polystachya (non-swamp generalist)",
               "S. ceylanica (swamp generalist)",
               "P. pinnata (swamp specialist)"
       ))
plot(usual$h.ppoS ~ c(1:length(usual$h.ppoS)), lwd=5, col=col.pal[2], type="l",
     ylab="Mean seedling height (log-transformed)", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=c(min(c(usual$h.ppoS, usual$h.sceS, usual$h.ppiS)), 
            max(c(usual$h.ppoS, usual$h.sceS, usual$h.ppiS))))
lines(usual$h.sceS ~ c(1:length(usual$h.sceS)), lwd=5, col=col.pal[1])
lines(usual$h.ppiS ~ c(1:length(usual$h.ppiS)), lwd=5, col=col.pal[3])
dev.off()

########################
# HABITAT ASSOCIATIONS #
########################

# Calculate SSI of all four species
ppo.ssi <- with(usual$ba.ppoT, swamp/landscape.swamp.prop / ((swamp/landscape.swamp.prop) + (nonswamp/(1-landscape.swamp.prop))))
sce.ssi <- with(usual$ba.sceT, swamp/landscape.swamp.prop / ((swamp/landscape.swamp.prop) + (nonswamp/(1-landscape.swamp.prop))))
ppi.ssi <- with(usual$ba.ppiT, swamp/landscape.swamp.prop / ((swamp/landscape.swamp.prop) + (nonswamp/(1-landscape.swamp.prop))))

#pdf(file = "D:\\Dropbox\\twn idiwn\\Post doc\\IBM temp\\ssi 3spp Aug21.pdf", width=6, height=8)
par(mfrow=c(1,1), mar=c(5,5,2,2))
plot(ppo.ssi ~ c(1:length(ppo.ssi)), lwd=5, col=col.pal[2], type="l",
     ylab="Swamp specialization index", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=c(0,1))
rect(0, 0.5, 25, 1.1, col = adjustcolor(col.t[5], alpha.f = 0.4), border=F)
abline(h=0.5, lty=2, lwd=3)
lines(sce.ssi ~ c(1:length(sce.ssi)), lwd=5, col=col.pal[1], type="l")
lines(ppi.ssi ~ c(1:length(ppi.ssi)), lwd=5, col=col.pal[3], type="l")
legend(x=1, y=0.8, col=col.pal[c(2,1,3)], lwd=3, lty=1, title="Species", bg="white",
       legend=c(
               "P. polystachya (non-swamp generalist)",
               "S. ceylanica (swamp generalist)",
               "P. pinnata (swamp specialist)"
       ))
dev.off()

######################
# SIZE DISTRIBUTIONS #
######################

# sample the height/dbh distributions instead of plotting all
n=200

# calculate average transition size of seedling (DBH=1cm)
tran.parm <- read.csv("data/transition params CAA Jun21.csv", header=T, row.names=1)
rec.parm <- t(read.csv("data/dispersal kernel parameters Apr21.csv", header=T, row.names=1))
tran.parm["transition.height",] <- (log(1) - tran.parm["tran.int",]) / tran.parm["tran.h",]

#pdf(file = "D:\\Dropbox\\twn idiwn\\Post doc\\IBM temp\\size distributions 3spp Aug21.pdf", width=16, height=9)
par(mfrow=c(1,2), mar=c(5,5,2,2))
hist(sample(usual$ppoS$logheight, n), col=col.t[2], main="", xlab="Seedling height (log-transformed)",
     xlim=c(2,6), breaks=seq(0,8,0.25))
hist(sample(usual$sceS$logheight, n), col=col.t[1], add=T, breaks=seq(0,8,0.25))
hist(sample(usual$ppiS$logheight, n), col=col.t[3], add=T, breaks=seq(0,8,0.25))
abline(v=tran.parm["transition.height", c("Gironniera.nervosa", "Prunus.polystachya", "Strombosia.ceylanica", "Pometia.pinnata")],
       col=col.pal[c(4,2,1,3)], lwd=3, lty=2)
abline(v=rec.parm["rcsz", c("Gironniera.nervosa", "Prunus.polystachya", "Strombosia.ceylanica", "Pometia.pinnata")],
       col=col.pal[c(4,2,1,3)], lwd=3)
legend('topleft', bty="n", lwd=3, lty=c(1,2), legend=c("average recruitment size", "average transition size"))
hist(sample(usual$ppiT$logdbh, n), col=col.t[3], main="", xlab="Adult DBH (log-transformed)", 
     xlim=c(0,5), breaks=seq(0,8,0.25))
hist(sample(usual$sceT$logdbh, n), col=col.t[1], add=T, breaks=seq(0,8,0.25))
hist(sample(usual$ppoT$logdbh, n), col=col.t[2], add=T, breaks=seq(0,8,0.25))
legend('topright', fill=col.t[c(2,1,3)], title="Species", bty="n",
       legend=c(
               "P. polystachya (non-swamp generalist)",
               "S. ceylanica (swamp generalist)",
               "P. pinnata (swamp specialist)"
       ))
dev.off()



