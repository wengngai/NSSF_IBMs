#################
# PLOT OUTPUTS  #
#################

usual <- readRDS(file = "out/sim_out3_usual.rds")
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


# Plot out 10% of population
#pdf(file = "D:\\Dropbox\\twn idiwn\\Post doc\\IBM temp\\distributions 3spp Aug21.pdf", width=16, height=9)
# Initial condition
par(mfrow=c(1,2), mar=c(5,5,3,2))
plot(nssf.usual[[1]], legend=F, col=col.pal[c(6,5)], main="2021")
# PPO
points(usual$ppoT.init[sample(nrow(usual$ppoT.init), size=nrow(usual$ppoT.init)/10),], cex=usual$ppoT.init$logdbh/2, col=col.t[1], pch=16)
# SCE
points(usual$sceT.init[sample(nrow(usual$sceT.init), size=nrow(usual$sceT.init)/10),], cex=usual$sceT.init$logdbh/2, col=col.t[2], pch=16)
# PPI
points(usual$ppiT.init[sample(nrow(usual$ppiT.init), size=nrow(usual$ppiT.init)/10),], cex=usual$ppiT.init$logdbh/2, col=col.t[3], pch=16)
# All other species (static--just for a background interspecific competition pressure)
#points(aosT, cex=aosT$logdbh/2,  pch=1, col="grey70")
legend('topleft', col=col.t[c(1,2,3)], pch=16, bty="n", pt.cex=2,
       legend=c(
               "P. polystachya (non-swamp generalist)",
               "S. ceylanica (swamp generalist)",
               "P. pinnata (swamp specialist)"
               ))

# Final condition
plot(nssf.usual[[22]], legend=F, col=col.pal[c(6,5)], main="2042")
# PPO
points(usual$ppoT[sample(nrow(usual$ppoT), size=nrow(usual$ppoT)/10),], cex=usual$ppoT$logdbh/2, col=col.t[1], pch=16)
# SCE
points(usual$sceT[sample(nrow(usual$sceT), size=nrow(usual$sceT)/10),], cex=usual$sceT$logdbh/2, col=col.t[2], pch=16)
# PPI
points(usual$ppiT[sample(nrow(usual$ppiT), size=nrow(usual$ppiT)/10),], cex=usual$ppiT$logdbh/2, col=col.t[3], pch=16)
# All other species (static--just for a background interspecific competition pressure)
#points(aosT, cex=aosT$logdbh/2,  pch=1, col="grey70")
scalebar(400, xy=c(367800, 151550), type="bar", lonlat=F, below="metres", divs=2)
dev.off()

#pdf(file = "D:\\Dropbox\\twn idiwn\\Post doc\\IBM temp\\populations 3spp Aug21.pdf", width=16, height=9)
par(mfrow=c(1,2), mar=c(5.5,5.5,2,2))
plot(log(usual$n.ppoT) ~ c(1:length(usual$n.ppoT)), lwd=5, col=col.pal[2], type="l",
     ylab="Adult population size (log-transformed)", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=log(c(min(c(usual$n.ppoT, usual$n.sceT, usual$n.ppiT)), 
            max(c(usual$n.ppoT, usual$n.sceT, usual$n.ppiT))))
     )
lines(log(usual$n.sceT) ~ c(1:length(usual$n.sceT)), lwd=5, col=col.pal[1])
lines(log(usual$n.ppiT) ~ c(1:length(usual$n.ppiT)), lwd=5, col=col.pal[3])
legend('topleft', col=col.pal[c(2,1,3)], bty="n", lwd=3, lty=1, title="Species",
       legend=c(
               "P. polystachya (non-swamp generalist)",
               "S. ceylanica (swamp generalist)",
               "P. pinnata (swamp specialist)"
       ))
plot(log(usual$n.ppoS) ~ c(1:length(usual$n.ppoS)), lwd=5, col=col.pal[2], type="l",
     ylab="Seedling population size (log-transformed)", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=log(c(min(c(usual$n.ppoS, usual$n.sceS, usual$n.ppiS)), 
            max(c(usual$n.ppoS, usual$n.sceS, usual$n.ppiS))))
     )
lines(log(usual$n.sceS) ~ c(1:length(usual$n.sceS)), lwd=5, col=col.pal[1])
lines(log(usual$n.ppiS) ~ c(1:length(usual$n.ppiS)), lwd=5, col=col.pal[3])
dev.off()

#pdf(file = "D:\\Dropbox\\twn idiwn\\Post doc\\IBM temp\\sizes 3spp Aug21.pdf", width=16, height=9)
par(mfrow=c(1,2), mar=c(5.5,5.5,2,2))
plot(usual$z.ppoT ~ c(1:length(usual$z.ppoT)), lwd=5, col=col.pal[2], type="l",
     ylab="Mean adult DBH (log-transformed)", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=c(min(c(usual$z.ppoT, usual$z.sceT, usual$z.ppiT)), 
            max(c(usual$z.ppoT, usual$z.sceT, usual$z.ppiT))))
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

# sample the height/dbh distributions instead of plotting all
n=500

# calculate average transition size of seedling (DBH=1cm)
tran.parm <- read.csv("data/transition params CAA Jun21.csv", header=T, row.names=1)
rec.parm <- t(read.csv("data/dispersal kernel parameters Apr21.csv", header=T, row.names=1))
tran.parm["transition.height",] <- (log(1) - tran.parm["tran.int",]) / tran.parm["tran.h",]

#pdf(file = "D:\\Dropbox\\twn idiwn\\Post doc\\IBM temp\\size distributions 3spp Aug21.pdf", width=16, height=9)
par(mfrow=c(1,2), mar=c(5,5,2,2))
hist(sample(usual$ppoS$logheight, n), col=col.t[2], main="", xlab="Seedling height (log-transformed)",
     xlim=c(0,6), breaks=seq(0,8,0.5))
hist(sample(usual$sceS$logheight, n), col=col.t[1], add=T, breaks=seq(0,8,0.5))
hist(sample(usual$ppiS$logheight, n), col=col.t[3], add=T, breaks=seq(0,8,0.5))
abline(v=tran.parm["transition.height", c("Gironniera.nervosa", "Prunus.polystachya", "Strombosia.ceylanica", "Pometia.pinnata")],
       col=col.pal[c(4,2,1,3)], lwd=3, lty=2)
abline(v=rec.parm["rcsz", c("Gironniera.nervosa", "Prunus.polystachya", "Strombosia.ceylanica", "Pometia.pinnata")],
       col=col.pal[c(4,2,1,3)], lwd=3)
legend('topleft', bty="n", lwd=3, lty=c(1,2), legend=c("average recruitment size", "average transition size"))
hist(sample(usual$ppiT$logdbh, n), col=col.t[3], main="", xlab="Adult DBH (log-transformed)", 
     xlim=c(0,6), breaks=seq(0,8,0.5))
hist(sample(usual$sceT$logdbh, n), col=col.t[1], add=T, breaks=seq(0,8,0.5))
hist(sample(usual$ppoT$logdbh, n), col=col.t[2], add=T, breaks=seq(0,8,0.5))
legend('topright', fill=col.t[c(4,2,1,3)], title="Species", bty="n",
       legend=c(
               "P. polystachya (non-swamp generalist)",
               "S. ceylanica (swamp generalist)",
               "P. pinnata (swamp specialist)"
       ))
dev.off()



