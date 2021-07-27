#################
# PLOT OUTPUTS  #
#################

usual <- readRDS(file = "out/sim_out.rds")
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
#pdf(file = "D:\\Dropbox\\twn idiwn\\Post doc\\IBM temp\\distributions Jul21.pdf", width=16, height=9)
# Initial condition
par(mfrow=c(1,2), mar=c(5,5,3,2))
plot(nssf.usual[[1]], legend=F, col=col.pal[c(6,5)], main="2021")
# GNE
points(usual$gneT.init[sample(nrow(usual$gneT.init), size=nrow(usual$gneT.init)/10),], cex=usual$gneT.init$logdbh/2, col=col.t[4], pch=16)
# PPO
points(usual$ppoT.init[sample(nrow(usual$ppoT.init), size=nrow(usual$ppoT.init)/10),], cex=usual$ppoT.init$logdbh/2, col=col.t[1], pch=16)
# SCE
points(usual$sceT.init[sample(nrow(usual$sceT.init), size=nrow(usual$sceT.init)/10),], cex=usual$sceT.init$logdbh/2, col=col.t[2], pch=16)
# PPI
points(usual$ppiT.init[sample(nrow(usual$ppiT.init), size=nrow(usual$ppiT.init)/10),], cex=usual$ppiT.init$logdbh/2, col=col.t[3], pch=16)
# All other species (static--just for a background interspecific competition pressure)
#points(aosT, cex=aosT$logdbh/2,  pch=1, col="grey70")
legend('topleft', col=col.t[c(4,1,2,3)], pch=16, bty="n", pt.cex=2,
       legend=c(
               "G. nervosa (non-swamp specialist)",
               "P. polystachya (non-swamp generalist)",
               "S. ceylanica (swamp generalist)",
               "P. pinnata (swamp specialist)"
               ))

# Final condition
plot(nssf.usual[[22]], legend=F, col=col.pal[c(6,5)], main="2042")
# GNE
points(usual$gneT[sample(nrow(usual$gneT), size=nrow(usual$gneT)/10),], cex=usual$gneT$logdbh/2, col=col.t[4], pch=16)
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

#pdf(file = "D:\\Dropbox\\twn idiwn\\Post doc\\IBM temp\\populations Jul21.pdf", width=16, height=9)
par(mfrow=c(1,2), mar=c(5.5,5.5,2,2))
plot(log(usual$n.ppoT) ~ c(1:length(usual$n.ppoT)), lwd=5, col=col.pal[2], type="l",
     ylab="Adult population size (log-transformed)", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=log(c(min(c(usual$n.ppoT, usual$n.sceT, usual$n.ppiT, usual$n.gneT)), 
            max(c(usual$n.ppoT, usual$n.sceT, usual$n.ppiT, usual$n.gneT))))
     )
lines(log(usual$n.sceT) ~ c(1:length(usual$n.sceT)), lwd=5, col=col.pal[1])
lines(log(usual$n.gneT) ~ c(1:length(usual$n.gneT)), lwd=5, col=col.pal[4])
lines(log(usual$n.ppiT) ~ c(1:length(usual$n.ppiT)), lwd=5, col=col.pal[3])
legend('topright', col=col.pal[c(4,2,1,3)], bty="n", lwd=3, lty=1, title="Species",
       legend=c(
               "G. nervosa (non-swamp specialist)",
               "P. polystachya (non-swamp generalist)",
               "S. ceylanica (swamp generalist)",
               "P. pinnata (swamp specialist)"
       ))
plot(log(usual$n.ppoS) ~ c(1:length(usual$n.ppoS)), lwd=5, col=col.pal[2], type="l",
     ylab="Seedling population size (log-transformed)", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=log(c(min(c(usual$n.ppoS, usual$n.sceS, usual$n.ppiS, usual$n.gneS)), 
            max(c(usual$n.ppoS, usual$n.sceS, usual$n.ppiS, usual$n.gneS))))
     )
lines(log(usual$n.sceS) ~ c(1:length(usual$n.sceS)), lwd=5, col=col.pal[1])
lines(log(usual$n.ppiS) ~ c(1:length(usual$n.ppiS)), lwd=5, col=col.pal[3])
lines(log(usual$n.gneS) ~ c(1:length(usual$n.gneS)), lwd=5, col=col.pal[4])
dev.off()

#pdf(file = "D:\\Dropbox\\twn idiwn\\Post doc\\IBM temp\\sizes Jul21.pdf", width=16, height=9)
par(mfrow=c(1,2), mar=c(5.5,5.5,2,2))
plot(usual$z.ppoT ~ c(1:length(usual$z.ppoT)), lwd=5, col=col.pal[2], type="l",
     ylab="Mean adult DBH (log-transformed)", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=c(min(c(usual$z.ppoT, usual$z.sceT, usual$z.ppiT, usual$z.gneT)), 
            max(c(usual$z.ppoT, usual$z.sceT, usual$z.ppiT, usual$z.gneT))))
lines(usual$z.sceT ~ c(1:length(usual$z.sceT)), lwd=5, col=col.pal[1])
lines(usual$z.ppiT ~ c(1:length(usual$z.ppiT)), lwd=5, col=col.pal[3])
lines(usual$z.gneT ~ c(1:length(usual$z.gneT)), lwd=5, col=col.pal[4])
legend('topleft', col=col.pal[c(4,2,1,3)], bty="n", lwd=3, lty=1, title="Species",
       legend=c(
               "G. nervosa (non-swamp specialist)",
               "P. polystachya (non-swamp generalist)",
               "S. ceylanica (swamp generalist)",
               "P. pinnata (swamp specialist)"
       ))
plot(usual$h.ppoS ~ c(1:length(usual$h.ppoS)), lwd=5, col=col.pal[2], type="l",
     ylab="Mean seedling height (log-transformed)", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=c(min(c(usual$h.ppoS, usual$h.sceS, usual$h.ppiS, usual$h.gneS)), 
            max(c(usual$h.ppoS, usual$h.sceS, usual$h.ppiS, usual$h.gneS))))
lines(usual$h.sceS ~ c(1:length(usual$h.sceS)), lwd=5, col=col.pal[1])
lines(usual$h.ppiS ~ c(1:length(usual$h.ppiS)), lwd=5, col=col.pal[3])
lines(usual$h.gneS ~ c(1:length(usual$h.gneS)), lwd=5, col=col.pal[4])
dev.off()


# Calculate SSI of all four species. Note that labels of the ba.sppT dataframes were wrongly named
gne.ssi <- 1-with(usual$ba.gneT, swamp*landscape.swamp.prop / ((swamp*landscape.swamp.prop) + (nonswamp*(1-landscape.swamp.prop))))
ppo.ssi <- 1-with(usual$ba.ppoT, swamp*landscape.swamp.prop / ((swamp*landscape.swamp.prop) + (nonswamp*(1-landscape.swamp.prop))))
sce.ssi <- 1-with(usual$ba.sceT, swamp*landscape.swamp.prop / ((swamp*landscape.swamp.prop) + (nonswamp*(1-landscape.swamp.prop))))
ppi.ssi <- 1-with(usual$ba.ppiT, swamp*landscape.swamp.prop / ((swamp*landscape.swamp.prop) + (nonswamp*(1-landscape.swamp.prop))))

#pdf(file = "D:\\Dropbox\\twn idiwn\\Post doc\\IBM temp\\ssi Jul21.pdf", width=6, height=8)
par(mfrow=c(1,1), mar=c(5,5,2,2))
plot(gne.ssi ~ c(1:length(gne.ssi)), lwd=5, col=col.pal[4], type="l",
     ylab="Swamp specialization index", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=c(0,1))
rect(0, 0.5, 25, 1.1, col = adjustcolor(col.t[5], alpha.f = 0.4), border=F)
abline(h=0.5, lty=2, lwd=3)
lines(ppo.ssi ~ c(1:length(ppo.ssi)), lwd=5, col=col.pal[2], type="l")
lines(sce.ssi ~ c(1:length(sce.ssi)), lwd=5, col=col.pal[1], type="l")
lines(ppi.ssi ~ c(1:length(ppi.ssi)), lwd=5, col=col.pal[3], type="l")
legend('bottomright', col=col.pal[c(4,2,1,3)], bty="n", lwd=3, lty=1, title="Species",
       legend=c(
               "G. nervosa (non-swamp specialist)",
               "P. polystachya (non-swamp generalist)",
               "S. ceylanica (swamp generalist)",
               "P. pinnata (swamp specialist)"
       ))
dev.off()







par(mfrow=c(1,2))
hist(usual$ppoS$logheight, col=col.t[2], main="", xlab="Seedling height (log-transformed)")
hist(usual$sceS$logheight, col=col.t[1], add=T)
hist(usual$gneS$logheight, col=col.t[4], add=T)
hist(usual$ppiS$logheight, col=col.t[3], add=T)
hist(usual$ppoT$logdbh, col=col.t[2], main="", xlab="Adult DBH (log-transformed)")
hist(usual$sceT$logdbh, col=col.t[1], add=T)
hist(usual$ppiT$logdbh, col=col.t[3], add=T)
hist(usual$gneT$logdbh, col=col.t[4], add=T)
legend('topright', fill=col.t[1:2], legend=c("Prunus", "Strombosia"), title="Species", bty="n")