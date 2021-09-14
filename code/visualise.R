#################
# PLOT OUTPUTS  #
#################

extreme <- readRDS(file = "out/sim_out3_extreme.rds")
summary(extreme)
usual <- readRDS(file = "out/sim_out3_usual_100yrs.rds")
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

#pdf(file = "D:\\Dropbox\\twn idiwn\\Post doc\\IBM temp\\distributions 3spp Aug21.pdf", width=14, height=9)
prop = 0.5

par(mfrow=c(1,2), mar=c(2,2,3,2))
## As usual scenario
# Initial condition
plot(nssf.usual[[1]], legend=F, col=col.pal[c(6,5)], yaxt="n", xaxt="n")
mtext(side=3, line=1, text="a) As usual scenario, year 2021", adj=0, cex=1.5)
# PPO
i.ppoT <- sample(1:nrow(usual$ppoT), prop*nrow(usual$ppoT))
points(usual$ppoT.init[i.ppoT,], cex=usual$ppoT.init$logdbh[i.ppoT]/2, col=col.t[1], pch=16)
# SCE
i.sceT <- sample(1:nrow(usual$sceT), prop*nrow(usual$sceT))
points(usual$sceT.init[i.sceT,], cex=usual$sceT.init$logdbh[i.sceT]/2, col=col.t[2], pch=16)
# PPI
i.ppiT <- sample(1:nrow(usual$ppiT), prop*nrow(usual$ppiT))
points(usual$ppiT.init[i.ppiT,], cex=usual$ppiT.init$logdbh[i.ppiT]/2, col=col.t[4], pch=16)
# All other species (static--just for a background iproperspecific competition pressure)
i.aosT <- sample(1:nrow(usual$aosT), prop*nrow(usual$aosT))
#points(aosT[i.aosT,], cex=aosT$logdbh[i.aosT]/2,  pch=1, col="grey70")
legend('topleft', col=col.t[c(1,2,4)], pch=16, pt.cex=3, y.intersp=1.5,
       legend=c(
               "P. polystachya\n(non-swamp generalist)",
               "S. ceylanica\n(swamp generalist)",
               "P. pinnata\n(swamp specialist)"
       ), bg="white", title="Species")

# Final condition
plot(nssf.usual[[22]], legend=F, col=col.pal[c(6,5)], xaxt="n", yaxt="n")
mtext(side=3, line=1, text="b) As usual scenario, year 2042", adj=0, cex=1.5)
# PPO
points(usual$ppoT[i.ppoT,], cex=usual$ppoT$logdbh[i.ppoT]/2, col=col.t[1], pch=16)
# SCE
points(usual$sceT[i.sceT,], cex=usual$sceT$logdbh[i.sceT]/2, col=col.t[2], pch=16)
# PPI
points(usual$ppiT[i.ppiT,], cex=usual$ppiT$logdbh[i.ppiT]/2, col=col.t[4], pch=16)
# All other species (static--just for a background interspecific competition pressure)
#points(aosT[i.aosT,], cex=aosT$logdbh[i.aosT]/2,  pch=1, col="grey70")

## Extreme drying scenario
# Initial condition
plot(nssf.extreme[[1]], legend=F, col=col.pal[c(6,5)], xaxt="n", yaxt="n")
mtext(side=3, line=1, text="c) Extreme drying scenario, year 2021", adj=0, cex=1.5)
# PPO
i.ppoT <- sample(1:nrow(extreme$ppoT), prop*nrow(extreme$ppoT))
points(extreme$ppoT.init[i.ppoT,], cex=extreme$ppoT.init$logdbh[i.ppoT]/2, col=col.t[1], pch=16)
# SCE
i.sceT <- sample(1:nrow(extreme$sceT), prop*nrow(extreme$sceT))
points(extreme$sceT.init[i.sceT,], cex=extreme$sceT.init$logdbh[i.sceT]/2, col=col.t[2], pch=16)
# PPI
i.ppiT <- sample(1:nrow(extreme$ppiT), prop*nrow(extreme$ppiT))
points(extreme$ppiT.init[i.ppiT,], cex=extreme$ppiT.init$logdbh[i.ppiT]/2, col=col.t[4], pch=16)
# All other species (static--just for a background iproperspecific competition pressure)
i.aosT <- sample(1:nrow(extreme$aosT), prop*nrow(extreme$aosT))
#points(aosT[i.aosT,], cex=aosT$logdbh[i.aosT]/2,  pch=1, col="grey70")


# Final condition
plot(nssf.extreme[[22]], legend=F, col=col.pal[c(6,5)], xaxt="n", yaxt="n")
mtext(side=3, line=1, text="d) Extreme drying scenario, year 2042", adj=0, cex=1.5)
# PPO
points(extreme$ppoT[i.ppoT,], cex=extreme$ppoT$logdbh[i.ppoT]/2, col=col.t[1], pch=16)
# SCE
points(extreme$sceT[i.sceT,], cex=extreme$sceT$logdbh[i.sceT]/2, col=col.t[2], pch=16)
# PPI
points(extreme$ppiT[i.ppiT,], cex=extreme$ppiT$logdbh[i.ppiT]/2, col=col.t[4], pch=16)
# All other species (static--just for a background interspecific competition pressure)
#points(aosT[i.aosT,], cex=aosT$logdbh[i.aosT]/2,  pch=1, col="grey70")
scalebar(100, xy=c(367525, 153025), type="bar", lonlat=F, below="metres", divs=2)

dev.off()


####################
# POPULATION SIZES #
####################

### BY SPECIES
#pdf(file = "D:\\Dropbox\\twn idiwn\\Post doc\\IBM temp\\populations 3spp Aug21.pdf", width=14, height=9)

## Adults
par(mfrow=c(2,3), mar=c(2,3.5,3.5,2), oma=c(3.5,6,1,1))

## PPO
plot(usual$n.ppoT ~ c(1:length(usual$n.ppoT)), lwd=5, col=col.pal[1], type="l", las=1,
     ylab="", xlab="", cex.lab=2, cex.axis=1.5,
     ylim=c(min(c(usual$n.ppoT, extreme$n.ppoT)), 
            max(c(usual$n.ppoT, extreme$n.ppoT)))
     )
lines(extreme$n.ppoT ~ c(1:length(extreme$n.ppoT)), lwd=5, col=col.pal[1], lty=2)
mtext(side=3, line=1, adj=0, text="a) Prunus polystachya (non-swamp generalist)", cex=1.2)
legend('topleft', bty="n", lwd=2, lty=c(1,2), title="Scenario", cex=1.5,
       legend=c("As usual", "Extreme drying"))

## SCE
plot(usual$n.sceT ~ c(1:length(usual$n.sceT)), lwd=5, col=col.pal[2], type="l", las=1,
     ylab="", xlab="", cex.lab=1.5, cex.axis=1.5,
     ylim=c(min(c(usual$n.sceT, extreme$n.sceT)), 
            max(c(usual$n.sceT, extreme$n.sceT)))
)
lines(extreme$n.sceT ~ c(1:length(extreme$n.sceT)), lwd=5, col=col.pal[2], lty=2)
mtext(side=3, line=1, adj=0, text="b) Strombosia ceylanica (swamp generalist)", cex=1.2)

## PPI
plot(usual$n.ppiT ~ c(1:length(usual$n.ppiT)), lwd=5, col=col.pal[4], type="l", las=1,
     ylab="", xlab="", cex.lab=1.5, cex.axis=1.5,
     ylim=c(min(c(usual$n.ppiT, extreme$n.ppiT)), 
            max(c(usual$n.ppiT, extreme$n.ppiT)))
)
lines(extreme$n.ppiT ~ c(1:length(extreme$n.ppiT)), lwd=5, col=col.pal[4], lty=2)
mtext(side=3, line=1, adj=0, text="c) Pometia pinnata (swamp specialist)", cex=1.2)
mtext(side=2, line=3.5, outer=T, adj=0.9, text="Adult population size", cex=1.5)

## AOS
#plot(usual$n.aosT ~ c(1:length(usual$n.aosT)), lwd=5, col="grey", type="l", las=1,
#     ylab="", xlab="", cex.lab=1.5, cex.axis=1.5,
#     ylim=c(min(c(usual$n.aosT, extreme$n.aosT)), 
#            max(c(usual$n.aosT, extreme$n.aosT)))
#)
#lines(extreme$n.aosT ~ c(1:length(extreme$n.aosT)), lwd=5, col="grey", lty=2)
#mtext(side=3, line=1, adj=0, text="All other species", cex=1.2)

## Seedlings

# PPO
plot(usual$n.ppoS ~ c(1:length(usual$n.ppoS)), lwd=5, col=col.pal[1], type="l", las=1,
     ylab="", xlab="", cex.lab=2, cex.axis=1.5,
     ylim=c(min(c(usual$n.ppoS, extreme$n.ppoS)), 
            max(c(usual$n.ppoS, extreme$n.ppoS)))
)
lines(extreme$n.ppoS ~ c(1:length(extreme$n.ppoS)), lwd=5, col=col.pal[1], lty=2)
mtext(side=3, line=1, adj=0, text="d) Prunus polystachya (non-swamp generalist)", cex=1.2)

# SCE
plot(usual$n.sceS ~ c(1:length(usual$n.sceS)), lwd=5, col=col.pal[2], type="l", las=1,
     ylab="", xlab="", cex.lab=1.5, cex.axis=1.5,
     ylim=c(min(c(usual$n.sceS, extreme$n.sceS)), 
            max(c(usual$n.sceS, extreme$n.sceS)))
)
lines(extreme$n.sceS ~ c(1:length(extreme$n.sceS)), lwd=5, col=col.pal[2], lty=2)
mtext(side=3, line=1, adj=0, text="e) Strombosia ceylanica (swamp generalist)", cex=1.2)

# PPI
plot(usual$n.ppiS ~ c(1:length(usual$n.ppiS)), lwd=5, col=col.pal[4], type="l", las=1,
     ylab="", xlab="", cex.lab=1.5, cex.axis=1.5,
     ylim=c(min(c(usual$n.ppiS, extreme$n.ppiS)), 
            max(c(usual$n.ppiS, extreme$n.ppiS)))
)
lines(extreme$n.ppiS ~ c(1:length(extreme$n.ppiS)), lwd=5, col=col.pal[4], lty=2)
mtext(side=3, line=1, adj=0, text="f) Pometia pinnata (swamp specialist)", cex=1.2)

mtext(side=2, line=3.5, outer=T, adj=0.15, text="Seedling population size", cex=1.5)
mtext(side=1, line=1.5, outer=T, adj=0.5, text="Time (years)", cex=1.5)

dev.off()

### BY SCENARIO
#pdf(file = "D:\\Dropbox\\twn idiwn\\Post doc\\IBM temp\\populations 3spp Aug21.pdf", width=12, height=9)
par(mfrow=c(2,2), mar=c(5.5,5.5,2,2))

## As usual scenario
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

## Extreme scenario
plot(log(extreme$n.ppoT) ~ c(1:length(extreme$n.ppoT)), lwd=5, col=col.pal[2], type="l",
     ylab="Adult population size (log-transformed)", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=log(c(min(c(extreme$n.ppoT, extreme$n.sceT, extreme$n.ppiT, extreme$n.aosT))+1, 
                max(c(extreme$n.ppoT, extreme$n.sceT, extreme$n.ppiT, extreme$n.aosT))))
)
lines(log(extreme$n.sceT) ~ c(1:length(extreme$n.sceT)), lwd=5, col=col.pal[1])
lines(log(extreme$n.ppiT) ~ c(1:length(extreme$n.ppiT)), lwd=5, col=col.pal[3])
lines(log(extreme$n.aosT) ~ c(1:length(extreme$n.aosT)), lwd=5, col="grey")
legend('topleft', col=c(col.pal[c(2,1,3)], "grey"), bty="n", lwd=3, lty=1, title="Species",
       legend=c(
               "P. polystachya (non-swamp generalist)",
               "S. ceylanica (swamp generalist)",
               "P. pinnata (swamp specialist)",
               "All other species"
       ))
plot(log(extreme$n.ppoS+1) ~ c(1:length(extreme$n.ppoS)), lwd=5, col=col.pal[2], type="l",
     ylab="Seedling population size (log-transformed)", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=log(c(min(c(extreme$n.ppoS, extreme$n.sceS, extreme$n.ppiS))+1, 
                max(c(extreme$n.ppoS, extreme$n.sceS, extreme$n.ppiS))))
)
lines(log(extreme$n.sceS+1) ~ c(1:length(extreme$n.sceS)), lwd=5, col=col.pal[1])
lines(log(extreme$n.ppiS+1) ~ c(1:length(extreme$n.ppiS)), lwd=5, col=col.pal[3])
dev.off()


#####################
# DEMOGRAPHIC RATES #
#####################

ppoT.mort.usual <- usual$deaths.ppoT / usual$n.ppoT[-length(usual$n.ppoT)]
ppoS.mort.usual <- usual$deaths.ppoS / usual$n.ppoS[-length(usual$n.ppoS)]
sceT.mort.usual <- usual$deaths.sceT / usual$n.sceT[-length(usual$n.sceT)]
sceS.mort.usual <- usual$deaths.sceS / usual$n.sceS[-length(usual$n.sceS)]
ppiT.mort.usual <- usual$deaths.ppiT / usual$n.ppiT[-length(usual$n.ppiT)]
ppiS.mort.usual <- usual$deaths.ppiS / usual$n.ppiS[-length(usual$n.ppiS)]
aosT.mort.usual <- usual$deaths.aosT / usual$n.aosT[-length(usual$n.aosT)]

ppoT.mort.extreme <- extreme$deaths.ppoT / extreme$n.ppoT[-length(extreme$n.ppoT)]
ppoS.mort.extreme <- extreme$deaths.ppoS / extreme$n.ppoS[-length(extreme$n.ppoS)]
sceT.mort.extreme <- extreme$deaths.sceT / extreme$n.sceT[-length(extreme$n.sceT)]
sceS.mort.extreme <- extreme$deaths.sceS / extreme$n.sceS[-length(extreme$n.sceS)]
ppiT.mort.extreme <- extreme$deaths.ppiT / extreme$n.ppiT[-length(extreme$n.ppiT)]
ppiS.mort.extreme <- extreme$deaths.ppiS / extreme$n.ppiS[-length(extreme$n.ppiS)]
aosT.mort.extreme <- extreme$deaths.aosT / extreme$n.aosT[-length(extreme$n.aosT)]


par(mfrow=c(2,2), mar=c(5.5,5.5,2,2))
plot(ppoT.mort.usual ~ c(1:length(ppoT.mort.usual)), lwd=5, col=col.pal[1], type="l",
     ylab="Annual mortality rate", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=c(0, max(c(ppoT.mort.usual, ppoS.mort.usual, sceT.mort.usual, sceS.mort.usual, ppiT.mort.usual, ppiS.mort.usual), na.rm=T)))
lines(sceT.mort.usual ~ c(1:length(sceT.mort.usual)), lwd=5, col=col.pal[2])
lines(ppiT.mort.usual ~ c(1:length(ppiT.mort.usual)), lwd=5, col=col.pal[4])
lines(aosT.mort.usual ~ c(1:length(aosT.mort.usual)), lwd=5, col="grey")

lines(ppoS.mort.usual ~ c(1:length(ppoS.mort.usual)), lwd=5, col=col.pal[1], lty=2)
lines(sceS.mort.usual ~ c(1:length(sceS.mort.usual)), lwd=5, col=col.pal[2], lty=2)
lines(ppiS.mort.usual ~ c(1:length(ppiS.mort.usual)), lwd=5, col=col.pal[4], lty=2)

legend('topleft', col=col.pal[c(2,1,3)], bty="n", lwd=3, lty=1, title="Species",
       legend=c(
               "P. polystachya (non-swamp generalist)",
               "S. ceylanica (swamp generalist)",
               "P. pinnata (swamp specialist)"
       ))

plot(usual$recs.ppoS ~ c(1:length(usual$recs.ppoS)), lwd=5, col=col.pal[1], type="l",
     ylab="Annual recruitment volume", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=c(0, 
            max(c(usual$recs.ppoS, usual$recs.sceS, usual$recs.ppiS))))
lines(usual$recs.sceS ~ c(1:length(usual$recs.sceS)), lwd=5, col=col.pal[2])
lines(usual$recs.ppiS ~ c(1:length(usual$recs.ppiS)), lwd=5, col=col.pal[4])

## Extreme drying scenario
plot(ppoT.mort.extreme ~ c(1:length(ppoT.mort.extreme)), lwd=5, col=col.pal[1], type="l",
     ylab="Annual mortality rate", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=c(0, max(c(ppoT.mort.extreme, ppoS.mort.extreme, sceT.mort.extreme, sceS.mort.extreme, ppiT.mort.extreme, ppiS.mort.extreme), na.rm=T)))
lines(sceT.mort.extreme ~ c(1:length(sceT.mort.extreme)), lwd=5, col=col.pal[2])
lines(ppiT.mort.extreme ~ c(1:length(ppiT.mort.extreme)), lwd=5, col=col.pal[4])
lines(aosT.mort.extreme ~ c(1:length(aosT.mort.extreme)), lwd=5, col="grey")

lines(ppoS.mort.extreme ~ c(1:length(ppoS.mort.extreme)), lwd=5, col=col.pal[1], lty=2)
lines(sceS.mort.extreme ~ c(1:length(sceS.mort.extreme)), lwd=5, col=col.pal[2], lty=2)
lines(ppiS.mort.extreme ~ c(1:length(ppiS.mort.extreme)), lwd=5, col=col.pal[4], lty=2)

legend('topleft', col=col.pal[c(2,1,3)], bty="n", lwd=3, lty=1, title="Species",
       legend=c(
               "P. polystachya (non-swamp generalist)",
               "S. ceylanica (swamp generalist)",
               "P. pinnata (swamp specialist)"
       ))

plot(extreme$recs.ppoS ~ c(1:length(extreme$recs.ppoS)), lwd=5, col=col.pal[1], type="l",
     ylab="Annual recruitment volume", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=c(0, 
            max(c(extreme$recs.ppoS, extreme$recs.sceS, extreme$recs.ppiS))))
lines(extreme$recs.sceS ~ c(1:length(extreme$recs.sceS)), lwd=5, col=col.pal[2])
lines(extreme$recs.ppiS ~ c(1:length(extreme$recs.ppiS)), lwd=5, col=col.pal[4])
dev.off()

####################
# INDIVIDUAL SIZES #
####################

### BY SPECIES
#pdf(file = "D:\\Dropbox\\twn idiwn\\Post doc\\IBM temp\\sizes 3spp Aug21.pdf", width=14, height=9)

## Adults
par(mfrow=c(2,3), mar=c(2,3.5,3.5,2), oma=c(3.5,6,1,1))

## PPO
plot(usual$z.ppoT ~ c(1:length(usual$z.ppoT)), lwd=5, col=col.pal[1], type="l", las=1,
     ylab="", xlab="", cex.lab=2, cex.axis=1.5,
     ylim=c(min(c(usual$z.ppoT, extreme$z.ppoT)), 
            max(c(usual$z.ppoT, extreme$z.ppoT)))
)
lines(extreme$z.ppoT ~ c(1:length(extreme$z.ppoT)), lwd=5, col=col.pal[1], lty=2)
mtext(side=3, line=1, adj=0, text="a) Prunus polystachya (non-swamp generalist)", cex=1.2)
legend('bottomleft', bty="n", lwd=2, lty=c(1,2), title="Scenario", cex=1.5,
       legend=c("As usual", "Extreme drying"))

## SCE
plot(usual$z.sceT ~ c(1:length(usual$z.sceT)), lwd=5, col=col.pal[2], type="l", las=1,
     ylab="", xlab="", cex.lab=1.5, cex.axis=1.5,
     ylim=c(min(c(usual$z.sceT, extreme$z.sceT)), 
            max(c(usual$z.sceT, extreme$z.sceT)))
)
lines(extreme$z.sceT ~ c(1:length(extreme$z.sceT)), lwd=5, col=col.pal[2], lty=2)
mtext(side=3, line=1, adj=0, text="b) Strombosia ceylanica (swamp generalist)", cex=1.2)

## PPI
plot(usual$z.ppiT ~ c(1:length(usual$z.ppiT)), lwd=5, col=col.pal[4], type="l", las=1,
     ylab="", xlab="", cex.lab=1.5, cex.axis=1.5,
     ylim=c(min(c(usual$z.ppiT, extreme$z.ppiT)), 
            max(c(usual$z.ppiT, extreme$z.ppiT)))
)
lines(extreme$z.ppiT ~ c(1:length(extreme$z.ppiT)), lwd=5, col=col.pal[4], lty=2)
mtext(side=3, line=1, adj=0, text="c) Pometia pinnata (swamp specialist)", cex=1.2)
mtext(side=2, line=1, outer=T, adj=0.8, text="Mean adult DBH\n(log-transformed)", cex=1.5)

## Seedlings

# PPO
plot(usual$h.ppoS ~ c(1:length(usual$h.ppoS)), lwd=5, col=col.pal[1], type="l", las=1,
     ylab="", xlab="", cex.lab=2, cex.axis=1.5,
     ylim=c(min(c(usual$h.ppoS, extreme$h.ppoS)), 
            max(c(usual$h.ppoS, extreme$h.ppoS)))
)
lines(extreme$h.ppoS ~ c(1:length(extreme$h.ppoS)), lwd=5, col=col.pal[1], lty=2)
mtext(side=3, line=1, adj=0, text="d) Prunus polystachya (non-swamp generalist)", cex=1.2)

# SCE
plot(usual$h.sceS ~ c(1:length(usual$h.sceS)), lwd=5, col=col.pal[2], type="l", las=1,
     ylab="", xlab="", cex.lab=1.5, cex.axis=1.5,
     ylim=c(min(c(usual$h.sceS, extreme$h.sceS)), 
            max(c(usual$h.sceS, extreme$h.sceS)))
)
lines(extreme$h.sceS ~ c(1:length(extreme$h.sceS)), lwd=5, col=col.pal[2], lty=2)
mtext(side=3, line=1, adj=0, text="e) Strombosia ceylanica (swamp generalist)", cex=1.2)

# PPI
plot(usual$h.ppiS ~ c(1:length(usual$h.ppiS)), lwd=5, col=col.pal[4], type="l", las=1,
     ylab="", xlab="", cex.lab=1.5, cex.axis=1.5,
     ylim=c(min(c(usual$h.ppiS, extreme$h.ppiS)), 
            max(c(usual$h.ppiS, extreme$h.ppiS)))
)
lines(extreme$h.ppiS ~ c(1:length(extreme$h.ppiS)), lwd=5, col=col.pal[4], lty=2)
mtext(side=3, line=1, adj=0, text="f) Pometia pinnata (swamp specialist)", cex=1.2)

mtext(side=2, line=1, outer=T, adj=0.15, text="Mean seedling height\n(log-transformed)", cex=1.5)
mtext(side=1, line=1.5, outer=T, adj=0.5, text="Time (years)", cex=1.5)

dev.off()

#pdf(file = "D:\\Dropbox\\twn idiwn\\Post doc\\IBM temp\\sizes 3spp Aug21.pdf", width=14, height=9)
par(mfrow=c(1,2), mar=c(5.5,5.5,2,2))
plot(usual$z.ppoT ~ c(1:length(usual$z.ppoT)), lwd=5, col=col.pal[1], type="l",
     ylab="Mean adult DBH (log-transformed)", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=c(min(c(usual$z.ppoT, usual$z.sceT, usual$z.ppiT), na.rm=T), 
            max(c(usual$z.ppoT, usual$z.sceT, usual$z.ppiT), na.rm=T)))
lines(usual$z.sceT ~ c(1:length(usual$z.sceT)), lwd=5, col=col.pal[2])
lines(usual$z.ppiT ~ c(1:length(usual$z.ppiT)), lwd=5, col=col.pal[4])
legend('bottomleft', col=col.pal[c(1,2,4)], bty="n", lwd=3, lty=1, title="Species",
       legend=c(
               "P. polystachya (non-swamp generalist)",
               "S. ceylanica (swamp generalist)",
               "P. pinnata (swamp specialist)"
       ))
plot(usual$h.ppoS ~ c(1:length(usual$h.ppoS)), lwd=5, col=col.pal[1], type="l",
     ylab="Mean seedling height (log-transformed)", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=c(min(c(usual$h.ppoS, usual$h.sceS, usual$h.ppiS)), 
            max(c(usual$h.ppoS, usual$h.sceS, usual$h.ppiS))))
lines(usual$h.sceS ~ c(1:length(usual$h.sceS)), lwd=5, col=col.pal[2])
lines(usual$h.ppiS ~ c(1:length(usual$h.ppiS)), lwd=5, col=col.pal[4])
dev.off()

########################
# HABITAT ASSOCIATIONS #
########################

# Calculate SSI of all four species
ppo.ssi.usual <- with(usual$ba.ppoT, swamp/landscape.swamp.prop / ((swamp/landscape.swamp.prop) + (nonswamp/(1-landscape.swamp.prop))))
sce.ssi.usual <- with(usual$ba.sceT, swamp/landscape.swamp.prop / ((swamp/landscape.swamp.prop) + (nonswamp/(1-landscape.swamp.prop))))
ppi.ssi.usual <- with(usual$ba.ppiT, swamp/landscape.swamp.prop / ((swamp/landscape.swamp.prop) + (nonswamp/(1-landscape.swamp.prop))))

ppo.ssi.extreme <- with(extreme$ba.ppoT, swamp/landscape.swamp.prop / ((swamp/landscape.swamp.prop) + (nonswamp/(1-landscape.swamp.prop))))
sce.ssi.extreme <- with(extreme$ba.sceT, swamp/landscape.swamp.prop / ((swamp/landscape.swamp.prop) + (nonswamp/(1-landscape.swamp.prop))))
ppi.ssi.extreme <- with(extreme$ba.ppiT, swamp/landscape.swamp.prop / ((swamp/landscape.swamp.prop) + (nonswamp/(1-landscape.swamp.prop))))

#pdf(file = "D:\\Dropbox\\twn idiwn\\Post doc\\IBM temp\\ssi 3spp Aug21.pdf", width=6, height=8)
par(mfrow=c(1,1), mar=c(5,5,2,2))
plot(ppo.ssi.usual ~ c(1:length(ppo.ssi.usual)), lwd=5, col=col.pal[1], type="l",
     ylab="Swamp specialization index", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
     ylim=c(0,1))
rect(0, 0.5, 25, 1.1, col = adjustcolor(col.t[5], alpha.f = 0.4), border=F)
abline(h=0.5, lty=2, lwd=3)
lines(sce.ssi.usual ~ c(1:length(sce.ssi.usual)), lwd=5, col=col.pal[2], type="l")
lines(ppi.ssi.usual ~ c(1:length(ppi.ssi.usual)), lwd=5, col=col.pal[4], type="l")

lines(ppo.ssi.extreme ~ c(1:length(ppo.ssi.extreme)), lwd=5, col=col.pal[1], type="l", lty=2)
lines(sce.ssi.extreme ~ c(1:length(sce.ssi.extreme)), lwd=5, col=col.pal[2], type="l", lty=2)
lines(ppi.ssi.extreme ~ c(1:length(ppi.ssi.extreme)), lwd=5, col=col.pal[4], type="l", lty=2)

legend(x=1, y=0.85, col=col.pal[c(1,2,4)], lwd=5, lty=1, title="Species", bg="white",
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
n=100

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



