#################
# PLOT OUTPUTS  #
#################

library(raster)

# read output
out <- readRDS("out/sim_out.rds")

# Define a colour palette
col.pal <- c("#E3C16F", "#946846", "#FAFF70", "#6D213C", "#65DEF1", "#F3F5F6", "#BAAB68", "#CBC5EA")
col.t <- adjustcolor(col.pal, alpha.f=0.6)

with(out, 
     { 
             par(mfrow=c(1,2), mar=c(4,4,2,2))
             # initial condition
             plot(nssf.m, legend=F, col=col.pal[c(6,5)], main="Initial condition")
             # PPO
             points(ppoT.init, cex=ppoT$logdbh, col=col.t[1], pch=16)
             points(data.frame(ppoS.init), col=col.t[1], pch=4, cex=0.1)
             # SCE
             points(sceT.init, cex=sceT$logdbh, col=col.t[2], pch=16)
             points(data.frame(sceS.init), col=col.t[2], pch=4, cex=0.1)
             legend('bottomleft', bg="white", legend=c("Adult tree", "Seedling"), 
                    pch=c(1,4), pt.cex=c(3,0.1), col=c("black","grey"))
             # All other spp
             points(aosT, cex=aosT$logdbh,  pch=1)
     })

# final condition
with(out, 
     {       
             plot(nssf.m, col=col.pal[c(6,5)], main=paste0("After ", time, " years"), legend=F)
             # PPO
             points(ppoT, cex=ppoT$logdbh, col=col.t[1], pch=16)
             points(data.frame(ppoS), col=col.t[1], pch=4, cex=0.1)
             # SCE
             points(sceT, cex=sceT$logdbh, col=col.t[2], pch=16)
             points(data.frame(sceS), col=col.t[2], pch=4, cex=0.1)
             scalebar(100, xy=c(367200, 152900), type="bar", lonlat=F, below="metres", divs=4)
             # All other spp
             points(aosT, cex=aosT$logdbh,  pch=1)
             legend('bottomleft', bg="white", legend=c("Prunus", "Strombosia", "All other spp."), 
                    pch=c(16,16,1), pt.cex=2, col=c(col.t[1:2],"black"))
     })

with(out, 
     {  
             par(mfrow=c(1,2), mar=c(5.5,5.5,2,2))
             plot(n.ppoT ~ c(1:length(n.ppoT)), lwd=5, col=col.pal[1], type="l",
                  ylab="Adult population size", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
                  ylim=c(min(c(n.ppoT,n.sceT)), max(c(n.ppoT,n.sceT))))
             lines(n.sceT ~ c(1:length(n.sceT)), lwd=5, col=col.pal[2])
             legend('topright', 
                    legend=c("Prunus", "Strombosia"), title="Species", 
                    col=col.pal[1:2], lwd=3, cex=1.5, bty="n")
             plot(n.ppoS ~ c(1:length(n.ppoS)), lwd=5, col=col.pal[1], type="l",
                  ylab="Seedling population size", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
                  ylim=c(min(c(n.ppoS,n.sceS)), max(c(n.ppoS,n.sceS))))
             lines(n.sceS ~ c(1:length(n.sceS)), lwd=5, col=col.pal[2])
     })

with(out, 
     {  
             par(mfrow=c(1,2), mar=c(5.5,5.5,2,2))
             plot(z.ppoT ~ c(1:length(z.ppoT)), lwd=5, col=col.pal[1], type="l",
                  ylab="Mean adult DBH (log-transformed)", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
                  ylim=c(min(c(z.ppoT,z.sceT)), max(c(z.ppoT,z.sceT))))
             lines(z.sceT ~ c(1:length(z.sceT)), lwd=5, col=col.pal[2])
             plot(h.ppoS ~ c(1:length(h.ppoS)), lwd=5, col=col.pal[1], type="l",
                  ylab="Mean seedling height (log-transformed)", xlab="Time (years)", cex.lab=1.5, cex.axis=1.5,
                  ylim=c(min(c(h.ppoS,h.sceS)), max(c(h.ppoS,h.sceS))))
             lines(h.sceS ~ c(1:length(h.sceS)), lwd=5, col=col.pal[2])
             legend('bottomright', 
                    legend=c("Prunus", "Strombosia"), title="Species", 
                    col=col.pal[1:2], lwd=3, cex=1.5, bty="n")
     })

par(mfrow=c(1,2))
hist(out$h.ppoS, col=col.t[1], main="", xlab="Seedling height (log-transformed)")
hist(out$h.sceS, col=col.t[2], add=T)
hist(out$z.ppoT, col=col.t[1], main="", xlab="Adult DBH (log-transformed)")
hist(out$z.sceT, col=col.t[2], add=T)
legend('topright', fill=col.t[1:2], legend=c("Prunus", "Strombosia"), title="Species", bty="n")

plot(out$ba.ppoT)
plot(out$ba.sceT)
