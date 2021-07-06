#################
# SETUP BASEMAP #
#################

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

#par(mfrow=c(4,6), mar=c(1.5,1,2,0))
#for(i in 1:22) plot(nssf.usual[[i]], main=2020+i, legend=F)
#par(mfrow=c(4,6), mar=c(1.5,1,2,0))
#for(i in 1:22) plot(nssf.extreme[[i]], main=2020+i, legend=F)
#par(mfrow=c(4,6), mar=c(1.5,1,2,0))
#for(i in 1:22) plot(nssf.low[[i]], main=2020+i, legend=F)
#par(mfrow=c(4,6), mar=c(1.5,1,2,0))
#for(i in 1:22) plot(nssf.high[[i]], main=2020+i, legend=F)
#par(mfrow=c(2,2), mar=c(1.5,1,2,0))
#plot(nssf.usual[[22]], main="As usual", legend=F)
#plot(nssf.extreme[[22]], main="Extreme drying", legend=F)
#plot(nssf.low[[22]], main="Low rainfall", legend=F)
#plot(nssf.high[[22]], main="High rainfall", legend=F)

## initialize PPO
init.ppo.n <- 4000

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
ppoT.ras <- sampleRandom(nssf.usual[[1]], size=length(ppoT.logdbh), na.rm=T, asRaster=T)
ppoT.ras[!is.na(ppoT.ras)] <- ppoT.logdbh
ppoT <- data.frame(rasterToPoints(ppoT.ras, spatial=F))
names(ppoT) <- c("x", "y", "logdbh")

# create random points for ppoS, assign them logheight values
ppoS.ras <- sampleRandom(nssf.usual[[1]], size=length(ppoS.logheight), na.rm=T, asRaster=T)
ppoS.ras[!is.na(ppoS.ras)] <- ppoS.logheight
ppoS <- data.frame(rasterToPoints(ppoS.ras, spatial=F))
names(ppoS) <- c("x", "y", "logheight")

## initialize SCE
init.sce.n <- 6000

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
sceT.ras <- sampleRandom(nssf.usual[[1]], size=length(sceT.logdbh), na.rm=T, asRaster=T)
sceT.ras[!is.na(sceT.ras)] <- sceT.logdbh
sceT <- data.frame(rasterToPoints(sceT.ras, spatial=F))
names(sceT) <- c("x", "y", "logdbh")

# create random points for sceS, assign them logheight values
sceS.ras <- sampleRandom(nssf.usual[[1]], size=length(sceS.logheight), na.rm=T, asRaster=T)
sceS.ras[!is.na(sceS.ras)] <- sceS.logheight
sceS <- data.frame(rasterToPoints(sceS.ras, spatial=F))
names(sceS) <- c("x", "y", "logheight")

## initialize "All other species" for background interspecific competition
init.aos.n <- 20000

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
aosT.ras <- sampleRandom(nssf.usual[[1]], size=length(aosT.logdbh), na.rm=T, asRaster=T)
aosT.ras[!is.na(aosT.ras)] <- aosT.logdbh
aosT <- data.frame(rasterToPoints(aosT.ras, spatial=F))
names(aosT) <- c("x", "y", "logdbh")

# create random points for aosS, assign them logheight values
aosS.ras <- sampleRandom(nssf.usual[[1]], size=length(aosS.logheight), na.rm=T, asRaster=T)
aosS.ras[!is.na(aosS.ras)] <- aosS.logheight
aosS <- data.frame(rasterToPoints(aosS.ras, spatial=F))
names(aosS) <- c("x", "y", "logheight")

# Set up grid for reducing memory consumption when computing inter-individual distances

# create low resolution raster (100 x 100m)
nssf.100 <- aggregate(nssf.usual[[1]], fact=5)
res(nssf.100)
# assign cell values = index
nssf.100 <- setValues(nssf.100, 1:length(nssf.100))
# create list of indexed grid neighbours (NOT including focal grid cell)
grid.neighbours <- list()
for(i in 1:length(nssf.100)) grid.neighbours[[i]] <- adjacent(nssf.100, i, directions=8, matrix=F)[,2]

# Add another attribute (column) to individual data
ppoT$grid <- extract(nssf.100, ppoT[,c("x","y")])
ppoS$grid <- extract(nssf.100, ppoS[,c("x","y")])

sceT$grid <- extract(nssf.100, sceT[,c("x","y")])
sceS$grid <- extract(nssf.100, sceS[,c("x","y")])

aosT$grid <- extract(nssf.100, aosT[,c("x","y")])
aosS$grid <- extract(nssf.100, aosS[,c("x","y")])


# Define a colour palette
col.pal <- c("#E3C16F", "#946846", "#FAFF70", "#6D213C", "#65DEF1", "#F3F5F6", "#BAAB68", "#CBC5EA")
col.t <- adjustcolor(col.pal, alpha.f=0.6)

plot(nssf.usual[[1]], legend=F, col=col.pal[c(6,5)])
# PPO
points(ppoT, cex=ppoT$logdbh, col=col.t[1], pch=16)
points(data.frame(ppoS), col=col.t[1], pch=4, cex=0.1)
# SCE
points(sceT, cex=sceT$logdbh, col=col.t[2], pch=16)
points(data.frame(sceS), col=col.t[2], pch=4, cex=0.1)
# All other species (static--just for a background interspecific competition pressure)
points(aosT, cex=aosT$logdbh,  pch=1, col="grey70")
scalebar(500, xy=c(368500, 152000), type="bar", lonlat=F, below="metres", divs=2)
