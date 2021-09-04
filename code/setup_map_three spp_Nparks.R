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

# to crop landscape to an area in the centre of x by y m

nssf.bbox <- bbox(nssf.usual[[1]])

if(is.null(crop.dim)) {
    crop.extent <- nssf.bbox
    crop.area <- prod(nssf.bbox[,2]-nssf.bbox[,1])
} else {
    crop.centre <- apply(nssf.bbox, 1, mean)
    crop.extent <- extent(
        crop.centre[1]-crop.dim[1]/2,
        crop.centre[1]+crop.dim[1]/2,
        crop.centre[2]-crop.dim[2]/2,
        crop.centre[2]+crop.dim[2]/2
    )
    crop.area <- crop.dim[1]*crop.dim[2]
}

#plot(nssf.usual[[1]])
#rect(crop.centre[1]-crop.dim[1]/2, crop.centre[2]-crop.dim[2]/2, crop.centre[1]+crop.dim[1]/2, crop.centre[2]+crop.dim[2]/2)

for(i in 1:22){
    nssf.usual[[i]] <- crop(nssf.usual[[i]], crop.extent)
    nssf.extreme[[i]] <- crop(nssf.extreme[[i]], crop.extent)
    nssf.low[[i]] <- crop(nssf.low[[i]], crop.extent)
    nssf.high[[i]] <- crop(nssf.high[[i]], crop.extent)
}

##################################
# Estimate baseline stem numbers #
#  for initializing populations  #
##################################

# read in survey stem densities of the four focal spp
distrib <- read.csv("data/Size distribution.csv")

# Among trees (>5cm DBH), the summed BA in swamp/non-swamp plots for each sp:
BA.plots <- matrix(c(0, 7373.44, 2408.03, 11833.72, 1147.71, 2552.13), 
       ncol=2, dimnames=list(c("PPI","PPO", "SCE"), c("non-swamp", "swamp")))

# And the BA relative to total BA across all plots of that type
BA.plots.prop <- BA.plots / matrix(rep(c(207239.22, 244171.15), each=3), 3, 2)

# proportion of the landscape that is swamp/non-swamp:
if(scenario=="usual")  prop.swamp <- sum(values(nssf.usual[[1]]), na.rm=T)/length(na.omit(values(nssf.usual[[1]])))
if(scenario=="extreme")  prop.swamp <- sum(values(nssf.extreme[[1]]), na.rm=T)/length(na.omit(values(nssf.extreme[[1]])))
prop.ns <- 1-prop.swamp

# density of stems across NSSF estimated from raw data (stems per square metres):
# trees: 0.1351282
# saplings: 0.4200
# seedlings: 4.425926
# seedlings are likely overestimated because seedling quadrats are not randomly placed (maybe 2x overestimated?)
# therefore true density probably about 2.2 stems/m2

n.trees <- 0.1351282 * crop.area / dilution
n.saplings <- 0.42 * crop.area / dilution
n.seedlings <- 2.2 * crop.area / dilution

n.init.trees <- round(BA.plots.prop * matrix(rep(c(prop.ns, prop.swamp), each=3), 3, 2) * n.trees, 0)
n.init.saplings <- round(BA.plots.prop * matrix(rep(c(prop.ns, prop.swamp), each=3), 3, 2) * n.saplings, 0)
n.init.seedlings <- round(BA.plots.prop * matrix(rep(c(prop.ns, prop.swamp), each=3), 3, 2) * n.seedlings, 0)
# cannot have zeroes
n.init.trees[1,1] <- 1 
n.init.saplings[1,1] <- 1 
n.init.seedlings[1,1] <- 1

# Generating random initial locations

# convert baselayer from raster to polygons for random point generation
if(scenario=="usual"){
    swamp.poly <- rasterToPolygons(nssf.usual[[1]], fun=function(x){x==1})
    nonswamp.poly <- rasterToPolygons(nssf.usual[[1]], fun=function(x){x==0})
}
if(scenario=="extreme"){
    swamp.poly <- rasterToPolygons(nssf.extreme[[1]], fun=function(x){x==1})
    nonswamp.poly <- rasterToPolygons(nssf.extreme[[1]], fun=function(x){x==0})
}

# abbreviated species names for for loop
tran.abbrev <- c("AAN", "AFR", "ALU", "ACL","BBR", "BPA", "CRU", "CSQ", "CBR", 
                   "EMA", "GPA", "GNE", "GWA", "GAX", "KMA", "MBA", "PAX", "PPI", "PPO", 
                   "PCO", "PEC", "SCE", "SPA", "TFL", "TWA", "XFL", "AOS")
for(i in 1:3){
    # abbreviated sp name
    sp <- rownames(n.init.trees)[i]
    
    n.swamp.adults <- round(n.init.trees[i,"swamp"], 0)
    n.swamp.saplings <- round(n.init.saplings[i,"swamp"], 0)
    n.swamp.seedlings <- round(n.init.seedlings[i,"swamp"], 0)
    n.nonswamp.adults <- round(n.init.trees[i,"non-swamp"], 0)
    n.nonswamp.saplings <- round(n.init.saplings[i,"non-swamp"], 0)
    n.nonswamp.seedlings <- round(n.init.seedlings[i,"non-swamp"], 0)
    
    stopifnot(
        "Some populations have zero count! spsample wil fail. Check if dilution is too high?" = all(
            n.swamp.adults > 0,
            n.swamp.saplings > 0,
            n.swamp.seedlings > 0,
            n.nonswamp.adults > 0,
            n.nonswamp.saplings > 0,
            n.nonswamp.seedlings > 0
        )
    )
    
    # sample swamp/non-swamp areas according to adults (adults + saplings) and seedling groups
    swamp.pts.T <- spsample(swamp.poly, n.swamp.adults+n.swamp.saplings, type = 'random')
    nonswamp.pts.T <- spsample(nonswamp.poly, n.nonswamp.adults+n.nonswamp.saplings, type = 'random')
    swamp.pts.S <- spsample(nonswamp.poly, n.swamp.seedlings, type = 'random')
    nonswamp.pts.S <- spsample(nonswamp.poly, n.nonswamp.seedlings, type = 'random')
    
    # combine random coordinates with sizes drawn from raw data
    assign(paste0(tolower(sp),"T"), rbind(
        data.frame(coordinates(swamp.pts.T), 
              logdbh=log(c(
                  sample(distrib$dbh[which(distrib$species==sp & distrib$cat=="tree")], size=n.swamp.adults, replace=T),
                  sample(distrib$dbh[which(distrib$species==sp & distrib$cat=="sapling")], size=n.swamp.saplings, replace=T)
                  ))
              ),
        data.frame(coordinates(nonswamp.pts.T),
              logdbh=log(c(
                  sample(distrib$dbh[which(distrib$species==sp & distrib$cat=="tree")], size=n.nonswamp.adults, replace=T),
                  sample(distrib$dbh[which(distrib$species==sp & distrib$cat=="sapling")], size=n.nonswamp.saplings, replace=T)
                  ))
              )
        ))
    
    # for seedlings, draw from raw data first, transform (DBH to height) later
    z <- log(c(
        sample(distrib$dbh[which(distrib$species==sp & distrib$cat=="seedling")], size=n.swamp.seedlings, replace=T),
        sample(distrib$dbh[which(distrib$species==sp & distrib$cat=="seedling")], size=n.nonswamp.seedlings, replace=T)
        ))
    
    # combine random coordinates with sizes drawn from raw data
    assign(paste0(tolower(sp),"S"),rbind(
        data.frame(
            rbind(coordinates(swamp.pts.S), coordinates(nonswamp.pts.S)),
            # transformation of logDBH to logheight
            logheight = (z - tran.parm["tran.int",which(tran.abbrev==sp)]) / tran.parm["tran.h",which(tran.abbrev==sp)]
        )
    ))
    
}

#par(mfrow=c(1,3)); hist(ppoT[,"logdbh"], main="PPO"); hist(sceT[,"logdbh"], main="SCE"); hist(ppiT[,"logdbh"], main="PPI")

## initialize "All other species" for background interspecific competition

# Approx. number of all other species stems (ADULTS ONLY--using 1:3:20 adults to saplings to seedlings ratio)
init.aos.n <- round(n.trees - sum(n.init.trees), 0)

# Use power law with 1cm min DBH for all other spp. Note that seedlings not modelled
init.aos.dbh <- rplcon(init.aos.n, 1, 2.2)
aosT.logdbh <- log(init.aos.dbh)
#hist(aosT.logdbh)

# create random points for aosT, assign them logdbh values
crop.poly <- as(nssf.usual[[1]], 'SpatialPolygons')
aos.loc <- spsample(crop.poly, init.aos.n, type = 'random')
aosT <- data.frame(coordinates(aos.loc), logdbh = aosT.logdbh)
aosS.loc <- spsample(crop.poly, 10, type = 'random')
aosS <- data.frame(coordinates(aosS.loc), logheight = log(20))


###################################
# Set up grid for reducing memory #
#   consumption when computing    #
#   inter-individual distances    #
###################################

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

ppiT$grid <- extract(nssf.100, ppiT[,c("x","y")])
ppiS$grid <- extract(nssf.100, ppiS[,c("x","y")])

aosT$grid <- extract(nssf.100, aosT[,c("x","y")])
aosS$grid <- extract(nssf.100, aosS[,c("x","y")])

#############
# Visualize #
#############

# # Define a colour palette
# col.pal <- c("#E3C16F", "#946846", "#FAFF70", "#6D213C", "#65DEF1", "#F3F5F6", "#BAAB68", "#CBC5EA")
# col.t <- adjustcolor(col.pal, alpha.f=0.6)

# pdf(file = "D:\\Dropbox\\twn idiwn\\Post doc\\IBM temp\\Fig 2-3.pdf", width=14, height=9)
 par(mfrow=c(1,2))
 plot(fullmap, legend=F, col=col.pal[c(6,5)])
 rect(367000, 153000, 367500, 153500)
 scalebar(500, xy=c(368000, 151650), type="bar", lonlat=F, below="metres", divs=2)
 mtext(side=3, adj=0, text="a) Full landscape", cex=1.5)
 
 plot(nssf.usual[[1]], legend=F, col=col.pal[c(6,5)])
 # All other species (static--just for a background interspecific competition pressure)
 points(aosT, cex=aosT$logdbh/2,  pch=1, col="grey90")
 # PPO
 points(ppoT, cex=ppoT$logdbh/2, col=col.t[1], pch=16)
# points(data.frame(ppoS), col=col.t[1], pch=4, cex=0.2)
 # SCE
 points(sceT, cex=sceT$logdbh/2, col=col.t[2], pch=16)
# points(data.frame(sceS), col=col.t[2], pch=4, cex=0.2)
 # PPI
 points(ppiT, cex=ppiT$logdbh/2, col=col.t[4], pch=16)
# points(data.frame(ppiS), col=col.t[4], pch=4, cex=0.2)
 scalebar(100, xy=c(367350, 152970), type="bar", lonlat=F, below="metres", divs=2)
 mtext(side=3, adj=0, text="b) Simulation extent", cex=1.5)
dev.off()
