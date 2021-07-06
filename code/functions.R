#################################
# DEFINING VITAL RATE FUNCTIONS #
#################################

# Note on notations:
# Seedling phase size units are log-transformed height, h
# Adult phase size units are log-transformed DBH, z
# At the end of each year interval, sizes are denoted by h1 and z1

# intraS, intraA denote symmetric and asymmetric intraspecific competition, respectively
# interS, interA, interspecific
# CAE (conspecific adult effect) and HAE (heterospecific) denote similar, but are used for seedling growth only
# also, for adult tree growth, inter/intra are calculated with distance from focal tree incorporated


## TREE Survival function, logistic regression
# note: need to add interS effect eventually

sT_z <- function(terrain, trees, sp, intraS, intraA, interS, interA)
{
    m.par <- surv.parm
    z <- trees[,"logdbh"]
    plot_type <- ifelse(extract(terrain, trees[,c("x","y")])==2, "wet", "dry")
    plot_type[is.na(plot_type)] <- "dry"
    
    # scale intraA, intraS, interS
    intraA.scaled <- (log(intraA + surv.parm.unscale["intraA.unscale.logoffset",]) -
                          surv.parm.unscale["intraA.unscale.mu",]) / surv.parm.unscale["intraA.unscale.sigma",]
    intraS.scaled <- (log(intraS + surv.parm.unscale["intraS.unscale.logoffset",]) -
                          surv.parm.unscale["intraS.unscale.mu",]) / surv.parm.unscale["intraS.unscale.sigma",]
    interS.scaled <- (log(interS) -
                          surv.parm.unscale["interS.unscale.mu",]) / surv.parm.unscale["interS.unscale.sigma",]
    # calculate interA*HM, then scale it
    interAHM <- interA * HM.parm[paste0(plot_type,".hm.stem"), sp]  # use HM.stem for this calculation for now
    interAHM.scaled <- (log(interAHM + surv.parm.unscale["interAHM.unscale.logoffset",]) 
                        - surv.parm.unscale["interAHM.unscale.mu",]) /  surv.parm.unscale["interAHM.unscale.sigma",]
    
    # calculate p1 and p2 from these
    p1 <- m.par["p1", sp] + (intraA.scaled * m.par["p1.intraA", sp])
    r1 <- m.par["r1", sp] + (interAHM.scaled * m.par["r1.interAHM", sp])
    p2 <- m.par["p2", sp] + (intraS.scaled * m.par["p2.intraS", sp]) + (interS.scaled * m.par["p2.interS", sp])
    
    surv.prob <- ifelse( z <= 1,
                         m.par["K",sp] / ( 1 + exp(-r1 * (z - p1)) ),
                         m.par["K",sp] / ( 1 + exp(-m.par["r2",sp] * (z - p2)) )
    )
    
    # return binary outcome
    surv <- rbinom(n = length(z), prob = surv.prob, size = 1)
    surv[is.na(surv)] <- 0
    return(surv)
}
#sT_z(nssf.m, ppoT, "Prunus.polystachya", 1, 1, 1, 1)
#sT_z(nssf.m, ppoT, "Prunus.polystachya", intra.calc(ppoT, ppoS, "Prunus.polystachya")[[1]])

## TREE Growth function
# this is the only function that uses dbh instead of logdbh
# note: need to add interS effect eventually

GT_z1z <- function(terrain, trees, sp, intraA, interS)
{
    m.par <- grow.T.parm
    z <- trees[,"logdbh"]
    dbh <- exp(z)
    plot_type <- ifelse(extract(terrain, trees[,c("x","y")])==2, "wet", "dry")
    plot_type[is.na(plot_type)] <- "dry"
    
    # scale intraA and interS
    intraA.scaled <- (log(intraA + grow.T.parm.unscale["grow.intraA.unscale.logoffset",]) -
                          grow.T.parm.unscale["grow.intraA.unscale.mu",]) / grow.T.parm.unscale["grow.intraA.unscale.sigma",]
    interS.scaled <- (log(interS) - grow.T.parm.unscale["grow.interS.unscale.mu",]) / 
        grow.T.parm.unscale["grow.interS.unscale.sigma",]
    interS.scaled[which(interS.scaled < -2)] <- -2 # in case too little interspecific competition causes infinite growth
    
    # calculate param b and c values
    b <- m.par["b", sp] + (intraA.scaled * m.par["b.intraA", sp]) + (interS.scaled * m.par["b.interS", sp])
    c <- m.par["c", sp] + (intraA.scaled * m.par["c.intraA", sp]) + (interS.scaled * m.par["c.interS", sp])
    
    loggrowth <- (
        m.par["a",sp] * 
            (dbh ^ b) *
            exp(-dbh * c)
    ) + rnorm(n = length(z), mean = 0, sd = m.par["sigma",sp]) # add the noise in (sigma)
    
    dbh1 <- dbh + (exp(loggrowth) - 1)
    dbh1[which(dbh1 < 1)] <- 1
    z1 <- log(dbh1)
    names(z1) <- NULL
    return(z1)
}
#plot(GT_z1z(nssf.m, ppoT, "Prunus.polystachya", intra.dist.calc(ppoT), inter.dist.calc(rbind(aosT,sceT), ppoT)) ~ ppoT[,"logdbh"], cex=3)
#abline(0,1, lty=2)

## FLOWERING function, logistic regression

p_bz <- function(terrain, trees, sp)
{
    m.par <- fruit.parm
    z <- trees[,"logdbh"]
    plot_type <- ifelse(extract(terrain, trees[,c("x","y")])==2, "wet", "dry")
    plot_type[is.na(plot_type)] <- "dry"
    
    linear.p <- m.par[paste0("fruit.int.",plot_type), sp] + m.par["fruit.z", sp] * z      # linear predictor
    
    # return binary outcome
    surv <- rbinom(n = length(linear.p), prob = 1/(1+exp(-linear.p)), size = 1)	# logistic transformation to probability
    surv[is.na(surv)] <- 0
    return(surv)
}
#p_bz(nssf.m, ppoT, "Prunus.polystachya")
#p_bz(nssf.m, sceT, "Strombosia.ceylanica")

## SEEDLING production function (this is seedling recruitment, not seed production)

b_z <- function(terrain, trees, sp)
{
    # first few lines can be modified for later expansion
    m.par <- fruit.parm
    z <- trees[,"logdbh"]
    plot_type <- ifelse(extract(terrain, trees[,c("x","y")])==2, "wet", "dry")
    plot_type[is.na(plot_type)] <- "dry"
    
    N <- m.par["rec.int",sp] + m.par["fruited.z",sp] * z 	# seed production of a size z plant
    N[is.na(N)] <- 0
    return(round(exp(N),0))						# convert log(recruitment) to actual recruitment no.
}
#b_z(nssf.m, ppoT, "Prunus.polystachya")
#p_bz(nssf.m, ppoT, "Prunus.polystachya")*b_z(nssf.m, ppoT, "Prunus.polystachya")

## SEEDLING recruit size pdf

c_0h1 <- function(n, sp)
{
    # first few lines can be modified for later expansion
    m.par <- rec.parm
    h1 <- rnorm(n, mean = m.par["rcsz",sp], sd = m.par["rcsd",sp])
    return(h1)
}
#c_0h1(100, "Prunus.polystachya")

## SEEDLING Growth function
# note: need to add HAE effect eventually

GS_h1h <- function(terrain, seedlings, sp, CAE, HAE)
{
    # first few lines can be modified for later expansion
    m.par <- grow.S.parm
    h <- seedlings[,"logheight"]
    CAE.scaled <- ( log(CAE+10) - m.par["grow.S.CAE.unscale.mu",sp] ) / m.par["grow.S.CAE.unscale.sigma",sp]
    HAE.scaled <- ( log(HAE) - m.par["grow.S.HAE.unscale.mu",sp] ) / m.par["grow.S.HAE.unscale.sigma",sp]
    HAE.scaled[which(HAE.scaled < -2)] <- -2 # replace HAE with a min value if too low
    
    mu <- m.par["grow.S.int",sp] + 					# intercept
        m.par["grow.S.CAE", sp] * CAE.scaled +			# CAE
        m.par["grow.S.HAE", sp] * HAE.scaled +			# HAE
        m.par["grow.S.h", sp] * h +        				# height
        m.par["grow.S.CAExh", sp] * h * CAE.scaled +			# CAE interaction
        m.par["grow.S.HAExh", sp] * h * HAE.scaled			# HAE interaction
    
    h1 <- mu + rnorm(n = length(mu), mean = 0, sd = m.par["grow.S.sigma",sp])
    
    names(h1) <- NULL
    return(h1)
}
#plot(GS_h1h(nssf.m, ppoS, "Prunus.polystachya", 10, 10) ~ ppoS[,"logheight"])
#points(GS_h1h(nssf.m, ppoS, "Prunus.polystachya", 1000) ~ ppoS[,"logheight"], col="red")
#abline(0,1,lty=2,lwd=2)

## SEEDLING Survival function, logistic regression

sS_h <- function(terrain, seedlings, sp, intraA, interA)
{
    m.par <- surv.parm
    h <- seedlings[,"logheight"]
    plot_type <- ifelse(extract(terrain, seedlings[,c("x","y")])==2, "wet", "dry")
    plot_type[is.na(plot_type)] <- "dry"
    
    # transform logheight to logdbh
    z <- tran.parm["tran.int", sp] + h*tran.parm["tran.h", sp]
    
    # scale intraA, intraS, interS
    intraA.scaled <- (log(intraA + surv.parm.unscale["intraA.unscale.logoffset",]) -
                          surv.parm.unscale["intraA.unscale.mu",]) / surv.parm.unscale["intraA.unscale.sigma",]
    # calculate interA*HM, then scale it
    interAHM <- interA * HM.parm[paste0(plot_type,".hm.stem"), sp]  # use HM.stem for this calculation for now
    interAHM.scaled <- (log(interAHM + surv.parm.unscale["interAHM.unscale.logoffset",]) 
                        - surv.parm.unscale["interAHM.unscale.mu",]) /  surv.parm.unscale["interAHM.unscale.sigma",]
    
    # calculate p1 and p2 from these
    p1 <- m.par["p1", sp] + (intraA.scaled * m.par["p1.intraA", sp])
    r1 <- m.par["r1", sp] + (interAHM.scaled * m.par["r1.interAHM", sp])
    
    surv.prob <- m.par["K", sp] / ( 1 + exp(-r1 * (z - p1)) )
    
    # return binary outcome
    surv <- rbinom(n = length(z), prob = surv.prob, size = 1)
    surv[is.na(surv)] <- 0
    return(surv)
}
#sS_h(nssf.m, ppoS, "Prunus.polystachya", 1, 1) 
#sS_h(nssf.m, ppoS, "Prunus.polystachya", intra.calc(ppoT, ppoS, "Prunus.polystachya")[[2]]) * GS_h1h(nssf.m, ppoS, "Prunus.polystachya", 100)
#plot(jitter(sS_h(nssf.m, ppoS, "Prunus.polystachya")) ~ sqrt(ppoS[,"logheight"]))

## SEEDLING-TO-SAPLING transition function

# takes h1 from seedlings, thus seedlings df must have already completed growth/mortality
T_z1h1 <- function(trees, seedlings, sp)
{
    m.par <- tran.parm
    h1 <- seedlings[,"logheight"]
    
    # compute logdbh from logheight
    z1 <- (m.par["tran.int", sp] + m.par["tran.h", sp] * h1)
    
    # only initiate transitions if there is at least one seedling larger than 1cm DBH
    if(sum(z1>0) > 0){
        tran.index <- which(z1 > 0) # z1 > 0 because log(1) = 0
        # add sapling to tree df
        newsaplings <- cbind(seedlings[tran.index, c("x","y")], z1[tran.index], seedlings[tran.index, "grid"])
        names(newsaplings) <- names(trees)
        trees <-  rbind(trees, newsaplings)
        # remove seedling from seedling stack
        seedlings <- seedlings[-tran.index,]
    }
    return(list(trees, seedlings))
}
# ppoTS <- T_z1h1(ppoT, ppoS, "Prunus.polystachya")
# nrow(ppoTS[[1]]); nrow(ppoTS[[2]])

## DISPERSAL KERNEL function

# define a 2DT pdf
twoDT.pdf <- function(dist, S, L) { S / 
        (pi * L * (1 + (dist^2/L))^(S+1))
}

# this function takes input of germinant number, generates vector of distances from mother tree 
# up to 100m from mother tree, and to a spatial resolution of 0.1m (because of discretization of pdf in "support")
twoDT.sample <- function(n, sp, max_dist = 100) {
    support <- seq(0.1, max_dist, 0.1)
    m.par <- rec.parm
    probs <- twoDT.pdf(support, m.par["S",sp], m.par["L",sp])
    probs <- probs/sum(probs) # normalize probability 
    sample(
        support, 
        n, # n draws
        TRUE,  # with replacement
        probs # using these probabilities
    )
}
# par(mfrow = c(1, 2))
# hist(twoDT.sample(200, "Gironniera.nervosa", max_dist = 100), breaks=100)
# hist(twoDT.sample(200, "Gironniera.nervosa", max_dist = 410), breaks=100)

##################################################
# FUNCTIONS FOR EXTRACTING COMPETITION VARIABLES #
##################################################

## Conspecific adult density around seedlings (CAE)
# CAE is sum BA of all adult stems in 20x20 plot
# we want to convert this 400m^2 into a circle about each focal seedling
# 400 = pi*r^2
#CAE.calc <- function(trees, seedlings, r = sqrt(400/pi)){
    # spDists returns distance matrix with rows as trees and cols as seedlings
#    Ts.in.r <- ifelse(spDists(as.matrix(trees), as.matrix(seedlings)) < r, 1, 0)
#    raw.CAE <- as.vector(pi*(exp(trees$logdbh)/2)^2) %*% Ts.in.r
#    return(as.vector(raw.CAE))
#}
#CAE.calc(ppoT, ppoS)

# new CAE.calc
CAE.calc <- function(trees, seedlings, r = sqrt(400/pi)){
    raw.CAE <- rep(NA, nrow(seedlings))
    for(i in 1:length(nssf.100)){
        # seedlings in grid i
        focal <- seedlings[seedlings$grid==i,]
        # trees in neighbouring grids
        neighbours <- rbind(
            trees[trees$grid==i,],
            trees[trees$grid %in% grid.neighbours[[i]],]
        )
        # compute distance matrix only if there are both seedlings and trees in the grids
        if(nrow(focal)>0 & nrow(neighbours)>0){
            dist.mat <- spDists(
                as.matrix(focal[,c("x","y")]), 
                as.matrix(neighbours[,c("x","y")])
            )
            Ts.in.r <- ifelse(dist.mat < r, 1, 0)
            raw.CAE[which(seedlings$grid==i)] <-  Ts.in.r %*% as.vector(pi*(exp(neighbours$logdbh)/2)^2)
            # if only have seedlings (but no trees), assign zero CAE. if not, ignore
        } else if (nrow(focal)>0 & nrow(neighbours)==0) {
            raw.CAE[which(seedlings$grid==i)] <- 0
        }
    }
    return(raw.CAE)
}

## Heterospecific adult density around seedlings (HAE)
# HAE is sum BA of all non focal species adult stems in 20x20 plot
# we want to convert this 400m^2 into a circle about each focal seedling
# 400 = pi*r^2
#HAE.calc <- function(trees.hetero, seedlings.con, r = sqrt(400/pi)){
    # spDists returns distance matrix with rows as trees and cols as seedlings
#    Ts.in.r <- ifelse(spDists(as.matrix(trees.hetero), as.matrix(seedlings.con)) < r, 1, 0)
#    raw.HAE <- as.vector(pi*(exp(trees.hetero$logdbh)/2)^2) %*% Ts.in.r
#    return(as.vector(raw.HAE))
#}
#HAE.calc(sceT, ppoS)

# new HAE.calc
HAE.calc <- function(trees.hetero, seedlings.con, r = sqrt(400/pi)){
    raw.HAE <- rep(NA, nrow(seedlings.con))
    for(i in 1:length(nssf.100)){
        # seedlings in grid i
        focal <- seedlings.con[seedlings.con$grid==i,]
        # trees in neighbouring grids
        neighbours <- rbind(
            trees.hetero[trees.hetero$grid==i,],
            trees.hetero[trees.hetero$grid %in% grid.neighbours[[i]],]
        )
        # compute distance matrix only if there are both seedlings and trees in the grids
        if(nrow(focal)>0 & nrow(neighbours)>0){
            dist.mat <- spDists(
                as.matrix(focal[,c("x","y")]), 
                as.matrix(neighbours[,c("x","y")])
            )
            Ts.in.r <- ifelse(dist.mat < r, 1, 0)
            raw.HAE[which(seedlings.con$grid==i)] <-  Ts.in.r %*% as.vector(pi*(exp(neighbours$logdbh)/2)^2)
            # if only have seedlings (but no trees), assign zero CAE. if not, ignore
        } else if (nrow(focal)>0 & nrow(neighbours)==0) {
            raw.HAE[which(seedlings.con$grid==i)] <- 0
        }
    }
    return(raw.HAE)
}
#HAE.calc(sceT, ppoS)

## Intraspecific competition (intra)
# almost same as above, except includes seedlings
# returns intraS (symmetric) for adults [[1]] and seedlings [[2]]
# and intraA (asymmetric) for adults [[3]] and seedlings [[4]]
#intra.calc <- function(trees, seedlings, sp, r = sqrt(400/pi)){
    # convert seedling heights to dbh
#    h <- seedlings$logheight
#    z <- tran.parm["tran.int",sp] + h*tran.parm["tran.h",sp]
    # concatenate adult and seedling zs
#    z <- c(trees$logdbh, z)
    
    # rbind coordinates of trees and seedlings
#    ALL <- as.matrix(rbind(trees[,1:2], seedlings[,1:2]))
    # calc distance and find all within 400m2 radius of focal
#    dists <- spDists(ALL)
#    Ts.in.r <- ifelse(dists < sqrt(400/pi), 1, 0)
#    diag(Ts.in.r) <- 0
    # for asymmetric competition, we only want pairs in which focal is the larger of the two, so need to compute another matrix
#    zz <- matrix(rep(z, each=length(z)), ncol=length(z))
#    big.mat <- ifelse(t(zz)>zz, 1, 0)
#    pairs.to.count <- Ts.in.r * big.mat
    # use matrix multiplication to obtain summed intra for each individual
#    raw.intraS <- pi*(exp(z)/2)^2 %*% Ts.in.r
#    raw.intraA <- pi*(exp(z)/2)^2 %*% pairs.to.count
#    return(list(
        # intraS for adults [[1]] and seedlings [[2]]
#        as.vector(raw.intraS)[1:nrow(trees)], as.vector(raw.intraS)[(nrow(trees)+1):(nrow(trees)+nrow(seedlings))],
        # intraA for adults [[3]] and seedlings [[4]]
#        as.vector(raw.intraA)[1:nrow(trees)], as.vector(raw.intraA)[(nrow(trees)+1):(nrow(trees)+nrow(seedlings))]
#    ))
#}

# New intra.calc
intra.calc <- function(trees, seedlings, sp, r = sqrt(400/pi)){
    # convert seedling heights to dbh
    seedlings$logdbh <- tran.parm["tran.int",sp] + seedlings$logheight*tran.parm["tran.h",sp]
    
    T.raw.intraS <- rep(NA, nrow(trees))
    T.raw.intraA <- rep(NA, nrow(trees))
    S.raw.intraS <- rep(NA, nrow(seedlings))
    S.raw.intraA <- rep(NA, nrow(seedlings))
    
    for(i in 1:length(nssf.100)){
        # need to separate focal trees from focal seedlings so that matching back is easier
        focalTs <- trees[trees$grid==i,]
        focalSs <- seedlings[seedlings$grid==i,]
        # neighbours can be rbind-ed
        neighbours <- rbind(
            trees[trees$grid %in% grid.neighbours[[i]], c("x","y","logdbh")],
            seedlings[seedlings$grid %in% grid.neighbours[[i]], c("x","y","logdbh")]
        )
        # compute the distance matrices
        if(nrow(focalTs > 0)){
            # neighbours to trees: first elements must be focal trees themselves (form the diagonal of the spDists matrix)
            # then, seedlings within the focal grid, then everything else in neighbouring grids
            Tneighbours <- rbind(focalTs[, c("x","y","logdbh")], focalSs[,c("x","y","logdbh")], neighbours)
            dist.Ts <- spDists(as.matrix(focalTs[,1:2]), as.matrix(Tneighbours[,1:2]))
            Ts.in.r <- ifelse(dist.Ts < sqrt(400/pi), 1, 0)
            diag(Ts.in.r) <- 0
            
            # for asymmetric competition, we only want pairs in which focal is the smaller of the two, so need to compute another matrix
            T.big.mat <- ifelse(
                matrix(rep(focalTs$logdbh, nrow(Tneighbours)), nrow=nrow(focalTs)) <
                    matrix(rep(Tneighbours$logdbh, each=nrow(focalTs)), nrow=nrow(focalTs))
                , 1, 0)
            T.pairs.to.count <- Ts.in.r * T.big.mat
            
            # use matrix multiplication to obtain summed intra for each individual
            T.raw.intraS[trees$grid==i] <- Ts.in.r %*% (pi*(exp(Tneighbours$logdbh)/2)^2)
            T.raw.intraA[trees$grid==i] <- T.pairs.to.count %*% (pi*(exp(Tneighbours$logdbh)/2)^2) 
        }
        if(nrow(focalSs > 0)){
            # Sneighbours to seedlings: first elements must be focal seedlings themselves (form the diagonal of the spDists matrix)
            # then, trees within the focal grid, then everything else in neighbouring grids
            Sneighbours <- rbind(focalSs[, c("x","y","logdbh")], focalTs[,c("x","y","logdbh")], neighbours)
            dist.Ss <- spDists(as.matrix(focalSs[,1:2]), as.matrix(Sneighbours[,1:2]))
            Ss.in.r <- ifelse(dist.Ss < sqrt(400/pi), 1, 0)
            diag(Ss.in.r) <- 0
            
            S.big.mat <- ifelse(
                matrix(rep(focalSs$logdbh, nrow(Sneighbours)), nrow=nrow(focalSs)) <
                    matrix(rep(Sneighbours$logdbh, each=nrow(focalSs)), nrow=nrow(focalSs))
                , 1, 0)
            S.pairs.to.count <- Ss.in.r * S.big.mat
            
            S.raw.intraS[seedlings$grid==i] <- Ss.in.r %*% (pi*(exp(Sneighbours$logdbh)/2)^2)
            S.raw.intraA[seedlings$grid==i] <- S.pairs.to.count %*% (pi*(exp(Sneighbours$logdbh)/2)^2)
        }
    }
    
    return(list(
        # intraS for adults [[1]] and seedlings [[2]]
        T.raw.intraS, S.raw.intraS,
        # intraA for adults [[3]] and seedlings [[4]]
        T.raw.intraA, S.raw.intraA
    ))
}

#intra.calc(ppoT, ppoS, sp="Prunus.polystachya")
#intraS <- intra.calc(ppoT, ppoS, sp="Prunus.polystachya")[[2]]
#intraA <- intra.calc(ppoT, ppoS, sp="Prunus.polystachya")[[4]]
#plot(log(intraA+0.001) ~ log(intraS+0.001))
#hist((log(intraA + surv.parm.unscale["surv.intraS.unscale.logoffset",]) -
#      surv.parm.unscale["surv.intraS.unscale.mu",]) / surv.parm.unscale["surv.intraS.unscale.sigma",])


## Interspecific competition (inter)
# almost same as above except competition from hetero- instead of conspecific individuals
# note: trees.hetero must include all non-focal species rbinded tgt in a 3 column matrix: [,1:2] = coordinates; [,3] = logdbh
# returns interS for adults [[1]] and seedlings [[2]]
# and interA for adults [[3]] and seedlings [[4]]
#inter.calc <- function(trees.hetero, seedlings.hetero, trees.con, seedlings.con, sp, r = sqrt(400/pi)){
    
    # convert conspecific (i.e., focal species) seedling heights to dbh, combine with adult logdbh
#    h <- seedlings.con$logheight
#    z <- c(
#        (tran.parm["tran.int",sp] + h*tran.parm["tran.h",sp]),
#        trees.con$logdbh)
    
    # convert heterospecific (i.e., competitor species) seedling heights to dbh, combine with adult logdbh
#    h.hetero <- seedlings.hetero$logheight
#    z.hetero <- c(
#        (tran.parm["tran.int","All.other.spp"] + h.hetero*tran.parm["tran.h","All.other.spp"]),
#        trees.hetero$logdbh)
    
    # rbind coordinates of conspecific (i.e., focal species) and heterospecific (i.e., competitor) trees and seedlings
#    con.ALL <- as.matrix(rbind(trees.con[,1:2], seedlings.con[,1:2]))
#    hetero.ALL <- as.matrix(rbind(trees.hetero[,1:2], seedlings.hetero[,1:2]))
    
    # calc distance between heterospecific competitors and focal species individuals, find all within 400m2 radius of focal
#    dists <- spDists(hetero.ALL, con.ALL)
#    Ts.in.r <- ifelse(dists < sqrt(400/pi), 1, 0)
    
    # for asymmetric competition, we only want pairs in which focal is the larger of the two, so need to compute another matrix
#    z.con.mat <- matrix(rep(z, each=length(z.hetero)), ncol=length(z))
#    z.hetero.mat <- matrix(rep(z.hetero, length(z)), ncol=length(z))
#    big.mat <- ifelse(z.hetero.mat > z.con.mat, 1, 0)
#    pairs.to.count <- Ts.in.r * big.mat
    # use matrix multiplication to obtain summed inter for each individual
#    raw.interS <- pi*(exp(z.hetero)/2)^2 %*% Ts.in.r
#    raw.interA <- pi*(exp(z.hetero)/2)^2 %*% pairs.to.count
    
#    return(list(
        # interS for adults [[1]] and seedlings [[2]]
#        as.vector(raw.interS)[1:nrow(trees.con)], as.vector(raw.interS)[(nrow(trees.con)+1):(nrow(trees.con)+nrow(seedlings.con))],
        # and interA for adults [[3]] and seedlings [[4]]
#        as.vector(raw.interA)[1:nrow(trees.con)], as.vector(raw.interA)[(nrow(trees.con)+1):(nrow(trees.con)+nrow(seedlings.con))]
#    ))
#}

# New inter.calc
inter.calc <- function(trees.hetero, seedlings.hetero, trees.con, seedlings.con, sp, r = sqrt(400/pi)){
    # convert seedling heights to dbh
    seedlings.con$logdbh <- tran.parm["tran.int",sp] + seedlings.con$logheight*tran.parm["tran.h",sp]
    seedlings.hetero$logdbh <- tran.parm["tran.int", "All.other.spp"] + seedlings.hetero$logheight*tran.parm["tran.h", "All.other.spp"]
    
    T.raw.interS <- rep(NA, nrow(trees.con))
    T.raw.interA <- rep(NA, nrow(trees.con))
    S.raw.interS <- rep(NA, nrow(seedlings.con))
    S.raw.interA <- rep(NA, nrow(seedlings.con))
    
    for(i in 1:length(nssf.100)){
        # need to separate focal trees from focal seedlings so that matching back is easier
        focalTs <- trees.con[trees.con$grid==i,]
        focalSs <- seedlings.con[seedlings.con$grid==i,]
        # neighbours are all heterospecific Ts and Ss from both focal grid and neighbouring grids
        neighbours <- rbind(
            trees.hetero[trees.hetero$grid %in% grid.neighbours[[i]], c("x","y","logdbh")],
            seedlings.hetero[seedlings.hetero$grid %in% grid.neighbours[[i]], c("x","y","logdbh")]
        )
        # compute the distance matrices
        if(nrow(focalTs > 0)){
            dist.Ts <- spDists(as.matrix(focalTs[,1:2]), as.matrix(neighbours[,1:2]))
            Ts.in.r <- ifelse(dist.Ts < sqrt(400/pi), 1, 0)
            
            # for asymmetric competition, we only want pairs in which focal is the smaller of the two, so need to compute another matrix
            T.big.mat <- ifelse(
                matrix(rep(focalTs$logdbh, nrow(neighbours)), nrow=nrow(focalTs)) <
                    matrix(rep(neighbours$logdbh, each=nrow(focalTs)), nrow=nrow(focalTs))
                , 1, 0)
            T.pairs.to.count <- Ts.in.r * T.big.mat
            
            # use matrix multiplication to obtain summed inter for each individual
            T.raw.interS[trees.con$grid==i] <- Ts.in.r %*% (pi*(exp(neighbours$logdbh)/2)^2)
            T.raw.interA[trees.con$grid==i] <- T.pairs.to.count %*% (pi*(exp(neighbours$logdbh)/2)^2) 
        }
        if(nrow(focalSs > 0)){
            # Sneighbours to seedlings: first elements must be focal seedlings themselves (form the diagonal of the spDists matrix)
            # then, trees within the focal grid, then everything else in neighbouring grids
            dist.Ss <- spDists(as.matrix(focalSs[,1:2]), as.matrix(neighbours[,1:2]))
            Ss.in.r <- ifelse(dist.Ss < sqrt(400/pi), 1, 0)
            
            S.big.mat <- ifelse(
                matrix(rep(focalSs$logdbh, nrow(neighbours)), nrow=nrow(focalSs)) <
                    matrix(rep(neighbours$logdbh, each=nrow(focalSs)), nrow=nrow(focalSs))
                , 1, 0)
            S.pairs.to.count <- Ss.in.r * S.big.mat
            
            S.raw.interS[seedlings.con$grid==i] <- Ss.in.r %*% (pi*(exp(neighbours$logdbh)/2)^2)
            S.raw.interA[seedlings.con$grid==i] <- S.pairs.to.count %*% (pi*(exp(neighbours$logdbh)/2)^2)
        }
    }
    
    return(list(
        # interS for adults [[1]] and seedlings [[2]]
        T.raw.interS, S.raw.interS,
        # interA for adults [[3]] and seedlings [[4]]
        T.raw.interA, S.raw.interA
    ))
}

#inter.calc(sceT, sceS, ppoT, ppoS, sp="Prunus.polystachya")
#interS <- inter.calc(sceT, sceS, ppoT, ppoS, sp="Prunus.polystachya")[[2]]
#interA <- inter.calc(sceT, sceS, ppoT, ppoS, sp="Prunus.polystachya")[[4]]
#plot(log(interS+0.01) ~ log(interA+0.01))
#hist((log(interS) -
#      surv.parm.unscale["surv.interS.unscale.mu",]) / surv.parm.unscale["surv.interS.unscale.sigma",])

# Below: Distance-explicit measures--only for adult growth models

## Distance-explicit asymmetric intraspecific competition (intraA.dist)
# almost same as above, except only for adult trees, and more spatially sensitive
#intra.dist.calc <- function(trees, r=sqrt(1600/pi)){
    # calc distance
#    dists <- spDists(trees)
#    dbh <- exp(trees$logdbh)
    # diagonal elements and trees which are > sqrt(400/pi) m from focal should not be computed (note: 1600m2 is the area of a 40x40m plot)
    # to get rid simply set to a very large value
#    dists[which(dists > r)] <- 1e10
#    diag(dists) <- 1e10
#    dists <- dists^-1
#    raw.intraS.dist <- (dists %*% dbh)
    # intraA is simply intraS divided by dbh
#    return(raw.intraS.dist/dbh)
#}

# new intra.dist.calc
intra.dist.calc <- function(trees, r=sqrt(1600/pi)){
    raw.intraS.dist <- rep(NA, nrow(trees))
    for(i in 1:length(nssf.100)){
        focal <- trees[trees$grid==i,]
        neighbours <- rbind(
            trees[trees$grid==i,],
            trees[trees$grid %in% grid.neighbours[[i]],]
        )
        # compute distance matrix only if there are both focals and neighbours in the grids
        if(nrow(focal)>0){
            dist.mat <- spDists(
                as.matrix(focal[,c("x","y")]), 
                as.matrix(neighbours[,c("x","y")])
            )
            # remove diagonals and distant competition by making very large value 
            dist.mat[which(dist.mat > r)] <- 1e10
            diag(dist.mat) <- 1e10
            # for overlapping stems, make distance = 10cm
            dist.mat[which(dist.mat==0)] <- 0.1
            dist.mat <- dist.mat^-1        
            raw.intraS.dist[which(trees$grid==i)] <- dist.mat %*% exp(neighbours$logdbh)
        } 
    }
    # intraA is simply intraS divided by dbh
    return(raw.intraS.dist/exp(trees$logdbh))
}
#intra.dist.calc(ppoT)
#hist(intraA <- intra.dist.calc(ppoT))
#hist((log(intraA + grow.T.parm.unscale["grow.intraA.unscale.logoffset",]) -
#      grow.T.parm.unscale["grow.intraA.unscale.mu",]) / grow.T.parm.unscale["grow.intraA.unscale.sigma",])

## Distance-explicit symmetric interspecific competition (interS.dist)
# almost same as above, except only for adult trees, and more spatially sensitive
#inter.dist.calc <- function(trees.hetero, trees.con, r=sqrt(1600/pi)){
    # calc distance
#    dists <- spDists(as.matrix(trees.hetero[,1:2]), as.matrix(trees.con[,1:2]))
#    dbh <- exp(trees.hetero$logdbh)
    # diagonal elements and trees which are > sqrt(400/pi) m from focal should not be computed (note: 1600m2 is the area of a 40x40m plot)
    # to get rid simply set to a very large value
#    dists[which(dists > r)] <- 1e10
#    diag(dists) <- 1e10
#    dists <- dists^-1
#    raw.interS.dist <- matrix(dbh, nrow=1) %*% dists
#    return(as.vector(raw.interS.dist))
#}

# new inter.dist.calc
inter.dist.calc <- function(trees.hetero, trees.con, r=sqrt(1600/pi)){
    raw.interS.dist <- rep(NA, nrow(trees.con))
    for(i in 1:length(nssf.100)){
        focal <- trees.con[trees.con$grid==i,]
        neighbours <- rbind(
            trees.hetero[trees.hetero$grid==i,],
            trees.hetero[trees.hetero$grid %in% grid.neighbours[[i]],]
        )
        if(nrow(focal)>0 & nrow(neighbours)>0){
            dist.mat <- spDists(
                as.matrix(focal[,c("x","y")]), 
                as.matrix(neighbours[,c("x","y")])
            )
            # remove distant competition by making very large value 
            dist.mat[which(dist.mat > r)] <- 1e10
            # for overlapping stems, make distance = 10cm
            dist.mat[which(dist.mat==0)] <- 0.1
            dist.mat <- dist.mat^-1       
            raw.interS.dist[which(trees.con$grid==i)] <- dist.mat %*% exp(neighbours$logdbh)
        } else if(nrow(focal)>0 & nrow(neighbours)==0) {
            raw.interS.dist[which(trees.con$grid==i)] <- 0
        } 
    }
    return(raw.interS.dist)
}
#inter.dist.calc(sceT, ppoT)
#hist(interS <- inter.dist.calc(sceT, ppoT))
#hist((log(interS) - grow.T.parm.unscale["grow.interS.unscale.mu",]) / grow.T.parm.unscale["grow.interS.unscale.sigma",])

###############################################
# FUNCTION FOR PREVENTING RECRUITMENT OVERLAP #
###############################################

# This function takes in locations (x, y, grid) of new recruits (that haven't entered landscape yet)
# compares them against locations of existing trees/seedlings
# and removes recruits which are recruiting in the physical space of existing trees/seedlings

rm.overlap <- function(recruits, trees, seedlings){
    # define a selection vector
    sel <- c()
    for(i in 1:nrow(recruits)){
        rgrid <- recruits[, "grid"][i]
        tree.neighbours <- trees[trees$grid %in% c(rgrid, grid.neighbours[[rgrid]]),]
        seedling.neighbours <- seedlings[seedlings$grid %in% c(rgrid, grid.neighbours[[rgrid]]),]
        tree.dists <- spDists(recruits[i, c("x","y"), drop = FALSE], as.matrix(tree.neighbours[,c("x","y")]))
        seedling.dists <- spDists(recruits[i, c("x","y"), drop = FALSE], as.matrix(seedling.neighbours[,c("x","y")]))
        if(
            # seedling is only allowed to recruit if it is not found within 10cm of existing seedlings
            min(seedling.dists) > 0.1 &
            # and not found within the DBH radii of existing trees
            sum(tree.dists < exp(tree.neighbours$logdbh)/2) == 0
        ){sel <- c(sel, i)}
    }
    return(recruits[sel, , drop = FALSE])
}
#rm.overlap(sceS[1:10,c("x","y","grid")], ppoT, ppoS)
