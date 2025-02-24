CA.update <- function(Develop, Water, DEM, centers, k, steps, aim_1, aim_2) {
  value_to_clusters <- function(x,centers) {
    values = seq(nrow(centers))*0
    for (i in seq(nrow(centers))) {
      loc <- dist(rbind(x,centers$cell[i,]))
      values[i] <-centers$last[i]-centers$last[i]*exp(-6+0.6*loc)/(1+exp(-6+0.6*loc))
    }
    return(sum(values))
  }
  n <- k$n
  # Prepare expanded Develop matrix
  x.expan <- base::matrix(NA, nrow = base::nrow(Develop)+n, ncol = base::ncol(Develop)+n)
  x.expan[(n/2+1):(n/2+nrow(Develop)),(n/2+1):(n/2+ncol(Develop))] <- Develop
  # Parallelize the outer loops using foreach
  x.develop <- foreach(i = 1:nrow(Develop), .combine = rbind, .packages = c('foreach')) %dopar% {
    row_result <- base::matrix(NA, nrow = 1, ncol = ncol(Develop))
    for (j in 1:ncol(Develop)) {
      if ((!is.na(Develop[i,j]))&(!is.na(Water[i,j]))) {
        if ((Develop[i,j]==0)&(Water[i,j]==1)&(DEM[i,j]==1)) {
          x.window <- x.expan[i:(i+n), j:(j+n)]
          row_result[j] <- sum(x.window[-(n/2*(n+2)+1)],na.rm=TRUE)/((n+1)^2-1)
        }
      }
    }
    row_result
  }
  # dist to zones_centers
  # Parallelize the outer loops using foreach
  x.dist <- foreach(i = 1:nrow(Water), .combine = rbind, .packages = c('foreach')) %dopar% {
    row_result <- base::matrix(NA, nrow = 1, ncol = ncol(Develop))
    for (j in 1:ncol(Water)) {
      if ((!is.na(Develop[i,j]))&(!is.na(Water[i,j]))) {
        if ((Develop[i,j]==0)&(Water[i,j]==1)&(DEM[i,j]==1)) {
          row_result[j] <- value_to_clusters(c(i,j),centers)
        }
      }
    }
    row_result
  }
  X = k$k0+k$k1*x.develop+k$k2*x.dist
  p.develop <- Water/(1+exp(-X))
  # Roulette
  x.new.develop <- foreach(i = 1:nrow(p.develop), .combine = rbind, .packages = c('foreach')) %dopar% {
    row_result <- numeric(ncol(p.develop))
    for (j in 1:ncol(p.develop)) {
      if (!is.na(p.develop[i, j])) {
        if ((Water[i,j]==1)&&(p.develop[i,j] >= runif(1))) {
          row_result[j] <- 1
        }
      }
    }
    row_result
  }
  df <- as.data.frame(which(x.new.develop == 1, arr.ind = TRUE))
  df$prob <- p.develop[cbind(df$row, df$col)]
  df_1 <- df[order(df$prob, decreasing = TRUE)[-c(1:round(nrow(df) * aim_1))], ]
  df_2 <- df[order(df$prob, decreasing = TRUE)[-c(1:aim_2)], ]
  x.new.develop_1 <- x.new.develop
  x.new.develop_1[cbind(df_1$row, df_1$col)] <- 0
  x.new.develop_2 <- x.new.develop
  x.new.develop_2[cbind(df_2$row, df_2$col)] <- 0
  return(list('self' = x.new.develop_1, 'plan' = x.new.develop_2))
}

library(doParallel)
library(foreach)
library(terra)
library(sf)
library(data.table)
library(vcd)	
# Set up parallel backend
num_cores <- detectCores() - 1  # Use one less than the number of available cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)
setwd("/mainfs/scratch/wz6g23/CA")
# read Nanjing shape file
Nanjing.shp <- terra::vect('Nanjing/Nanjing.shp')
# import GHSL data
# Total
GHSL.Total <- terra::rast('settlement-GHSL/Total/GHS_BUILT_S_E1975_GLOBE_R2023A_54009_1000_V1_0_R6_C29.tif')
for (i in 1:11) {
  terra::add(GHSL.Total) <- terra::rast(paste0('settlement-GHSL/Total/GHS_BUILT_S_E',1975+5*i,'_GLOBE_R2023A_54009_1000_V1_0_R6_C29.tif'))
}
names(GHSL.Total) <- seq(1975,2030,5)
GHSL.Total <- terra::project(GHSL.Total,Nanjing.shp)
GHSL.Total <- terra::crop(GHSL.Total,Nanjing.shp, mask = TRUE)
# import DEM
DEM.1 <- terra::rast('DEM/ASTGTMV003_N31E118_dem.tif')
DEM.2 <- terra::rast('DEM/ASTGTMV003_N31E119_dem.tif')
DEM.3 <- terra::rast('DEM/ASTGTMV003_N32E118_dem.tif')
DEM.4 <- terra::rast('DEM/ASTGTMV003_N32E119_dem.tif')
DEM <- terra::merge(DEM.1,DEM.2,DEM.3,DEM.4)
DEM <- terra::project(DEM,GHSL.Total,'average')
DEM <- terra::crop(DEM, Nanjing.shp, mask = TRUE)
DEM <- 1-base::as.numeric(DEM>80)
GHSL.Total[['1975']] <- GHSL.Total[['1975']]*DEM
GHSL.Total[['2020']] <- GHSL.Total[['2020']]*DEM
DEM <- terra::as.matrix(DEM, wide=TRUE)
# import Water
Water.1 <- terra::rast('Water/ASTWBDV001_N31E118_att.tif')
Water.2 <- terra::rast('Water/ASTWBDV001_N31E119_att.tif')
Water.3 <- terra::rast('Water/ASTWBDV001_N32E118_att.tif')
Water.4 <- terra::rast('Water/ASTWBDV001_N32E119_att.tif')
Water <- terra::merge(Water.1,Water.2,Water.3,Water.4)
Water <- base::as.numeric(Water>0)
Water <- terra::project(Water,GHSL.Total,'average')
Water <- terra::crop(Water,Nanjing.shp, mask = TRUE)
Water <- 1-base::as.numeric(Water>0.5)# more than half area was water was defined as a water cell which cannot be developed
GHSL.Total[['1975']] <- GHSL.Total[['1975']]*Water
GHSL.Total[['2020']] <- GHSL.Total[['2020']]*Water
Water <- terra::as.matrix(Water, wide=TRUE)
# transfer to bool
GHSL.Total.bool <- base::as.numeric(GHSL.Total>0.5e5)
# identify centers
InitialPattern <- GHSL.Total.bool['1975']
InitialPattern.polygons <- terra::as.polygons(InitialPattern,values=T,na.rm=TRUE)
InitialPattern.settlement <- InitialPattern.polygons[terra::values(InitialPattern.polygons)==1]
InitialPattern.settlement <- sf::st_as_sf(InitialPattern.settlement)
InitialPattern.settlement <- sf::st_cast(InitialPattern.settlement, "POLYGON")
InitialPattern.settlement <- terra::vect(InitialPattern.settlement)
InitialPattern.settlement.center <- terra::centroids(InitialPattern.settlement,inside=T)
terra::values(InitialPattern.settlement.center) <- terra::expanse(InitialPattern.settlement)
InitialPattern.center <- terra::rasterize(InitialPattern.settlement.center,InitialPattern,field='value')
InitialPattern.center[InitialPattern.center<1e6] <- NA
centers <- terra::as.data.frame(InitialPattern.center,cells=T)
centers$cell <- terra::rowColFromCell(InitialPattern.center, centers$cell)
centers$last <- centers$last/max(centers$last)
#calculation of development probability
D.Total <- c()
df.tbl <- c()
kappa.tbl <- c()
Moran.df <- c()
k <- data.frame(n  = 2,
                k0 = 0,
                k1 = 1,
                k2 = 1)
for (lower in seq(0,1,by=0.1)) {
  for (upper in seq(1,10,by=0.1)) {
    Develop_1 <- list(terra::as.matrix(GHSL.Total.bool[[1]], wide=TRUE))
    Develop_2 <- list(terra::as.matrix(GHSL.Total.bool[[1]], wide=TRUE))
    #
    Steps <- 41
    aim.plan <- rep(75,Steps)
    aim.self <- 0.02
    Gap <- rep(0,Steps)
    lever.p <- rep(1,Steps)
    lever.s <- rep(1,Steps)
    for (i in 1:Steps) {
      if (((i>2)&(i<=12))|((i>14)&(i<=22))|((i>24)&(i<=32))|(i>34)) {
        if ((Gap[i-1]<0)&(abs(Gap[i-2])<abs(Gap[i-1]))&(Gap[i-2]!=0)) {lever.p[i] <- abs(Gap[i-1]/Gap[i-2])}
        if ((Gap[i-1]>0)&(abs(Gap[i-2])<abs(Gap[i-1]))) {lever.p[i] <- abs(Gap[i-2]/Gap[i-1])}
        if ((Gap[i-2]==Gap[i-1])&(Gap[i-2]>0)) {lever.p[i] <- 1/Gap[i-1]}
        if ((Gap[i-2]==Gap[i-1])&(Gap[i-2]<0)) {lever.p[i] <- 1-1/Gap[i-1]}
        if ((Gap[i-2]==0)&(Gap[i-1]==0)) {lever.p[i] <- 1}
      }
      if (lever.p[i]>upper) {lever.p[i]=upper}
      if (lever.p[i]<lower) {lever.p[i]=lower}
      if (((i>4)&(i<=12))|((i>16)&(i<=22))|((i>26)&(i<=32))|(i>36)) {
        if ((lever.p[i-2]>=1)&(lever.p[i-1]>=1)) {lever.s[i] <- lever.p[i-2]/(lever.p[i-1]+lever.p[i-2])}
        if ((lever.p[i-2]>1)&(lever.p[i-1]==lever.p[i-2])) {lever.s[i] <- 1/lever.p[i-1]}
        if ((lever.p[i-2]<=1)&(lever.p[i-1]<=1)) {lever.s[i] <- 1+(lever.p[i-2]/(lever.p[i-1]+lever.p[i-2]))}
        if ((lever.p[i-1]<1)&(lever.p[i-1]==lever.p[i-2])) {lever.s[i] <- 1+lever.p[i-1]}
      }
      New.developed <- CA.update(Develop_1[[i]],
                                 Water,
                                 DEM,
                                 centers,
                                 k,
                                 i,
                                 aim.self*lever.p[i],
                                 aim.plan[i]*lever.s[i])
      Develop_1[[i+1]] <- Develop_1[[i]] + New.developed[['self']]
      Gap[i] <- sum(New.developed[['self']])-sum(New.developed[['plan']])
      #plan
      aim.plan[i+1] = aim.plan[i]*lever.s[i]-5
      if (aim.plan[i+1]<5) {
        aim.plan[i+1] = 5
      }
      if (i==12) {
        aim.plan[i+1] = 125
      }
      if (i==22) {
        aim.plan[i+1] = 225
      }
      if (i==32) {
        aim.plan[i+1] = 175
      }
    }
    df <- table(as.data.frame(na.omit(cbind(c(t(Develop_1[[i+1]])),terra::values(GHSL.Total.bool[['2020']])))))
    kappa.results <- vcd::Kappa(df)
    df.tbl <- rbind(df.tbl,c(4.1,lower,upper,sum(Develop_1[[i+1]],na.rm=T) - sum(terra::values(GHSL.Total.bool[['2020']]),na.rm=T), kappa.results$Unweighted[1]))
    #
    aim.plan <- rep(75,Steps)
    aim.self <- 0.02
    Gap <- rep(0,Steps)
    lever.p <- rep(1,Steps)
    lever.s <- rep(1,Steps)
    for (i in 1:Steps) {
      if (i>2) {
        if ((Gap[i-1]<0)&(abs(Gap[i-2])<abs(Gap[i-1]))&(Gap[i-2]!=0)) {lever.p[i] <- abs(Gap[i-1]/Gap[i-2])}
        if ((Gap[i-1]>0)&(abs(Gap[i-2])<abs(Gap[i-1]))) {lever.p[i] <- abs(Gap[i-2]/Gap[i-1])}
        if ((Gap[i-2]==Gap[i-1])&(Gap[i-2]>0)) {lever.p[i] <- 1/Gap[i-1]}
        if ((Gap[i-2]==Gap[i-1])&(Gap[i-2]<0)) {lever.p[i] <- 1-1/Gap[i-1]}
        if ((Gap[i-2]==0)&(Gap[i-1]==0)) {lever.p[i] <- 1}
      } 
      if (lever.p[i]>upper) {lever.p[i]=upper}
      if (lever.p[i]<lower) {lever.p[i]=lower}
      if (i>4) {
        if ((lever.p[i-2]>=1)&(lever.p[i-1]>=1)) {lever.s[i] <- lever.p[i-2]/(lever.p[i-1]+lever.p[i-2])}
        if ((lever.p[i-2]>1)&(lever.p[i-1]==lever.p[i-2])) {lever.s[i] <- 1/lever.p[i-1]}
        if ((lever.p[i-2]<=1)&(lever.p[i-1]<=1)) {lever.s[i] <- 1+(lever.p[i-2]/(lever.p[i-1]+lever.p[i-2]))}
        if ((lever.p[i-1]<1)&(lever.p[i-1]==lever.p[i-2])) {lever.s[i] <- 1+lever.p[i-1]}
      }
      New.developed <- CA.update(Develop_2[[i]],
                                 Water,
                                 DEM,
                                 centers,
                                 k,
                                 i,
                                 aim.self*lever.p[i],
                                 aim.plan[i]*lever.s[i])
      Develop_2[[i+1]] <- Develop_2[[i]] + New.developed[['self']]
      Gap[i] <- sum(New.developed[['self']])-sum(New.developed[['plan']])
      # plan
      aim.plan[i+1] = aim.plan[i]*lever.s[i]-5
      if (aim.plan[i+1]<5) {
        aim.plan[i+1] = 5
      }
      if (i%%10==0) {
        aim.plan[i+1] = aim.plan[i-9] * mean(lever.s[(i-9):i])+25
      }
    }
    df <- table(as.data.frame(na.omit(cbind(c(t(Develop_2[[i+1]])),terra::values(GHSL.Total.bool[['2020']])))))
    kappa.results <- vcd::Kappa(df)
    df.tbl <- rbind(df.tbl,c(4.3,lower,upper,sum(Develop_2[[i+1]],na.rm=T) - sum(terra::values(GHSL.Total.bool[['2020']]),na.rm=T), kappa.results$Unweighted[1]))
  }
}
stopCluster(cl)
write.csv(df.tbl,'result1.csv')