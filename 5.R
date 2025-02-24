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
  r <- sum(Develop==2,na.rm=T)/sum(Develop==1,na.rm=T)
  Develop.raw <- Develop
  Develop[Develop==2] <- 1
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
  # Parallelize the outer loops using foreach
  x.dist <- foreach(i = 1:nrow(Water), .combine = rbind, .packages = c('foreach')) %dopar% {
    row_result <- base::matrix(NA, nrow = 1, ncol = ncol(Develop))
    for (j in 1:ncol(Water)) {
      if ((!is.na(Develop[i,j]))&(!is.na(Water[i,j]))) {
        if ((Develop[i,j]==0)&(Water[i,j]==1)&(DEM[i,j]==1)) {
          row_result[j] <- value_to_clusters(c(i,j),centers)        }
      }
    }
    row_result
  }
  if (steps<=10) {
    X = k$k0+k$k1*x.develop+k$k2*x.dist
  }
  if ((10<steps)&(steps<=20)) {
    X = k$k0+k$k1*x.develop+k$k2*x.dist
  }
  if ((20<steps)&(steps<=28)) {
    X = k$k0+k$k1*x.develop+k$k2*x.dist
  }
  if (steps>28) {
    X = k$k0+k$k1*x.develop+k$k2*x.dist
  }
  p.develop <- Water/(1+exp(-X))
  # Roulette
  x.expan <- base::matrix(NA, nrow = base::nrow(Develop)+n, ncol = base::ncol(Develop)+n)
  x.expan[(n/2+1):(n/2+nrow(Develop)),(n/2+1):(n/2+ncol(Develop))] <- Develop.raw
  x.new.develop <- foreach(i = 1:nrow(p.develop), .combine = rbind, .packages = c('foreach')) %dopar% {
    row_result <- numeric(ncol(p.develop))
    for (j in 1:ncol(p.develop)) {
      if (!is.na(p.develop[i, j])) {
        if ((Water[i,j]==1)&&(p.develop[i,j] >= runif(1))) {
          x.window <- x.expan[i:(i+n), j:(j+n)]
          x.window <- sum(x.window==2,na.rm=T)/8
          if (r<0.5) {
            if (((exp(-2.3+2.3*x.window)-exp(-2.3))/(1-exp(-2.3))) >=7/8) {
              row_result[j] <- 1
            }
            if ((((exp(-2.3+2.3*x.window)-exp(-2.3))/(1-exp(-2.3))) <7/8) & 
                (((exp(-2.3+2.3*x.window)-exp(-2.3))/(1-exp(-2.3))) >4/8)) {
              row_result[j] <- 1
            }
            if (((exp(-2.3+2.3*x.window)-exp(-2.3))/(1-exp(-2.3))) ==4/8)  {
              row_result[j] <- 2
            }
            if ((((exp(-2.3+2.3*x.window)-exp(-2.3))/(1-exp(-2.3))) <4/8) & 
                (((exp(-2.3+2.3*x.window)-exp(-2.3))/(1-exp(-2.3))) >=2/8)) {
              row_result[j] <- 2
            }
            if (((exp(-2.3+2.3*x.window)-exp(-2.3))/(1-exp(-2.3))) < 2/8) {
              row_result[j] <- 2
            }
          }
          if ((r>=0.5)&(r<=0.6)) {
            if (x.window>=7/8) {
              row_result[j] <- 1
            }
            if ((x.window<7/8)&(x.window>4/8)) {
              row_result[j] <- 2
            }
            # if (x.window==4/8) {
            #   row_result[j] <- sample(1:2,1)
            # }
            if (x.window==4/8) {
              row_result[j] <- 1
            }
            if ((x.window<4/8)&(x.window>2/8)) {
              row_result[j] <- 1
            }
            if (x.window<=2/8) {
              row_result[j] <- 2
            }
          }
          if (r>0.6) {
            if (((1-exp(-2.44*x.window))/(1-exp(-2.44))) >= 7/8) {
              row_result[j] <- 2
            }
            if (((1-exp(-2.44*x.window))/(1-exp(-2.44))) < 7/8 & 
                ((1-exp(-2.44*x.window))/(1-exp(-2.44))) > 4/8) {
              row_result[j] <- 2
            }
            if ((1-exp(-2.44*x.window))/(1-exp(-2.44)) ==4/8)  {
              row_result[j] <- 1
            }
            if (((1-exp(-2.44*x.window))/(1-exp(-2.44))) < 4/8 & 
                ((1-exp(-2.44*x.window))/(1-exp(-2.44))) >= 2/8) {
              row_result[j] <- 1
            }
            if (((1-exp(-2.44*x.window))/(1-exp(-2.44))) < 2/8) {
              row_result[j] <- 1
            }
          }
        }
      }
    }
    row_result
  }
  df <- as.data.frame(which(x.new.develop >= 1, arr.ind = TRUE))
  aim_1 <- round(nrow(df) * aim_1)
  df$zones <- zones[cbind(df$row,df$col)]
  df <- df[df$zones!=7,]
  df$prob <- p.develop[cbind(df$row, df$col)]
  df_1 <- df[order(df$prob, decreasing = TRUE)[c(1:aim_1)], ]
  df_2 <- df[order(df$prob, decreasing = TRUE)[c(1:aim_2)], ]
  x.new.develop_1 <- base::matrix(0, nrow = base::nrow(p.develop), ncol = base::ncol(p.develop))
  x.new.develop_1[cbind(df_1$row, df_1$col)] <- x.new.develop[cbind(df_1$row, df_1$col)]
  x.new.develop_2 <- base::matrix(0, nrow = base::nrow(p.develop), ncol = base::ncol(p.develop))
  x.new.develop_2[cbind(df_2$row, df_2$col)] <- x.new.develop[cbind(df_2$row, df_2$col)]
  return(list('self' = x.new.develop_1, 'plan' = x.new.develop_2))
}

library(doParallel)
library(foreach)
library(terra)
library(sf)
library(data.table)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
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
GHSL.NRES <- terra::rast('settlement-GHSL/NRES/GHS_BUILT_S_NRES_E1975_GLOBE_R2023A_54009_1000_V1_0_R6_C29.tif')
for (i in 1:11) {
  terra::add(GHSL.NRES) <- terra::rast(paste0('settlement-GHSL/NRES/GHS_BUILT_S_NRES_E',1975+5*i,'_GLOBE_R2023A_54009_1000_V1_0_R6_C29.tif'))
}
names(GHSL.NRES) <- seq(1975,2030,5)
GHSL.NRES <- terra::project(GHSL.NRES,Nanjing.shp)
GHSL.NRES <- terra::crop(GHSL.NRES,Nanjing.shp, mask = TRUE)
GHSL.RES <- GHSL.Total - GHSL.NRES
names(GHSL.RES) <- seq(1975,2030,5)
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
Water <- 1-base::as.numeric(Water>0.5)
GHSL.Total[['1975']] <- GHSL.Total[['1975']]*Water
GHSL.Total[['2020']] <- GHSL.Total[['2020']]*Water
Water <- terra::as.matrix(Water, wide=TRUE)
# transfer to bool
GHSL.Total.bool <- base::as.numeric(GHSL.Total>0.5e5)
GHSL.NRES.bool <- base::as.numeric(GHSL.NRES>0.15e5&GHSL.Total>0.5e5)
GHSL.Total.bool <- GHSL.Total.bool+GHSL.NRES.bool # 1 residential, 2 non-residential
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
InitialPattern.center[InitialPattern.center<5e6] <- NA
centers <- terra::as.data.frame(InitialPattern.center,cells=T)
centers$cell <- terra::rowColFromCell(InitialPattern.center, centers$cell)
centers$last <- centers$last/max(centers$last)
#
# Set up parallel backend
num_cores <- detectCores() - 1  # Use one less than the number of available cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)
df.tbl <- c()
k <- data.frame(n  = 2,
                k0 = 0,
                k1 = 0.75,
                k2 = 0.5,
                k3 = 0.25)
Develop <- list(terra::as.matrix(GHSL.Total.bool[[1]], wide=TRUE))
Steps <- 41
aim.plan <- rep(75,Steps)
aim.self <- 0.02
Gap <- rep(0,Steps)
lever.p <- rep(1,Steps)
lever.s <- rep(1,Steps)
for (i in 1:Steps) {
  if (i <= 10) {
    lower <- 0.7
    upper <- 1.5
  }
  if (10 < i & i <= 20) {
    lower <- 0.7
    upper <- 2.5
  }
  if (20 < i & i <= 30) {
    lower <- 0.5
    upper <- 3.5
  }
  if (30 < i & i <= 36) {
    lower <- 0.7
    upper <- 2.5
  }
  if (36 < i) {
    lower <- 0.9
    upper <- 1.5
  }
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
  New.developed <- CA.update(Develop[[i]],
                             Water,
                             DEM,
                             centers,
                             k,
                             i,
                             aim.self*lever.p[i],
                             aim.plan[i]*lever.s[i])
  Develop[[i+1]] <- Develop[[i]] + New.developed[['self']]
  Gap[i] <- sum(New.developed[['self']]>0,na.rm=T)-sum(New.developed[['plan']]>0,na.rm=T)
  # plan
  aim.plan[i+1] = aim.plan[i]*lever.s[i]-5
  if (aim.plan[i+1]<5) {
    aim.plan[i+1] = 5
  }
  if (i%%10==0) {
    aim.plan[i+1] = aim.plan[i-9] * mean(lever.s[(i-9):i])+25
  }
}

df <- table(as.data.frame(na.omit(cbind(c(t(Develop[[i+1]])),terra::values(GHSL.Total.bool[['2020']])))))
kappa.results <- vcd::Kappa(df)
df.tbl <- rbind(df.tbl,c(sum(Develop[[i+1]]>0,na.rm=T) - sum(terra::values(GHSL.Total.bool[['2020']])>0,na.rm=T),
                         kappa.results$Unweighted[1],raster::Moran(raster::raster(Develop[[i+1]]))))
df
df.tbl
terra::plot(c(terra::rast(Develop[[i+1]],crs = crs(GHSL.Total.bool[['2020']]),extent=ext(GHSL.Total.bool[['2020']])),
              GHSL.Total.bool[['2020']]))
terra::writeRaster(terra::rast(Develop[[i+1]],crs = crs(GHSL.Total.bool[['2020']]),extent=ext(GHSL.Total.bool[['2020']])),
                   "final_final_final.tif",overwrite=T)
terra::writeRaster(GHSL.Total.bool[['2020']], "2020_compete.tif",overwrite=T)
terra::writeRaster(GHSL.Total.bool[['1975']], "1975_compete.tif",overwrite=T)
sum(Develop[[i+1]]==2,na.rm = T)/sum(Develop[[i+1]]==1,na.rm = T)
sum(terra::values(GHSL.Total.bool[['2020']])==2,na.rm=T)/sum(terra::values(GHSL.Total.bool[['2020']])==1,na.rm=T)
sum(terra::values(GHSL.Total.bool[['1975']])==2,na.rm=T)/sum(terra::values(GHSL.Total.bool[['1975']])==1,na.rm=T)
stopCluster(cl)

raster <- terra::rast("final_final_final.tif")  
raster2 <-terra::rast("2020_compete.tif")  
df <- table(as.data.frame(na.omit(cbind(terra::values(raster),terra::values(raster2)))))
kappa.results <- vcd::Kappa(df)
c(sum(terra::values(raster),na.rm=T) - sum(terra::values(raster2),na.rm=T),
  kappa.results$Unweighted[1],raster::Moran(raster::raster(raster)))

