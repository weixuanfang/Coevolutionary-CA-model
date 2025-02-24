
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
#terra::writeRaster(InitialPattern.center,'center.tif')
centers <- terra::as.data.frame(InitialPattern.center,cells=T)
centers$cell <- terra::rowColFromCell(InitialPattern.center, centers$cell)
centers$last <- centers$last/max(centers$last)
# calculation of development probability
source("CA.update.self.R")
source("CA.update.plan.R")
k <- data.frame(n  = 2,#n
                k0 = 0,#k0
                k1 = 1,#DEVELOP
                k2 = 1)#Dist to center
Develop <- list(terra::as.matrix(GHSL.Total.bool[[1]], wide=TRUE))
#
New.developed.self <- list(Water*0)
New.developed.plan <- list(Water*0)
x.self <- GHSL.Total.bool[[1]]
x.plan <- GHSL.Total.bool[[1]]
x.self.new <- GHSL.Total.bool[[1]]
x.plan.new <- GHSL.Total.bool[[1]]
Steps <- 41
aim.plan <- rep(50,Steps)
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
  if (lever.p[i]>1.4) {lever.p[i]=1.4}
  if (lever.p[i]<0.6) {lever.p[i]=0.6}
  if (i>4) {
    if ((lever.p[i-2]>=1)&(lever.p[i-1]>=1)) {lever.s[i] <- lever.p[i-2]/(lever.p[i-1]+lever.p[i-2])}
    if ((lever.p[i-2]>1)&(lever.p[i-1]==lever.p[i-2])) {lever.s[i] <- 1/lever.p[i-1]}
    if ((lever.p[i-2]<=1)&(lever.p[i-1]<=1)) {lever.s[i] <- 1+(lever.p[i-2]/(lever.p[i-1]+lever.p[i-2]))}
    if ((lever.p[i-1]<1)&(lever.p[i-1]==lever.p[i-2])) {lever.s[i] <- 1+lever.p[i-1]}
  }
  New.developed.self[[i]] <- CA.update.self(Develop[[i]],
                                            Water,
                                            DEM,
                                            centers, 
                                            k,
                                            aim.self*lever.p[i])
  terra::values(x.self.new) <- Develop[[i]]+2*New.developed.self[[i]]
  x.self <- c(x.self,x.self.new)
  New.developed.plan[[i]] <- CA.update.plan(Develop[[i]],
                                            Water,
                                            DEM,
                                            centers, 
                                            k,
                                            aim.plan[i]*lever.s[i])
  terra::values(x.plan.new) <- Develop[[i]]+2*New.developed.plan[[i]]
  x.plan <- c(x.plan,x.plan.new)
  Gap[i] <- sum(New.developed.self[[i]],na.rm=T)-sum(New.developed.plan[[i]],na.rm=T)
  Develop[[i+1]] <-Develop[[i]]+New.developed.self[[i]]
  aim.plan[i+1] = aim.plan[i]*lever.s[i]
  if (i%%10==0) {
    aim.plan[i+1] = aim.plan[i-9] * mean(lever.s[(i-9):i])
  }
  # if (i==12) {
  #   aim.plan[i+1] = aim.plan[i-11] * mean(lever.s[(i-11):i])
  # }
  # if (i==22) {
  #   aim.plan[i+1] = aim.plan[i-9] * mean(lever.s[(i-9):i])
  # }
  # if (i==32) {
  #   aim.plan[i+1] = aim.plan[i-9] * mean(lever.s[(i-9):i])
  # }
}
D.Total <- unlist(lapply(Develop,FUN = sum, na.rm=T))
N.Total <- unlist(lapply(New.developed.self,FUN = sum, na.rm=T))
## compare
Ref <- GHSL.Total.bool[["2020"]]
terra::values(x.self.new) <- Develop[[Steps+1]]
Ref.num <- sum(terra::values(Ref),na.rm=T)
model.num <- sum(terra::values(x.self.new),na.rm=T)
# kappa
df <- table(as.data.frame(na.omit(cbind(terra::values(Ref),terra::values(x.self.new)))))
kappa.results <- vcd::Kappa(df)
# Morans'I
raster::Moran(raster::raster(Ref))
raster::Moran(raster::raster(x.self.new))
#
D.Total
D.Total[2:42]-D.Total[1:41]
df
kappa.results
aim.plan

