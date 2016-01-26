#OOOOOOOOOOOOOOOOOOOOOOOOOOOO
##ANALISIS DE DIVERSIDAD ALFA 
#OOOOOOOOOOOOOOOOOOOOOOOOOOOO

library(vegan)
library(fossil)
library(maptools)
library(data.table)
library(sp)
library(raster)
library(rgdal)
library(plyr)

#Leer tabla con datos de especies, latitud y longitud
data.registros <- read.csv("Ruta Tabla con registros")
#data.registros$X <- NULL
#data.registros$optional <- NULL

coordinates(data.registros) <- c("lon", "lat")
registros <- SpatialPoints(data.registros) #pasando de SpatialPointsDataFrame a SpatialPoints
proj4string(registros) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
writeOGR(registros,"Ruta de salida")

######################
##UNIDADES DE ANALISIS

##subir shape de unidades de analisis (shp.ua)
shp.ua <- readOGR("Ruta carpeta shape","nombre del shape") #Este shapefile debe tener un campo llamado "ID" con un identificador unico para cada poligono
shp.ua <- spTransform(shp.ua, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#Tamano de celda (tc)
tc <- readGDAL("Raster con Tamaño de celda")
tc <- raster(tc)
crs(tc) <-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#Rasterizar shp.ua con tc
raster.ua <- rasterize(shp.ua,tc,field="ID")
crs(raster.ua) <-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

####################################################
##NUMERO DE REGISTROS POR UNIDAD DE ANALISIS (nr.ua)

nr <- extract(raster.ua,data.registros)
nr.table <- as.data.frame(table(nr))
names(nr.table)[1]<-paste("ID")
nr.ua <- merge(shp.ua,nr.table, by="ID")
writeSpatialShape(nr.ua, "Ruta de salida")

############################################################################
##NUMERO DE ESPECIES (A PARTIR DE REGISTROS) POR UNIDAD DE ANALISIS (nsp.ua)

nsp <- data.frame(registros,data.registros$'nombre del campo que tiene los nombres de las especies',nr)
nsp <- na.omit(nsp)
nsp <- nsp[,3:4]
names(nsp) <- c('especies','nr')

nsp.table <- with(nsp,tapply(especies,nr,function(x) length(unique(x))))
nsp.data.frame <- data.frame(nsp.table)

setDT(nsp.data.frame, keep.rownames = TRUE)[]
names(nsp.data.frame)[1]<-paste("ID")
nsp.ua <- merge(shp.ua,nsp.data.frame, by="ID")
writeSpatialShape(nsp.ua, "Ruta de salida")

###################################################################################
##MATRIZ DE UNIDADES DE ANALISIS VS ESPECIES (A PARTIR DE REGISTROS) (matrix.sp_ua)

tabla.sp_ua <- table(as.data.frame(nsp))
df.sp_ua <- data.frame(tabla.sp_ua)
matrix.sp_ua <- with(df.sp_ua,tapply(Freq,list(nr,especies),max))
df.spVsua <- as.data.frame(matrix.sp_ua)
setDT(df.spVsua, keep.rownames = TRUE)[]
names(df.spVsua)[1]<-paste("ID")
write.csv(matrix.sp_ua, "Ruta salida")
#####################################################################################################################################
##ESTIMADORES DE RIQUEZA NO PARAMETRICOS

#Opcion 1: continuar con la matriz creada
SpeciesMatrix <- matrix.sp_ua

#Opcion 2: cargar la matriz desde un archivo
SpeciesMatrix <- read.csv("Ruta de la matriz", row.names=1, head=T)

#######
#Chao 1

chaoRichness <- apply(SpeciesMatrix,1,function(x){chao1(x,taxa.row=F)})
chaoRichness1 <- as.data.frame(chaoRichness)
setDT(chaoRichness1, keep.rownames = TRUE)[]
names(chaoRichness1)[1]<-paste("ID")
ChaoRichness <- merge(shp.ua,chaoRichness1)
writeSpatialShape(ChaoRichness, "Ruta salida")

##########
#Bootstrap

bootRichness = apply(SpeciesMatrix,1,function(i){bootstrap(i,taxa.row=F, samples=nrow(SpeciesMatrix))})
bootRichness1 <- as.data.frame(bootRichness)
setDT(bootRichness1, keep.rownames = TRUE)[]
names(bootRichness1)[1]<-paste("ID")
BootRichness <- merge(shp.ua,bootRichness1)
writeSpatialShape(BootRichness, "Ruta salida")

############
#Jackknife 2

jackRichness = apply(SpeciesMatrix,1,function(x){jack2(x,taxa.row=F,abund=F)})
jackRichness1 <- as.data.frame(jackRichness)
setDT(jackRichness1, keep.rownames = TRUE)[]
names(jackRichness1)[1]<-paste("ID")
JackRichness <- merge(shp.ua,jackRichness1)
writeSpatialShape(JackRichness, "Ruta salida")

####
#ACE

aceRichness <- apply(SpeciesMatrix,1,function(x){ACE(x,taxa.row=F)})
aceRichness1 <- as.data.frame(aceRichness)
setDT(aceRichness1, keep.rownames = TRUE)[]
names(aceRichness1)[1]<-paste("ID")
ACERichness <- merge(shp.ua,aceRichness1)
writeSpatialShape(ACERichness, "Ruta salida")

####
#ICE

iceRichness <- apply(SpeciesMatrix,1,function(x){ICE(x,taxa.row=F)})
iceRichness1 <- as.data.frame(iceRichness)
setDT(iceRichness1, keep.rownames = TRUE)[]
names(iceRichness1)[1]<-paste("ID")
ICERichness <- merge(shp.ua,iceRichness1)
writeSpatialShape(ICERichness, "Ruta salida")

############
#Rarefaction 

#Elegir tres diferentes valores de numeros de registros para hacer
#comparaciones

# excluir celdas con menos de 500 registros
index1 = which(rowSums(SpeciesMatrix)>500) #Cambiar el numero de registros elegido
rarMatrix1 = SpeciesMatrix[index1,]
rarRichness1 = apply(rarMatrix1,1,function(x){rarefy(x,500,se=FALSE,MARGIN=1)})
rarRichness1.table <- as.data.frame(rarRichness1)
setDT(rarRichness1.table, keep.rownames = TRUE)[]
names(rarRichness1.table)[1]<-paste("ID")
RarRichness1 <- merge(shp.ua,rarRichness1.table)
writeSpatialShape(RarRichness1, "Ruta salida")

# excluir celdas con menos de 5.000 registros
index2 = which(rowSums(SpeciesMatrix)>5000) #Cambiar el numero de registros elegido
rarMatrix2 = SpeciesMatrix[index2,]
rarRichness2 = apply(rarMatrix2,1,function(x){rarefy(x,5000,se=FALSE,MARGIN=1)})
rarRichness2.table <- as.data.frame(rarRichness2)
setDT(rarRichness2.table, keep.rownames = TRUE)[]
names(rarRichness2.table)[1]<-paste("ID")
RarRichness2 <- merge(shp.ua,rarRichness2.table)
writeSpatialShape(RarRichness2, "Ruta salida")

# excluir celdas con menos de 50.000 registros
index3 = which(rowSums(SpeciesMatrix)>50000) #Cambiar el numero de registros elegido
rarMatrix3 = SpeciesMatrix[index3,]
rarRichness3 = apply(rarMatrix3,1,function(x){rarefy(x,50000,se=FALSE,MARGIN=1)})
rarRichness3.table <- as.data.frame(rarRichness3)
setDT(rarRichness3.table, keep.rownames = TRUE)[]
names(rarRichness3.table)[1]<-paste("ID")
RarRichness3 <- merge(shp.ua,rarRichness3.table)
writeSpatialShape(RarRichness3, "Ruta salida")

#RarCurve <- rarecurve(SpeciesMatrix.rar,step=1,sample=raremax, col = "blue", cex = 0.6, ylab='# especies', xlab='# de muestras')

###############################
##RIQUEZA A PATIR DE BIOMODELOS

alfaBiomodelos <- read.csv('Ruta tabla alfa Biomodelos')
BiomodRichness <- merge(shp.ua,alfaBiomodelos)
writeSpatialShape(BiomodRichness, "Ruta salida")

##################################################################
##MATRIZ DE CORRELACION DE ESTIMADORES DE LAS ESTIMACIONES DE ALFA

alfaBiomodelos$ID <- rownames(alfaBiomodelos)
nsp.data.frame$ID <- rownames(nsp.data.frame)
chaoRichness1$ID <- rownames(chaoRichness1)
bootRichness1$ID <- rownames(bootRichness1)
jackRichness1$ID <- rownames(jackRichness1)
aceRichness1$ID <- rownames(aceRichness1)
iceRichness1$ID <- rownames(iceRichness1)
rarRichness1.table$ID <- rownames(rarRichness1.table)
rarRichness2.table$ID <- rownames(rarRichness2.table)
rarRichness3.table$ID <- rownames(rarRichness3.table)

dataMatrix <- join_all(list(alfaBiomodelos,nsp.data.frame,chaoRichness1,bootRichness1,jackRichness1,aceRichness1,iceRichness1,rarRichness1.table,rarRichness2.table,rarRichness3.table),by="ID")
dataMatrix[is.na(dataMatrix)] <- 0
write.csv(dataMatrix,"Ruta salida")
dataMatrix$ID <- NULL 
matrixCorr = cor(dataMatrix,use="complete.obs",method="pearson")

write.csv(matrixCorr,"Ruta salida")
