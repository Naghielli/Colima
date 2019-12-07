####################################################################################
################################  Proyección Colima, 2016 -2040 ###############################
#Alvarez Chombo Naghielli
#Arrieta Arrieta Alí 


#Instalación de la paquería a utilizar 
rm(list=ls())
gc()
require(forecast)
require(ggplot2)
require(mvtnorm)
require(simPop)

#### Bases de datos para mortalidad ###

mxCOL=read.csv("mx.csv",header=T)
mxCOL<-mxCOL[,-c(1,2)]


### Bases de datos migración internacional ###

InmxCOL=read.csv("Imx5_Colima.csv",header=T)
EmigxCOL=read.csv("Emx5_Colima.csv",header=T)
NxCOL=read.csv("Nx.csv",header=T)

InmxCOL.F5<-data.frame(InmxCOL[InmxCOL$Sexo=="Mujeres",-1])
InmxCOL.H5<-data.frame(InmxCOL[InmxCOL$Sexo=="Hombres",-1])
EmigxCOL.F5<-data.frame(EmigxCOL[EmigxCOL$Sexo=="Mujeres",-1])
EmigxCOL.H5<-data.frame(EmigxCOL[EmigxCOL$Sexo=="Hombres",-1])

InmxCOL.F=matrix(0,81,46)
InmxCOL.H=matrix(0,81,46)
EmigxCOL.F=matrix(0,81,46)
EmigxCOL.H=matrix(0,81,46)

for(i in 0:45){
  InmxCOL.F[,i+1]=as.matrix(sprague(InmxCOL.F5[,i+2]))
  InmxCOL.H[,i+1]=as.matrix(sprague(InmxCOL.H5[,i+2]))
  EmigxCOL.F[,i+1]=as.matrix(sprague(EmigxCOL.F5[,i+2]))
  EmigxCOL.H[,i+1]=as.matrix(sprague(EmigxCOL.H5[,i+2]))
}

maNxM=as.matrix(NxCOL[NxCOL$Sexo=="Hombre",-c(1:2)])
maNxF=as.matrix(NxCOL[NxCOL$Sexo=="Mujer",-c(1:2)])

#Calculo de las tasas
ixt.MCOL<-abs(InmxCOL.H)/maNxM
ixt.FCOL<-abs(InmxCOL.F)/maNxF
ext.MCOL<-abs(EmigxCOL.H)/maNxM
ext.FCOL<-abs(EmigxCOL.F)/maNxF

#Seleccionamos solo la información a partir de 1990 (el ajuste es mejor)
ixt.FCOL<-ixt.FCOL[,c(21:46)]
ixt.MCOL<-ixt.MCOL[,c(21:46)]
ixt.MCOL<-ixt.MCOL[-81,]
ext.FCOL<-ext.FCOL[,c(21:46)]
ext.MCOL<-ext.MCOL[,c(21:46)]
ext.MCOL<-ext.MCOL[-81,]


### Bases de datos Fecundidad ###

fx5<- read.csv("Fx5_II.CSV",
               header=TRUE)

fx5q<- read.csv("Fx_quinquenios.CSV",
                header=TRUE)

mx<- read.csv("mx.CSV",
              header=TRUE)

asfr=function(fx5,year){
  
  fx0=fx5[fx5$Year==year,"fx"]
  
  Fx=5*cumsum(fx0)
  TGF=Fx[7]
  FxF=Fx/TGF
  
  x5=seq(17.5,47.5,5)
  Yx= log(-log(FxF))
  
  
  Yx.lm=lm(Yx[-7] ~ x5[-7])
  a=Yx.lm$coefficients[1]
  b=Yx.lm$coefficients[2]
  
  A=exp(-exp(a))
  B=exp(b)
  x1=c(15:50)
  
  (Fx.estim=TGF*A^(B^(x1)))
  
  fx=Fx.estim[2:36]-Fx.estim[1:35]
  fx=c(Fx.estim[1],fx)
  
  return(fx)
}

fx1=data.frame(matrix(0,36,46))
row.names(fx1)=c(15:50)
names(fx1)=c(1970:2015)


for(i in 1:46){
  
  fx1[,i]=asfr(fx5q,1969+i)
  
}

#Seleccionamos los años a partir de 1990, el ajuste es mejor. 
fx <-fx1[1:35,20:45]  


### Definición de la información para la proyección ###

edades <- dim(mxCOL)[1] 
tiempo.mort <- dim(mxCOL)[2]
añoini.mort <- 1990
añobase <- 2015
horizonte <- 25
añofin <- añobase+horizonte
tiempo.tot <- tiempo.mort+horizonte
edades.fec <-dim(fx)[1]
añoini.fec <- 1990
tiempo.fec <-dim(fx)[2]
edades.migCOL <- dim(ixt.FCOL)[1]
tiempo.migCOL <- dim(ixt.FCOL)[2]
añoini.migCOL <- 1990

### Inicia la estimacion con el método de Lee Carter (1992)


lc.svd <- function(m,edades,tiempo1, tiempo2,ln){
  if(ln==TRUE){
    lm <- log(m)
  }else{
    lm<-
      m
  }
  
  ax <- rowMeans(lm[,tiempo1:tiempo2])
  lm_a <- lm -ax
  d <- matrix(0,nr = min(edades,tiempo2),
              nc = min(edades,tiempo2))
  
  diag(d) <- svd(lm_a)$d
  
  kt <-(d%*%t(-svd(lm_a)$v))
  bx <- -svd(lm_a)$u
  
  lc.svd=list(ax = ax, bx = bx, kt = kt, D = d)
  
  
}
tabmort <- function(m,edades,sex){
  
  mx <- m
  
  nax <- matrix(0.5,dim(mx)[1],dim(mx)[2])
  ## 1 MUJERES 2 HOMBRES
  if(sex==1){
    for(i in 1:dim(mx)[2]){
      if(mx[1,i]<0.01724){
        nax[1,i] <- 0.14903-2.05527*mx[1,i]
      }else if(mx[1,i]>=0.01724 & mx[1,i]<0.06891){
        nax[1,i] <- 0.04667+3.88089*mx[1,i]
      }else{nax[1,i] <- 0.31411}
    }
  }else{
    for(i in 1:dim(mx)[2]){
      if(mx[1,i]<0.023){
        nax[1,i] <- 0.14929-1.99545*mx[1,i]
      }else if(mx[1,i]>=0.023 & mx[1,i]<0.08307){
        nax[1,i] <- 0.02832+3.26021*mx[1,i]
      }else{nax[1,i] <- 0.29915}
    }
  }
  
  
  nax[edades,] <- 1/mx[edades,]
  
  qx<-matrix(1,dim(mx)[1],dim(mx)[2])
  
  for(i in 1:(dim(mx)[1])){
    qx[i,]<-mx[i,]/(1+(1-nax[i,])*mx[i,])
  }
  
  px <- 1-qx
  
  lx<-matrix(1,dim(mx)[1],dim(mx)[2])
  
  for(i in 2:dim(mx)[1]){
    lx[i,] <- lx[i-1,]*px[i-1,]
  }
  
  dx <- matrix(0,dim(mx)[1],dim(mx)[2])
  dx[dim(mx)[1],] <- lx[dim(mx)[1],]
  for(i in 1:(dim(mx)[1]-1)){
    dx[i,]<-lx[i,]-lx[i+1,]
  }
  
  
  Lx<-matrix(0,dim(mx)[1],dim(mx)[2])
  Lx[1,] <- dx[1,]/mx[1,]
  Lx[edades,] <- dx[edades,]/mx[edades,]
  for(i in 2:(edades-1)){
    Lx[i,]<-(lx[i,]+lx[i+1,])/2
  }
  
  Tx<-matrix(0,dim(mx)[1],dim(mx)[2])
  Tx[edades,]<-Lx[edades,]
  for(i in (edades-1):1){
    Tx[i,]<-Lx[i,]+Tx[i+1,]
  }
  
  ex <- Tx/lx
  
  Sx<-matrix(0,(dim(mx)[1]+1),dim(mx)[2])
  Sx[1,]<-Lx[1,]/lx[1,]
  Sx[(edades+1),] <- Tx[edades,]/Tx[(edades-1),]
  for(i in 2:edades){
    Sx[i,]<-Lx[i,]/Lx[i-1,]
  }
  
  tabmort <- list(Edad=c(0:(edades-1)),mx=mx, nax=nax, qx=qx, 
                  px=px, lx=lx, dx=dx, Lx=Lx, Tx=Tx, ex=ex, Sx=Sx)
}

lc.mort <-lc.svd(mxCOL, edades, tiempo1 = 41, 
                 tiempo2 = tiempo.mort, 
                 ln=TRUE)

lc.fec <-lc.svd(fx, edades=edades.fec, 
                tiempo1 = 10, 
                tiempo2 = tiempo.fec,
                ln=TRUE)

lc.inmFCOL <-lc.svd(ixt.FCOL, edades=edades.migCOL, 
                   tiempo1 = 1, 
                   tiempo2 = tiempo.migCOL,
                   ln=TRUE)

lc.inmMCOL <-lc.svd(ixt.MCOL, edades=edades.migCOL, 
                   tiempo1 = 1, 
                   tiempo2 = tiempo.migCOL,
                   ln=TRUE)

lc.emigFCOL <-lc.svd(ext.FCOL, edades=edades.migCOL, 
                    tiempo1 = 1, 
                    tiempo2 = tiempo.migCOL,
                    ln=TRUE)

lc.emigMCOL <-lc.svd(ext.MCOL, edades=edades.migCOL, 
                    tiempo1 = 1, 
                    tiempo2 = tiempo.migCOL,
                    ln=TRUE)

kt1.fit <-auto.arima(lc.mort$kt[1,], trace=T, d=1)

ft1.fit <-auto.arima(lc.fec$kt[1,], trace=T, d=1)

it1F.fit <- auto.arima(lc.inmFCOL$kt[1,], trace=T, allowdrift = F,d=1)

it1M.fit <- auto.arima(lc.inmMCOL$kt[1,], trace=T, allowdrift = F)

et1F.fit <- auto.arima(lc.emigFCOL$kt[1,], trace=T, allowdrift = F)

et1M.fit <- auto.arima(lc.emigMCOL$kt[1,], trace=T, allowdrift = F)


# Inicia la proyeccion para un horizonte de 25 años

h <- 25

kt.for <- forecast(kt1.fit, h = h, c(95))

ft.for <- forecast(ft1.fit, h = h, c(95))

itF.for <- forecast(it1F.fit, h = h, c(95))
itM.for <- forecast(it1M.fit, h = h, c(95))

etF.for <- forecast(et1F.fit, h = h, c(95))
etM.for <- forecast(et1M.fit, h = h, c(95))


#Tasas centrales de mortalidad
mx.for <- exp(lc.mort$ax + lc.mort$bx[,1]%*%t(kt.for$mean))

#Superviviencia
SxF.for <- tabmort(mx.for[111:220,], edades = 110, sex=1)$Sx
SxM.for <- tabmort(mx.for[1:110,], edades = 110, sex=2)$Sx

#Tasas especifícas de fecundidad 
fx.for<-exp(lc.fec$ax + lc.fec$bx[,1]%*%t(ft.for$mean))
TGF=colSums(fx.for)           

#Tasas especificas de migración
ixF.For <- rbind(exp(lc.inmFCOL$ax + lc.inmFCOL$bx[,1]%*%t(itF.for$mean)),
                 matrix(0,0,25))

ixM.For <- rbind(exp(lc.inmMCOL$ax + lc.inmMCOL$bx[,1]%*%t(itM.for$mean)),
                 matrix(0,0,25))

exF.For <- rbind(exp(lc.emigFCOL$ax + lc.emigFCOL$bx[,1]%*%t(etF.for$mean)),
                 matrix(0,0,25))

exM.For <- rbind(exp(lc.emigMCOL$ax + lc.emigMCOL$bx[,1]%*%t(etM.for$mean)),
                 matrix(0,0,25))

#Cohorte Componente (Población)

Px=read.csv("Px.csv",header=T)

PxFCOL <- Px[Px$Sexo=="Mujeres", -c(1,2)]
PxMCOL <- Px[Px$Sexo=="Hombres", -c(1,2)]

NxFCOL <- Px[NxCOL$Sexo=="Mujer", -c(1,2)]
NxMCOL <- Px[NxCOL$Sexo=="Hombre", -c(1,2)]

PxF.for <- matrix(0,81,26)
PxM.for <- matrix(0,81,26)

NxF.for <- matrix(0,81,26)
NxM.for <- matrix(0,81,26)

PxF.for[,1] <- PxFCOL[,"X2016"]
PxM.for[,1] <- PxMCOL[,"X2016"]

NxF.for[,1] <- NxFCOL[,"X2015"]
NxM.for[,1] <- NxMCOL[,"X2015"]

Bx <- matrix(0,35,26)
BF <- vector(length=26)
BM <- vector(length=26)


## Primero para las mujeres, pues el modelo es femenino dominante 
for(i in 2:26){
  
  #Con pob a mitad de año 1 a 108
  PxF.for[2:80,i] <- (PxF.for[1:79,i-1] + 
                        0.5*NxF.for[1:79,i-1]*ixF.For[1:79,i-1])* SxF.for[1:79,i-1]+
    NxF.for[2:80,i-1]*0.5*ixF.For[2:80,i-1]-
    NxF.for[1:79,i-1]*exF.For[1:79,i-1]
  
  #Para el último grupo de edad 
  PxF.for[81,i]<-(PxF.for[80,i-1] + 
                    0.5*NxF.for[80,i-1]*ixF.For[80,i-1])*SxF.for[80,i-1] -
    NxF.for[80, i-1]*exF.For[80,i-1] +
    (PxF.for[81,i-1] + 
       NxF.for[81,i-1]*0.5*ixF.For[81,i-1])*SxF.for[81,i-1]+
    NxF.for[81,i-1]*0.5*ixF.For[81,i-1]-
    NxF.for[81,i-1]*0.5*exF.For[81,i-1]
  
  
  #Nacimientos 
  Bx[,i-1] <- fx.for[,i-1]*(PxF.for[16:50,i-1]+
                              0.5*NxF.for[16:50,i-1]*ixF.For[16:50,i-1] +
                              PxF.for[16:50,i]) * 0.5
  
  BF[i-1] <- (1/2.05)*sum(Bx[,i-1])
  
  
  #Primer grupo de edad
  
  PxF.for[1,i] <- BF[1]*SxF.for[1,i-1] + 
    NxF.for[1,i-1]*0.5*ixF.For[1,i-1]+
    NxF.for[1,i-1]*exF.For[1,i-1]
  
  
  #POb mitad de año
  
  NxF.for[,i] <- 0.5*(PxF.for[,i-1] + PxF.for[,i])
}

matplot(NxF.for[-81,], type="l")

### Hombres #####
for(i in 2:26){
  
  #Con pob a mitad de año 1 a 108
  PxM.for[2:80,i] <- (PxM.for[1:79,i-1] + 
                        0.5*NxM.for[1:79,i-1]*ixM.For[1:79,i-1])* SxM.for[1:79,i-1]+
    NxM.for[2:80,i-1]*0.5*ixM.For[2:80,i-1]-
    NxM.for[1:79,i-1]*exM.For[1:79,i-1]
  
  #ultimo grupo 109 y mas
  PxM.for[81,i]<-(PxM.for[80,i-1] + 
                    0.5*NxM.for[80,i-1]*ixM.For[80,i-1])*SxM.for[80,i-1] -
    NxM.for[80, i-1]*exM.For[80,i-1] +
    (PxM.for[81,i-1] + 
       NxM.for[81,i-1]*0.5*ixM.For[80,i-1])*SxM.for[81,i-1]+
    NxM.for[81,i-1]*0.5*ixM.For[80,i-1]-
    NxM.for[81,i-1]*0.5*exM.For[80,i-1]
  
  
  #Nacimientos 
  
  BM[i-1] <- (1.05/2.05)*sum(Bx[,i-1])
  
  
  #primer grupo de edad
  
  PxM.for[1,i] <- BM[1]*SxM.for[1,i-1] + 
    NxM.for[1,i-1]*0.5*ixM.For[1,i-1]+
    NxM.for[1,i-1]*exM.For[1,i-1]
  
  
  #POb mitad de año
  
  NxM.for[,i] <- 0.5*(PxM.for[,i-1] + PxM.for[,i])
}

matplot(NxM.for[-81,], type="l")

POBM =colSums(PxM.for)
POBF=colSums(PxF.for)       
POB <-   colSums(PxF.for) + colSums(PxM.for)          
POB


PxF.for <-round(PxF.for,0)
PxM.for <-round(PxM.for,0)

matplot(PxF.for[-81,],type="l")
colSums(PxF.for[-81,])       


### Graficas

### Fecundidad 

par(mfrow=c(1,1))

fx.forup<-exp(lc.fec$ax + lc.fec$bx[,1]%*%t(ft.for$upper))
TGFup=colSums(fx.forup)           

fx.forlow<-exp(lc.fec$ax + lc.fec$bx[,1]%*%t(ft.for$lower))
TGFlowe=colSums(fx.forlow)           


matplot(seq(2016,2040,by=1), cbind(TGF,TGFup,TGFlowe),type="l",col=c("slateblue1", "blue","turquoise3"),lwd = 3,
        xlab = "Años",
ylab= "TGF", 
main= " Proyección de la Tasa Global de Fecundidad, 2016 - 2040",
sub="Estado: Colima",
bg=2)
labels <- c("Medio","Alto","Bajo")
legend("topright",inset =c(0.8 ,0.75), labels,lwd =2, cex=0.6, lty=1:2, col=c("slateblue1", "blue","turquoise3"))
grid(col="lightgray")

edades=seq(15,49,1)
matplot(edades,cbind(fx.for[,25],fx.forup[,25],fx.forlow[,25],fx.for[,1],fx.forup[,1],fx.forlow[,1]),
        type="l",col=c("slateblue1", "blue","turquoise3", "olivedrab3",  "olivedrab", "palegreen"),lwd = 3,
        xlab = "Edades",
        ylab= "fx", 
        main= "Proyección de las Tasas específicas de fecundidad, 2016 - 2040",
        sub="Estado: Colima",
        bg=2)
labels <- c("Medio - 2040","Alto - 2040","Bajo - 2040", "Medio - 2016","Alto - 2016","Bajo - 2016")

legend("topright",inset =c(0.1 ,0.2), labels,lwd =2, cex=0.6, lty=1:2, col=c("slateblue1", "blue","turquoise3",  "olivedrab3","olivedrab", "palegreen"))
grid(col="lightgray")

matplot(fx.for, type="l")

# Mortalidad
require(faraway)
par(mfrow=c(1,1))
matplot(cbind(log(mx.for[1:110,]),log(mx.for[111:220,])),col=c("lightseagreen", "firebrick1"), type="l", lwd=2, 
        xlab = "Años",
        ylab= "ln (Tasas centrales de la mortalidad)", 
        main= "Tasas Centrales de la mortalidad, 2016 -2040",
        sub="Estado: Colima",)
labels <- c("Hombres","Mujeres")
legend("topright",inset =c(0.1 ,0.7), labels,lwd =2, cex=0.6, lty=1:2, col=c("lightseagreen", "firebrick1"))
grid(col="lightgray")


par(mfrow=c(1,2))

hopeM=tabmort(mx.for[1:110,], edades = 110, sex=2)$ex
hopeF=tabmort(mx.for[111:220,], edades = 110, sex=1)$ex
matplot(hopeF,type="l")
matplot(hopeM,type="l")

eF0=hopeF[1,]
eM0=hopeM[1,]
parmfr
plot(seq(2016,2040,by=1),eM0,type="l")

matplot(seq(2016,2040,by=1),eF0,type="l", lty = 2,lwd = 3,
        col="firebrick1", main= "Esperanza de Vida al Nacimiento Mujeres, 2016 -2040", sub="Estado: Colima",
        xlab = "Años",
        ylab= "Esperanza de vida al nacimiento",
        bg=2)

grid(col="lightgray")

matplot(seq(2016,2040,by=1),eM0,type="l", lty = 2,lwd = 3,
        col="slateblue3", main= "Esperanza de Vida al Nacimiento Hombres, 2016 -2040", sub="Estado: Colima",
        xlab = "Años",
        ylab= "Esperanza de vida al nacimiento",
        bg=2)

grid(col="lightgray")


### Migración ####

par(mfrow=c(1,2))

ixF.For1 =ixF.For * NxF.for[,-26]
exF.For1 = exF.For * NxF.for[,-26]
SNMF = -exF.For1 + ixF.For1

ixM.For1 =ixM.For * NxM.for[1:80,-26]
exM.For1 = exM.For * NxM.for[1:80,-26]
SNMM = ixM.For1 - exM.For1

SNMT= SNMF[1:80] + SNMM 
snmtota = colSums(SNMT)

par(mfrow=c(1,2))


matplot(cbind(SNMF[,1], SNMF[,25]),type="l", lty = 2,lwd = 3, col=c("lightseagreen", "firebrick1"),
  main= "Saldo Neto Migratorio, mujeres", sub="Estado: Colima",
        xlab = "Edades",
        ylab= "SNM",
        bg=2)

labels <- c("2016","2040")
legend("topright",inset =c(0.1 ,0.2), labels,lwd =2, cex=0.6, lty=1:2, col=c("lightseagreen", "firebrick1"))
grid(col="lightgray")

matplot(cbind(SNMM[,1], SNMM[,25]),type="l", lty = 2,lwd = 3, col=c("lightseagreen", "firebrick1"),
        main= "Saldo Neto Migratorio, hombres", sub="Estado: Colima",
        xlab = "Edades",
        ylab= "SNM",
        bg=2)

grid(col="lightgray")

labels <- c("2016","2040")
legend("topright",inset =c(0.1 ,0.1), labels,lwd =2, cex=0.6, lty=1:2, col=c("lightseagreen", "firebrick1"))

par(mfrow=c(1,2))

matplot(cbind(SNMT[,1], SNMT[,25]),type="l", lty = 2,lwd = 3, col=c("lightseagreen", "firebrick1"),
        main= "Saldo Neto Migratorio Total", sub="Estado: Colima",
        xlab = "Edades",
        ylab= "Personas",
        bg=2)

grid(col="lightgray")

labels <- c("2016","2040")
legend("topright",inset =c(0.1 ,0.6), labels,lwd =2, cex=0.6, lty=1:2, col=c("lightseagreen", "firebrick1"))



matplot(seq(2016,2040,by=1),snmtota, type = "l", col="firebrick1", lwd=2, lty=1,
        main= "Saldo Neto Migratorio Total, 2016 -2040", sub="Estado: Colima",
        xlab = "Años",
        ylab= "Personas",)
grid(col="lightgray")

## Población 

par(mfrow=c(1,2))

matplot(seq(2015,2040,by=1),cbind(POB), type = "l", col=c("firebrick1"), lwd=2, lty=1,
        main= "Población proyectada total, 2016 - 2040", sub="Estado: Colima",
        xlab = "Años",
        ylab= "Personas")
grid(col="lightgray")


matplot(seq(2015,2040,by=1),cbind(POBM, POBF), type = "l", col=c("lightseagreen", "firebrick1"), lwd=2, lty=1,
        main= "Población proyectada total por sexo, 2016 - 2040", sub="Estado: Colima",
        xlab = "Años",
        ylab= "Personas")
grid(col="lightgray")
labels <- c("Hombres","Mujeres")
legend("topright",inset =c(0.1 ,0.7), labels,lwd =2, cex=0.6, lty=1:2, col=c("lightseagreen", "firebrick1"))
grid(col="lightgray")





par(mfrow=c(1,2))

matplot(cbind(PxF.for[-81,]),type = "l", lwd=1, lty=1,
        main= "Población proyectada total mujeres, 2016 - 2040", sub="Estado: Colima",
        xlab = "Edades",
        ylab= "Personas")
grid(col="lightgray")

matplot(cbind(PxM.for[-81,]),type = "l", lwd=1, lty=1,
        main= "Población proyectada hombres, 2016 - 2040", sub="Estado: Colima",
        xlab = "Edades",
        ylab= "Personas")
grid(col="lightgray")

