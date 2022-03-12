## Recombiation map
## We use a parameter to sample recombination breakpoints
## We use a uniform distribution to sample the location of breakpoints
L         <- 35.9  # Size of the genome in Morgan -- Table 2 (Weir & Hill, 2011)
m         <- 1500  # Number of markers
dp        <- 1/m
grid      <- seq(dp,1-dp,by=dp)

recombine <- function(x,y){ # Haldane's map function
  nbSegments <- rpois(1,lambda = L)
  if(nbSegments>0){
    breakPoints <- sort(c(0,sample(grid,nbSegments-1),1))
    segments    <- cut(0:m,
                       breaks = unique(breakPoints) * m,
                       include.lowest = TRUE)
    idSegments  <- as.numeric(segments)
    if(rbinom(1,1,0.5)==1){
      z <- ifelse(idSegments%%2,y,x)
    }else{
      z <- ifelse(idSegments%%2,x,y)
    }
  }else{
    if(rbinom(1,1,0.5)==1){
      z <- x
    }else{
      z <- y
    }
  }
  return(z)
} 

#L <- 35.9
#m <- 1
#l <- L/m
#sqrt( ( 2*L - m + m * exp(-2*l) ) / (32 * L^2) )

## Pedigree (Fig. 1, W&H 2011)
##     A ---------  B ------------ C
##          |               |
##          |         -------------
##          |         |           |
## s(D) --- D     M---E           F---N
##       |          |               |
##       G---s(G)   H---s(H)        I---s(I)
##         |          |               |
##         J          K               L
##
#############################################

assessRelatedness <- function(x,y){
  cond1a <- (x[,1]==y[,1] & x[,2]!=y[,2]) | (x[,1]!=y[,1]  & x[,2]==y[,2])
  cond1b <- (x[,1]==y[,2] & x[,2]!=y[,1]) | (x[,1]!=y[,2]  & x[,2]==y[,1])
  
  IBD1  <- mean(cond1a | cond1b)
  IBD2  <- mean( (x[,1]==y[,1] & x[,2]==y[,2]) | (x[,1]==y[,2] & x[,2]==y[,1]) )
  0.5*IBD1 + IBD2
}

simulatePedigree <- function(){
  ## Simulate A, B, C, D, E and F
  A   <- cbind(rep("A1",m+1),rep("A2",m+1))
  B   <- cbind(rep("B1",m+1),rep("B2",m+1))
  C   <- cbind(rep("C1",m+1),rep("C2",m+1))
  D   <- cbind(recombine(A[,1],A[,2]),recombine(B[,1],B[,2]))
  E   <- cbind(recombine(C[,1],C[,2]),recombine(B[,1],B[,2]))
  F   <- cbind(recombine(C[,1],C[,2]),recombine(B[,1],B[,2]))
  
  ## Simulate M and N
  M   <- cbind(rep("M1",m+1),rep("M2",m+1))
  N   <- cbind(rep("N1",m+1),rep("N2",m+1))
  
  ## Simulate spouse of D
  sD   <- cbind(rep("sD1",m+1),rep("sD2",m+1))
  
  ## Simulate G 
  G    <- cbind(recombine(D[,1],D[,2]),recombine(sD[,1],sD[,2]))
  
  ## Simulate spouses of G, H and I
  sG   <- cbind(rep("sG1",m+1),rep("sG2",m+1))
  sH   <- cbind(rep("sH1",m+1),rep("sH2",m+1))
  sI   <- cbind(rep("sI1",m+1),rep("sI2",m+1))
  
  ## Simulate H
  H     <- cbind(recombine(M[,1],M[,2]),recombine(E[,1],E[,2]))
  
  ## Simulate I
  I     <- cbind(recombine(N[,1],N[,2]),recombine(F[,1],F[,2]))
  
  ## Simulate K
  K     <- cbind(recombine(H[,1],H[,2]),recombine(sH[,1],sH[,2]))
  
  ## Simulate J
  J     <- cbind(recombine(G[,1],G[,2]),recombine(sG[,1],sG[,2]))
  
  ## Simulate L
  L     <- cbind(recombine(I[,1],I[,2]),recombine(sI[,1],sI[,2]))
  return(list(A=A,B=B,C=C,D=D,E=E,F=F,G=G,H=H,I=I,J=J,K=K,L=L))
}

simulateRelDist <- function(id){
  ped     <- simulatePedigree()
  idInds  <- names(ped)
  nInds   <- length(ped)
  nPairs  <- nInds*(nInds-1)/2
  pairs   <- rep("",nPairs)
  R       <- rep(NA,nPairs)
  k <- 0
  for(i in 1:(nInds-1)){
    for(j in (i+1):nInds){
      k <- k + 1
      pairs[k] <- paste0(idInds[i],idInds[j])
      R[k] <- assessRelatedness(ped[[i]],ped[[j]])
    }
  }
  names(R) <- pairs
  return(R)
}

## Legends of Relationships
## -- Lineal relatives
## AD: Parent-Offspring (PO)
## AG: Grandparent - grandoffspring (GP-GO)
## AJ: Greatgrandparent - Greatgrandoffspring (GGP-GGO)
## DE: Halfsib (HS)
## DH: Half uncle–nephew (HUN)
## GH: Half cousins (HC)
## -- Full sibs and their descendants
## EF: Full sibs (FS)
## EI: Uncle–nephew (UN)
## EL: Great uncle - Great nephew (GUGN)
## HI: First cousins (C)
## HL: Cousins once removed (C1R)
## KL: Second cousins (2C)

expectedMean <- c(AD=1/2,AG=1/4,AJ=1/8,DE=1/4,DH=1/8,GH=1/16,
                  EF=1/2,EI=1/4,EL=1/8,HI=1/8,HL=1/16,KL=1/32)
expectedSD  <- c(AD=0,AG=0.0362,AJ=0.0291,DE=0.0277,DH=0.0256,GH=0.0188,
                 EF=0.0392,EI=0.0256,EL=0.0247,HI=0.0218,HL=0.0170,KL=0.0120) ## Table 2.

nSim   <- 10000
resSim <- do.call("rbind",lapply(1:nSim,simulateRelDist))
observedMean <- apply(resSim[,names(expectedMean)],2,mean)
observedSD <- apply(resSim[,names(expectedMean)],2,sd)

#boxplot(resSim[,names(expectedMean)])
#points(expectedMean,pch=19,col=2)
# op <- par(mfrow=c(1,2))
# par(mar=c(5,5,3,2))
# plot(expectedMean,observedMean,pch=19);abline(a=0,b=1,col=2)
# abline(v=0.5**c(1:5),col="grey",lty=4)
# par(mar=c(5,5,3,2))
# plot(expectedSD,observedSD,pch=19);abline(a=0,b=1,col=2)
# par(op)

rbind(expectedSD,observedSD)

## selected relationship for Figure 5
## FS (EF, 1/5), UN (EI, 1/4), C(HI, 1/8), C1R(HL, 1/16), 2C(KL, 1/32)
dt <- unique(resSim[,c("EF","EI","HI","HL","KL")])
colnames(dt) <- c("FS","AV","C","C1R","2C")

png("~/Desktop/Papers/TributeToBillHill2022/Fig5_revisited_12Mar2022_LY.png",width=2500,height=2000,res=300)
Cols <- c("coral1","dodgerblue","goldenrod","seagreen3","purple")
nc   <- 100
ymax <- 600
par(mar=c(5,5,3,2))
hist(dt[,1],breaks=nc,xlim=c(0,0.7),
     ylim=c(0,ymax),
     col=Cols[1],border=0,cex.lab=1.2,
     main="",density = 50,
     xlab="Realised relatedness",
     ylab="Frequency")

for(k in 2:5){
  hist(dt[,k],breaks=nc,add=TRUE,col=Cols[k],border=Cols[k],density = 50)
}
legend(0.0,ymax,legend=c("Full Siblings","Avuncular","First Cousin","First Cousin Once Removed","Second Cousin"),
       box.lty=0,fill=Cols,border=0,cex=1.1)
dev.off()


#abline(v=0.5**c(1:5),col=Cols,lty=2)
