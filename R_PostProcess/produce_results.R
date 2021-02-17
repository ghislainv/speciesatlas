
# ==================
# Creating Maps
# ==================

############
#== Init ==#
############

library(rasterVis)
library(shape)
library(viridis)
load("C:/Users/chocho/Desktop/Rapport/analysis/models.rda")
### Looping through :
fut.var <- list(c("cc","he","gs"),c("00","45","85"),c("2080"))
def <- c("_defo","")
deforested <- c("Yes","No")
dispersal <- c("full","zero") # In that order only

setwd("C:/Users/chocho/Desktop/Ddata/stageBio/Données/Madagascar/rastersv2")

model.var <-  c("temp","tseas","prec","pseas","cwd","foret")

environ <- stack(brick("current.tif"), brick("environ.tif"))
names(environ) <- c(paste0("current",1:36),"temp",paste0("current",38:39),"tseas",paste0("current",41:47),
                    "prec","bioclim1",
                    "bioclim2","pseas",paste0("current",52:67),"pet","cwd","ndm","alt",paste0("environ",2:7),"watershed","foret")
environ <- stack(environ) # Transform back from RasterBrick to rastertack

environFut <- stack(brick("cc_85_2080.tif"), brick("environ_2080.tif"))
names(environFut) <- c(paste0("current",1:36),"temp",paste0("current",38:39),"tseas",paste0("current",41:47),
                    "prec","bioclim1",
                    "bioclim2","pseas",paste0("current",52:67),"pet","cwd","ndm","alt",paste0("environ",2:7),"watershed","foret")
environFut <- stack(environFut) # Transform back from RasterBrick to rastertack

s <- stack(environ[[model.var]])
setwd("C:/Users/chocho/Desktop/Rapport/analysis/")

### Data structure
rangeShift <- data.frame(species=character(),x1=numeric(),
                         y1=numeric(),x2=numeric(),y2=numeric(),distance=numeric(),direction=numeric())
out.SDA.zero <- data.frame(species=character(),area.pres=numeric(),
                           area.fut=numeric(),perc.chang=numeric())
out.SDA.full <- data.frame(species=character(),area.pres=numeric(),
                           area.fut=numeric(),perc.chang=numeric())

####################################################
#== Gathering present niche & computing models performance ==#
####################################################
taxon.first <- TRUE
list.current.niche <- list()
for(i in 1:length(enough)){
  if(enough[i]==T){
    ### Loading Data
    load(paste0("figures/",taxon.names,"/",sp.dir[i],"/plotting.rda"))
    
    ### Species Statistical Model Perfomance
    Perf.mods <- get.Model.Perf(Perf.mods)
    
    ### Current species distribution
    #Number of models computed :
    nbCategory <- dim(stack(paste0("BIOMOD/",sp.dir[i],"/proj_current/proj_current_",sp.dir[i],".grd")))[[3]]+1##Adding +1 for zeros
    treshold <- seq(0,1000,length.out=nbCategory)[floor(nbCategory/2)+1]
    #Models outputs :
    ca <- stack(paste0("BIOMOD/",sp.dir[i],"/proj_current/proj_current_",sp.dir[i],"_ensemble.grd"))[[1]]
    ca[values(ca)<treshold] <- 0
    ca[values(ca)>=treshold] <- 1
    
    #Irreplaceability
    
    
    list.current.niche[[i]] <- ca
    if(taxon.first==TRUE){
      Perf.Ens.Tot <- Perf.ca
      Perf.mods.TSS.Tot <- Perf.mods[[1]]
      Perf.mods.ROC.Tot <- Perf.mods[[2]]
      taxon.first <- FALSE
    }else{
      Perf.Ens.Tot <- rbind(Perf.Ens.Tot,Perf.ca)
      Perf.mods.TSS.Tot <- rbind(Perf.mods.TSS.Tot,Perf.mods[[1]])
      Perf.mods.ROC.Tot <- rbind(Perf.mods.ROC.Tot,Perf.mods[[2]])
    }
  }else{
    list.current.niche[i] <- NULL
  }
}

### Computing models perfomance
out.TSS <- compute.Model.Perf.Tot(Perf.mods.TSS.Tot)
out.ROC <- compute.Model.Perf.Tot(Perf.mods.ROC.Tot)
out.Ens <- cbind(ROC=c(Min=min(Perf.Ens.Tot$ROC),Max=max(Perf.Ens.Tot$ROC),Median=median(Perf.Ens.Tot$ROC)),
                 OA=c(Min=min(Perf.Ens.Tot$OA),Max=max(Perf.Ens.Tot$OA),Median=median(Perf.Ens.Tot$OA)),
                 TSS=c(Min=min(Perf.Ens.Tot$TSS),Max=max(Perf.Ens.Tot$TSS),Median=median(Perf.Ens.Tot$TSS)))
### Writing models performance
write.table(out.TSS,paste0("figures/",taxon.names,"/","overall/","Perf.mods.TSS.txt"),sep=";",row.names=F)
write.table(out.ROC,paste0("figures/",taxon.names,"/","overall/","Perf.mods.ROC.txt"),sep=";",row.names=F)
write.table(out.Ens,paste0("figures/",taxon.names,"/","overall/","Perf.Ens.Tot.txt"),sep=";",row.names=T)


####################################################
#== Gathering future niche                       ==#
####################################################

list.futur.niche <- list()
for (m in 1:length(def)) { #Deforestation or not
  list.futur.niche[[m]] <- list()
  for (j in 1:length(fut.var[[2]])) { #Rcp: 4.5 or 8.5
    list.futur.niche[[m]][[j]] <- list()
    for (l in 1:length(fut.var[[3]])) { #Time period: 2085
      list.futur.niche[[m]][[j]][[l]] <- list()
      #For every species :
      print(paste0(fut.var[[2]][j],"_",fut.var[[3]][l],def[m]))
      for(i in 1:length(enough)){
        if(enough[i]==T & !(def[m]!="_defo"&fut.var[[2]][j]=="00")){
          caS <- stack()
          for (mc in 1:length(fut.var[[1]])) {
            name.f <- paste0("/proj_",fut.var[[1]][mc],"_",fut.var[[2]][j],"_",fut.var[[3]][l],def[m])
            pred.f <- stack(paste0("BIOMOD/",sp.dir[i],name.f,name.f,"_",sp.dir[i],"_ensemble.grd"))
            caS <- addLayer(caS,pred.f[[1]])
          }
          # Sum and binarize
          caFut <- stack()
          caFut <- sum(caS)
          #Number of models computed :
          nbCategoryFuture <- (nbCategory-1)*length(fut.var[[1]])+1
          tresholdFuture <- seq(0,1000*length(fut.var[[1]]),length.out=nbCategoryFuture)[floor(nbCategoryFuture/2)+1]
          #Models outputs :
          caFut[values(caFut)<tresholdFuture] <- 0
          caFut[values(caFut)>=tresholdFuture] <- 1
          
          list.futur.niche[[m]][[j]][[l]][[i]] <- caFut
        }
      }
    }
  }
}


#########################################################
#== Computing gain, lost, SDA etc.. for every models  ==#
##############################7###########################

for(d in 1:length(dispersal)){
  for (m in 1:length(def)) { #Deforestation or not
    for (j in 1:length(fut.var[[2]])) { #Rcp: 4.5 or 8.5
      for (l in 1:length(fut.var[[3]])) { #Time period: 2085
        print(paste0(fut.var[[2]][j],"_",fut.var[[3]][l],def[m],"_",dispersal[d]))
        taxon.first <- TRUE
        for(i in 1:length(enough)){
          if(enough[[i]]==T & !(def[m]!="_defo"&fut.var[[2]][j]=="00")){
            ### Loading SDA
            load(paste0("figures/",taxon.names,"/",sp.dir[i],"/plotting.rda"))
            SDA.fut <- SDA.fut[SDA.fut$rcp==fut.var[[2]][j]&SDA.fut$yr==fut.var[[3]][l]&
                               SDA.fut$disp==dispersal[d]&SDA.fut$defor==deforested[m],]
            SDA.fut <- unique(SDA.fut)
            ### Selecting current niche
            ca <- list.current.niche[[i]]
            ### Selecting futur niche
            caFut <- list.futur.niche[[m]][[j]][[l]][[i]]

            #Removing colonized area for zero dispersal
            if(dispersal[d]=="zero"){
              values(caFut)[values(ca)==0] <- 0
            }

            #Computing range lost for each species
            perte <- ca-caFut == 1
            gain <- ca-caFut == -1
            stable <- (ca+caFut) - (perte+gain)

            #Irreplaceability
            irr.pres <- ca*(1/SDA.fut$area.pres)
            if(SDA.fut$area.fut==0) irr.fut <- caFut*1
            else irr.fut <- caFut*(1/SDA.fut$area.fut) #i=7 ou 15

            #Computing range shift for each species
            centroidCa <- colMeans(xyFromCell(ca, which(ca[]==1)))
            centroidCaFut <- colMeans(xyFromCell(caFut, which(caFut[]==1)))
            vectCentroids <- c((centroidCaFut[1]-centroidCa[1])*(centroidCaFut[1]-centroidCa[1]),
                               (centroidCaFut[2]-centroidCa[2])*(centroidCaFut[2]-centroidCa[2]))
            direction <- vectorDirection(vectCentroids)
            distCentroids <- sqrt(abs(vectCentroids[1]+vectCentroids[2]))

            if(taxon.first==TRUE){
              rangeShift <- data.frame(species=sp.dir[i],x1=centroidCa[1],
                                       y1=centroidCa[2],x2=centroidCaFut[1],y2=centroidCaFut[2],distance=distCentroids,direction=direction)
              out.caTot <- ca
              out.caFutTot<-caFut
              out.perteTot <- perte
              out.gainTot <- gain
              out.stableTot <- stable
              out.irrPres <- irr.pres
              out.irrFut <- irr.fut
              SDATot <- data.frame(species=sp.dir[i],area.pres=SDA.fut$area.pres,
              area.fut=SDA.fut$area.fut,perc.chang=SDA.fut$perc.change)
              taxon.first <- FALSE
            }else {
              rangeShift <- rbind(rangeShift,data.frame(species=sp.dir[i],x1=centroidCa[1],
                                                        y1=centroidCa[2],x2=centroidCaFut[1],y2=centroidCaFut[2],distance=distCentroids,direction=direction))
              out.caTot <- ca + out.caTot
              out.caFutTot<-caFut + out.caFutTot
              out.perteTot <- out.perteTot + perte
              out.gainTot <- out.gainTot + gain
              out.stableTot <- out.stableTot + stable
              out.irrPres <- out.irrPres + irr.pres
              out.irrFut <- out.irrFut + irr.fut
              SDATot <- rbind(SDATot,data.frame(species=sp.dir[i],area.pres=SDA.fut$area.pres,
                                                            area.fut=SDA.fut$area.fut,perc.chang=SDA.fut$perc.change))
            }
          }
        }
        if(!(def[m]!="_defo"&fut.var[[2]][j]=="00")){
          ### Computing results ###
          #SDA
          out.SDA <- compute.sum.SDA(SDATot)

          #Lost : removing where there is no species
          values(out.perteTot)[values(out.caTot)==0]<-NA

          #Gain : removing where there is no species in the future
          values(out.gainTot)[values(out.caFutTot)==0]<-NA

          #Calculating turn-over
          out.turnover <- (out.perteTot+out.gainTot)/(out.caTot+out.gainTot)

          #Calculating areas of changes
          out.changes <- compute.areas.changes(out.caTot,out.perteTot)

          #Calculating irreplaceability difference
          out.irrFut[out.caFutTot==0]<-NA
          # out.irrDiff <- out.irrFut-out.irrPres
          # out.irrDiff[out.irrDiff[]== 0] <- NA #Areas = 0 are the one with no species
          # out.irrDiffScaled <- out.irrDiff
          # out.irrDiffScaled[] <- (out.irrDiff[]-min(out.irrDiff[],na.rm=T))/(max(out.irrDiff[],na.rm=T)-min(out.irrDiff[],na.rm=T))

          # #Normalizing present irreplaceability
          # out.irrPresScaled <- out.irrPres
          # out.irrPresScaled[] <- (out.irrPres[]-min(out.irrPres[],na.rm=T))/(max(out.irrPres[],na.rm=T)-min(out.irrPres[],na.rm=T))

          # #Computing refuges areas
          # out.refuges <- raster(s)
          # out.refuges[!is.na(s$temp)] <- 0
          # out.refuges[out.irrFut[]>quantile(out.irrFut[values(out.irrFut)!=0],0.9,na.rm=T)&out.turnover[]<0.5] <- 1

          #Stats about refuges areas
          # paMada <- raster::shapefile("data/shape/protected.areas.v2.shp")
          # out.refuges.stats <- get.refuges.stats(enough,list.current.niche,list.futur.niche,out.refuges,paMada)

          #### Writing outputs ###
          path <- paste0(fut.var[[2]][j],"_",fut.var[[3]][l],def[m],"_",dispersal[d])
          
          #Refuges
          # write.table(out.refuges.stats,paste0("figures/",taxon.names,"/","overall/refuges.",path,".txt"),sep=";",row.names=F)
          
          #Rastef
          # results <- stack(out.caFutTot,out.caTot,out.gainTot,out.perteTot,out.stableTot,out.turnover,out.changes,
          #                  out.refuges,out.irrFut)
          # writeRaster(results,paste0("figures/",taxon.names,"/","overall/",path,".tiff"),overwrite=T)
          # 
          # #Range shift
          # write.table(rangeShift,paste0("figures/",taxon.names,"/","overall/rangeShift.",path,".txt"),sep=";",row.names=F)
          # 
          # #SDA
          # write.table(out.SDA,paste0("figures/",taxon.names,"/","overall/SDA.",path,".txt"),sep=";",row.names=F)
          # write.table(SDATot,paste0("figures/",taxon.names,"/","overall/SDAeach",path,".txt"),sep=";",row.names=F)
          # 
          # ### Ploting ###
          # colors <- c("#EDEDED",colorRampPalette(c("khaki2","orange","red","black"))(max(values(out.caTot),na.rm=T)))
          # breakpoints <- -0.5:(max(values(out.caTot),na.rm=T)+0.5)
          # a.arg <- list(at=c(0,max(values(out.caTot),na.rm=T)), labels=as.character(c(0,max(values(out.caTot),na.rm=T))),cex.axis=2.5)
          # l.arg <- list(text="Number of species (Current)",side=2, line=0.5, cex=3)
          # 
          # #Current
          # png(paste0("figures/",taxon.names,"/","overall/richness_current",def[m],"_",dispersal[d],".png"),width=650,height=1000)
          # plot(out.caTot,col=colors,breaks=breakpoints,ext=ext,
          #      legend.width=2,legend.shrink=0.9,legend.mar=7,
          #      axis.args=a.arg,legend.arg=l.arg,main="Currently Predicted Richness",cex.main=2,
          #      axes=FALSE, box=FALSE, zlim=c(0,max(values(out.caTot),na.rm=T)))
          # dev.off()
          # 
          # #Future
          # l.arg <- list(text="Number of species (Future)",side=2, line=0.5, cex=3)
          # png(paste0("figures/",taxon.names,"/","overall/richness.",path,".png"),width=650,height=1000)
          # plot(out.caFutTot,col=colors,breaks=breakpoints,ext=ext,
          #      legend.width=2,legend.shrink=0.9,legend.mar=7,cex.main=2,
          #      axis.args=a.arg,legend.arg=l.arg,main=paste0("Predicted Richness for the year ",fut.var[[3]][l],"\n(rcp ",fut.var[[2]][j],", dispersal ",dispersal[d],", deforestation: ",deforested[m],")"),
          #      axes=FALSE, box=FALSE, zlim=c(0,max(values(out.caTot),na.rm=T)))
          # dev.off()
          # 
          # #Turnover
          # colors <- rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(100))
          # l.arg <- list(text="Species Turnover (%)",side=2, line=0.5, cex=2.5)
          # a.arg <- list(cex.axis=3)
          # png(paste0("figures/",taxon.names,"/","overall/turnover.",path,".png"),width=650,height=1000)
          # plot(out.turnover*100,col=colors,ext=ext,cex.main=1.5,
          #      legend.width=2,legend.shrink=0.9,legend.mar=7,cex.axis=3,
          #      main=paste0("Species Turn-Over*100 for the year ",fut.var[[3]][l],"\n(rcp ",fut.var[[2]][j],", dispersal ",dispersal[d],", deforestation: ",deforested[m],")"),
          #      axes=FALSE, box=FALSE,legend.arg=l.arg,axis.args=a.arg)
          # dev.off()
          # 
          # #Turnover Filtered
          # colors <- c("#EDEDED",colorRampPalette(c("khaki2","orange","red","black"))(100))
          # l.arg <- list(text="Species Turn-Over * 100",side=2, line=0.5, cex=1.3)
          # png(paste0("figures/",taxon.names,"/","overall/turnover.filtered.",path,".png"),width=650,height=1000)
          # plot(out.turnoverFiltered*100,col=colors,ext=ext,cex.main=1.5,
          #      main=paste0("Species Turn-Over*100 with more than 5 species for the year ",fut.var[[3]][l],"\n(rcp ",fut.var[[2]][j],", dispersal ",dispersal[d],", deforestation: ",deforested[m],")"),
          #      legend.width=1.5,legend.shrink=0.6,legend.mar=7,
          #      axes=FALSE, box=FALSE,legend.arg=l.arg)
          # dev.off()
          # 
          #Species lost
          # colors <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(max(values(out.perteTot),na.rm=T)+1))
          # breakpoints <- c(-0.5,0.5:max(values(out.perteTot)+1,na.rm=T))
          # a.arg <- list(at=c(0,max(values(out.perteTot),na.rm=T)), labels=as.character(c(0,max(values(out.perteTot),na.rm=T))),cex.axis=2.5)
          # l.arg <- list(text="Number of species lost",side=2, line=0.5, cex=3)
          # png(paste0("figures/",taxon.names,"/","overall/species.lost.",path,".png"),width=650,height=1000)
          # plot(out.perteTot,col=colors,breaks=breakpoints,cex.main=2.5,
          #      main=paste0("Species lost from present to ",fut.var[[3]][l],"\n(rcp ",fut.var[[2]][j],", dispersal ",dispersal[d],", deforestation: ",deforested[m],")"),
          #      legend.width=2,legend.shrink=0.9,legend.mar=7,
          #      axis.args=a.arg,legend.arg=l.arg,
          #      axes=FALSE, box=FALSE)
          # dev.off()
          # 
          # Irr futur
          # colors <- viridis(9)
          # breakpoints <- quantile(out.irrFut[out.irrFut[]!=0],c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),na.rm=T)
          # a.arg <- list(seq(10,100,10),labels=as.character(seq(10,100,10)))
          # png(paste0("figures/",taxon.names,"/","overall/irr.fut.",path,".png"),width=650,height=1000)
          # plot(out.irrFut,col=colors,breaks=breakpoints,legend=F,
          #      axis.args=a.arg,
          #      axes=FALSE, box=FALSE)
          # dev.off()
          #
          # #Species gain
          # colors <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(max(values(out.gainTot),na.rm=T)))
          # breakpoints <- c(-0.5,0.5:max(values(out.gainTot)+1,na.rm=T))
          # a.arg <- list(at=c(0,max(values(out.gainTot),na.rm=T)), labels=as.character(c(0,max(values(out.gainTot),na.rm=T))),cex.axis=2.5)
          # l.arg <- list(text="Number of species gained",side=2, line=0.5, cex=3)
          # png(paste0("figures/",taxon.names,"/","overall/species.gained.",path,".png"),width=650,height=1000)
          # plot(out.gainTot,cex.axis=3,col=colors,breaks=breakpoints,ext=ext,cex.main=2.5,
          #      legend.width=2,legend.shrink=0.9,legend.mar=7,
          #      axis.args=a.arg,legend.arg=l.arg,
          #      axes=FALSE, box=FALSE)
          # dev.off()
          # 
          # #Irremp
          # colors <- c("#EDEDED","#057400")
          # breakpoints <- c(-0.01,0,0.01)
          # a.arg <- list(at=c(0), labels=as.character(c(0)),cex.axis=2.5)
          # l.arg <- list(text="Delta Irremplaceability (Future-Current)",side=2, line=0.5, cex=3)
          # png(paste0("figures/",taxon.names,"/","overall/delta.irr.",path,".png"),width=650,height=1000)
          # plot(out.irrDiff,col=colors,breaks=breakpoints,ext=ext,
          #      legend.width=2,legend.shrink=0.9,legend.mar=7,
          #      axis.args=a.arg,legend.arg=l.arg,
          #      axes=FALSE, box=FALSE)
          # dev.off()
          # 
          # #Changes areas
          # out.changes <- ratify(out.changes)      #Adding RAT table for qualitative variable
          # ratTable <- c("no event","Complete Exctinction","Persistance and Extinction","Persistance","Persistance and Newly Colonized","Newly Colonized")[levels(out.changes)[[1]]$ID+1]
          # levels(out.changes)[[1]]$changes <- ratTable
          # color = c("#EDEDED","#B01717","#DC650C","#DFD108","#80E204","#12910B")[c(levels(out.changes)[[1]]$ID+1)]
          # png(paste0("figures/",taxon.names,"/","overall/changes.",path,".png"),width=1300,height=2000)
          # print(levelplot(out.changes, att='changes', col.regions=color,cex=2.5,maxpixels = ncell(out.changes),main=paste0("Biodiversity patterns from present to",fut.var[[3]][l],"\n(rcp ",fut.var[[2]][j],", dispersal ",dispersal[d],", deforestation: ",deforested[m],")")))
          # dev.off()
          # 
          # #Range shift
          # png(paste0("figures/",taxon.names,"/","overall/rangeShift.",path,".png"),width=1300,height=2000)
          # plot(environ$alt,main=paste0('Shift from present to ',fut.var[[3]][l],"\n(rcp ",fut.var[[2]][j],", dispersal ",dispersal[d],", deforestation: ",deforested[m],")"),cex.main=2.5)
          # apply(as.matrix(rangeShift), MARGIN=1, FUN=function(x){Arrows(x0=as.numeric(x[2]),y0=as.numeric(x[3]),
          #                                                             x1=as.numeric(x[4]),y1=as.numeric(x[5]),lwd=0.1,arr.adj = 1,arr.type="simple")})
          # dev.off()
          # 
          # #Refuges
          # png(paste0("figures/",taxon.names,"/","overall/refuges2.",path,".png"),width=650,height=1000)
          # plot(out.refuges,axes=FALSE, box=FALSE)
          # plot(paMada,add=T)
          # dev.off()
        }
      }
    }
  }
}

  get.Model.Perf <- function(Perf.mods){
  names(Perf.mods) <- c("wIndex","Index","Model","Run","PA","Value")
  Perf.mods.test <- Perf.mods[Perf.mods$Index=="Testing.data"&Perf.mods$wIndex!="KAPPA",c(1,3,4,6)]
  Perf.mods.TSS <- reshape(Perf.mods.test[Perf.mods.test$wIndex=="TSS",-1], 
                           timevar="Run", idvar="Model", direction="wide")
  Perf.mods.ROC <- reshape(Perf.mods.test[Perf.mods.test$wIndex=="ROC",-1], 
                           timevar="Run", idvar="Model", direction="wide")
  return(list(Perf.mods.TSS,Perf.mods.ROC))
}

compute.Model.Perf.Tot <- function(Perf.mods.Tot){
  out <- data.frame(Model=character(),MinCalib=numeric(),MaxCalib=numeric(),MedianCalib=numeric(),
                       MinFull=numeric(),MaxFull=numeric(),MedianFull=numeric())
  out.Min <- aggregate(Perf.mods.Tot[,c(2,3)],list(Perf.mods.Tot$Model),FUN=min,na.rm = TRUE)
  out.Max <- aggregate(Perf.mods.Tot[,c(2,3)],list(Perf.mods.Tot$Model),FUN=max,na.rm = TRUE)
  out.Median <- aggregate(Perf.mods.Tot[,c(2,3)],list(Perf.mods.Tot$Model),FUN=median,na.rm = TRUE)
  out<- rbind(out,cbind(as.character(out.Min[,1]),out.Min[,2],out.Max[,2],
                        out.Median[,2],out.Min[,3],out.Max[,3],out.Median[,3]))
  colnames(out) <- c("Model","MinCalib","MaxCalib","MedianCalib","MinFull","MaxFull","MedianFull")
  return(out)
}

compute.sum.SDA <- function(SDA){
  SDATot <- data.frame(Defor=character(),RCP=numeric(),Year=numeric(),
                       Dispersal=character(),Extinct=character(),
                       DiminSup90=numeric(),DiminSUp50=numeric(),
                       DiminSup0=numeric(),IncreaSup0=numeric())
  SDATot <- rbind(SDATot,cbind(
    Defor=SDA$defor,RCP=SDA$rcp,Year=SDA$yr,Dispersal=SDA$disp,
    Extinct=length(which(SDA$area.fut==0)),
    DiminSup90=length(which(SDA$perc.chang< -90&SDA$area.fut!=0)),
    DiminSup50=length(which(SDA$perc.chang>= -90&SDA$perc.chang< -50)),
    DiminSup0=length(which(SDA$perc.chang>= -50&SDA$perc.chang< -0)),
    IncreaSup0=length(which(SDA$perc.chang>=0))))
  return(SDATot)
}

compute.areas.changes <- function(out.caTot,out.perteTot){
  #Selecting areas with more than n species
  n <- 0
  p5 <- out.perteTot ; p5[values(p5<n)] <- 0
  s5 <- out.stableTot ; s5[values(s5<n)] <- 0
  g5 <- out.gainTot ; g5[values(g5<n)] <- 0
  ps5 <- (p5!=0 & s5!=0)
  gs5 <- (g5!=0 & s5!=0)
  
  #Aggregating results
  noice <- raster(out.perteTot)
  values(noice)[!is.na(values(out.caTot))] <-0
  values(noice)[which(s5[]!=0)] <- 3
  values(noice)[which(g5[]!=0)] <- 5
  values(noice)[which(p5[]!=0&p5[]>g5[])] <- 1
  values(noice)[which(gs5[]!=0)] <- 4
  values(noice)[which(ps5[]!=0&p5[]>g5[])] <- 2
  
  return(noice)
}

vectorDirection <- function(vector){
  if(is.nan(vector[1])){
    return(NA)
  }
  angle <- atan2(vector[2],vector[1])
  direction <- c("W","SW","S","SE","E","NE","N","NW","W") #Since tan2 goes from ]-pi;pi], easy solution is to put "W" at both end
  
  increment <- (2*pi)/8
  initPi <- -pi-(increment/2)
  
  i <- 1
  while(angle>initPi){
    i <- i+1
    initPi <- initPi+(increment)
  }
  return(direction[i-1])
}

get.refuges.stats <- function(enough,list.current.niche,list.futur.niche,out.refuges,paMada){
  insidePres <- vector()
  insideFut <- vector()
  for(i in 1:length(enough)){
    if(enough[[i]]==T){
      ### Selecting current niche
      ca <- list.current.niche[[i]]
      
      ### Selecting futur niche
      caFut <- list.futur.niche[[m]][[j]][[l]][[i]]
      
      ca[ca==0]<-NA
      caFut[caFut==0]<-NA
      
      if(maxValue(ca==out.refuges)==1) insidePres <- c(insidePres,length(which((ca==out.refuges)[]==1))/length(which(ca[]==1)))
      else insidePres <- c(insidePres,F)
      if(is.na(maxValue(caFut==out.refuges))) insideFut <- c(insideFut,NA)
      else if(maxValue(caFut==out.refuges)==1) insideFut <- c(insideFut,length(which((caFut==out.refuges)[]==1))/length(which(caFut[]==1)))
      else insideFut <- c(insideFut,F)
    }
  }
  refuge <- out.refuges
  refuge[refuge!=1] <- NA 
  presenceCurrent <- out.caTot
  presenceCurrent[presenceCurrent>0]<-1
  refuge <- length(which(refuge[]==presenceCurrent[]))
  presenceCurrent <- length(which((presenceCurrent[]==1)))
  
  pourcTotHab <- refuge/presenceCurrent

  ex <- extract(out.refuges, paMada, fun=sum, na.rm=TRUE, df=TRUE)
  
  results <- data.frame(RefugesArea=refuge,PourcToTtHab=pourcTotHab,PourcInPA=sum(ex$layer)/refuge,nbSpProtected=length(which(insidePres!=0)))
  return(results)  
}

