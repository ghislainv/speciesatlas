############
#== Init ==#
############

library(reshape2)
library(ggplot2)
library(rasterVis)
library(shape)
library(ade4)
library(dendextend)
library(dplyr)
library(randomcoloR)

# setwd("C:/Users/etcla/Documents/ANALYSE/Atlas_Corentin/atlas-run/")

### Looping through :
fut.var <- list(c("cc","he","gs"),c("00","45","85"),c("2080"))
def <- c("_defo","")
deforested <- c("Yes","No")
dispersal <- c("full","zero") # In that order only


#################
#== Hist Loss ==#
#################

list.futur.niche<- list()

d <- 1
i <- 1
for (m in 1:length(def)) { #Deforestation or not
  for (j in 1:length(fut.var[[2]])) { #Rcp: 4.5 or 8.5
    for (l in 1:length(fut.var[[3]])) { #Time period: 2085
      if(!(def[m]!="_defo"&fut.var[[2]][j]=="00")){
        path=paste0("figures/lemurs/overall/tif/",
                    fut.var[[2]][[j]],"_",fut.var[[3]][[l]],def[[m]],"_",dispersal[d],".tif")
        print(path)
        list.futur.niche[[paste0(fut.var[[2]][[j]],"_",fut.var[[3]][[l]],def[[m]])]] <- as.data.frame(table(stack(path)[[4]][]))$Freq
        i <- i+1
      }
    }
  }
}

list.futur.niche.full.df <- rowr::cbind.fill(0:22,list.futur.niche[[1]],list.futur.niche[[2]]
                                             ,list.futur.niche[[3]],list.futur.niche[[4]],list.futur.niche[[5]],fill=0)
colnames(list.futur.niche.full.df) <- c("Nb Individus","No Climate + Deforest","RCP 4.5 Climate + Deforet","RCP8.5 Climate +Deforest",
                                        "RCP4.5 Climate","RCP8.5 Climate")

dfp1 <- melt(list.futur.niche.full.df,id.vars="Nb Individus")
colnames(dfp1) <- c("Species richness","2080 Scenarios","Surface (km²)")
ggplot(dfp1, aes(x = `Species richness`, y= `Surface (km²)`, fill = `2080 Scenarios`)) +
  geom_bar(stat="identity", position = "dodge",width=0.9) +
  scale_fill_brewer(palette="Set1") +
  ylab("Cumulative surface area (km²)")+
  xlab("Number of species lost")+
  labs(fill="Scenario for the year 2080")+
  theme_bw()+
 theme(legend.position = c(0.725,0.8),
       legend.background = element_rect(fill="white",size=0.5, linetype="solid",colour ="darkgrey"),
       legend.text = element_text(size=rel(1)),
       legend.title = element_text(size=rel(1.1),face="bold"),
       panel.grid.minor = element_blank(),
       panel.grid.major = element_blank(),
       axis.title = element_text(size=rel(1.2)),
       axis.text = element_text(size=rel(1.1))
       
       )

##################
#== SDA Change ==#
##################

#### PCA ####
# sda<-read.table("figures/lemurs/overall/SDA/sda.each.brute.txt",sep="\t",header=T,colClasses = c("factor",rep("numeric",11),rep("factor"),rep("numeric",2)),dec=",")
sda<-read.table("figures/lemurs/overall/SDA/sda.each.perc.txt",sep="\t",header=T,colClasses = c("factor",rep("numeric",11),rep("factor",0)))

# sda <- sda[-which(sda$species== "Eulemur.cinereiceps"
                  # |sda$class=="1" 
                  # |sda$class=="2"
                  # |sda$class=="3"
                  # ),
                  # -c(3,4,5,6)
                  # -c(7,8,9,10)
          # ]
z<- dudi.pca(df = log(sda[,-c(1,2,3)]+101), center = T, scale = T, scannf = FALSE, nf = 2)
# z<- dudi.pca(df = sda[,-c(1,2,7,8)], center = T, scale = T, scannf = FALSE, nf = 2)

s.corcircle(z$co,clabel=0.5)
s.class(z$li,sda[,1],cellipse=T,clabel=0.5)
s.class(z$li,sda[,length(sda)-1],cellipse=T,clabel=0.5,col=c("red","blue"))
s.class(z$li,sda[,length(sda)-2],cellipse=T,clabel=0.5,col=c("darkgreen","blue","#FFBE38","#FF6928","#FF1419","grey"))
s.class(z$li,as.factor(length(sda)),cellipse=T,clabel=0.5,col=distinctColorPalette(unique(length(sda[,length(sda)]))))

#### CAH ####
sda.cr <- scale(log(sda[,-c(1,2,3)]+101),center=T,scale=T)
rownames(sda.cr) <- paste0(" ",sda$species)
clust <- sda.cr %>% dist %>%
  hclust(method="ward.D")

inertie <- sort(clust$height, decreasing = TRUE)
plot(inertie[1:20], type = "s", xlab = "Number of classes", ylab = "Inertia",cex.lab=1.5,cex.axis=1.5,lwd=2)
points(c(4,5), inertie[c(4,5)], col = c("red3","blue3"), cex = 2, lwd = 3)

sda$cahLogScale4 <- as.character(cutree(clust,k=4))
sda$cahLogScale5 <- as.character(cutree(clust,k=5))
clust<-prune(clust,c(" Eulemur.cinereiceps"))
#### Â¨Plotting hclust
k=5
dend <- clust %>% as.dendrogram %>%
  set("branches_k_color",k=k) %>% set("branches_lwd", 0.7) %>%
  set("labels_cex", 0.6) %>% set("labels_colors",k=k) %>%
  set("leaves_pch", 19) %>% set("leaves_cex", 0.5) 
ggd1 <- as.ggdend(dend)

## Resize width and height plotting area
ggplot(ggd1, horiz = TRUE)

#### Prepare data ####
n <- length(sda$sp)
sda.new <- data.frame(species=rep(sda$sp,each=12),current=rep(sda$current.area,each=12),
                      rcp=rep(c(rep("0.0",4),rep("4.5",4),rep("8.5",4)),n),
                      dispersal=rep(rep(c("zero","full"),6),n),
                      deforestation=rep(rep(c("Yes","Yes","No","No"),3),n),pourcChange=NA)
sda.new <- sda.new[order(sda.new$species),]

sda.new$pourcChange[which(sda.new$rcp=="0.0"&sda.new$dispersal=="full"&sda.new$deforestation=="Yes")] <- sda$Cli_Def_00_full
sda.new$pourcChange[which(sda.new$rcp=="4.5"&sda.new$dispersal=="full"&sda.new$deforestation=="Yes")] <- sda$Cli_Def_4.5_full
sda.new$pourcChange[which(sda.new$rcp=="8.5"&sda.new$dispersal=="full"&sda.new$deforestation=="Yes")] <- sda$Cli_Def_8.5_full
sda.new$pourcChange[which(sda.new$rcp=="0.0"&sda.new$dispersal=="full"&sda.new$deforestation=="No")] <- sda$NoCli_NoDef
sda.new$pourcChange[which(sda.new$rcp=="4.5"&sda.new$dispersal=="full"&sda.new$deforestation=="No")] <- sda$Cli_4.5_full
sda.new$pourcChange[which(sda.new$rcp=="8.5"&sda.new$dispersal=="full"&sda.new$deforestation=="No")] <- sda$Cli_8.5_full

sda.new$pourcChange[which(sda.new$rcp=="0.0"&sda.new$dispersal=="zero"&sda.new$deforestation=="Yes")] <- sda$Cli_Def_00_zero
sda.new$pourcChange[which(sda.new$rcp=="4.5"&sda.new$dispersal=="zero"&sda.new$deforestation=="Yes")] <- sda$Cli_Def_4.5_zero
sda.new$pourcChange[which(sda.new$rcp=="8.5"&sda.new$dispersal=="zero"&sda.new$deforestation=="Yes")] <- sda$Cli_Def_8.5_zero
sda.new$pourcChange[which(sda.new$rcp=="0.0"&sda.new$dispersal=="zero"&sda.new$deforestation=="No")] <- sda$NoCli_NoDef
sda.new$pourcChange[which(sda.new$rcp=="4.5"&sda.new$dispersal=="zero"&sda.new$deforestation=="No")] <- sda$Cli_4.5_zero
sda.new$pourcChange[which(sda.new$rcp=="8.5"&sda.new$dispersal=="zero"&sda.new$deforestation=="No")] <- sda$Cli_8.5_zero

sda.new <- sda.new[-which(sda.new$rcp=="0.0"&sda.new$deforestation=="No"),]
rownames(sda.new)<-1:(n*10)

#### Modelisation ####
dotchart(sda.new$pourcChange)
hist(sda.new$pourcChange)

hist(log(sda.new$pourcChange+101))

# mod2 <- lm(log(pourcChange+101)~rcp+deforestation+dispersal+current,data=sda.new)
# summary(mod2);AIC(mod2)

mod3<-lme4::lmer(log(pourcChange+101)~rcp+deforestation+ (1|species),data=sda.new)
summary(mod3);AIC(mod3)
# 
# mod3<-gls(log(pourcChange+101)~rcp+deforestation+dispersal+current,data=sda.new) # generalized least square regression
# summary(mod3)    # equivalent to a GLM
# mod4<-lme(log(pourcChange+101)~rcp+deforestation+dispersal+current,
#           random=~1|species, method='REML', data=sda.new) # mixed-effects model with random intercept for nest ID
# # using REML method (with ML, you cannot perform the LRT below)
# summary(mod4)
# anova(mod3,mod4)
# 
# mod4<-lme(log(pourcChange+101)~rcp+deforestation+dispersal,
#           random=~1|species, method='REML', data=sda.new) # mixed-effects model with random intercept for nest ID
# summary(mod4)
# 
# hist(resid(mod4),breaks=c(-4.5:4.5))
# par(mfrow=c(2,2))

# plot(mod2)
ano <- anova(mod3)

ano$percSqTot <- (ano$`Sum Sq`/sum(ano$`Sum Sq`))*100
ano$Df <- c("RCP","Deforestation","Residuals")
# 
# mod2 <- lm(log(pourcChange+101)~deforestation+rcp+dispersal+current,data=sda.new)
# summary(mod2);AIC(mod2)

hist(resid(mod4),breaks=c(-4.5:4.5))
par(mfrow=c(2,2))
plot(mod3)


##################
#== SDA Change ==#
##################

library(spatstat)
library(data.table)
file="rangeShift.85_2080_defo_full"
rangeShift <- read.table(paste0("figures/",taxon.names,"/","overall/range shift/",file,".txt"),sep=";",header=T,
                         colClasses=c("character",rep("numeric",5),"character"))
rangeShift <- na.omit(rangeShift)
rangeShift$id <- rownames(rangeShift)
listSP <- unique(rangeShift$species)

for(i in 1:length(listSP)){
  rangeShift <- rangeShift[rangeShift$species==listSP[i],]
  
  ####Method 1
  dt <- as.data.table(rangeShift[,c(1,2,3,4,5,8)])
  
  sf <- dt[
    , {
      geometry <- sf::st_linestring(x = matrix(c(x1,x2, y1, y2), ncol = 2))
      geometry <- sf::st_sfc(geometry)
      geometry <- sf::st_sf(geometry = geometry,crs=as.character(crs(s)))
    }
    , by = id
    ]
  
  test <- sf::st_as_sf(sf)
  sf::st_write(test,paste0("figures/",taxon.names,"/","overall/range shift/shp/",file,".shp"))
  
}

  Sldf <- rgdal::readOGR(paste0("figures/",taxon.names,"/","overall/range shift/shp/rangeShift.85_2080_defo_full.shp"))

#### METHOD 2
rangeShift$Line <- apply(rangeShift[,],1, FUN = function (x){sp::Line(cbind(c(as.numeric(x[2]),as.numeric(x[3])),c(as.numeric(x[4]),as.numeric(x[5]))))})
rangeShift.Lines <- vector("list")
for(i in 1:length(rangeShift$Line)){
  rangeShift.Lines[[i]] <- sp::Lines(rangeShift$Line[i], ID = rangeShift$species[i])
}

Sl <- SpatialLines(rangeShift.Lines,proj4string = crs(s$temp))
rownames(rangeShift) <- rangeShift$species
Sldf <- SpatialLinesDataFrame(Sl, data = rangeShift[,c(1,6,7)])

#####AFTER

full <- raster("E:/actual_results/rangeShift/full.tif")
# colors <- rev(colorRampPalette(brewer.pal(11, "YlGnBu"))(max(values(full),na.rm=T)+1))
colors <- viridis(max(values(full),na.rm=T)+1)
# # breakpoints <- c(-0.5,0.5:max(values(out.perteTot)+1,na.rm=T))
a.arg <- list(at=quantile(seq(0,max(full[],na.rm=T),1),c(0.1,0.9)), labels=c("Low","High"),cex.axis=4)
l.arg <- list(text="Line density",side=2, line=1, cex=4)
png(paste0("E:/actual_results/rangeShift/full.png"),width=1300,height=2000)
plot(full,col=colors,cex.main=2.5,
     legend.width=2,legend.shrink=0.9,legend.mar=13,
     axis.args=a.arg,legend.arg=l.arg,
     axes=FALSE, box=FALSE)
dev.off()

arr <- read.table(paste0("figures/",taxon.names,"/","overall/range shift/rangeShift.85_2080_defo_full.txt"),header=T,sep=";",colClasses=c("character",rep("numeric",5),"character"))
colors <- colorRampPalette(brewer.pal(11, "YlGn"))(max(values(s$foret),na.rm=T)+1)
a.arg <- list(at=quantile(seq(0,max(s$foret[],na.rm=T),1),c(0.1,0.9)), labels=c("10","90"),cex.axis=4)
l.arg <- list(text="Forest cover (%)",side=2, line=1, cex=4)
png(paste0("E:/actual_results/rangeShift/full_vector.png"),width=1300,height=2000)
plot(s$foret,col=colors,cex.main=2.5,
     legend.width=2,legend.shrink=0.9,legend.mar=10,
     axis.args=a.arg,legend.arg=l.arg,
     axes=FALSE, box=FALSE)
apply(as.matrix(arr), MARGIN=1, FUN=function(x){Arrows(x0=as.numeric(x[2]),y0=as.numeric(x[3]),
                                                              x1=as.numeric(x[4]),y1=as.numeric(x[5]),col="red",lwd=2,arr.adj = 1,arr.type="simple")})
dev.off()

#############################
#==Range Shift Aggregatedf==#
#############################

### For every species
df.rangeShift <- data.frame(dispersal=factor(),defor=factor(),rcp=factor(),meanDist=numeric(),sdDist=numeric(),N=character(),NE=character(),E=character(),
                            SE=character(),S=character(),SW=character(),W=character(),NW=character(),Nsp=numeric())
first = T
deforested <- c("Deforested","NotDeforested")
for(d in 1:length(dispersal)){
  for (m in 1:length(def)) { #Deforestation or not
    for (j in 1:length(fut.var[[2]])) { #Rcp: 4.5 or 8.5
      for (l in 1:length(fut.var[[3]])) { #Time period: 2085
        if(!(def[m]!="_defo"&fut.var[[2]][j]=="00")){
          path=paste0("figures/lemurs/overall/range shift/rangeShift.",
                      fut.var[[2]][[j]],"_",fut.var[[3]][[l]],def[[m]],"_",dispersal[d],".txt")
          print(path)
          arr <- read.table(path,header=T,sep=";",colClasses=c("character",rep("numeric",5),"character"))
          arr<-na.omit(arr)
          if(first==T){
            df.rangeShift <-data.frame(dispersal=dispersal[d],defor=deforested[[m]],rcp=fut.var[[2]][[j]],
                                       meanDist=mean(arr$distance),sdDist=sd(arr$distance),N=length(which(arr$direction=="N")),NE=length(which(arr$direction=="NE")),
                                       E=length(which(arr$direction=="E")),SE=length(which(arr$direction=="SE")),S=length(which(arr$direction=="S")),
                                       SW=length(which(arr$direction=="SW")),W=length(which(arr$direction=="W")),
                                       NW=length(which(arr$direction=="NW")),Nsp=dim(arr)[1])
            first=F
          }else{
            df.rangeShift <-rbind.data.frame(df.rangeShift,data.frame(dispersal=dispersal[d],defor=deforested[[m]],rcp=fut.var[[2]][[j]],
                                                                      meanDist=mean(arr$distance),sdDist=sd(arr$distance),N=length(which(arr$direction=="N")),NE=length(which(arr$direction=="NE")),
                                                                      E=length(which(arr$direction=="E")),SE=length(which(arr$direction=="SE")),S=length(which(arr$direction=="S")),
                                                                      SW=length(which(arr$direction=="SW")),W=length(which(arr$direction=="W")),
                                                                      NW=length(which(arr$direction=="NW")),Nsp=dim(arr)[1]))
          }
        }
      }
    }
  }
}
write.table(df.rangeShift,"figures/lemurs/overall/range shift/aggregated.full.txt",row.names=F)

###For specific cases
path=paste0("~/ANALYSE/Atlas_Corentin/atlas-run/figures/lemurs/overall/range shift/",
            c("deforNo","deforYes"),".txt")
print(path)
first=T
for(i in 1:length(path)){
  arr <- read.table(path[i],header=T,sep=";",colClasses=c("character",rep("numeric",5),"character"))
  arr<-na.omit(arr)
  if(first==T){
    df.rangeShift <-data.frame(file=path[i],meanDist=mean(arr$distance),sdDist=sd(arr$distance),N=length(which(arr$direction=="N")),NE=length(which(arr$direction=="NE")),
                               E=length(which(arr$direction=="E")),SE=length(which(arr$direction=="SE")),S=length(which(arr$direction=="S")),
                               SW=length(which(arr$direction=="SW")),W=length(which(arr$direction=="W")),
                               NW=length(which(arr$direction=="NW")),Nsp=dim(arr)[1])
    first=F
  }else{
    df.rangeShift <-rbind.data.frame(df.rangeShift,data.frame(file=path[i],meanDist=mean(arr$distance),sdDist=sd(arr$distance),N=length(which(arr$direction=="N")),NE=length(which(arr$direction=="NE")),
                                                              E=length(which(arr$direction=="E")),SE=length(which(arr$direction=="SE")),S=length(which(arr$direction=="S")),
                                                              SW=length(which(arr$direction=="SW")),W=length(which(arr$direction=="W")),
                                                              NW=length(which(arr$direction=="NW")),Nsp=dim(arr)[1]))
  }
}
write.table(df.rangeShift,"~/ANALYSE/Atlas_Corentin/atlas-run/figures/lemurs/overall/range shift/aggregated.defor.txt",row.names=F)

###For every species
first = T
deforested <- c("Deforested","NotDeforested")
for(d in 1:length(dispersal)){
  for (m in 1:length(def)) { #Deforestation or not
    for (j in 1:length(fut.var[[2]])) { #Rcp: 4.5 or 8.5
      for (l in 1:length(fut.var[[3]])) { #Time period: 2085
        if(!(def[m]!="_defo"&fut.var[[2]][j]=="00")){
          path=paste0("figures/lemurs/overall/range shift/rangeShift.",
                      fut.var[[2]][[j]],"_",fut.var[[3]][[l]],def[[m]],"_",dispersal[d],".txt")
          print(path)
          arr <- read.table(path,header=T,sep=";",colClasses=c("character",rep("numeric",5),"character"))
          if(first==T){
            df.rangeShift <-cbind.data.frame(expand.grid(dispersal=dispersal[d],defor=deforested[[m]],rcp=fut.var[[2]][[j]],species=arr$species),
                                             arr[,-1])
            first=F
          }else{
            df.rangeShift <-rbind.data.frame(df.rangeShift,cbind.data.frame(expand.grid(dispersal=dispersal[d],defor=deforested[[m]],rcp=fut.var[[2]][[j]],species=arr$species),
                                                                            arr[,-1]))
          }
        }
      }
    }
  }
}
df.rangeShift <- na.omit(df.rangeShift)
dotchart(df.rangeShift$distance)

dotchart(log(df.rangeShift$distance))

dotchart(sqrt(df.rangeShift$distance))
df.rangeShift$defor <- ordered(df.rangeShift$defor, levels = c("NotDeforested", "Deforested"))

mod1<-lm(sqrt(distance)~defor+rcp,data=df.rangeShift)
summary(mod1);AIC(mod1)
hist(resid(mod1))
par(mfrow=c(2,2))
plot(mod1)

mod3<-gls(sqrt(distance)~rcp+defor,data=df.rangeShift) # generalized least square regression
summary(mod3)    # equivalent to a GLM
mod4<-lme(sqrt(distance)~rcp+defor,
           random=~1|species, method='REML', data=df.rangeShift) # mixed-effects model with random intercept for nest ID
summary(mod4)

mod3<-lme4::lmer(sqrt(distance)~rcp+defor+ (1|species),data=df.rangeShift)
summary(mod3);AIC(mod3)

ano <- anova(mod3)
ano$percSqTot <- (ano$`Sum Sq`/sum(ano$`Sum Sq`))*100

r.squaredGLMM(mod3)

plot(distance~defor,data=df.rangeShift)
write.table(summary(mod1)$coefficient,"~/ANALYSE/Atlas_Corentin/atlas-run/figures/lemurs/overall/range shift/model.txt")
#############
#==Sum tif==#
#############
filesall <- dir("figures/lemurs/overall/range shift/tif_species_masked")
files <- filesall[grepl(".tif",files)]
sumRaster <- raster("C:/Users/etcla/Desktop/test_masked.tif")
sumRaster <- resample(sumRaster,s)
sumRaster[!is.na(s$temp[])&is.na(sumRaster[])]<-0
sumRaster[is.na(s$temp[])]<-NA
for(i in 1:length(files)){
  print(files[i])
  r <- raster(paste0("figures/lemurs/overall/range shift/tif_species_masked/",files[i]))
  print(max(r[],na.rm=T))
  # if(yres(r)!=1000|xres(r)!=1000){
  #   print("resampling")
  #   r <- resample(r,s)
  # }
  # r[is.na(r[])&!is.na(s$temp[])] <- 0
  # sumRaster[] <- sumRaster[]-r[]
  # plot(stack(r,sumRaster))
}

writeRaster(sumRaster,paste0("figures/lemurs/overall/range shift/lineDensities.tif"))

sumBis <- sumRaster
sumBis[sumRaster[]<423] <- 0
plot(sumBis)

########################
#==Richness explained==#
########################
list.richness <- data.frame(richness=numeric(),dispersal=factor(),defor=factor(),rcp=factor(),
                            temp=numeric(),prec=numeric(),tseas=numeric(),pseas=numeric(),cwd=numeric(),foret=numeric())
first = T
deforested <- c("Deforested","NotDeforested")
for(d in 1:length(dispersal)){
  for (m in 1:length(def)) { #Deforestation or not
    for (j in 1:length(fut.var[[2]])) { #Rcp: 4.5 or 8.5
      for (l in 1:length(fut.var[[3]])) { #Time period: 2085
        if(!(def[m]!="_defo"&fut.var[[2]][j]=="00")){
          path=paste0("figures/lemurs/overall/tif/",
                      fut.var[[2]][[j]],"_",fut.var[[3]][[l]],def[[m]],"_",dispersal[d],".tif")
          print(path)
          if(first==T){
            list.richness <-data.frame(richness=stack(path)[[1]][],dispersal=dispersal[d],defor=deforested[[m]],rcp=fut.var[[2]][[j]],
                                       temp=as.numeric(s$temp[]),prec=as.numeric(s$prec[]),tseas=as.numeric(s$tseas[]),pseas=as.numeric(s$pseas[]),
                                       cwd=as.numeric(s$cwd[]),foret=as.numeric(s$foret[]))
            first=F
          }else{
            list.richness <-rbind.data.frame(list.richness,cbind.data.frame(richness=stack(path)[[1]][],dispersal=dispersal[d],defor=deforested[[m]],rcp=fut.var[[2]][[j]],
                                                                            temp=as.numeric(s$temp[]),prec=as.numeric(s$prec[]),tseas=as.numeric(s$tseas[]),pseas=as.numeric(s$pseas[]),
                                                                            cwd=as.numeric(s$cwd[]),foret=as.numeric(s$foret[])))
          } 
        }
      }
    }
  }
}
str(list.richness)
list.richness <- na.omit(list.richness)
list.richness$defor <- ordered(list.richness$defor, levels = c("NotDeforested", "Deforested"))
# samples.list.richness <- list.richness[sample(dim(list.richness)[1],1000000),]
mod1 <- lm(log(1+richness)~rcp+defor, data=list.richness)
summary(mod1);AIC(mod1)
aa<-anova(mod1)
aa$perVar <- aa$`Sum Sq`/sum(aa$`Sum Sq`)
hist(resid(mod1))

var(log(1+list.richness$richness))/mean(log(1+list.richness$richness))
mod2 <- glm(log(1+richness)~rcp+defor,family="poisson",data=list.richness)
mod2<-glm.nb(richness~rcp+defor, data=list.richness[list.richness$richness>5,])
summary(mod2)

write.table(summary(mod1)$coefficient,"~/ANALYSE/Atlas_Corentin/atlas-run/figures/lemurs/overall/richness/model.txt")
write.table(aa,"~/ANALYSE/Atlas_Corentin/atlas-run/figures/lemurs/overall/richness/anova.txt")

########################
#==Turnover explained==#
########################
list.turnover <- data.frame(richness=numeric(),dispersal=factor(),defor=factor(),rcp=factor(),
                            temp=numeric(),prec=numeric(),tseas=numeric(),pseas=numeric(),cwd=numeric(),foret=numeric())
first = T
deforested <- c("Deforested","NotDeforested")
for(d in 1:length(dispersal)){
  for (m in 1:length(def)) { #Deforestation or not
    for (j in 1:length(fut.var[[2]])) { #Rcp: 4.5 or 8.5
      for (l in 1:length(fut.var[[3]])) { #Time period: 2085
        if(!(def[m]!="_defo"&fut.var[[2]][j]=="00")){
          path=paste0("figures/lemurs/overall/tif/",
                      fut.var[[2]][[j]],"_",fut.var[[3]][[l]],def[[m]],"_",dispersal[d],".tif")
          print(path)
          if(first==T){
            list.turnover <-data.frame(turnover=stack(path)[[6]][],dispersal=dispersal[d],defor=deforested[[m]],rcp=fut.var[[2]][[j]],
                                       temp=as.numeric(s$temp[]),prec=as.numeric(s$prec[]),tseas=as.numeric(s$tseas[]),pseas=as.numeric(s$pseas[]),
                                       cwd=as.numeric(s$cwd[]),foret=as.numeric(s$foret[]))
            first=F
          }else{
            list.turnover <-rbind.data.frame(list.richness,cbind.data.frame(richness=stack(path)[[6]][],dispersal=dispersal[d],defor=deforested[[m]],rcp=fut.var[[2]][[j]],
                                                                            temp=as.numeric(s$temp[]),prec=as.numeric(s$prec[]),tseas=as.numeric(s$tseas[]),pseas=as.numeric(s$pseas[]),
                                                                            cwd=as.numeric(s$cwd[]),foret=as.numeric(s$foret[])))
          } 
        }
      }
    }
  }
}
str(list.turnover)
list.turnover <- na.omit(list.turnover)
list.turnover$defor <- ordered(list.turnover$defor, levels = c("NotDeforested", "Deforested"))
# samples.list.richness <- list.richness[sample(dim(list.richness)[1],1000000),]
mod1 <- lm(turnover~rcp+defor, data=list.turnover)
summary(mod1);AIC(mod1)
aa<-anova(mod1)
aa$perVar <- aa$`Sum Sq`/sum(aa$`Sum Sq`)
hist(resid(mod1))

var(list.turnover$turnover)/mean(list.turnover$turnover)
mod2<-glm.nb(turnover~rcp+defor, data=list.turnover)
summary(mod2)

write.table(summary(mod1)$coefficient,"~/ANALYSE/Atlas_Corentin/atlas-run/figures/lemurs/overall/richness/model.txt")
write.table(aa,"~/ANALYSE/Atlas_Corentin/atlas-run/figures/lemurs/overall/richness/anova.txt")

############################
#== Deforestation effect ==#
############################

for(d in 1:length(dispersal)){
    for (j in 1:length(fut.var[[2]])) { #Rcp: 4.5 or 8.5
      for (l in 1:length(fut.var[[3]])) { #Time period: 2085
        if(!(fut.var[[2]][j]=="00")){
          path=paste0("figures/lemurs/overall/tif/",
                      fut.var[[2]][[j]],"_",fut.var[[3]][[l]],c("_defo",""),"_",dispersal[d],".tif")
          print(path)
          refugeDef <- stack(path[1])[[8]]
          refugeNoDef <- stack(path[2])[[8]]
          aresOfInterest <- raster(refugeDef)
          aresOfInterest[which(refugeDef[]==1&refugeNoDef[]==1)] <- 3 #Refuges & Nothing
          aresOfInterest[which(refugeDef[]==1&refugeNoDef[]==0)] <-1 #Increasing interest
          aresOfInterest[which(refugeDef[]==0&refugeNoDef[]==1)] <-2 #Higly threatened
          aresOfInterest[which(refugeNoDef[]==0&refugeDef[]==0)] <-1 #No interest
          
          print(table(aresOfInterest[])[3]/(table(aresOfInterest[])[3]+table(aresOfInterest[])[2]))
          # writeRaster(aresOfInterest,paste0("figures/lemurs/overall/tif/refuges.defor.zero.tiff"),overwrite=T)
          
          ###PLotting
          # aresOfInterest <- ratify(aresOfInterest)      #Adding RAT table for qualitative variable
          # ratTable <- c("","Highly Threatened Refuges","Refuges")[levels(aresOfInterest)[[1]]$ID]
          # levels(aresOfInterest)[[1]]$conservation <- ratTable
          # color = c("#EDEDED","#B01717","#80E204")[c(levels(aresOfInterest)[[1]]$ID)]#FFD259
          # png(paste0("figures/lemurs/","overall/defor.refuges",fut.var[[2]][[j]],"_",fut.var[[3]][[l]],"_",dispersal[d],".png"),width=1300,height=2000)
          # print(levelplot(aresOfInterest, att='conservation', col.regions=color,cex.legends=2.5,
          #                 maxpixels = ncell(aresOfInterest), par.settings = list(axis.line = list(col = "transparent")), 
          #                 scales = list(draw=FALSE)))
          # dev.off()
        }
      }
  }
}

############################
#==  RICHNESS EVOLUTION  ==#
############################
list.richness <- list()

for(d in 1:length(dispersal)){
  for(m in 1:length(def)){
    for (j in 1:length(fut.var[[2]])) { #Rcp: 4.5 or 8.5
      for (l in 1:length(fut.var[[3]])) { #Time period: 2085
        path=paste0("figures/lemurs/overall/tif/",
                    fut.var[[2]][[j]],"_",fut.var[[3]][[l]],def[m],"_",dispersal[d],".tif")
        print(path)
        refugeDef <- stack(path)
        names(refugeDef) <- c("out.caFutTot","out.caTot","out.gainTot","out.perteTot","out.stableTot","out.turnover","out.changes","out.refuges")
        list.richness[path] <-refugeDef
      }
    }
  }
}
  
############################
#== Is turn(over)Â² rated ?=#
############################

tif <- stack("figures/lemurs/overall/tif/85_2080_zero.tif")
turnover <- stack("figures/lemurs/overall/tif/85_2080_zero.tif")[[6]]
perclost <- raster(turnover)
perclost[] <- (stack("figures/lemurs/overall/tif/85_2080_zero.tif")[[4]][]/stack("figures/lemurs/overall/tif/85_2080_zero.tif")[[2]][])

test <- out.turnover
test[out.irrFut[]<=3.260742e-04]<-NA
out.irrFut[]>3.260742e-04&out.turnover[]<0.5
