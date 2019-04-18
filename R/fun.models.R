#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent, Aurelien Colas
# email           :ghislain.vieilledent@cirad.fr, aurelien.colas@insa-lyon.fr
# license         :GPLv3
# ==============================================================================

# ==================
# Modelling
# ==================

fun.models.run <- function(name,spdir,p,s,spname,model.var,future,fut.var,maxent.path){
  
  if (length(p)<30){
    models_used <- c("GAM","RF","MAXENT.Phillips","ANN")
  }else{
    models_used <- c("GLM","GAM","RF","MAXENT.Phillips","ANN")
  }
  
  ## BIOMOD_FormatingData
  set.seed(1234) ## Reproducible pseudo-absences
  BiomodData <- BIOMOD_FormatingData(resp.var=p,
                                     expl.var=s,
                                     resp.name=spname,
                                     PA.nb.rep=1,
                                     PA.nb.absences=10000,
                                     PA.strategy="random",
                                     na.rm=TRUE)
  saveRDS(BiomodData,paste0(spdir,"/BiomodData.rds"))

  ## BIOMOD_ModelingOptions
  BiomodOptions <- BIOMOD_ModelingOptions(GLM=list(type="quadratic",interaction.level=0,myFormula=NULL,family=binomial(link="logit"),
                                                    test="AIC",  control = glm.control(epsilon = 1e-08, maxit = 100, trace = F)),
                                          GAM=list(algo="GAM_mgcv",type="s_smoother",k=4,interaction.level=0,
                                                   myFormula=NULL,
                                                   family=binomial(link="logit")),
                                          RF=list(do.classif=TRUE, ntree=500),
                                          MAXENT.Phillips=list(path_to_maxent.jar=maxent.path,
                                                               visible=FALSE,maximumiterations=500,
                                                               memory_allocated=512,
                                                               # To avoid overparametrization (Merow  et al.  2013)
                                                               product=FALSE, threshold=FALSE, hinge=FALSE),
                                          ANN=list(NbCV=5,size=NULL,decay=NULL,rang=0.1,maxit=200))
  saveRDS(BiomodOptions,paste0(spdir,"/BiomodOptions.rds"))

  ## BIOMOD_Modeling
  set.seed(1234) ## Reproducible results
  BiomodModel <- BIOMOD_Modeling(BiomodData,
                                 models=models_used,
                                 models.options=BiomodOptions,
                                 NbRunEval=1,
                                 DataSplit=70,
                                 VarImport=3,
                                 models.eval.meth=c("KAPPA","TSS","ROC"),
                                 rescal.all.models=TRUE,
                                 do.full.models=TRUE,
                                 modeling.id="5mod") ## 5 statistical models
  
  # ## Removing models not converging 
  # good_models <- models_used
  # warn = grep("not converge|fit", capture.output(warnings()))
  # if (length(warn)>0){
  #   if (grepl("glm", warn)==T){good_models<-good_models[!grepl("GLM", good_models)]}
  #   if (grepl("gam", warn)==T){good_models<-good_models[!grepl("GAM", good_models)]}
  # }
  # 
  # ## If some models didn't worked we try again without them
  # if (length(good_models)!= length(models_used)){
  #   BiomodModel <- BIOMOD_Modeling(BiomodData,
  #                                  models=good_models,
  #                                  models.options=BiomodOptions,
  #                                  NbRunEval=1,
  #                                  DataSplit=70,
  #                                  VarImport=3,
  #                                  models.eval.meth=c("KAPPA","TSS","ROC"),
  #                                  rescal.all.models=TRUE,
  #                                  do.full.models=TRUE,
  #                                  modeling.id="5mod") ## 5 statistical models
  # }
  
  ## Building ensemble-models
  BiomodEM <- BIOMOD_EnsembleModeling(modeling.output=BiomodModel,
                                      chosen.models= grep("_Full_",
                                                          get_built_models(BiomodModel),
                                                          value=TRUE),
                                      em.by="all",
                                      eval.metric=c("TSS"),
                                      eval.metric.quality.threshold=0.5,
                                      prob.mean=FALSE,
                                      prob.ci=FALSE,
                                      prob.ci.alpha=0.05,
                                      committee.averaging=TRUE,
                                      prob.mean.weight=TRUE,
                                      prob.mean.weight.decay="proportional")

  ## BIOMOD_Projection == PRESENT == ## Individual model projection
  BiomodProj <- BIOMOD_Projection(modeling.output=BiomodModel,
                                  new.env=s,
                                  proj.name="current",
                                  selected.models=grep("_Full_",
                                                       get_built_models(BiomodModel),
                                                       value=TRUE), ## Full models only
                                  binary.meth=c("TSS"),
                                  filtered.meth=c("TSS"),
                                  compress=TRUE,
                                  build.clamping.mask=FALSE, ## No need for present
                                  omi.na=TRUE,
                                  on_0_1000=TRUE,
                                  output.format=".grd")

  ## BIOMOD_EnsembleForecasting == PRESENT == ## Ensemble forecasting
  BiomodEF <- BIOMOD_EnsembleForecasting(EM.output=BiomodEM, ## Rules for assembling
                                         projection.output=BiomodProj, ## Individual model projection
                                         binary.meth=c("TSS"),
                                         filtered.meth=c("TSS"),
                                         compress=TRUE,
                                         on_0_1000=TRUE)

  ## Future distribution

  i.mod <- 1

  for (j in 1:length(fut.var[[2]])) {
    for (l in 1:length(fut.var[[3]])) {
      for (mc in 1:length(fut.var[[1]])) {
        ## Projections by model
        BiomodProjFuture <- BIOMOD_Projection(modeling.output=BiomodModel,
                                              new.env=future[[i.mod]][[mc]][[model.var]],
                                              proj.name=paste0(fut.var[[1]][mc],"_",fut.var[[2]][j],"_",fut.var[[3]][l]),
                                              selected.models=grep("_Full_",
                                                                   get_built_models(BiomodModel),
                                                                   value=TRUE), ## Full models only
                                              binary.meth=c("TSS"),
                                              filtered.meth=c("TSS"),
                                              compress=TRUE,
                                              build.clamping.mask=FALSE, ## TRUE?
                                              omi.na=TRUE,
                                              on_0_1000=TRUE,
                                              output.format=".grd")

        ## BIOMOD_EnsembleForecasting == FUTURE == ## Ensemble forecasting
        BiomodEF_Future <- BIOMOD_EnsembleForecasting(EM.output=BiomodEM, ## Rules for assembling
                                                      projection.output=BiomodProjFuture, ## Individual model projection
                                                      binary.meth=c("TSS"),
                                                      filtered.meth=c("TSS"),
                                                      compress=TRUE,
                                                      on_0_1000=TRUE)
      }
      i.mod <- i.mod+1
    }
  }

  return(list(BiomodData,BiomodModel))
}

fun.models.no.run <- function(name,spdir){

  ## BIOMOD_FormatingData
  BiomodData <- readRDS(paste0(spdir,"/BiomodData.rds"))

  ## BIOMOD_ModelingOptions
  BiomodOptions <- readRDS(paste0(spdir,"/BiomodOptions.rds"))

  ## BIOMOD_Modeling
  BiomodModel <- get(load(paste0(spdir,"/",spdir,".5mod.models.out")))

  ## Building ensemble-models
  BiomodEM <- get(load(paste0(spdir,"/",spdir,".5modensemble.models.out")))

  ## BIOMOD_Projection == PRESENT == ## Individual model projection
  BiomodProj <- get(load(paste0(spdir,"/proj_current/",spdir,".current.projection.out")))

  ## BIOMOD_EnsembleForecasting == PRESENT == ## Ensemble forecasting
  BiomodEF <- get(load(paste0(spdir,"/proj_current/",spdir,".current.ensemble.projection.out")))

  return(list(BiomodData,BiomodModel))
}
