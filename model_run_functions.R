library(rjags)
library(modeest)

#install.packages("rjags", lib="C:/Program Files/R/R-3.5.1/library")
#library(rjags, lib.loc = "C:/Program Files/R/R-3.5.1/library")

# ALPHA TO COLOR - cosmetic function
# function to add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#colmax and colmeans functions
colMax <- function (colData) {
  apply(colData, MARGIN=c(2), max)
}

colMin <- function (colData) {
  apply(colData, MARGIN=c(2), min)
}

#red the priors for the constants from previous calibration
constants<-read.csv("priors.matrix.csv")

#main function, model for the simulation
ICBM_Agro<-function(SOC_data, sim_length){
  
  sim_length=180
  
  #declaring the parameters for JAGS
  N.ADAPTS=1000
  N.RUNS=5000
  sampling.nr=500

  #simulation
  #create the simulation time series from the available data and the required length
  sim_years=seq(from=SOC_data[1,1], to=SOC_data[1,1]+sim_length)
  SOC_long<-mat.or.vec(2, length(sim_years))
  SOC_long[1,]<-seq(from=1, to=sim_length+1)
  SOC_data_points<-(SOC_data[,1]-SOC_data[1,1])+1
  SOC_long[2,]<-NA
  SOC_long[2,SOC_data_points]<-SOC_data[,2]

  model.ICBM_agro<- jags.model('./ICBM_Agro_0_24_SS_forward_exp_0.1.R',
                               data = list(
                                 'SOC' = SOC_long[2,],
                                 'Init' = SOC_long[2,1],
                                 'error_SOC' = 3.4,
                                 'y_elapsed' = sim_length+1,
                                 'aboveground_start'= 20,
                                 'max_biomass'= 100,
                                 'max_biomass_error'= 20,                         
                                 'annual_increment'= 1.2,                         
                                 'annual_increment_error'= 0.5,                         
                                 'trees.nontrees'= 0.6,
                                 'trees.nontrees_err'= 0.2,
                                 'amendments'= 0,
                                 'amendments_err'= 0.1,
                                 'S.R.plants'= 0.26,
                                 'S.R.plants_err'= 0.2,
                                 'S.R.trees'= 0.26,
                                 'S.R.trees_err'= 0.2,
                                 'D.R.plants'= 0.5, # TO KEEP FIXED
                                 'D.R.plants_err'= 0.1, # TO KEEP FIXED
                                 'D.R.trees'= 0.5, # TO KEEP FIXED
                                 'D.R.trees_err'= 0.1, # TO KEEP FIXED
                                 'Woody_root_mortality_y'= 20,
                                 'Woody_root_mortality_y_err'= 2,
                                 'Climatic.factor' = 5,
                                 'error_re' = 3,
                                 'k1_prior'=constants[,2],
                                 'k2_prior'=constants[,3],
                                 'hr_prior'=constants[,4],
                                 'hw_prior'=constants[,5],
                                 'hs_prior'=constants[,6],
                                 'ha_prior'=constants[,7],
                                 'prior_size'=dim(constants)[1]),
                               n.chains = 1,
                               n.adapt = N.ADAPTS)

  update(model.ICBM_agro, n.iter=N.RUNS, n.burnin=N.BURNIN)

  parameter_list<-c("k1",
                    "k2",
                    "Tot",
                    "h_r",
                    "h_w",
                    "h_s",
                    "h_a",
                    "I_r",
                    "I_w",
                    "I_s",
                    "I_a",
                    "AG_annual",
                    "AG_plants",
                    "AG_trees_annual",
                    "re")
  

  mcmc.array.ICBM.agro<-(jags.samples(model.ICBM_agro,parameter_list, sampling.nr))

  return(list(mcmc.array.ICBM.agro, SOC_long))

}



#### plotting

#total SOC
plot.SOC<-function(model_out){
  model_out1<-model_out[[1]]
  SOC_long<-model_out[[2]]
mcmc.list.SOC<-as.mcmc.list(model_out1$Tot, chains=F)
mcmc.unlist.SOC<-mcmc.list.SOC[[1]]
prediction_array.SOC.min<-colMin(mcmc.unlist.SOC)
prediction_array.SOC.max<-colMax(mcmc.unlist.SOC)
prediction_array.SOC.mean<-colMeans(mcmc.unlist.SOC)
yearseq<-seq(from=0, to=(length(prediction_array.SOC.mean)-1))

plot(yearseq,prediction_array.SOC.mean, type="l", col="darkred",
     main="SOC",
     xlab="Elapsed years", ylab=expression(paste("C (Mg ha"^-1,")")), ylim=c(0,350))
polygon(c(yearseq,rev(yearseq)),
        c(prediction_array.SOC.max,rev(prediction_array.SOC.min)),
        col=add.alpha("darkred",0.65),border=add.alpha("darkred",0.40))
points( SOC_long[1,], SOC_long[2,],
        pch=8)
}

plot.AG<-function(model_out){
  model_out1<-model_out[[1]]
#AG_annual
mcmc.list.AG_annual<-as.mcmc.list(model_out1$AG_annual, chains=F)
mcmc.unlist.AG_annual<-mcmc.list.AG_annual[[1]]
prediction_array.AG_annual.min<-colMin(mcmc.unlist.AG_annual)
prediction_array.AG_annual.max<-colMax(mcmc.unlist.AG_annual)
prediction_array.AG_annual.mean<-colMeans(mcmc.unlist.AG_annual)

plot(yearseq,prediction_array.AG_annual.mean, type="l", col="darkgreen",
     main="Abovegorund annual production",
     xlab="Elapsed years", ylab=expression(paste("C (Mg ha"^-1,")")), ylim=c(0,350))
polygon(c(yearseq,rev(yearseq)),
        c(prediction_array.AG_annual.max,rev(prediction_array.AG_annual.min)),
        col=add.alpha("darkgreen",0.65),border=add.alpha("darkgreen",0.40))
}



plot.inputs<-function(model_out){

  model_out1<-model_out[[1]]
  
  #inputs

#I_r
mcmc.list.I_r<-as.mcmc.list(model_out1$I_r, chains=F)
mcmc.unlist.I_r<-mcmc.list.I_r[[1]]
prediction_array.I_r.min<-colMin(mcmc.unlist.I_r)
prediction_array.I_r.max<-colMax(mcmc.unlist.I_r)
prediction_array.I_r.mean<-colMeans(mcmc.unlist.I_r)
#total I_w
mcmc.list.I_w<-as.mcmc.list(model_out1$Tot, chains=F)
mcmc.unlist.I_w<-mcmc.list.I_w[[1]]
prediction_array.I_w.min<-colMin(mcmc.unlist.I_w)
prediction_array.I_w.max<-colMax(mcmc.unlist.I_w)
prediction_array.I_w.mean<-colMeans(mcmc.unlist.I_w)
#total I_a
mcmc.list.I_a<-as.mcmc.list(model_out1$Tot, chains=F)
mcmc.unlist.I_a<-mcmc.list.I_a[[1]]
prediction_array.I_a.min<-colMin(mcmc.unlist.I_a)
prediction_array.I_a.max<-colMax(mcmc.unlist.I_a)
prediction_array.I_a.mean<-colMeans(mcmc.unlist.I_a)
#total I_s
mcmc.list.I_s<-as.mcmc.list(model_out1$Tot, chains=F)
mcmc.unlist.I_s<-mcmc.list.I_s[[1]]
prediction_array.I_s.min<-colMin(mcmc.unlist.I_s)
prediction_array.I_s.max<-colMax(mcmc.unlist.I_s)
prediction_array.I_s.mean<-colMeans(mcmc.unlist.I_s)

par(mfrow=c(2,2))
Input_palette<-c("chocolate3", "burlywood4", "brown3", "cadetblue")
plot(yearseq,prediction_array.I_r.mean, type="l", col=Input_palette[1],
     main="Annual fine roots inputs",
     xlab="Elapsed years", ylab=expression(paste("C (Mg ha"^-1,")")), ylim=c(0,350))
polygon(c(yearseq,rev(yearseq)),
        c(prediction_array.I_r.max,rev(prediction_array.I_r.min)),
        col=add.alpha(Input_palette[1],0.65),border=add.alpha(Input_palette[1],0.40))

plot(yearseq,prediction_array.I_w.mean, type="l", col=Input_palette[2],
     main="Annual coarse roots inputs",
     xlab="Elapsed years", ylab=expression(paste("C (Mg ha"^-1,")")), ylim=c(0,350))
polygon(c(yearseq,rev(yearseq)),
        c(prediction_array.I_w.max,rev(prediction_array.I_w.min)),
        col=add.alpha(Input_palette[2],0.65),border=add.alpha(Input_palette[2],0.40))

plot(yearseq,prediction_array.I_s.mean, type="l", col=Input_palette[3],
     main="Annual shoots inputs",
     xlab="Elapsed years", ylab=expression(paste("C (Mg ha"^-1,")")), ylim=c(0,350))
polygon(c(yearseq,rev(yearseq)),
        c(prediction_array.I_s.max,rev(prediction_array.I_s.min)),
        col=add.alpha(Input_palette[3],0.65),border=add.alpha(Input_palette[3],0.40))

plot(yearseq,prediction_array.I_a.mean, type="l", col=Input_palette[4],
     main="Annual amendments inputs",
     xlab="Elapsed years", ylab=expression(paste("C (Mg ha"^-1,")")), ylim=c(0,350))
polygon(c(yearseq,rev(yearseq)),
        c(prediction_array.I_a.max,rev(prediction_array.I_a.min)),
        col=add.alpha(Input_palette[4],0.65),border=add.alpha(Input_palette[4],0.40))
}
        
        
SS<-function(model_out){
  model_out1<-model_out[[1]]
  
mcmc.list.k1<-unlist(as.mcmc.list(model_out1$k1, chains=F))
mcmc.list.k2<-unlist(as.mcmc.list(model_out1$k2, chains=F))
mcmc.list.hr<-unlist(as.mcmc.list(model_out1$h_r, chains=F))
mcmc.list.hw<-unlist(as.mcmc.list(model_out1$h_w, chains=F))
mcmc.list.hs<-unlist(as.mcmc.list(model_out1$h_s, chains=F))
mcmc.list.ha<-unlist(as.mcmc.list(model_out1$h_a, chains=F))
mcmc.list.re<-unlist(as.mcmc.list(model_out1$re, chains=F))
mcmc.list.I_r<-unlist(as.mcmc.list(model_out1$I_r, chains=F))
mcmc.list.I_w<-unlist(as.mcmc.list(model_out1$I_w, chains=F))
mcmc.list.I_a<-unlist(as.mcmc.list(model_out1$I_a, chains=F))
mcmc.list.I_s<-unlist(as.mcmc.list(model_out1$I_s, chains=F))

YrSS_vector<-mcmc.list.I_r/((mcmc.list.k1)*(mcmc.list.re))
YsSS_vector<-mcmc.list.I_s/((mcmc.list.k1)*(mcmc.list.re))
YwSS_vector<-mcmc.list.I_w/((mcmc.list.k1)*(mcmc.list.re))
YaSS_vector<-mcmc.list.I_a/((mcmc.list.k1)*(mcmc.list.re))
OSS_vector<-(mcmc.list.hr*mcmc.list.I_r+mcmc.list.hs*mcmc.list.I_s+
                                    mcmc.list.ha*mcmc.list.I_a+mcmc.list.hw*mcmc.list.I_w)/
                           (mcmc.list.k2*mcmc.list.re)

SS<-(OSS_vector+
       YrSS_vector+
       YsSS_vector+
       YaSS_vector+
       YwSS_vector)

return(SS)
}


plot.SS<-function(SS){
dev.off()
plot(density(SS), main = "Steady state SOC", xlab=expression(paste("C (Mg ha"^-1,")")))
polygon(density(SS), col=add.alpha("darkblue", 0.4))
abline(v=mlv(SS, method = "naive"), lty=2)
abline(v=mean(SS))
legend("topright", c("Mean", "Mode"), lty=c(1, 2), bty="n")
}



