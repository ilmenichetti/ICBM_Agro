model{
  
    #pool initializatyion
    Y_r[1]<-(Init*(1-Init_ratio))*0.5
    Y_w[1]<-(Init*(1-Init_ratio))*0.3
    Y_s[1]<-(Init*(1-Init_ratio))*0.2
    Y_a[1]<-0

    O[1]  <-Init*Init_ratio
    
    
    # interpolation of the aboveground woody biomass
    AG_annual_change ~ dnorm(annual_increment,1/annual_increment_error)  T(0.001, annual_increment*2)
    max_biomass_stochastic ~ dnorm(max_biomass,1/max_biomass_error) T(0.001, max_biomass*2)
   
      
    #### Loop for time
    for (i in 1:(y_elapsed)){
      
      #resample k1 from the distribution coming from the calibration
      k1_sampling[i]~dunif(1, (prior_size-1))
      k1[i]<-k1_prior[abs(round(k1_sampling[i]))]
      k2[i]<-k2_prior[abs(round(k1_sampling[i]))]
      h_r[i]<-hr_prior[abs(round(k1_sampling[i]))]
      h_w[i]<-hw_prior[abs(round(k1_sampling[i]))]
      h_s[i]<-hs_prior[abs(round(k1_sampling[i]))]
      h_a[i]<-ha_prior[abs(round(k1_sampling[i]))]
     
      ### Input calculation
      #calculation of the annual aboveground
      AG_annual[i]        <- ifelse((aboveground_start + AG_annual_change*i)<max_biomass_stochastic, 
                                      aboveground_start + AG_annual_change*i, 
                                      max_biomass_stochastic)
      AG_trees_annual[i]  <- AG_annual[i]*J.trees.nontrees
      AG_plants[i]        <- AG_annual[i]*(1-J.trees.nontrees)
      
      #calculate litter production according to Mohan Kumar B. Litter Dynamics in Plantation and Agroforestry Systems of the Tropics-A Review of Observations and Methods. Ecological Basis of Agroforestry. 2010. doi:10.1201/9781420043365.ch10
      Litter.production_t_ha_y[i]<-0.0325*AG_annual[i]+2.328

      # calculation of the woody roots, gradual release
      I_w[i]              <- (AG_trees_annual[i]*J.S.R.trees)*1/J.Woody_root_mortality
      
      # calculation of the fine roots and rhizodeposition
      I_r_plants[i]       <- AG_plants[i]*J.S.R.plants
      I_r_rhizodep[i]     <- I_r_plants[i]*J.D.R.plants
      I_w_rhizodep[i]     <- I_w[i]*J.D.R.trees
      I_r[i]              <- I_r_plants[i]+I_r_rhizodep[i]+I_w_rhizodep[i]
      
      #calculation of aboveground inputs
      I_s_trees[i]        ~ dnorm(Litter.production_t_ha_y[i],1/(Litter.production_t_ha_y[i]*0.17)) T(0.001,300)
      #I_s_trees[i]        <- Litter.production_t_ha_y[i]
      I_s_plants[i]       <- AG_plants[i]*Stubble_proportion
      I_s[i]              <-I_s_trees[i]+I_s_plants[i]
      
      #eventual_amendments
      I_a[i]              <- J.amendments
      

      ### SOC simulation
      # young pool
      Y_r[i+1] 		<-  (I_r[i]+Y_r[i])*exp(-k1[i]*re)
      Y_w[i+1] 		<-  (I_w[i]+Y_w[i])*exp(-k1[i]*re)
      Y_s[i+1] 		<-  (I_s[i]+Y_s[i])*exp(-k1[i]*re)
      Y_a[i+1] 		<-  (I_a[i]+Y_a[i])*exp(-k1[i]*re)
      
      # fluxes
      flux_r[i] 		  <-  h_r[i]*(k1[i]*(Y_r[i]+I_r[i])/(k2[i]-k1[i]))
      flux_w[i] 		  <-  h_w[i]*(k1[i]*(Y_w[i]+I_w[i])/(k2[i]-k1[i]))
      flux_s[i] 		  <-  h_s[i]*(k1[i]*(Y_s[i]+I_s[i])/(k2[i]-k1[i]))
      flux_a[i] 		  <-  h_a[i]*(k1[i]*(Y_a[i]+I_a[i])/(k2[i]-k1[i]))
      
      # old pool
      O[i+1]   	<-  (O[i]-
                       flux_r[i]-
                       flux_w[i]-
                       flux_s[i]-
                       flux_r[i])*exp(-k2[i]*re) +
                       flux_r[i]*exp(-k1[i]*re) +
                       flux_w[i]*exp(-k1[i]*re) +
                       flux_s[i]*exp(-k1[i]*re) +
                       flux_a[i]*exp(-k1[i]*re)
      #O[i]   	<-  1
      
      # total C
      Tot[i] <- Y_r[i] + 
                  Y_w[i] + 
                  Y_s[i] + 
                  Y_a[i] + 
                  O[i]

      #Error of the measurement (assumed proportional to the measurement) 
      SOC[i]  ~ dnorm(Tot[i],1/(error_SOC))
    
      }
    

    #climatic factro plus error
    re      ~ dunif(Climatic.factor-Climatic.factor*error_re,Climatic.factor+Climatic.factor*error_re)
    
    #priors for calculation of inputs
    J.trees.nontrees ~ dnorm(trees.nontrees,1/trees.nontrees_err) T(0.001,0.8)
    J.S.R.plants ~ dnorm(S.R.plants,1/S.R.plants_err) T(0.001,100)
    J.S.R.trees ~ dnorm(S.R.trees,1/S.R.trees_err) T(0.001,100)
    J.D.R.plants ~ dnorm(D.R.plants,1/D.R.plants_err) T(0.001,100)
    J.D.R.trees ~ dnorm(D.R.trees,1/D.R.trees_err) T(0.001,100)
    J.Woody_root_mortality ~ dnorm(Woody_root_mortality_y,1/Woody_root_mortality_y_err) T(0.001,20)
    J.amendments  ~ dnorm(amendments,amendments_err) T(0.001,100)
    
  
  # getting the values valid for all sites
  
  ##Parameters

  Stubble_proportion ~ dunif(0.1, 0.3)
  Init_ratio ~ dnorm(0.9291667,1/0.01) T(0.8,0.98)

  #error terms for all sites
  error_k<-0.05
  limits_k<-0.1
  error_h<-0.05
  limits_h<-0.1

}
