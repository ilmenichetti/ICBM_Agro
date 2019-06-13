

##Running the functions
#read the time series of SOC
SOC_data<-read.csv("SOC_fake_data.csv")
#run the sims
model_out<-ICBM_Agro(SOC_data = SOC_data , sim_length = 180)
SS<-SS(model_out)
#plotting results
plot.SOC(model_out)
plot.AG(model_out)
plot.inputs(model_out)
plot.SS(SS)
