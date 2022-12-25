#####  Calling STAN Model for the radius behaviour


load(file=paste(".ROSData", sep=""))



index.angular.change<-index.angular.change.cluster.1



################################################################


library(matlib)

library(sde)

library("rstan")

library("MotilityLab")



radius_cluster<-radius.complete.full.data[index.angular.change]


##########################################################################################################################

#####  stan file
  
                      
                  
  stanmodelcoderadius="data {
	
  int<lower=1> N;        

  vector[N] radius;    


      
}


parameters {

     
     
    real<lower=-5, upper=5> mean_radius;

    real<lower=0, upper=5> sigma_radius;

  


    
}




model {


##   Introducing no priors means we consider flat prior for the parameter. Note that the domain of the flat prior comes from the bound we have already defined for parameters.
   

  



###  the persistent and non_persistent might be with different length


   for(n in 1:N){
   	
   	           	target += lognormal_lpdf(radius[n] | mean_radius, sigma_radius); 
   	          	 	          	
   	            }  	 	          	
   	          	 	          	
   


                       
  
 }"
                                
                                
##########################################################################################################################################  STAN fit



                     data_radius <- list(N = length(radius_cluster), radius= radius_cluster)






                    save(list = ls(all=TRUE), file =".ROSData")


                    date() 
   

    
                    rstan_options(auto_write = TRUE)
                    
                    
                    options(mc.cores = parallel::detectCores())


                    
                    fit_radius <- stan(model_code = stanmodelcoderadius, model_name = "Radius",
                data = data_radius, iter = 1000, chains = 8, verbose = TRUE, control = list(adapt_delta = 0.9999, stepsize = 0.001, max_treedepth =15), cores=8)
                
                
                    fit_radius_ss <- summary(fit_radius) # fit_ss is a list


                    print(fit_radius_ss$summary)
                    
                    
                    save(list = ls(all=TRUE), file = ".ROSData")


                    date() 

                                
                                
                                
# # # #                       mean      se_mean          sd          2.5%           25%
# # mean_radius   5.238565e-01 0.0001634102 0.008176240  5.077565e-01  5.183323e-01
# # sigma_radius  6.748772e-01 0.0001443586 0.005770411  6.636358e-01  6.708982e-01
# # lp__         -1.058751e+04 0.0257870504 0.991169665 -1.059027e+04 -1.058789e+04
                       # # 50%           75%         97.5%    n_eff     Rhat
# # mean_radius   5.238702e-01  5.293540e-01  5.399454e-01 2503.506 1.000840
# # sigma_radius  6.747841e-01  6.788388e-01  6.862465e-01 1597.820 1.002051
# # lp__         -1.058719e+04 -1.058681e+04 -1.058653e+04 1477.382 1.001514
# # >                     

