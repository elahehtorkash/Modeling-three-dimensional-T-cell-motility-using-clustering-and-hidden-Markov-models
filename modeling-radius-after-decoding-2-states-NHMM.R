#####  Calling STAN Model for the radius behaviour


load(file=paste(".RData", sep=""))


index.angular.change<-index.angular.change.cluster.1
################################################################


library(matlib)

library(sde)

library("rstan")

library("MotilityLab")



P_1<-matrix(c(fit_ss$summary[3:6]), nrow=2, ncol=2, byrow=T)

P_2<-matrix(c(fit_ss$summary[7:10]), nrow=2, ncol=2, byrow=T)

prior_pi<-c(fit_ss$summary[1], fit_ss$summary[2])



beta_persistent<-fit_ss$summary[12]

alpha_persistent<-fit_ss$summary[11]

beta_non_persistent<-fit_ss$summary[14]

alpha_non_persistent<-fit_ss$summary[13]



markov_probability_cluster<-matrix(0, nrow= number.of.positions.total, ncol=2)

for(i in 1: number.of.cells){
					 
					 index.i<-which(data_cell_motility[, 8]==id[i])
					 
					 markov_probability_cluster[index.i[1], ]<-prior_pi
					 
					 for(l in 2:(length(index.i))){
					 	
					 	if(number.of.neighbors.full.data[index.i[l-1]]==0){
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_1
					 	}
					 	
					 	
					 	if(number.of.neighbors.full.data[index.i[l-1]]>0){
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_2
					 	}
					 	
					 } 
					 
				}
				
				
######  Defining the ratio for the purpose of classification

t<-which(angular.change.complete.full.data==0)

angular.change.complete.full.data[t]<-0.00001

classification.matrix<-matrix(c(angular.change.complete.full.data, rep(0, length(angular.change.complete.full.data))), nrow=length(angular.change.complete.full.data), ncol=2)

group.probability<-matrix(0, nrow=length(angular.change.complete.full.data), ncol=2)


group.probability[, 1]<-markov_probability_cluster[, 1]*dbeta((angular.change.complete.full.data/pi), alpha_non_persistent, beta_non_persistent)/(markov_probability_cluster[, 1]*dbeta((angular.change.complete.full.data/pi), alpha_non_persistent, beta_non_persistent) +markov_probability_cluster[, 2]*dbeta((angular.change.complete.full.data/pi), alpha_persistent, beta_persistent))


group.probability[, 2]<-markov_probability_cluster[, 2]*dbeta((angular.change.complete.full.data/pi), alpha_persistent, beta_persistent)/(markov_probability_cluster[, 1]*dbeta((angular.change.complete.full.data/pi), alpha_non_persistent, beta_non_persistent) +markov_probability_cluster[, 2]*dbeta((angular.change.complete.full.data/pi), alpha_persistent, beta_persistent))

sum(group.probability[index.angular.change, 1]==group.probability[index.angular.change, 2])

for(i in 1:length(index.angular.change)){
	
	if(group.probability[index.angular.change[i], 1]>=group.probability[index.angular.change[i], 2]){classification.matrix[index.angular.change[i], 2]<-1}
	
	
	if(group.probability[index.angular.change[i], 1]<group.probability[index.angular.change[i], 2]){classification.matrix[index.angular.change[i], 2]<-2}
	
	
}



weight<-NULL

weight[1]<-mean(classification.matrix[index.angular.change, 2]==1)

weight[2]<-mean(classification.matrix[index.angular.change, 2]==2)



### 1:black <- non-persistent

### 2:darkred <- persistent


index.persistent<-which(classification.matrix[index.angular.change, 2]==2)

index.non.persistent<-which(classification.matrix[index.angular.change, 2]==1)



persist<-radius.complete.full.data[index.angular.change[index.persistent]]

non.persist<-radius.complete.full.data[index.angular.change[index.non.persistent]]



##########################################################################################################################

#####  stan file
  
                      
                  
  stanmodelcoderadius="data {
	
  int<lower=1> N_persist;        

  vector[N_persist] radius_persist;    


  int<lower=1> N_non_persist;       

  vector[N_non_persist] radius_non_persist;    

   
  real mu_0;
	  
  real<lower=0> sigma_0;

      
}


parameters {

     
     
    real mean_non_per_radius;

    real<lower=0> sigma_non_per_radius;

    real mean_per_radius;

    real<lower=0> sigma_per_radius;
   
    
    real<lower=0> alpha_inv;
	  
	real<lower=0> beta_inv;
	  


    
}




model {


##   Introducing no priors means we consider flat prior for the parameter. Note that the domain of the flat prior comes from the bound we have already defined for parameters.
   
  
           target+=inv_gamma_lpdf(sigma_non_per_radius|alpha_inv, beta_inv);

           target+=inv_gamma_lpdf(sigma_per_radius|alpha_inv, beta_inv);

           # target+=inv_gamma_lpdf(beta_hesi_radius|alpha_inv, beta_inv);

     

           target+=normal_lpdf(mean_non_per_radius| mu_0, sigma_0);

           target+=normal_lpdf(mean_per_radius| mu_0, sigma_0);

  



###  the persistent and non_persistent might be with different length


   for(n in 1:N_persist){
   	
   	           	target += lognormal_lpdf(radius_persist[n] | mean_per_radius, sigma_per_radius); 
   	          	 	          	
   	            }  	 	          	
   	          	 	          	

   for(n in 1:N_non_persist){
   	
   	           	target += lognormal_lpdf(radius_non_persist[n] | mean_non_per_radius, sigma_non_per_radius); 
   	          	 	          	
   	            }  	 	          	

   


                       
  
 }"
                                
                                
##########################################################################################################################################  STAN fit


                    data_radius <- list(N_persist = length(persist), radius_persist= persist, N_non_persist = length(non.persist), radius_non_persist= non.persist, mu_0=0, sigma_0=1.5)







                    save(list = ls(all=TRUE), file =".RData")


                    date() 
   

    
                    rstan_options(auto_write = TRUE)
                    
                    
                    options(mc.cores = parallel::detectCores())


                    
                    fit_radius <- stan(model_code = stanmodelcoderadius, model_name = "Radius",
                data = data_radius, iter = 2000, chains = 8, verbose = TRUE, control = list(adapt_delta = 0.9999, stepsize = 0.001, max_treedepth =15), cores=8)
                
                
                    fit_radius_ss <- summary(fit_radius) # fit_ss is a list


                    print(fit_radius_ss$summary)
                    
                    
                    save(list = ls(all=TRUE), file = ".RData")


                    date() 

                                
                                
                                

