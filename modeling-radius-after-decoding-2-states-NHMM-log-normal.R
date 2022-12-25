#####  Calling STAN Model for the radius behaviour


load(file=paste(".RData", sep=""))

index.angular.change<-index.angular.change.cluster.1


cluster.1<-c(1,   2,   6,   8,   9,  12,  17,  18,  19,  20,  21,  22,  46,  52,  54,  57,  60,  62,  68,  69,  77,  83,  84,  85,  86,  90,  91, 104, 105, 113, 124, 127, 128, 132, 135, 150, 151, 152, 153, 157, 160, 162, 163, 166, 169, 171, 173, 178, 183, 187, 191, 192, 193, 199, 206, 207, 216, 218, 234, 238, 241, 244, 247, 248, 252, 253, 256, 257, 258, 260, 266, 267, 273, 276, 277, 283, 289, 293, 296, 298, 302, 303, 304, 310, 313, 321, 322, 324, 325, 332, 341, 345, 346)
			         
			                                     
cluster.2<-c( 3,   7,  10,  13,  15,  16,  23,  24,  28,  32,  35,  36,  37,  39,  45,  53,  55,  71,  73,  75,  78,  92,  93, 103, 108, 109, 110, 111, 130, 133, 134, 137, 139, 143, 144, 146, 148, 161, 165, 172, 174, 197, 201, 202, 204, 210, 211, 217, 221, 228, 230, 231, 240, 242, 243, 251, 254, 255, 263, 264, 265, 269, 271, 281, 287, 294, 295, 299, 300, 315, 318, 323, 327, 329, 330, 331, 335, 336, 338, 339, 342, 343, 349, 352)


cluster.3<-c( 4,   5,  11,  14,  25,  26,  29,  30,  31,  33,  34,  38,  40,  41,  47 , 48 , 50 , 51 , 58,  59,  61,  63,  64,  65,  66,  67,  70,  72,  74,  79,  80,  81,  82,  87,  88,  94,  95,  96 , 97,  98, 100, 101, 102, 106, 107, 114, 115, 116, 118, 120 ,121, 122, 123, 125, 126, 129, 136, 140, 141, 142, 145, 147, 149, 154, 155, 156, 158, 164, 167, 168, 170, 175, 179, 180 ,181, 182, 184, 185, 186, 188 ,189, 190 ,194, 200, 203, 205, 208, 209, 213, 214, 219, 220, 222, 223, 224, 225 ,227, 233, 236, 237, 239 ,250, 259, 261, 262, 268, 270, 272, 274, 275, 278, 279, 280, 282, 284, 286, 288, 290 ,292, 297, 305, 307 ,308, 311, 312, 316, 320, 326, 328, 333, 334, 337, 340, 344, 348, 350, 351)


cluster.4<-c(27,  42,  43,  44,  49,  56 , 76,  89,  99, 112, 117, 119, 131, 138, 159 ,176 ,177 ,195 ,196 ,198 ,212, 215, 226, 229, 232, 235, 245, 246, 249, 285, 291, 301, 306, 309, 314, 317, 319, 347)



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

   
  # real mu_0;
	  
  # real<lower=0> sigma_0;

      
}


parameters {

     
     
    real<lower=-5, upper=5> mean_non_per_radius;

    real<lower=0, upper=5> sigma_non_per_radius;

    real<lower=-5, upper=5> mean_per_radius;

    real<lower=0, upper=5> sigma_per_radius;
   
    
    # real<lower=0> alpha_inv;
	  
	# real<lower=0> beta_inv;
	  


    
}




model {


##   Introducing no priors means we consider flat prior for the parameter. Note that the domain of the flat prior comes from the bound we have already defined for parameters.
   
  
# # # #            target+=inv_gamma_lpdf(sigma_non_per_radius|alpha_inv, beta_inv);

           # # target+=inv_gamma_lpdf(sigma_per_radius|alpha_inv, beta_inv);

           # # # target+=inv_gamma_lpdf(beta_hesi_radius|alpha_inv, beta_inv);

     

           # # target+=normal_lpdf(mean_non_per_radius| mu_0, sigma_0);

           # # target+=normal_lpdf(mean_per_radius| mu_0, sigma_0);

  



###  the persistent and non_persistent might be with different length


   for(n in 1:N_persist){
   	
   	           	target += lognormal_lpdf(radius_persist[n] | mean_per_radius, sigma_per_radius); 
   	          	 	          	
   	            }  	 	          	
   	          	 	          	

   for(n in 1:N_non_persist){
   	
   	           	target += lognormal_lpdf(radius_non_persist[n] | mean_non_per_radius, sigma_non_per_radius); 
   	          	 	          	
   	            }  	 	          	

   


                       
  
 }"
                                
                                
##########################################################################################################################################  STAN fit


                    # data_radius <- list(N_persist = length(persist), radius_persist= persist, N_non_persist = length(non.persist), radius_non_persist= non.persist, mu_0=0, sigma_0=1.5)



                     data_radius <- list(N_persist = length(persist), radius_persist= persist, N_non_persist = length(non.persist), radius_non_persist= non.persist)






                    save(list = ls(all=TRUE), file =".RData")


                    date() 
   

    
                    rstan_options(auto_write = TRUE)
                    
                    
                    options(mc.cores = parallel::detectCores())


                    
                    fit_radius <- stan(model_code = stanmodelcoderadius, model_name = "Radius",
                data = data_radius, iter = 1000, chains = 8, verbose = TRUE, control = list(adapt_delta = 0.9999, stepsize = 0.001, max_treedepth =15), cores=8)
                
                
                    fit_radius_ss <- summary(fit_radius) # fit_ss is a list


                    print(fit_radius_ss$summary)
                    
                    
                    save(list = ls(all=TRUE), file = ".RData")


                    date() 

                                
                                
                                

