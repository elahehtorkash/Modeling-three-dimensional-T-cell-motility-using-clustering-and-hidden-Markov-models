library(rstan)

remove(list=ls())


aic.clusters.4.cluster<-matrix(0, nrow=4)


aic.clusters.2.cluster<-matrix(0, nrow=2)



setwd("/global/scratch/torkashe/Clustering-Results/cluster-1-3-states-HMM/")

load(".RData")



cluster.c<-1



###################################################


log.likelihood.f<-function(parameter){
	
	
	
	
	
	prior_pi <-parameter[1:3]
	
	P_1 <- matrix(parameter[4:12], nrow=3, ncol=3, byrow=T)
	
	# P_2<-matrix(parameter[13:21], nrow=3, ncol=3, byrow=T)
	
	alpha_persistent <- parameter[13]
	
	beta_persistent <- parameter[14]
	
	alpha_non_persistent <- parameter[15]
	
	beta_non_persistent <- parameter[16]
	
	alpha_transition<- parameter[17]
	
	beta_transition<- parameter[18]
	
	counter.cluster<-parameter[19]
	
	

markov_probability_cluster<-matrix(0, nrow= number.of.positions.total, ncol=3)


likelihood_cell<-matrix(1, nrow= length(cluster[[counter.cluster]]), ncol=1)



for(i in 1: length(cluster[[counter.cluster]])){

					 index.i<-which(data_cell_motility[, 8]==id[cluster[[counter.cluster]][i]])
					 
				 
				     n.i<-length(index.i)
				 
					 
					 		
					 markov_probability_cluster[index.i[1], ]= prior_pi %*% P_1
					 	
					  	
					 
					 
					 for(l in 2:(length(index.i))){
					 	
					 	
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_1 %*% diag(c(dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_non_persistent, beta_non_persistent), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_transition, beta_transition), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_persistent, beta_persistent)))
					 	
					 	
					 	
					 	
					 	
					 } 
					 
					 
					 likelihood_cell[i]<-sum((markov_probability_cluster[index.i[n.i-1], ]))
					 
					 # print(likelihood_cell[i])
					 
					 
					}


              t<-sum(log(likelihood_cell[1:counter.cluster]))
              
              # t<-likelihood_cell[1:counter.cluster]
              
              return(t)
	
	
	
}


###################################################




P_1<-matrix(c(fit_ss$summary[4:12]), nrow=3, ncol=3, byrow=T)

# P_2<-matrix(c(fit_ss$summary[4:12]), nrow=3, ncol=3, byrow=T)

prior_pi<-c(fit_ss$summary[1:3])

beta_non_persistent<-fit_ss$summary[16]

alpha_non_persistent<-fit_ss$summary[15]

beta_persistent<-fit_ss$summary[14]

alpha_persistent<-fit_ss$summary[13]

alpha_transition<-beta_transition<-fit_ss$summary[17]

parameter.bayes<-cbind(prior_pi[1], prior_pi[2], prior_pi[3], P_1[1, 1], P_1[1, 2], P_1[1, 3], P_1[2, 1], P_1[2, 2], P_1[2, 3], P_1[3, 1], P_1[3, 2], P_1[3, 3], alpha_persistent, beta_persistent, alpha_non_persistent, beta_non_persistent, alpha_transition, beta_transition, cluster.c)


##########################################################################


aic.clusters.4.cluster[cluster.c]<-2*(length(parameter.bayes)-1)-2*log.likelihood.f(parameter.bayes)



print("aic.clusters.4.cluster")

print(aic.clusters.4.cluster)


###########################################################################################################################################################################################



setwd("/global/scratch/torkashe/Clustering-Results/cluster-2-3-states-HMM/")

load(".RData")



cluster.c<-2



###################################################


log.likelihood.f<-function(parameter){
	
	
	
	
	
	prior_pi <-parameter[1:3]
	
	P_1<-matrix(parameter[4:12], nrow=3, ncol=3, byrow=T)
	
	# P_2<-matrix(parameter[13:21], nrow=3, ncol=3, byrow=T)
	
	alpha_persistent <-parameter[13]
	
	beta_persistent <-parameter[14]
	
	alpha_non_persistent <-parameter[15]
	
	beta_non_persistent <-parameter[16]
	
	alpha_transition<-parameter[17]
	
	beta_transition<-parameter[18]
	
	counter.cluster<-parameter[19]
	
	

markov_probability_cluster<-matrix(0, nrow= number.of.positions.total, ncol=3)


likelihood_cell<-matrix(1, nrow= length(cluster[[counter.cluster]]), ncol=1)



for(i in 1: length(cluster[[counter.cluster]])){

					 index.i<-which(data_cell_motility[, 8]==id[cluster[[counter.cluster]][i]])
					 
				 
				     n.i<-length(index.i)
				 
					 
					 		
					 markov_probability_cluster[index.i[1], ]= prior_pi %*% P_1
					 	
					  	
					 
					 
					 for(l in 2:(length(index.i))){
					 	
					 	
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_1 %*% diag(c(dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_non_persistent, beta_non_persistent), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_transition, beta_transition), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_persistent, beta_persistent)))
					 	
					 	
					 	
					 	
					 	
					 } 
					 
					 
					 likelihood_cell[i]<-sum((markov_probability_cluster[index.i[n.i-1], ]))
					 
					 # print(likelihood_cell[i])
					 
					 
					}


              t<-sum(log(likelihood_cell[1:counter.cluster]))
              
              # t<-likelihood_cell[1:counter.cluster]
              
              return(t)
	
	
	
}


###################################################




P_1<-matrix(c(fit_ss$summary[4:12]), nrow=3, ncol=3, byrow=T)

# P_2<-matrix(c(fit_ss$summary[4:12]), nrow=3, ncol=3, byrow=T)

prior_pi<-c(fit_ss$summary[1:3])

beta_non_persistent<-fit_ss$summary[16]

alpha_non_persistent<-fit_ss$summary[15]

beta_persistent<-fit_ss$summary[14]

alpha_persistent<-fit_ss$summary[13]

alpha_transition<-beta_transition<-fit_ss$summary[17]

parameter.bayes<-cbind(prior_pi[1], prior_pi[2], prior_pi[3], P_1[1, 1], P_1[1, 2], P_1[1, 3], P_1[2, 1], P_1[2, 2], P_1[2, 3], P_1[3, 1], P_1[3, 2], P_1[3, 3], alpha_persistent, beta_persistent, alpha_non_persistent, beta_non_persistent, alpha_transition, beta_transition, cluster.c)


##########################################################################


aic.clusters.4.cluster[cluster.c]<-2*(length(parameter.bayes)-1)-2*log.likelihood.f(parameter.bayes)



print("aic.clusters.4.cluster")

print(aic.clusters.4.cluster)


###########################################################################################################################################################################################



setwd("/global/scratch/torkashe/Clustering-Results/cluster-3-3-states-HMM/")

load(".RData")



cluster.c<-3



###################################################


log.likelihood.f<-function(parameter){
	
	
	
	
	
	prior_pi <-parameter[1:3]
	
	P_1<-matrix(parameter[4:12], nrow=3, ncol=3, byrow=T)
	
	# P_2<-matrix(parameter[13:21], nrow=3, ncol=3, byrow=T)
	
	alpha_persistent <-parameter[13]
	
	beta_persistent <-parameter[14]
	
	alpha_non_persistent <-parameter[15]
	
	beta_non_persistent <-parameter[16]
	
	alpha_transition<-parameter[17]
	
	beta_transition<-parameter[18]
	
	counter.cluster<-parameter[19]
	
	

markov_probability_cluster<-matrix(0, nrow= number.of.positions.total, ncol=3)


likelihood_cell<-matrix(1, nrow= length(cluster[[counter.cluster]]), ncol=1)



for(i in 1: length(cluster[[counter.cluster]])){

					 index.i<-which(data_cell_motility[, 8]==id[cluster[[counter.cluster]][i]])
					 
				 
				     n.i<-length(index.i)
				 
					 
					 		
					 markov_probability_cluster[index.i[1], ]= prior_pi %*% P_1
					 	
					  	
					 
					 
					 for(l in 2:(length(index.i))){
					 	
					 	
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_1 %*% diag(c(dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_non_persistent, beta_non_persistent), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_transition, beta_transition), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_persistent, beta_persistent)))
					 	
					 	
					 	
					 	
					 	
					 } 
					 
					 
					 likelihood_cell[i]<-sum((markov_probability_cluster[index.i[n.i-1], ]))
					 
					 # print(likelihood_cell[i])
					 
					 
					}


              t<-sum(log(likelihood_cell[1:counter.cluster]))
              
              # t<-likelihood_cell[1:counter.cluster]
              
              return(t)
	
	
	
}


###################################################




P_1<-matrix(c(fit_ss$summary[4:12]), nrow=3, ncol=3, byrow=T)

# P_2<-matrix(c(fit_ss$summary[4:12]), nrow=3, ncol=3, byrow=T)

prior_pi<-c(fit_ss$summary[1:3])

beta_non_persistent<-fit_ss$summary[16]

alpha_non_persistent<-fit_ss$summary[15]

beta_persistent<-fit_ss$summary[14]

alpha_persistent<-fit_ss$summary[13]

alpha_transition<-beta_transition<-fit_ss$summary[17]

parameter.bayes<-cbind(prior_pi[1], prior_pi[2], prior_pi[3], P_1[1, 1], P_1[1, 2], P_1[1, 3], P_1[2, 1], P_1[2, 2], P_1[2, 3], P_1[3, 1], P_1[3, 2], P_1[3, 3], alpha_persistent, beta_persistent, alpha_non_persistent, beta_non_persistent, alpha_transition, beta_transition, cluster.c)


##########################################################################


aic.clusters.4.cluster[cluster.c]<-2*(length(parameter.bayes)-1)-2*log.likelihood.f(parameter.bayes)



print("aic.clusters.4.cluster")

print(aic.clusters.4.cluster)


###########################################################################################################################################################################################



setwd("/global/scratch/torkashe/Clustering-Results/cluster-4-3-states-HMM/")

load(".RData")



cluster.c<-4



###################################################


log.likelihood.f<-function(parameter){
	
	
	
	
	
	prior_pi <-parameter[1:3]
	
	P_1<-matrix(parameter[4:12], nrow=3, ncol=3, byrow=T)
	
	# P_2<-matrix(parameter[13:21], nrow=3, ncol=3, byrow=T)
	
	alpha_persistent <-parameter[13]
	
	beta_persistent <-parameter[14]
	
	alpha_non_persistent <-parameter[15]
	
	beta_non_persistent <-parameter[16]
	
	alpha_transition<-parameter[17]
	
	beta_transition<-parameter[18]
	
	counter.cluster<-parameter[19]
	
	

markov_probability_cluster<-matrix(0, nrow= number.of.positions.total, ncol=3)


likelihood_cell<-matrix(1, nrow= length(cluster[[counter.cluster]]), ncol=1)



for(i in 1: length(cluster[[counter.cluster]])){

					 index.i<-which(data_cell_motility[, 8]==id[cluster[[counter.cluster]][i]])
					 
				 
				     n.i<-length(index.i)
				 
					 
					 		
					 markov_probability_cluster[index.i[1], ]= prior_pi %*% P_1
					 	
					  	
					 
					 
					 for(l in 2:(length(index.i))){
					 	
					 	
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_1 %*% diag(c(dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_non_persistent, beta_non_persistent), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_transition, beta_transition), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_persistent, beta_persistent)))
					 	
					 	
					 	
					 	
					 	
					 } 
					 
					 
					 likelihood_cell[i]<-sum((markov_probability_cluster[index.i[n.i-1], ]))
					 
					 # print(likelihood_cell[i])
					 
					 
					}


              t<-sum(log(likelihood_cell[1:counter.cluster]))
              
              # t<-likelihood_cell[1:counter.cluster]
              
              return(t)
	
	
	
}


###################################################




P_1<-matrix(c(fit_ss$summary[4:12]), nrow=3, ncol=3, byrow=T)

# P_2<-matrix(c(fit_ss$summary[4:12]), nrow=3, ncol=3, byrow=T)

prior_pi<-c(fit_ss$summary[1:3])

beta_non_persistent<-fit_ss$summary[16]

alpha_non_persistent<-fit_ss$summary[15]

beta_persistent<-fit_ss$summary[14]

alpha_persistent<-fit_ss$summary[13]

alpha_transition<-beta_transition<-fit_ss$summary[17]

parameter.bayes<-cbind(prior_pi[1], prior_pi[2], prior_pi[3], P_1[1, 1], P_1[1, 2], P_1[1, 3], P_1[2, 1], P_1[2, 2], P_1[2, 3], P_1[3, 1], P_1[3, 2], P_1[3, 3], alpha_persistent, beta_persistent, alpha_non_persistent, beta_non_persistent, alpha_transition, beta_transition, cluster.c)


##########################################################################


aic.clusters.4.cluster[cluster.c]<-2*(length(parameter.bayes)-1)-2*log.likelihood.f(parameter.bayes)



print("aic.clusters.4.cluster")

print(aic.clusters.4.cluster)

sum(aic.clusters.4.cluster)


###########################################################################################################################################################################################


library(rstan)



# # # aic.clusters.4.cluster<-matrix(0, nrow=4)


aic.clusters.2.cluster<-matrix(0, nrow=2)



setwd("/global/scratch/torkashe/Clustering-Results/2-clusters-3-states-HMM-1/")

load(".RData")



cluster.c<-1



###################################################


log.likelihood.f<-function(parameter){
	
	
	
	
	
	prior_pi <-parameter[1:3]
	
	P_1<-matrix(parameter[4:12], nrow=3, ncol=3, byrow=T)
	
	# P_2<-matrix(parameter[13:21], nrow=3, ncol=3, byrow=T)
	
	alpha_persistent <-parameter[13]
	
	beta_persistent <-parameter[14]
	
	alpha_non_persistent <-parameter[15]
	
	beta_non_persistent <-parameter[16]
	
	alpha_transition<-parameter[17]
	
	beta_transition<-parameter[18]
	
	counter.cluster<-parameter[19]
	
	

markov_probability_cluster<-matrix(0, nrow= number.of.positions.total, ncol=3)


likelihood_cell<-matrix(1, nrow= length(cluster[[counter.cluster]]), ncol=1)



for(i in 1: length(cluster[[counter.cluster]])){

					 index.i<-which(data_cell_motility[, 8]==id[cluster[[counter.cluster]][i]])
					 
				 
				     n.i<-length(index.i)
				 
					 
					 		
					 markov_probability_cluster[index.i[1], ]= prior_pi %*% P_1
					 	
					  	
					 
					 
					 for(l in 2:(length(index.i))){
					 	
					 	
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_1 %*% diag(c(dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_non_persistent, beta_non_persistent), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_transition, beta_transition), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_persistent, beta_persistent)))
					 	
					 	
					 	
					 	
					 	
					 } 
					 
					 
					 likelihood_cell[i]<-sum((markov_probability_cluster[index.i[n.i-1], ]))
					 
					 # print(likelihood_cell[i])
					 
					 
					}


              t<-sum(log(likelihood_cell[1:counter.cluster]))
              
              # t<-likelihood_cell[1:counter.cluster]
              
              return(t)
	
	
	
}


###################################################




P_1<-matrix(c(fit_ss$summary[4:12]), nrow=3, ncol=3, byrow=T)

# P_2<-matrix(c(fit_ss$summary[4:12]), nrow=3, ncol=3, byrow=T)

prior_pi<-c(fit_ss$summary[1:3])

beta_non_persistent<-fit_ss$summary[16]

alpha_non_persistent<-fit_ss$summary[15]

beta_persistent<-fit_ss$summary[14]

alpha_persistent<-fit_ss$summary[13]

alpha_transition<-beta_transition<-fit_ss$summary[17]

parameter.bayes<-cbind(prior_pi[1], prior_pi[2], prior_pi[3], P_1[1, 1], P_1[1, 2], P_1[1, 3], P_1[2, 1], P_1[2, 2], P_1[2, 3], P_1[3, 1], P_1[3, 2], P_1[3, 3], alpha_persistent, beta_persistent, alpha_non_persistent, beta_non_persistent, alpha_transition, beta_transition, cluster.c)


##########################################################################


aic.clusters.2.cluster[cluster.c]<-2*(length(parameter.bayes)-1)-2*log.likelihood.f(parameter.bayes)



print("aic.clusters.2.cluster")

print(aic.clusters.2.cluster)


###########################################################################################################################################################################################



setwd("/global/scratch/torkashe/Clustering-Results/2-clusters-3-states-HMM-2/")

load(".RData")



cluster.c<-2



###################################################


log.likelihood.f<-function(parameter){
	
	
	
	
	
	prior_pi <-parameter[1:3]
	
	P_1<-matrix(parameter[4:12], nrow=3, ncol=3, byrow=T)
	
	# P_2<-matrix(parameter[13:21], nrow=3, ncol=3, byrow=T)
	
	alpha_persistent <-parameter[13]
	
	beta_persistent <-parameter[14]
	
	alpha_non_persistent <-parameter[15]
	
	beta_non_persistent <-parameter[16]
	
	alpha_transition<-parameter[17]
	
	beta_transition<-parameter[18]
	
	counter.cluster<-parameter[19]
	
	

markov_probability_cluster<-matrix(0, nrow= number.of.positions.total, ncol=3)


likelihood_cell<-matrix(1, nrow= length(cluster[[counter.cluster]]), ncol=1)



for(i in 1: length(cluster[[counter.cluster]])){

					 index.i<-which(data_cell_motility[, 8]==id[cluster[[counter.cluster]][i]])
					 
				 
				     n.i<-length(index.i)
				 
					 
					 		
					 markov_probability_cluster[index.i[1], ]= prior_pi %*% P_1
					 	
					  	
					 
					 
					 for(l in 2:(length(index.i))){
					 	
					 	
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_1 %*% diag(c(dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_non_persistent, beta_non_persistent), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_transition, beta_transition), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_persistent, beta_persistent)))
					 	
					 	
					 	
					 	
					 	
					 } 
					 
					 
					 likelihood_cell[i]<-sum((markov_probability_cluster[index.i[n.i-1], ]))
					 
					 # print(likelihood_cell[i])
					 
					 
					}


              t<-sum(log(likelihood_cell[1:counter.cluster]))
              
              # t<-likelihood_cell[1:counter.cluster]
              
              return(t)
	
	
	
}


###################################################




P_1<-matrix(c(fit_ss$summary[4:12]), nrow=3, ncol=3, byrow=T)

# P_2<-matrix(c(fit_ss$summary[4:12]), nrow=3, ncol=3, byrow=T)

prior_pi<-c(fit_ss$summary[1:3])

beta_non_persistent<-fit_ss$summary[16]

alpha_non_persistent<-fit_ss$summary[15]

beta_persistent<-fit_ss$summary[14]

alpha_persistent<-fit_ss$summary[13]

alpha_transition<-beta_transition<-fit_ss$summary[17]

parameter.bayes<-cbind(prior_pi[1], prior_pi[2], prior_pi[3], P_1[1, 1], P_1[1, 2], P_1[1, 3], P_1[2, 1], P_1[2, 2], P_1[2, 3], P_1[3, 1], P_1[3, 2], P_1[3, 3], alpha_persistent, beta_persistent, alpha_non_persistent, beta_non_persistent, alpha_transition, beta_transition, cluster.c)


##########################################################################


aic.clusters.2.cluster[cluster.c]<-2*(length(parameter.bayes)-1)-2*log.likelihood.f(parameter.bayes)



print("aic.clusters.2.cluster")

print(aic.clusters.2.cluster)


###########################################################################################################################################################################################


remove(list=ls())



setwd("/global/scratch/torkashe/Clustering-Results/Position-3-states-HMM/")

load(".RData")







###################################################


log.likelihood.f<-function(parameter){
	
	
	
	
	
	prior_pi <-parameter[1:3]
	
	P_1<-matrix(parameter[4:12], nrow=3, ncol=3, byrow=T)
	
	# P_2<-matrix(parameter[13:21], nrow=3, ncol=3, byrow=T)
	
	alpha_persistent <-parameter[13]
	
	beta_persistent <-parameter[14]
	
	alpha_non_persistent <-parameter[15]
	
	beta_non_persistent <-parameter[16]
	
	alpha_transition<-parameter[17]
	
	beta_transition<-parameter[18]
	
	# counter.cluster<-parameter[19]
	
	

markov_probability_cluster<-matrix(0, nrow= number.of.positions.total, ncol=3)


likelihood_cell<-matrix(1, nrow= number.of.cells, ncol=1)



for(i in 1: number.of.cells ){

					 index.i<-which(data_cell_motility[, 8]==id[i])
					 
				 
				     n.i<-length(index.i)
				 
					 
					 		
					 markov_probability_cluster[index.i[1], ]= prior_pi %*% P_1
					 	
					  	
					 
					 
					 for(l in 2:(length(index.i))){
					 	
					 	
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_1 %*% diag(c(dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_non_persistent, beta_non_persistent), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_transition, beta_transition), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_persistent, beta_persistent)))
					 	
					 	
					 	
					 	
					 	
					 } 
					 
					 
					 likelihood_cell[i]<-sum((markov_probability_cluster[index.i[n.i-1], ]))
					 
					 # print(likelihood_cell[i])
					 
					 
					}


              t<-sum(log(likelihood_cell[1:counter.cluster]))
              
              # t<-likelihood_cell[1:counter.cluster]
              
              return(t)
	
	
	
}


###################################################




P_1<-matrix(c(fit_ss$summary[4:12]), nrow=3, ncol=3, byrow=T)

# P_2<-matrix(c(fit_ss$summary[4:12]), nrow=3, ncol=3, byrow=T)

prior_pi<-c(fit_ss$summary[1:3])

beta_non_persistent<-fit_ss$summary[16]

alpha_non_persistent<-fit_ss$summary[15]

beta_persistent<-fit_ss$summary[14]

alpha_persistent<-fit_ss$summary[13]

alpha_transition<-beta_transition<-fit_ss$summary[17]

parameter.bayes<-cbind(prior_pi[1], prior_pi[2], prior_pi[3], P_1[1, 1], P_1[1, 2], P_1[1, 3], P_1[2, 1], P_1[2, 2], P_1[2, 3], P_1[3, 1], P_1[3, 2], P_1[3, 3], alpha_persistent, beta_persistent, alpha_non_persistent, beta_non_persistent, alpha_transition, beta_transition)


##########################################################################


aic.no.cluster<-2*(length(parameter.bayes))-2*log.likelihood.f(parameter.bayes)



print("aic.no.cluster")

print(aic.no.cluster)


###########################################################################################################################################################################################


