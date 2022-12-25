remove(list=ls())



library(rstan)



aic.clusters.4.cluster<-matrix(0, nrow=4)


aic.clusters.2.cluster<-matrix(0, nrow=2)



setwd("/global/scratch/torkashe/Clustering-Results/201-hyper-inf-2-states/")

load(".RData-22")



cluster.c<-1



###################################################


log.likelihood.f<-function(parameter){
	
	
	
	
	
	prior_pi <-parameter[1:2]
	
	P_1<-matrix(parameter[3:6], nrow=2, ncol=2, byrow=T)
	
	P_2<-matrix(parameter[7:10], nrow=2, ncol=2, byrow=T)
	
	alpha_persistent <-parameter[11]
	
	beta_persistent <-parameter[12]
	
	alpha_non_persistent <-parameter[13]
	
	beta_non_persistent <-parameter[14]
	
	
	counter.cluster<-parameter[15]
	
	

markov_probability_cluster<-matrix(0, nrow= number.of.positions.total, ncol=2)


likelihood_cell<-matrix(1, nrow= length(cluster[[counter.cluster]]), ncol=1)



for(i in 1: length(cluster[[counter.cluster]])){

					 index.i<-which(data_cell_motility[, 8]==id[cluster[[counter.cluster]][i]])
					 
				 
				     n.i<-length(index.i)
				 
					 
					 if(number.of.neighbors.full.data[index.i[1]]==0){
					 		
					 		markov_probability_cluster[index.i[1], ]= prior_pi %*% P_1
					 	}
					 	
					 	
					  if(number.of.neighbors.full.data[index.i[1]]>0){
					 		
					 		markov_probability_cluster[index.i[1], ]= prior_pi %*% P_2
					 		
					 	}
					 	
					  	
					 
					 
					 for(l in 2:(length(index.i))){
					 	
					 	
					 	if(number.of.neighbors.full.data[index.i[l]]==0){
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_1 %*% diag(c(dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_non_persistent, beta_non_persistent), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_persistent, beta_persistent)))
					 	}
					 	
					 	
					 	if(number.of.neighbors.full.data[index.i[l]]>0){
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_2 %*% diag(c(dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_non_persistent, beta_non_persistent), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_persistent, beta_persistent)))
					 	}
					 	
					 	
					 	
					 	
					 } 
					 
					 
					 
					 likelihood_cell[i]<-sum((markov_probability_cluster[index.i[n.i-1], ]))
					 
					 
					 # print(likelihood_cell[i])
					 
					 
					}


              t<-sum(log(likelihood_cell[1:counter.cluster]))
              
              # t<-likelihood_cell[1:counter.cluster]
              
              return(t)
	
	
	
}


###################################################




P_1<-matrix(c(fit_ss$summary[3:6]), nrow=2, ncol=2, byrow=T)

P_2<-matrix(c(fit_ss$summary[7:10]), nrow=2, ncol=2, byrow=T)

prior_pi<-c(fit_ss$summary[1], fit_ss$summary[2])

beta_non_persistent<-fit_ss$summary[14]

alpha_non_persistent<-fit_ss$summary[13]

beta_persistent<-fit_ss$summary[12]

alpha_persistent<-fit_ss$summary[11]



parameter.bayes<-cbind(prior_pi, P_1[1, 1], P_1[1, 2], P_1[2, 1], P_1[2, 2], P_2[1, 1], P_2[1, 2], P_2[2, 1], P_2[2, 2], alpha_persistent, beta_persistent, alpha_non_persistent, beta_non_persistent, cluster.c)


##########################################################################


aic.clusters.4.cluster[cluster.c] <- 2*(length(parameter.bayes)-1) - 2*log.likelihood.f(parameter.bayes)



print("aic.clusters.4.cluster")

print(aic.clusters.4.cluster)


###########################################################################################################################################################################################

library(rstan)


setwd("/global/scratch/torkashe/Clustering-Results/202-hyper-inf-2-states/")

load(".RData-22")



cluster.c<-2

###################################################


log.likelihood.f<-function(parameter){
	
	
	
	
	
	prior_pi <-parameter[1:2]
	
	P_1<-matrix(parameter[3:6], nrow=2, ncol=2, byrow=T)
	
	P_2<-matrix(parameter[7:10], nrow=2, ncol=2, byrow=T)
	
	alpha_persistent <-parameter[11]
	
	beta_persistent <-parameter[12]
	
	alpha_non_persistent <-parameter[13]
	
	beta_non_persistent <-parameter[14]
	
	
	counter.cluster<-parameter[15]
	
	

markov_probability_cluster<-matrix(0, nrow= number.of.positions.total, ncol=2)


likelihood_cell<-matrix(1, nrow= length(cluster[[counter.cluster]]), ncol=1)



for(i in 1: length(cluster[[counter.cluster]])){

					 index.i<-which(data_cell_motility[, 8]==id[cluster[[counter.cluster]][i]])
					 
				 
				     n.i<-length(index.i)
				 
					 
					 if(number.of.neighbors.full.data[index.i[1]]==0){
					 		
					 		markov_probability_cluster[index.i[1], ]= prior_pi %*% P_1
					 	}
					 	
					 	
					  if(number.of.neighbors.full.data[index.i[1]]>0){
					 		
					 		markov_probability_cluster[index.i[1], ]= prior_pi %*% P_2
					 		
					 	}
					 	
					  	
					 
					 
					 for(l in 2:(length(index.i))){
					 	
					 	
					 	if(number.of.neighbors.full.data[index.i[l]]==0){
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_1 %*% diag(c(dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_non_persistent, beta_non_persistent), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_persistent, beta_persistent)))
					 	}
					 	
					 	
					 	if(number.of.neighbors.full.data[index.i[l]]>0){
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_2 %*% diag(c(dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_non_persistent, beta_non_persistent), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_persistent, beta_persistent)))
					 	}
					 	
					 	
					 	
					 	
					 } 
					 
					 
					 likelihood_cell[i]<-sum((markov_probability_cluster[index.i[n.i-1], ]))
					 
					 # print(likelihood_cell[i])
					 
					 
					}


              t<-sum(log(likelihood_cell[1:counter.cluster]))
              
              # t<-likelihood_cell[1:counter.cluster]
              
              return(t)
	
	
	
}


###################################################




P_1<-matrix(c(fit_ss$summary[3:6]), nrow=2, ncol=2, byrow=T)

P_2<-matrix(c(fit_ss$summary[7:10]), nrow=2, ncol=2, byrow=T)

prior_pi<-c(fit_ss$summary[1], fit_ss$summary[2])

beta_non_persistent<-fit_ss$summary[14]

alpha_non_persistent<-fit_ss$summary[13]

beta_persistent<-fit_ss$summary[12]

alpha_persistent<-fit_ss$summary[11]



parameter.bayes<-cbind(prior_pi, P_1[1, 1], P_1[1, 2], P_1[2, 1], P_1[2, 2], P_2[1, 1], P_2[1, 2], P_2[2, 1], P_2[2, 2], alpha_persistent, beta_persistent, alpha_non_persistent, beta_non_persistent, cluster.c)


##########################################################################



aic.clusters.4.cluster[cluster.c]<-2*(length(parameter.bayes)-1)-2*log.likelihood.f(parameter.bayes)



print("aic.clusters.4.cluster")

print(aic.clusters.4.cluster)


###########################################################################################################################################################################################


library(rstan)


setwd("/global/scratch/torkashe/Clustering-Results/203-hyper-inf-2-states/")

load(".RData-22")



cluster.c<-3

###################################################


log.likelihood.f<-function(parameter){
	
	
	
	
	
	prior_pi <-parameter[1:2]
	
	P_1<-matrix(parameter[3:6], nrow=2, ncol=2, byrow=T)
	
	P_2<-matrix(parameter[7:10], nrow=2, ncol=2, byrow=T)
	
	alpha_persistent <-parameter[11]
	
	beta_persistent <-parameter[12]
	
	alpha_non_persistent <-parameter[13]
	
	beta_non_persistent <-parameter[14]
	
	
	counter.cluster<-parameter[15]
	
	

markov_probability_cluster<-matrix(0, nrow= number.of.positions.total, ncol=2)


likelihood_cell<-matrix(1, nrow= length(cluster[[counter.cluster]]), ncol=1)



for(i in 1: length(cluster[[counter.cluster]])){

					 index.i<-which(data_cell_motility[, 8]==id[cluster[[counter.cluster]][i]])
					 
				 
				     n.i<-length(index.i)
				 
					 
					 if(number.of.neighbors.full.data[index.i[1]]==0){
					 		
					 		markov_probability_cluster[index.i[1], ]= prior_pi %*% P_1
					 	}
					 	
					 	
					  if(number.of.neighbors.full.data[index.i[1]]>0){
					 		
					 		markov_probability_cluster[index.i[1], ]= prior_pi %*% P_2
					 		
					 	}
					 	
					  	
					 
					 
					 for(l in 2:(length(index.i))){
					 	
					 	
					 	if(number.of.neighbors.full.data[index.i[l]]==0){
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_1 %*% diag(c(dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_non_persistent, beta_non_persistent), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_persistent, beta_persistent)))
					 	}
					 	
					 	
					 	if(number.of.neighbors.full.data[index.i[l]]>0){
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_2 %*% diag(c(dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_non_persistent, beta_non_persistent), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_persistent, beta_persistent)))
					 	}
					 	
					 	
					 	
					 	
					 } 
					 
					 
					 likelihood_cell[i]<-sum((markov_probability_cluster[index.i[n.i-1], ]))
					 
					 # print(likelihood_cell[i])
					 
					 
					}


              t<-sum(log(likelihood_cell[1:counter.cluster]))
              
              # t<-likelihood_cell[1:counter.cluster]
              
              return(t)
	
	
	
}


###################################################




P_1<-matrix(c(fit_ss$summary[3:6]), nrow=2, ncol=2, byrow=T)

P_2<-matrix(c(fit_ss$summary[7:10]), nrow=2, ncol=2, byrow=T)

prior_pi<-c(fit_ss$summary[1], fit_ss$summary[2])

beta_non_persistent<-fit_ss$summary[14]

alpha_non_persistent<-fit_ss$summary[13]

beta_persistent<-fit_ss$summary[12]

alpha_persistent<-fit_ss$summary[11]



parameter.bayes<-cbind(prior_pi, P_1[1, 1], P_1[1, 2], P_1[2, 1], P_1[2, 2], P_2[1, 1], P_2[1, 2], P_2[2, 1], P_2[2, 2], alpha_persistent, beta_persistent, alpha_non_persistent, beta_non_persistent, cluster.c)


##########################################################################



aic.clusters.4.cluster[cluster.c]<-2*(length(parameter.bayes)-1)-2*log.likelihood.f(parameter.bayes)



print("aic.clusters.4.cluster")

print(aic.clusters.4.cluster)



###########################################################################################################################################################################################


library(rstan)



setwd("/global/scratch/torkashe/Clustering-Results/204-hyper-inf-2-states/")

load(".RData-22")



cluster.c<-4

###################################################


log.likelihood.f<-function(parameter){
	
	
	
	
	
	prior_pi <-parameter[1:2]
	
	P_1<-matrix(parameter[3:6], nrow=2, ncol=2, byrow=T)
	
	P_2<-matrix(parameter[7:10], nrow=2, ncol=2, byrow=T)
	
	alpha_persistent <-parameter[11]
	
	beta_persistent <-parameter[12]
	
	alpha_non_persistent <-parameter[13]
	
	beta_non_persistent <-parameter[14]
	
	
	counter.cluster<-parameter[15]
	
	

markov_probability_cluster<-matrix(0, nrow= number.of.positions.total, ncol=2)


likelihood_cell<-matrix(1, nrow= length(cluster[[counter.cluster]]), ncol=1)



for(i in 1: length(cluster[[counter.cluster]])){

					 index.i<-which(data_cell_motility[, 8]==id[cluster[[counter.cluster]][i]])
					 
				 
				     n.i<-length(index.i)
				 
					 
					 if(number.of.neighbors.full.data[index.i[1]]==0){
					 		
					 		markov_probability_cluster[index.i[1], ]= prior_pi %*% P_1
					 	}
					 	
					 	
					  if(number.of.neighbors.full.data[index.i[1]]>0){
					 		
					 		markov_probability_cluster[index.i[1], ]= prior_pi %*% P_2
					 		
					 	}
					 	
					  	
					 
					 
					 for(l in 2:(length(index.i))){
					 	
					 	
					 	if(number.of.neighbors.full.data[index.i[l]]==0){
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_1 %*% diag(c(dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_non_persistent, beta_non_persistent), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_persistent, beta_persistent)))
					 	}
					 	
					 	
					 	if(number.of.neighbors.full.data[index.i[l]]>0){
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_2 %*% diag(c(dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_non_persistent, beta_non_persistent), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_persistent, beta_persistent)))
					 	}
					 	
					 	
					 	
					 	
					 } 
					 
					 
					 likelihood_cell[i]<-sum((markov_probability_cluster[index.i[n.i-1], ]))
					 
					 # print(likelihood_cell[i])
					 
					 
					}


              t<-sum(log(likelihood_cell[1:counter.cluster]))
              
              # t<-likelihood_cell[1:counter.cluster]
              
              return(t)
	
	
	
}


###################################################




P_1<-matrix(c(fit_ss$summary[3:6]), nrow=2, ncol=2, byrow=T)

P_2<-matrix(c(fit_ss$summary[7:10]), nrow=2, ncol=2, byrow=T)

prior_pi<-c(fit_ss$summary[1], fit_ss$summary[2])

beta_non_persistent<-fit_ss$summary[14]

alpha_non_persistent<-fit_ss$summary[13]

beta_persistent<-fit_ss$summary[12]

alpha_persistent<-fit_ss$summary[11]



parameter.bayes<-cbind(prior_pi, P_1[1, 1], P_1[1, 2], P_1[2, 1], P_1[2, 2], P_2[1, 1], P_2[1, 2], P_2[2, 1], P_2[2, 2], alpha_persistent, beta_persistent, alpha_non_persistent, beta_non_persistent, cluster.c)


##########################################################################



aic.clusters.4.cluster[cluster.c]<-2*(length(parameter.bayes)-1)-2*log.likelihood.f(parameter.bayes)



print("aic.clusters.4.cluster")

print(aic.clusters.4.cluster)

sum(aic.clusters.4.cluster)


###########################################################################################################################################################################################


library(rstan)



setwd("/global/scratch/torkashe/Clustering-Results/Position-2-states-NHMM/")

load(".RData-22")




###################################################


log.likelihood.f<-function(parameter){
	
	
	
	
	
	prior_pi <-parameter[1:2]
	
	P_1<-matrix(parameter[3:6], nrow=2, ncol=2, byrow=T)
	
	P_2<-matrix(parameter[7:10], nrow=2, ncol=2, byrow=T)
	
	alpha_persistent <-parameter[11]
	
	beta_persistent <-parameter[12]
	
	alpha_non_persistent <-parameter[13]
	
	beta_non_persistent <-parameter[14]
	
	
	# counter.cluster<-parameter[15]
	
	

markov_probability_cluster<-matrix(0, nrow= number.of.positions.total, ncol=2)


likelihood_cell<-matrix(1, nrow= number.of.cells, ncol=1)



for(i in 1: number.of.cells){

					 index.i<-which(data_cell_motility[, 8]==id[i])
					 
				 
				     n.i<-length(index.i)
				 
					 
					 if(number.of.neighbors.full.data[index.i[1]]==0){
					 		
					 		markov_probability_cluster[index.i[1], ]= prior_pi %*% P_1
					 	}
					 	
					 	
					  if(number.of.neighbors.full.data[index.i[1]]>0){
					 		
					 		markov_probability_cluster[index.i[1], ]= prior_pi %*% P_2
					 		
					 	}
					 	
					  	
					 
					 
					 for(l in 2:(length(index.i))){
					 	
					 	
					 	if(number.of.neighbors.full.data[index.i[l]]==0){
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_1 %*% diag(c(dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_non_persistent, beta_non_persistent), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_persistent, beta_persistent)))
					 	}
					 	
					 	
					 	if(number.of.neighbors.full.data[index.i[l]]>0){
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_2 %*% diag(c(dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_non_persistent, beta_non_persistent), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_persistent, beta_persistent)))
					 	}
					 	
					 	
					 	
					 	
					 } 
					 
					 
					 likelihood_cell[i]<-sum((markov_probability_cluster[index.i[n.i-1], ]))
					 
					 # print(likelihood_cell[i])
					 
					 
					}


              t<-sum(log(likelihood_cell[1: number.of.cells]))
              
              # t<-likelihood_cell[1:counter.cluster]
              
              return(t)
	
	
	
}


###################################################




P_1<-matrix(c(fit_ss$summary[3:6]), nrow=2, ncol=2, byrow=T)

P_2<-matrix(c(fit_ss$summary[7:10]), nrow=2, ncol=2, byrow=T)

prior_pi<-c(fit_ss$summary[1], fit_ss$summary[2])

beta_non_persistent<-fit_ss$summary[14]

alpha_non_persistent<-fit_ss$summary[13]

beta_persistent<-fit_ss$summary[12]

alpha_persistent<-fit_ss$summary[11]



parameter.bayes<-cbind(prior_pi, P_1[1, 1], P_1[1, 2], P_1[2, 1], P_1[2, 2], P_2[1, 1], P_2[1, 2], P_2[2, 1], P_2[2, 2], alpha_persistent, beta_persistent, alpha_non_persistent, beta_non_persistent)


##########################################################################

aic.no.cluster<-2*(length(parameter.bayes))-2*log.likelihood.f(parameter.bayes)


print("aic.no.cluster")

print(aic.no.cluster)


###########################################################################################################################################################################################





setwd("/global/scratch/torkashe/Clustering-Results/2-clusters-NHMM-1/")

load(".RData-22")



cluster.c<-1


###################################################


log.likelihood.f<-function(parameter){
	
	
	
	
	
	prior_pi <-parameter[1:2]
	
	P_1<-matrix(parameter[3:6], nrow=2, ncol=2, byrow=T)
	
	P_2<-matrix(parameter[7:10], nrow=2, ncol=2, byrow=T)
	
	alpha_persistent <-parameter[11]
	
	beta_persistent <-parameter[12]
	
	alpha_non_persistent <-parameter[13]
	
	beta_non_persistent <-parameter[14]
	
	
	counter.cluster<-parameter[15]
	
	

markov_probability_cluster<-matrix(0, nrow= number.of.positions.total, ncol=2)


likelihood_cell<-matrix(1, nrow= length(cluster[[counter.cluster]]), ncol=1)



for(i in 1: length(cluster[[counter.cluster]])){

					 index.i<-which(data_cell_motility[, 8]==id[cluster[[counter.cluster]][i]])
					 
				 
				     n.i<-length(index.i)
				 
					 
					 if(number.of.neighbors.full.data[index.i[1]]==0){
					 		
					 		markov_probability_cluster[index.i[1], ]= prior_pi %*% P_1
					 	}
					 	
					 	
					  if(number.of.neighbors.full.data[index.i[1]]>0){
					 		
					 		markov_probability_cluster[index.i[1], ]= prior_pi %*% P_2
					 		
					 	}
					 	
					  	
					 
					 
					 for(l in 2:(length(index.i))){
					 	
					 	
					 	if(number.of.neighbors.full.data[index.i[l]]==0){
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_1 %*% diag(c(dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_non_persistent, beta_non_persistent), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_persistent, beta_persistent)))
					 	}
					 	
					 	
					 	if(number.of.neighbors.full.data[index.i[l]]>0){
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_2 %*% diag(c(dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_non_persistent, beta_non_persistent), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_persistent, beta_persistent)))
					 	}
					 	
					 	
					 	
					 	
					 } 
					 
					 
					 likelihood_cell[i]<-sum((markov_probability_cluster[index.i[n.i-1], ]))
					 
					 # print(likelihood_cell[i])
					 
					 
					}


              t<-sum(log(likelihood_cell[1:counter.cluster]))
              
              # t<-likelihood_cell[1:counter.cluster]
              
              return(t)
	
	
	
}


###################################################




P_1<-matrix(c(fit_ss$summary[3:6]), nrow=2, ncol=2, byrow=T)

P_2<-matrix(c(fit_ss$summary[7:10]), nrow=2, ncol=2, byrow=T)

prior_pi<-c(fit_ss$summary[1], fit_ss$summary[2])

beta_non_persistent<-fit_ss$summary[14]

alpha_non_persistent<-fit_ss$summary[13]

beta_persistent<-fit_ss$summary[12]

alpha_persistent<-fit_ss$summary[11]



parameter.bayes<-cbind(prior_pi, P_1[1, 1], P_1[1, 2], P_1[2, 1], P_1[2, 2], P_2[1, 1], P_2[1, 2], P_2[2, 1], P_2[2, 2], alpha_persistent, beta_persistent, alpha_non_persistent, beta_non_persistent, cluster.c)


##########################################################################



aic.clusters.2.cluster[cluster.c]<-2*(length(parameter.bayes)-1)-2*log.likelihood.f(parameter.bayes)



print("aic.clusters.2.cluster")

print(aic.clusters.2.cluster)



###########################################################################################################################################################################################





setwd("/global/scratch/torkashe/Clustering-Results/2-clusters-NHMM-2/")

load(".RData-22")



cluster.c<-2



###################################################


log.likelihood.f<-function(parameter){
	
	
	
	
	
	prior_pi <-parameter[1:2]
	
	P_1<-matrix(parameter[3:6], nrow=2, ncol=2, byrow=T)
	
	P_2<-matrix(parameter[7:10], nrow=2, ncol=2, byrow=T)
	
	alpha_persistent <-parameter[11]
	
	beta_persistent <-parameter[12]
	
	alpha_non_persistent <-parameter[13]
	
	beta_non_persistent <-parameter[14]
	
	
	counter.cluster<-parameter[15]
	
	

markov_probability_cluster<-matrix(0, nrow= number.of.positions.total, ncol=2)


likelihood_cell<-matrix(1, nrow= length(cluster[[counter.cluster]]), ncol=1)



for(i in 1: length(cluster[[counter.cluster]])){

					 index.i<-which(data_cell_motility[, 8]==id[cluster[[counter.cluster]][i]])
					 
				 
				     n.i<-length(index.i)
				 
					 
					 if(number.of.neighbors.full.data[index.i[1]]==0){
					 		
					 		markov_probability_cluster[index.i[1], ]= prior_pi %*% P_1
					 	}
					 	
					 	
					  if(number.of.neighbors.full.data[index.i[1]]>0){
					 		
					 		markov_probability_cluster[index.i[1], ]= prior_pi %*% P_2
					 		
					 	}
					 	
					  	
					 
					 
					 for(l in 2:(length(index.i))){
					 	
					 	
					 	if(number.of.neighbors.full.data[index.i[l]]==0){
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_1 %*% diag(c(dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_non_persistent, beta_non_persistent), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_persistent, beta_persistent)))
					 	}
					 	
					 	
					 	if(number.of.neighbors.full.data[index.i[l]]>0){
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_2 %*% diag(c(dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_non_persistent, beta_non_persistent), dbeta((angular.change.complete.full.data[index.i[l]]/pi), alpha_persistent, beta_persistent)))
					 	}
					 	
					 	
					 	
					 	
					 } 
					 
					 
					 likelihood_cell[i]<-sum((markov_probability_cluster[index.i[n.i-1], ]))
					 
					 # print(likelihood_cell[i])
					 
					 
					}


              t<-sum(log(likelihood_cell[1:counter.cluster]))
              
              # t<-likelihood_cell[1:counter.cluster]
              
              return(t)
	
	
	
}


###################################################




P_1<-matrix(c(fit_ss$summary[3:6]), nrow=2, ncol=2, byrow=T)

P_2<-matrix(c(fit_ss$summary[7:10]), nrow=2, ncol=2, byrow=T)

prior_pi<-c(fit_ss$summary[1], fit_ss$summary[2])

beta_non_persistent<-fit_ss$summary[14]

alpha_non_persistent<-fit_ss$summary[13]

beta_persistent<-fit_ss$summary[12]

alpha_persistent<-fit_ss$summary[11]



parameter.bayes<-cbind(prior_pi, P_1[1, 1], P_1[1, 2], P_1[2, 1], P_1[2, 2], P_2[1, 1], P_2[1, 2], P_2[2, 1], P_2[2, 2], alpha_persistent, beta_persistent, alpha_non_persistent, beta_non_persistent, cluster.c)


##########################################################################



aic.clusters.2.cluster[cluster.c]<-2*(length(parameter.bayes)-1)-2*log.likelihood.f(parameter.bayes)



print("aic.clusters.2.cluster")

print(aic.clusters.2.cluster)

sum(aic.clusters.2.cluster)