remove(list=ls())

### Find the derived code

### The average computed log pointwise predictive density for a dataset

### BIC

### Checking residual fit for radius

### Checking markovian assumption

 
#####  The equilibirium property should be calculated at fixed time

equilibrium.assumption.function<-function(state.i, state.j, P.matrix.1, P.matrix.2, state.i.and.neighbor.zero, state.j.and.neighbor.zero, state.i.and.neighbor.at.least.one, state.j.and.neighbor.at.least.one){
	
	
	# # # p.i.given.j[i-1]<-(P.matrix.1[state.j, state.i]*state.j.and.neighbor.zero+P.matrix.2[state.j, state.i]*state.j.and.neighbor.at.least.one)/(states.probability[i, state.j])
	
	# # # p.j.given.i[i-1]<-(P.matrix.1[state.i, state.j]*state.i.and.neighbor.zero+P.matrix.2[state.i, state.j]*state.i.and.neighbor.at.least.one)/(states.probability[i, state.i])
	
	p.i.given.j<-(P.matrix.1[state.j, state.i]*state.j.and.neighbor.zero+P.matrix.2[state.j, state.i]*state.j.and.neighbor.at.least.one)
	
	p.j.given.i<-(P.matrix.1[state.i, state.j]*state.i.and.neighbor.zero+P.matrix.2[state.i, state.j]*state.i.and.neighbor.at.least.one)
	
	equilibrium.difference<-p.i.given.j-p.j.given.i
	
		
	return(equilibrium.difference)
	
}





equilibrium.of.states<-function(P.matrix.1, P.matrix.2, ratio.of.neighbors.per.states.matrix, number.of.states, states.probability){
	
	t<-matrix(0, nrow= number.of.states, ncol= number.of.states)
	###  ratio.of.neighbors.per.states.matrix shows the ratio without and with neighbors in each state.
	for(i in 1: number.of.states){
		
		for(j in 1: number.of.states){
			
			if(i!=j){
				###  For now, I run it in a loop to check the performance. Later we can reduce the loop as some parts are repititions.
				t[i, j]<-markovian.assumption.function(state.i=i, state.j=j, P.matrix.1, P.matrix.2, state.i.and.neighbor.zero=ratio.of.neighbors.per.states.matrix[i, 1], state.j.and.neighbor.zero=ratio.of.neighbors.per.states.matrix[j, 1], state.i.and.neighbor.at.least.one=ratio.of.neighbors.per.states.matrix[i, 2], state.j.and.neighbor.at.least.one= ratio.of.neighbors.per.states.matrix[j, 2], states.probability)
				
			}
			
		}
		
	}
	
	
} 


### Checking the predictive density assumption

predictive.density.new.data.point<-function(number.of.states, weight, alpha.parameter, beta.parameter, mu.parameter, variance.parameter, new.data.angle, new.data.radius){
	
	model.component<-NULL

	for(k in 1: number.of.states){
		
	     model.component[k]<- weight[k]*dbeta((new.data.angle/pi), alpha.parameter[k], beta.parameter[k])*dnorm(log(new.data.radius), mu.parameter[k], variance.parameter[k])	
		
		
	}
	
	model<-sum(model.component)
	
	return(model)
	
	
}



############################################## First Position-NHMM


setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/Position-2-states-HMM/")


load(file=paste(".RData", sep=""))

library(rstan)


fit_simulated <- extract(fit, permuted = TRUE)


R<-(length(fit_simulated[[1]])/2)




P_1<-matrix(c(fit_ss$summary[3:6]), nrow=2, ncol=2, byrow=T)

prior_pi<-c(fit_ss$summary[1], fit_ss$summary[2])



beta_persistent<-fit_ss$summary[8]

alpha_persistent<-fit_ss$summary[7]

beta_non_persistent<-fit_ss$summary[10]

alpha_non_persistent<-fit_ss$summary[9]



markov_probability_cluster<-matrix(0, nrow= number.of.positions.total, ncol=2)

for(i in 1: number.of.cells){

					 index.i<-which(data_cell_motility[, 8]==id[i])
					 
					 markov_probability_cluster[index.i[1], ]<-prior_pi
					 
					 for(l in 2:(length(index.i))){
					 	
					 	
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_1
					 	
					 	
					 	
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

for(i in 1: number.of.positions.total){
	
	
	if(group.probability[i, 1]>=group.probability[i, 2]){classification.matrix[i, 2]<-1}
	
	
	if(group.probability[i, 1]<group.probability[i, 2]){classification.matrix[i, 2]<-2}
	
	
}

weight<-NULL

weight[1]<-mean(classification.matrix[index.angular.change, 2]==1)

weight[2]<-mean(classification.matrix[index.angular.change, 2]==2)

### 1:black <- non-persistent

### 2:darkred <- persistent


index.persistent<-which(classification.matrix[index.angular.change, 2]==2)

index.non.persistent<-which(classification.matrix[index.angular.change, 2]==1)


### For now, as the program includes the radius devided by max for others.

# # # setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/cluster-2-hyper-informative-second-prior-2-states/")
# # # load(file=paste(".RData", sep=""))



persist<-radius.complete.full.data[index.angular.change[index.persistent]]

non.persist<-radius.complete.full.data[index.angular.change[index.non.persistent]]


######   checking the markovian assumption-non-persistent mode

min.track.length<-min(track.length)


###  note that we do not have angular change for the first position



probability.of.states<-matrix(0, nrow=(min.track.length-1), ncol=2)

probability.of.states[1, ]<-prior_pi




for(t in 2:(min.track.length-1)){
	
	
	probability.of.states[t, ]<-probability.of.states[(t-1), ] %*% P_1
	
}




p_non_per<-probability.of.states[2:(min.track.length-1), 1]


p_per<-probability.of.states[2:(min.track.length-1), 2]




### if we consider the same probability for neighbors and states through the time

### What I obtain here is for all cells. The geometric probability should be modified to show the same concept.



repitition<-21

empirical.probabilities.non.persistent<-matrix(0,  nrow= repitition)

waiting.time.matrix.non.persistent <-matrix(0, nrow= repitition)

for(i in 1:number.of.cells){
	
	index.i<-which(data_cell_motility[, 8]==id[i])
	
	for(k in 1: repitition){
		
	  for(j in 1:(length(index.i)-k)){

               if(j>=2){ ####  because we already removed the first position.

               if((sum(classification.matrix[index.i[j:(j+k-1)], 2]==2)==k) & classification.matrix[index.i[j-1], 2]!=2 & classification.matrix[index.i[(j+k)], 2]==1){
              
                  waiting.time.matrix.non.persistent[k]<-waiting.time.matrix.non.persistent[k]+1
                  
}

}

}
 
} 
 
}





waiting.time.matrix.non.persistent.zero<-0

for(i in 1:number.of.cells){
	
	index.i<-which(data_cell_motility[, 8]==id[i])
	
	for(j in 2:length(index.i)){

               ####  because we already removed the first position.

               if((classification.matrix[index.i[j], 2]==1) & classification.matrix[index.i[j-1], 2]!=2){
              
                  waiting.time.matrix.non.persistent.zero<-waiting.time.matrix.non.persistent.zero+1
                  
}


 
} 
 
}




freq.non.persistent<-c(waiting.time.matrix.non.persistent.zero, waiting.time.matrix.non.persistent)


for(k in 1: (repitition+1)){
	
	 empirical.probabilities.non.persistent[k]<-(sum(freq.non.persistent[1:k]))          
	
}


empirical.probabilities.non.persistent<-empirical.probabilities.non.persistent/sum(freq.non.persistent)


empirical.plot<-c(0, 0, empirical.probabilities.non.persistent)


x<-seq(1, 21, 1)

u<-matrix(0, nrow=number.of.cells, ncol=length(x))

u.mean<-NULL


par(mfrow=c(1, 2), lwd=3)



xlab.values.axes<-seq(-1, (repitition), 1)


xlab.values<-c(-1, 0, seq(0, (repitition), 1))



plot(xlab.values , empirical.plot, xlab="Trials",  main="(a)",  ylab="CDF", ylim=c(0, 1), type="l", xlim=c(-1, 21), lwd=0.7, xaxt = "n", cex.lab=1.4, cex.main=2, cex.axis=1.5)


axis(1, at=0:(repitition), labels= xlab.values.axes[2:(repitition+2)])


# # # # legend("bottomright", 
  # # # # legend = c("No Cluster", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"), 
  # # # # col = c("black", "darkgreen", "darkblue", "darkred", "darkorange"), 
  # # # # pch = c(20, 13, 17, 18, 23), 
  # # # # pt.cex = 2, 
  # # # # cex = 1.2, 
  # # # # text.col = "black", 
  # # # # horiz = F , 
  # # # # inset = c(0.1, 0.1))

  # # # # # lty=c(1,1, 1, 1, 1),




 
v<-NULL

v[1]<-mean(p_non_per)

for(t in 1: (repitition+2)){
	
	v[t]<-((1-v[1])^(t-1))*v[1]
  
 }
 
 
prob.estimated.non<-matrix(0, nrow=(min.track.length-1), ncol=1)

for(t in 1: (repitition+2)){
	
	
	prob.estimated.non[t]<-sum(v[1:t])
	
	
}
 


xlab.values<-seq(0, (repitition), 1)



points(xlab.values, prob.estimated.non[1:22], type="p", pch=20)





setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/cluster-1-2-states-HMM/")


load(file=paste(".RData", sep=""))

library(rstan)


# # # # fit_simulated <- extract(fit, permuted = TRUE)


# # R<-(length(fit_simulated[[1]])/2)




P_1<-matrix(c(fit_ss$summary[3:6]), nrow=2, ncol=2, byrow=T)

prior_pi<-c(fit_ss$summary[1], fit_ss$summary[2])



beta_persistent<-fit_ss$summary[8]

alpha_persistent<-fit_ss$summary[7]

beta_non_persistent<-fit_ss$summary[10]

alpha_non_persistent<-fit_ss$summary[9]



markov_probability_cluster<-matrix(0, nrow= number.of.positions.total, ncol=2)


for(i in 1: number.of.cells){

					 index.i<-which(data_cell_motility[, 8]==id[i])
					 
					 markov_probability_cluster[index.i[1], ]<-prior_pi
					 
					 for(l in 2:(length(index.i))){
					 	
					 	
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_1
					 	
					 	
					 	
					 } 
				}
######  Defining the ratio for the purpose of classification

t<-which(angular.change.complete.full.data==0)

angular.change.complete.full.data[t]<-0.00001

classification.matrix<-matrix(c(angular.change.complete.full.data, rep(0, length(angular.change.complete.full.data))), nrow=length(angular.change.complete.full.data), ncol=2)

group.probability<-matrix(0, nrow=length(angular.change.complete.full.data), ncol=2)


group.probability[, 1]<-markov_probability_cluster[, 1]*dbeta((angular.change.complete.full.data/pi), alpha_non_persistent, beta_non_persistent)/(markov_probability_cluster[, 1]*dbeta((angular.change.complete.full.data/pi), alpha_non_persistent, beta_non_persistent) +markov_probability_cluster[, 2]*dbeta((angular.change.complete.full.data/pi), alpha_persistent, beta_persistent))


group.probability[, 2]<-markov_probability_cluster[, 2]*dbeta((angular.change.complete.full.data/pi), alpha_persistent, beta_persistent)/(markov_probability_cluster[, 1]*dbeta((angular.change.complete.full.data/pi), alpha_non_persistent, beta_non_persistent) +markov_probability_cluster[, 2]*dbeta((angular.change.complete.full.data/pi), alpha_persistent, beta_persistent))



sum(group.probability[index.angular.change.cluster.1, 1]==group.probability[index.angular.change.cluster.1, 2])

for(i in 1: number.of.positions.total){
	
	
	if(group.probability[i, 1]>=group.probability[i, 2]){classification.matrix[i, 2]<-1}
	
	
	if(group.probability[i, 1]<group.probability[i, 2]){classification.matrix[i, 2]<-2}
	
	
}

weight<-NULL

weight[1]<-mean(classification.matrix[index.angular.change.cluster.1, 2]==1)

weight[2]<-mean(classification.matrix[index.angular.change.cluster.1, 2]==2)

### 1:black <- non-persistent

### 2:darkred <- persistent


index.persistent<-which(classification.matrix[index.angular.change.cluster.1, 2]==2)

index.non.persistent<-which(classification.matrix[index.angular.change.cluster.1, 2]==1)


### For now, as the program includes the radius devided by max for others.

# # # setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/cluster-2-hyper-informative-second-prior-2-states/")
# # # load(file=paste(".RData", sep=""))



persist<-radius.complete.full.data[index.angular.change.cluster.1[index.persistent]]

non.persist<-radius.complete.full.data[index.angular.change.cluster.1[index.non.persistent]]


######   checking the markovian assumption-non-persistent mode

min.track.length<-min(track.length)




####  note that we do not have angular change for the first position



probability.of.states<-matrix(0, nrow=(min.track.length-1), ncol=2)

probability.of.states[1, ]<-prior_pi




for(t in 2:(min.track.length-1)){
	
	
	probability.of.states[t, ]<-probability.of.states[(t-1), ] %*% P_1
	
}




p_non_per<-probability.of.states[2:(min.track.length-1), 1]


p_per<-probability.of.states[2:(min.track.length-1), 2]




### if we consider the same probability for neighbors and states through the time

### What I obtain here is for all cells. The geometric probability should be modified to show the same concept.



repitition<-21

empirical.probabilities.non.persistent<-matrix(0,  nrow= repitition)

waiting.time.matrix.non.persistent <-matrix(0, nrow= repitition)

	

for(i in 1:length(cluster.1)){
	
	
	index.i<-which(data_cell_motility[, 8]==id[cluster.1[i]])
	
	
	for(k in 1: repitition){
		
	  for(j in 1:(length(index.i)-k)){

               if(j>=2){ ####  because we already removed the first position.

               if((sum(classification.matrix[index.i[j:(j+k-1)], 2]==2)==k) & classification.matrix[index.i[j-1], 2]!=2 & classification.matrix[index.i[(j+k)], 2]==1){
              
                  waiting.time.matrix.non.persistent[k]<-waiting.time.matrix.non.persistent[k]+1
                  
}

}

}
 
} 
 
}





waiting.time.matrix.non.persistent.zero<-0

	

for(i in 1:length(cluster.1)){
	
    index.i<-which(data_cell_motility[, 8]==id[cluster.1[i]])

	
	for(j in 2:length(index.i)){

               ####  because we already removed the first position.

               if((classification.matrix[index.i[j], 2]==1) & classification.matrix[index.i[j-1], 2]!=2){
              
                  waiting.time.matrix.non.persistent.zero<-waiting.time.matrix.non.persistent.zero+1
                  
}


 
} 
 
}




freq.non.persistent<-c(waiting.time.matrix.non.persistent.zero, waiting.time.matrix.non.persistent)


for(k in 1: (repitition+1)){
	
	 empirical.probabilities.non.persistent[k]<-(sum(freq.non.persistent[1:k]))          
	
}


empirical.probabilities.non.persistent<-empirical.probabilities.non.persistent/sum(freq.non.persistent)



x<-seq(1, 21, 1)

u<-matrix(0, nrow=number.of.cells, ncol=length(x))

u.mean<-NULL



empirical.plot<-c(0, 0, empirical.probabilities.non.persistent)


xlab.values.axes<-seq(-1, (repitition), 1)


xlab.values<-c(-1, 0, seq(0, (repitition), 1))



points(xlab.values , empirical.plot, type="l", pch=20, col="darkgreen")

axis(1, at=-1:(repitition), labels= xlab.values.axes)



 
v<-NULL

v[1]<-mean(p_non_per)

for(t in 1: (repitition+2)){
	
	v[t]<-((1-v[1])^(t-1))*v[1]
  
 }
 
 
prob.estimated.non<-matrix(0, nrow=(min.track.length-1), ncol=1)

for(t in 1: (repitition+2)){
	
	
	prob.estimated.non[t]<-sum(v[1:t])
	
	
}
 

xlab.values<-seq(0, (repitition), 1)



points(xlab.values, prob.estimated.non[1:22], type="p", pch=13, col="darkgreen")






###################################################################################################################################################################################################################################################




setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/cluster-2-2-states-HMM/")


load(file=paste(".RData", sep=""))

library(rstan)




P_1<-matrix(c(fit_ss$summary[3:6]), nrow=2, ncol=2, byrow=T)

prior_pi<-c(fit_ss$summary[1], fit_ss$summary[2])



beta_persistent<-fit_ss$summary[8]

alpha_persistent<-fit_ss$summary[7]

beta_non_persistent<-fit_ss$summary[10]

alpha_non_persistent<-fit_ss$summary[9]



markov_probability_cluster<-matrix(0, nrow= number.of.positions.total, ncol=2)



for(i in 1: number.of.cells){

					 index.i<-which(data_cell_motility[, 8]==id[i])
					 
					 markov_probability_cluster[index.i[1], ]<-prior_pi
					 
					 for(l in 2:(length(index.i))){
					 	
					 	
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_1
					 	
					 	
					 	
					 } 
				}
######  Defining the ratio for the purpose of classification

t<-which(angular.change.complete.full.data==0)

angular.change.complete.full.data[t]<-0.00001

classification.matrix<-matrix(c(angular.change.complete.full.data, rep(0, length(angular.change.complete.full.data))), nrow=length(angular.change.complete.full.data), ncol=2)

group.probability<-matrix(0, nrow=length(angular.change.complete.full.data), ncol=2)


group.probability[, 1]<-markov_probability_cluster[, 1]*dbeta((angular.change.complete.full.data/pi), alpha_non_persistent, beta_non_persistent)/(markov_probability_cluster[, 1]*dbeta((angular.change.complete.full.data/pi), alpha_non_persistent, beta_non_persistent) +markov_probability_cluster[, 2]*dbeta((angular.change.complete.full.data/pi), alpha_persistent, beta_persistent))


group.probability[, 2]<-markov_probability_cluster[, 2]*dbeta((angular.change.complete.full.data/pi), alpha_persistent, beta_persistent)/(markov_probability_cluster[, 1]*dbeta((angular.change.complete.full.data/pi), alpha_non_persistent, beta_non_persistent) +markov_probability_cluster[, 2]*dbeta((angular.change.complete.full.data/pi), alpha_persistent, beta_persistent))



sum(group.probability[index.angular.change.cluster.2, 1]==group.probability[index.angular.change.cluster.2, 2])

for(i in 1: number.of.positions.total){
	
	
	if(group.probability[i, 1]>=group.probability[i, 2]){classification.matrix[i, 2]<-1}
	
	
	if(group.probability[i, 1]<group.probability[i, 2]){classification.matrix[i, 2]<-2}
	
	
}

weight<-NULL

weight[1]<-mean(classification.matrix[index.angular.change.cluster.2, 2]==1)

weight[2]<-mean(classification.matrix[index.angular.change.cluster.2, 2]==2)

### 1:black <- non-persistent

### 2:darkred <- persistent


index.persistent<-which(classification.matrix[index.angular.change.cluster.2, 2]==2)

index.non.persistent<-which(classification.matrix[index.angular.change.cluster.2, 2]==1)


### For now, as the program includes the radius devided by max for others.

# # # setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/cluster-2-hyper-informative-second-prior-2-states/")
# # # load(file=paste(".RData", sep=""))



persist<-radius.complete.full.data[index.angular.change.cluster.2[index.persistent]]

non.persist<-radius.complete.full.data[index.angular.change.cluster.2[index.non.persistent]]


######   checking the markovian assumption-non-persistent mode

min.track.length<-min(track.length)




####  note that we do not have angular change for the first position



probability.of.states<-matrix(0, nrow=(min.track.length-1), ncol=2)

probability.of.states[1, ]<-prior_pi




for(t in 2:(min.track.length-1)){
	
	
	probability.of.states[t, ]<-probability.of.states[(t-1), ] %*% P_1
	
}




p_non_per<-probability.of.states[2:(min.track.length-1), 1]


p_per<-probability.of.states[2:(min.track.length-1), 2]



### if we consider the same probability for neighbors and states through the time

### What I obtain here is for all cells. The geometric probability should be modified to show the same concept.



repitition<-21

empirical.probabilities.non.persistent<-matrix(0,  nrow= repitition)

waiting.time.matrix.non.persistent <-matrix(0, nrow= repitition)



	

for(i in 1:length(cluster.2)){
	
	index.i<-which(data_cell_motility[, 8]==id[cluster.2[i]])
	
	for(k in 1: repitition){
		
	  for(j in 1:(length(index.i)-k)){

               if(j>=2){ ####  because we already removed the first position.

               if((sum(classification.matrix[index.i[j:(j+k-1)], 2]==2)==k) & classification.matrix[index.i[j-1], 2]!=2 & classification.matrix[index.i[(j+k)], 2]==1){
              
                  waiting.time.matrix.non.persistent[k]<-waiting.time.matrix.non.persistent[k]+1
                  
}

}

}
 
} 
 
}





waiting.time.matrix.non.persistent.zero<-0

	

for(i in 1:length(cluster.2)){
	
	index.i<-which(data_cell_motility[, 8]==id[cluster.2[i]])

	
	for(j in 2:length(index.i)){

               ####  because we already removed the first position.

               if((classification.matrix[index.i[j], 2]==1) & classification.matrix[index.i[j-1], 2]!=2){
              
                  waiting.time.matrix.non.persistent.zero<-waiting.time.matrix.non.persistent.zero+1
                  
}


 
} 
 
}





freq.non.persistent<-c(waiting.time.matrix.non.persistent.zero, waiting.time.matrix.non.persistent)


for(k in 1: (repitition+1)){
	
	 empirical.probabilities.non.persistent[k]<-(sum(freq.non.persistent[1:k]))          
	
}


empirical.probabilities.non.persistent<-empirical.probabilities.non.persistent/sum(freq.non.persistent)



x<-seq(1, 21, 1)

u<-matrix(0, nrow=number.of.cells, ncol=length(x))

u.mean<-NULL



empirical.plot<-c(0, 0, empirical.probabilities.non.persistent)


xlab.values.axes<-seq(-1, (repitition), 1)


xlab.values<-c(-1, 0, seq(0, (repitition), 1))



points(xlab.values , empirical.plot, type="l", pch=20, col="darkblue")

axis(1, at=-1:(repitition), labels= xlab.values.axes)



 
v<-NULL

v[1]<-mean(p_non_per)

for(t in 1: (repitition+2)){
	
	v[t]<-((1-v[1])^(t-1))*v[1]
  
 }
 
 
prob.estimated.non<-matrix(0, nrow=(min.track.length-1), ncol=1)

for(t in 1: (repitition+2)){
	
	
	prob.estimated.non[t]<-sum(v[1:t])
	
	
}
 

xlab.values<-seq(0, (repitition), 1)



points(xlab.values, prob.estimated.non[1:22], type="p", pch=17, col="darkblue")






###################################################################################################################################################################################################################################################




setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/cluster-3-2-states-HMM/")


load(file=paste(".RData", sep=""))

library(rstan)




P_1<-matrix(c(fit_ss$summary[3:6]), nrow=2, ncol=2, byrow=T)

prior_pi<-c(fit_ss$summary[1], fit_ss$summary[2])



beta_persistent<-fit_ss$summary[8]

alpha_persistent<-fit_ss$summary[7]

beta_non_persistent<-fit_ss$summary[10]

alpha_non_persistent<-fit_ss$summary[9]



markov_probability_cluster<-matrix(0, nrow= number.of.positions.total, ncol=2)


for(i in 1: number.of.cells){

					 index.i<-which(data_cell_motility[, 8]==id[i])
					 
					 markov_probability_cluster[index.i[1], ]<-prior_pi
					 
					 for(l in 2:(length(index.i))){
					 	
					 	
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_1
					 	
					 	
					 	
					 } 
				}
######  Defining the ratio for the purpose of classification

t<-which(angular.change.complete.full.data==0)

angular.change.complete.full.data[t]<-0.00001

classification.matrix<-matrix(c(angular.change.complete.full.data, rep(0, length(angular.change.complete.full.data))), nrow=length(angular.change.complete.full.data), ncol=2)

group.probability<-matrix(0, nrow=length(angular.change.complete.full.data), ncol=2)


group.probability[, 1]<-markov_probability_cluster[, 1]*dbeta((angular.change.complete.full.data/pi), alpha_non_persistent, beta_non_persistent)/(markov_probability_cluster[, 1]*dbeta((angular.change.complete.full.data/pi), alpha_non_persistent, beta_non_persistent) +markov_probability_cluster[, 2]*dbeta((angular.change.complete.full.data/pi), alpha_persistent, beta_persistent))


group.probability[, 2]<-markov_probability_cluster[, 2]*dbeta((angular.change.complete.full.data/pi), alpha_persistent, beta_persistent)/(markov_probability_cluster[, 1]*dbeta((angular.change.complete.full.data/pi), alpha_non_persistent, beta_non_persistent) +markov_probability_cluster[, 2]*dbeta((angular.change.complete.full.data/pi), alpha_persistent, beta_persistent))



sum(group.probability[index.angular.change.cluster.3, 1]==group.probability[index.angular.change.cluster.3, 2])

for(i in 1: number.of.positions.total){
	
	
	if(group.probability[i, 1]>=group.probability[i, 2]){classification.matrix[i, 2]<-1}
	
	
	if(group.probability[i, 1]<group.probability[i, 2]){classification.matrix[i, 2]<-2}
	
	
}

weight<-NULL

weight[1]<-mean(classification.matrix[index.angular.change.cluster.3, 2]==1)

weight[2]<-mean(classification.matrix[index.angular.change.cluster.3, 2]==2)

### 1:black <- non-persistent

### 2:darkred <- persistent


index.persistent<-which(classification.matrix[index.angular.change.cluster.3, 2]==2)

index.non.persistent<-which(classification.matrix[index.angular.change.cluster.3, 2]==1)


### For now, as the program includes the radius devided by max for others.

# # # setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/cluster-2-hyper-informative-second-prior-2-states/")
# # # load(file=paste(".RData", sep=""))



persist<-radius.complete.full.data[index.angular.change.cluster.3[index.persistent]]

non.persist<-radius.complete.full.data[index.angular.change.cluster.3[index.non.persistent]]


######   checking the markovian assumption-non-persistent mode

min.track.length<-min(track.length)




####  note that we do not have angular change for the first position



probability.of.states<-matrix(0, nrow=(min.track.length-1), ncol=2)

probability.of.states[1, ]<-prior_pi




for(t in 2:(min.track.length-1)){
	
	
	probability.of.states[t, ]<-probability.of.states[(t-1), ] %*% P_1
	
}




p_non_per<-probability.of.states[2:(min.track.length-1), 1]


p_per<-probability.of.states[2:(min.track.length-1), 2]




### if we consider the same probability for neighbors and states through the time

### What I obtain here is for all cells. The geometric probability should be modified to show the same concept.


repitition<-21

empirical.probabilities.non.persistent<-matrix(0,  nrow= repitition)

waiting.time.matrix.non.persistent <-matrix(0, nrow= repitition)

	
	for(i in 1: length(cluster.3)){
		
	index.i<-which(data_cell_motility[, 8]==id[cluster.3[i]])
	
	for(k in 1: repitition){
		
	  for(j in 1:(length(index.i)-k)){

               if(j>=2){ ####  because we already removed the first position.

               if((sum(classification.matrix[index.i[j:(j+k-1)], 2]==2)==k) & classification.matrix[index.i[j-1], 2]!=2 & classification.matrix[index.i[(j+k)], 2]==1){
              
                  waiting.time.matrix.non.persistent[k]<-waiting.time.matrix.non.persistent[k]+1
                  
}

}

}
 
} 
 
}





waiting.time.matrix.non.persistent.zero<-0

	
	for(i in 1: length(cluster.3)){
		
	    index.i<-which(data_cell_motility[, 8]==id[cluster.3[i]])
	
	    for(j in 2:length(index.i)){

               ####  because we already removed the first position.

               if((classification.matrix[index.i[j], 2]==1) & classification.matrix[index.i[j-1], 2]!=2){
              
                  waiting.time.matrix.non.persistent.zero<-waiting.time.matrix.non.persistent.zero+1
                  
}


 
} 
 
}




freq.non.persistent<-c(waiting.time.matrix.non.persistent.zero, waiting.time.matrix.non.persistent)


for(k in 1: (repitition+1)){
	
	 empirical.probabilities.non.persistent[k]<-(sum(freq.non.persistent[1:k]))          
	
}


empirical.probabilities.non.persistent<-empirical.probabilities.non.persistent/sum(freq.non.persistent)



x<-seq(1, 21, 1)

u<-matrix(0, nrow=number.of.cells, ncol=length(x))

u.mean<-NULL



empirical.plot<-c(0, 0, empirical.probabilities.non.persistent)


xlab.values.axes<-seq(-1, (repitition), 1)


xlab.values<-c(-1, 0, seq(0, (repitition), 1))



points(xlab.values , empirical.plot, type="l", pch=20, col="darkred")

axis(1, at=-1:(repitition), labels= xlab.values.axes)



 
v<-NULL

v[1]<-mean(p_non_per)

for(t in 1: (repitition+2)){
	
	v[t]<-((1-v[1])^(t-1))*v[1]
  
 }
 
 
prob.estimated.non<-matrix(0, nrow=(min.track.length-1), ncol=1)

for(t in 1: (repitition+2)){
	
	
	prob.estimated.non[t]<-sum(v[1:t])
	
	
}
 

xlab.values<-seq(0, (repitition), 1)



points(xlab.values, prob.estimated.non[1:22], type="p", pch=18, col="darkred")





###################################################################################################################################################################################################################################################




setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/cluster-4-2-states-HMM/")


load(file=paste(".RData", sep=""))

library(rstan)




P_1<-matrix(c(fit_ss$summary[3:6]), nrow=2, ncol=2, byrow=T)

prior_pi<-c(fit_ss$summary[1], fit_ss$summary[2])



beta_persistent<-fit_ss$summary[8]

alpha_persistent<-fit_ss$summary[7]

beta_non_persistent<-fit_ss$summary[10]

alpha_non_persistent<-fit_ss$summary[9]



markov_probability_cluster<-matrix(0, nrow= number.of.positions.total, ncol=2)


for(i in 1: number.of.cells){

					 index.i<-which(data_cell_motility[, 8]==id[i])
					 
					 markov_probability_cluster[index.i[1], ]<-prior_pi
					 
					 for(l in 2:(length(index.i))){
					 	
					 	
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_1
					 	
					 	
					 	
					 } 
				}
######  Defining the ratio for the purpose of classification

t<-which(angular.change.complete.full.data==0)

angular.change.complete.full.data[t]<-0.00001

classification.matrix<-matrix(c(angular.change.complete.full.data, rep(0, length(angular.change.complete.full.data))), nrow=length(angular.change.complete.full.data), ncol=2)

group.probability<-matrix(0, nrow=length(angular.change.complete.full.data), ncol=2)


group.probability[, 1]<-markov_probability_cluster[, 1]*dbeta((angular.change.complete.full.data/pi), alpha_non_persistent, beta_non_persistent)/(markov_probability_cluster[, 1]*dbeta((angular.change.complete.full.data/pi), alpha_non_persistent, beta_non_persistent) +markov_probability_cluster[, 2]*dbeta((angular.change.complete.full.data/pi), alpha_persistent, beta_persistent))


group.probability[, 2]<-markov_probability_cluster[, 2]*dbeta((angular.change.complete.full.data/pi), alpha_persistent, beta_persistent)/(markov_probability_cluster[, 1]*dbeta((angular.change.complete.full.data/pi), alpha_non_persistent, beta_non_persistent) +markov_probability_cluster[, 2]*dbeta((angular.change.complete.full.data/pi), alpha_persistent, beta_persistent))



sum(group.probability[index.angular.change.cluster.4, 1]==group.probability[index.angular.change.cluster.4, 2])

for(i in 1: number.of.positions.total){
	
	
	if(group.probability[i, 1]>=group.probability[i, 2]){classification.matrix[i, 2]<-1}
	
	
	if(group.probability[i, 1]<group.probability[i, 2]){classification.matrix[i, 2]<-2}
	
	
}

weight<-NULL

weight[1]<-mean(classification.matrix[index.angular.change.cluster.4, 2]==1)

weight[2]<-mean(classification.matrix[index.angular.change.cluster.4, 2]==2)

### 1:black <- non-persistent

### 2:darkred <- persistent


index.persistent<-which(classification.matrix[index.angular.change.cluster.4, 2]==2)

index.non.persistent<-which(classification.matrix[index.angular.change.cluster.4, 2]==1)


### For now, as the program includes the radius devided by max for others.

# # # setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/cluster-2-hyper-informative-second-prior-2-states/")
# # # load(file=paste(".RData", sep=""))



persist<-radius.complete.full.data[index.angular.change.cluster.4[index.persistent]]

non.persist<-radius.complete.full.data[index.angular.change.cluster.4[index.non.persistent]]


######   checking the markovian assumption-non-persistent mode

min.track.length<-min(track.length)




####  note that we do not have angular change for the first position



probability.of.states<-matrix(0, nrow=(min.track.length-1), ncol=2)

probability.of.states[1, ]<-prior_pi




for(t in 2:(min.track.length-1)){
	
	
	probability.of.states[t, ]<-probability.of.states[(t-1), ] %*% P_1
	
}




p_non_per<-probability.of.states[2:(min.track.length-1), 1]


p_per<-probability.of.states[2:(min.track.length-1), 2]




### if we consider the same probability for neighbors and states through the time

### What I obtain here is for all cells. The geometric probability should be modified to show the same concept.


repitition<-21

empirical.probabilities.non.persistent<-matrix(0,  nrow= repitition)

waiting.time.matrix.non.persistent <-matrix(0, nrow= repitition)

	
	for(i in 1: length(cluster.4)){
		
	index.i<-which(data_cell_motility[, 8]==id[cluster.4[i]])
	
	for(k in 1: repitition){
		
	  for(j in 1:(length(index.i)-k)){

               if(j>=2){ ####  because we already removed the first position.

               if((sum(classification.matrix[index.i[j:(j+k-1)], 2]==2)==k) & classification.matrix[index.i[j-1], 2]!=2 & classification.matrix[index.i[(j+k)], 2]==1){
              
                  waiting.time.matrix.non.persistent[k]<-waiting.time.matrix.non.persistent[k]+1
                  
}

}

}
 
} 
 
}





waiting.time.matrix.non.persistent.zero<-0

	
	for(i in 1: length(cluster.4)){
		
	    index.i<-which(data_cell_motility[, 8]==id[cluster.4[i]])
	
	for(j in 2:length(index.i)){

               ####  because we already removed the first position.

               if((classification.matrix[index.i[j], 2]==1) & classification.matrix[index.i[j-1], 2]!=2){
              
                  waiting.time.matrix.non.persistent.zero<-waiting.time.matrix.non.persistent.zero+1
                  
}


 
} 
 
}




freq.non.persistent<-c(waiting.time.matrix.non.persistent.zero, waiting.time.matrix.non.persistent)


for(k in 1: (repitition+1)){
	
	 empirical.probabilities.non.persistent[k]<-(sum(freq.non.persistent[1:k]))          
	
}


empirical.probabilities.non.persistent<-empirical.probabilities.non.persistent/sum(freq.non.persistent)



x<-seq(1, 21, 1)

u<-matrix(0, nrow=number.of.cells, ncol=length(x))

u.mean<-NULL



empirical.plot<-c(0, 0, empirical.probabilities.non.persistent)


xlab.values.axes<-seq(-1, (repitition), 1)


xlab.values<-c(-1, 0, seq(0, (repitition), 1))



points(xlab.values , empirical.plot, type="l", pch=20, col="darkorange")

axis(1, at=-1:(repitition), labels= xlab.values.axes)



 
v<-NULL

v[1]<-mean(p_non_per)

for(t in 1: (repitition+2)){
	
	v[t]<-((1-v[1])^(t-1))*v[1]
  
 }
 
 
prob.estimated.non<-matrix(0, nrow=(min.track.length-1), ncol=1)

for(t in 1: (repitition+2)){
	
	
	prob.estimated.non[t]<-sum(v[1:t])
	
	
}
 

xlab.values<-seq(0, (repitition), 1)



points(xlab.values, prob.estimated.non[1:22], type="p", pch=23, col="darkorange")


###################################################################################################################################################################################################################################################



setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/Position-2-states-HMM/")


load(file=paste(".RData", sep=""))

library(rstan)


fit_simulated <- extract(fit, permuted = TRUE)


R<-(length(fit_simulated[[1]])/2)





P_1<-matrix(c(fit_ss$summary[3:6]), nrow=2, ncol=2, byrow=T)

prior_pi<-c(fit_ss$summary[1], fit_ss$summary[2])



beta_persistent<-fit_ss$summary[8]

alpha_persistent<-fit_ss$summary[7]

beta_non_persistent<-fit_ss$summary[10]

alpha_non_persistent<-fit_ss$summary[9]



markov_probability_cluster<-matrix(0, nrow= number.of.positions.total, ncol=2)


for(i in 1: number.of.cells){

					 index.i<-which(data_cell_motility[, 8]==id[i])
					 
					 markov_probability_cluster[index.i[1], ]<-prior_pi
					 
					 for(l in 2:(length(index.i))){
					 	
					 	
					 		
					 		markov_probability_cluster[index.i[l], ]=markov_probability_cluster[index.i[l-1], ] %*% P_1
					 	
					 	
					 	
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

for(i in 1: number.of.positions.total){
	
	
	if(group.probability[i, 1]>=group.probability[i, 2]){classification.matrix[i, 2]<-1}
	
	
	if(group.probability[i, 1]<group.probability[i, 2]){classification.matrix[i, 2]<-2}
	
	
}

weight<-NULL

weight[1]<-mean(classification.matrix[index.angular.change, 2]==1)

weight[2]<-mean(classification.matrix[index.angular.change, 2]==2)

### 1:black <- non-persistent

### 2:darkred <- persistent


index.persistent<-which(classification.matrix[index.angular.change, 2]==2)

index.non.persistent<-which(classification.matrix[index.angular.change, 2]==1)


### For now, as the program includes the radius devided by max for others.

# # # setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/cluster-2-hyper-informative-second-prior-2-states/")
# # # load(file=paste(".RData", sep=""))



persist<-radius.complete.full.data[index.angular.change[index.persistent]]

non.persist<-radius.complete.full.data[index.angular.change[index.non.persistent]]

######   checking the markovian assumption-non-persistent mode

min.track.length<-min(track.length)




####  note that we do not have angular change for the first position



probability.of.states<-matrix(0, nrow=(min.track.length-1), ncol=2)

probability.of.states[1, ]<-prior_pi




for(t in 2:(min.track.length-1)){
	
	
	probability.of.states[t, ]<-probability.of.states[(t-1), ] %*% P_1
	
}




p_non_per<-probability.of.states[2:(min.track.length-1), 1]


p_per<-probability.of.states[2:(min.track.length-1), 2]




### if we consider the same probability for neighbors and states through the time

### What I obtain here is for all cells. The geometric probability should be modified to show the same concept.


repitition<-21

empirical.probabilities.non.persistent<-matrix(0,  nrow= repitition)

waiting.time.matrix.non.persistent <-matrix(0, nrow= repitition)

for(i in 1:number.of.cells){
	
	index.i<-which(data_cell_motility[, 8]==id[i])
	
	for(k in 1: repitition){
		
	  for(j in 1:(length(index.i)-k)){

               if(j>=2){ ####  because we already removed the first position.

               if((sum(classification.matrix[index.i[j:(j+k-1)], 2]==2)==k) & classification.matrix[index.i[j-1], 2]!=2 & classification.matrix[index.i[(j+k)], 2]==1){
              
                  waiting.time.matrix.non.persistent[k]<-waiting.time.matrix.non.persistent[k]+1
                  
}

}

}
 
} 
 
}





waiting.time.matrix.non.persistent.zero<-0

for(i in 1:number.of.cells){
	
	index.i<-which(data_cell_motility[, 8]==id[i])
	
	for(j in 2:length(index.i)){

               ####  because we already removed the first position.

               if((classification.matrix[index.i[j], 2]==1) & classification.matrix[index.i[j-1], 2]!=2){
              
                  waiting.time.matrix.non.persistent.zero<-waiting.time.matrix.non.persistent.zero+1
                  
}


 
} 
 
}




freq.non.persistent<-c(waiting.time.matrix.non.persistent.zero, waiting.time.matrix.non.persistent)


for(k in 1: (repitition+1)){
	
	 empirical.probabilities.non.persistent[k]<-(sum(freq.non.persistent[1:k]))          
	
}


empirical.probabilities.non.persistent<-empirical.probabilities.non.persistent/sum(freq.non.persistent)


empirical.plot<-c(0, 0, empirical.probabilities.non.persistent)


x<-seq(1, 21, 1)

u<-matrix(0, nrow=number.of.cells, ncol=length(x))

u.mean<-NULL


########################################################################
####################################################  frequencies



xlab.values.axes<-seq(-1, (repitition), 1)


xlab.values<-c(-1, 0, seq(0, (repitition), 1))


xlab.values.axes<-seq(-1, (repitition), 1)


xlab.values<-c(-1, 0, seq(0, (repitition), 1))



plot(xlab.values[3:(repitition+2)] , log(freq.non.persistent[1: repitition]), xlab="Trials",  main="(b)",  type="p", xlim=c(0, repitition), lwd=0.7, xaxt = "n", pch=20, ylab="log(Frequency)", cex.lab=1.4, cex.main=2, cex.axis=1.5)


axis(1, at=0:(repitition), labels= xlab.values.axes[2:(repitition+2)])



regression.line<-lm(log(freq.non.persistent[1: repitition])~xlab.values[3:(repitition+2)])

abline(b=regression.line[[1]][2], a=regression.line[[1]][1])






