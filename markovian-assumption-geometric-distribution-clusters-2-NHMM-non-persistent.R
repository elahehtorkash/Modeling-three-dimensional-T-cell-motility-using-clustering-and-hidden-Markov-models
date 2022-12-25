remove(list=ls())


############################################## First Position-NHMM


setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/Position-2-states-NHMM/")


load(file=paste(".RData", sep=""))

library(rstan)


fit_simulated <- extract(fit, permuted = TRUE)


R<-(length(fit_simulated[[1]])/2)



prior_pi<-c(fit_ss$summary[1], fit_ss$summary[2])

P_1<-matrix(c(fit_ss$summary[3:6]), nrow=2, ncol=2, byrow=T)

P_2<-matrix(c(fit_ss$summary[7:10]), nrow=2, ncol=2, byrow=T)


beta_non_persistent<-fit_ss$summary[14]

alpha_non_persistent<-fit_ss$summary[13]

beta_persistent<-fit_ss$summary[12]

alpha_persistent<-fit_ss$summary[11]


# # # # # # beta_persistent<-fit_ss$summary[16]

# # # alpha_persistent<-fit_ss$summary[15]

# # # beta_non_persistent<-fit_ss$summary[18]

# # # alpha_non_persistent<-fit_ss$summary[17]



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



ratio.states.neighbor<-array(0, dim=c(2, 2, t.record=(min.track.length-1)))



counter<-array(0, dim=c(2, 2, t.record=(min.track.length-1)))


ratio.states<-matrix(0, nrow=(min.track.length-1), ncol=2)


for(t in 2:(min.track.length-1)){
	
	####  the minimum length of trajectories is (min.track.length-1). Note that we remove the last position which means index.i[(min.track.length-1)] 
	
	
	for(i in 1: number.of.cells){
		
	    index.i<-which(data_cell_motility[, 8]==id[i])
	    
	    if(classification.matrix[index.i[t], 2]==1 & number.of.neighbors.full.data[index.i[t]]==0) counter[1, 1, t]<-counter[1, 1, t]+1 
	    
	    if(classification.matrix[index.i[t], 2]==1 & number.of.neighbors.full.data[index.i[t]]>0) counter[1, 2, t]<-counter[1, 2, t]+1	
		
			    
	    if(classification.matrix[index.i[t], 2]==2 & number.of.neighbors.full.data[index.i[t]]==0) counter[2, 1, t]<-counter[2, 1, t]+1 
	    
	    if(classification.matrix[index.i[t], 2]==2 & number.of.neighbors.full.data[index.i[t]]>0) counter[2, 2, t]<-counter[2, 2, t]+1	
		

		
		
	}
	
		
}


for (t in 2:(min.track.length-1)){
	
	
	total.time.t<-sum(counter[, , t])
	
	
	for(r in 1:2){ ### states
				
		sum.row<-sum(counter[r, , t])
		
		
		ratio.states[t, r]<-sum.row/total.time.t
		
		for(c in 1:2){  ### whether there exists neighbor or not
			
			ratio.states.neighbor[r, c, t]<-counter[r, c, t]/total.time.t
			
			
		}
	}
	
}


P_independent_of_neighbors<-array(0, dim=c(2, 2, t.record=(min.track.length-1)))


for(t in 1:(min.track.length-1)){
	
	for(r in 1:2){
		
		for(c in 1:2){
			
			P_independent_of_neighbors[r, c, t]<-(P_1[r, c]* ratio.states.neighbor[r, 1, t]+P_2[r, c]* ratio.states.neighbor[r, 2, t])/ratio.states[t, r]
			
		}
		
	}
	
	
	
}



####  note that we do not have angular change for the first position



probability.of.states<-matrix(0, nrow=(min.track.length-1), ncol=2)

probability.of.states[1, ]<-prior_pi


probability.of.states[2, ]<-ratio.states[2, ]


for(t in 3:(min.track.length-1)){
	
	
	probability.of.states[t, ]<-probability.of.states[(t-1), ] %*% P_independent_of_neighbors[, , (t-1)]
	
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



     
# # # #       quartz()

      # # pWidth = 10

      # # pHeight = 8

      # # plot.window(c(0,pWidth),
             # # c(0,pHeight))     
             
             
                   
setwd("/Users/e2torkas/Desktop/Clustering-independent-of-direction/Revision-Aug-7/")

  
     
     # i<-2

     # m <- rbind(cbind(c(1, 1, 1), c(2, 3, 4)), c(5))

     # m <- cbind(rbind(c(1, 2), c(3, 4), c(5, 6), c(7, 8)))
     
     m <- rbind(c(1, 2))
     
     # # # # par(mar=c(4,7,2,1)) 
     
     # # # # m <- cbind(c(1, 1, 1), c(2, 3, 4))
     
     # # # # m <- rbind(c(5), cbind(c(1, 1, 1), c(2, 3, 4)))

     layout(m)
     

     # par(lwd=2, mar=c(5, 5, 4, 2))




par(lwd=2, cex.axis=2, mar=c(5, 5, 4, 2) + 0.1)


# # # pdf("Rplot-non-persistent-2-states-NHMM-Frequency.pdf", width=10, height=5, paper='special') 
       
 


xlab.values.axes<-seq(-1, (repitition), 1)


xlab.values<-c(-1, 0, seq(0, (repitition), 1))



plot(xlab.values , empirical.plot, xlab="Trials",  main="(a)",  ylab="CDF", ylim=c(0, 1), type="l", xlim=c(-1, 21), lwd=0.7, xaxt = "n", cex.lab=1.5, cex.main=1.5, cex.axis=1.5)


axis(1, at=0:(repitition), labels= xlab.values.axes[2:(repitition+2)], cex.lab=1.5)



  legend("bottomright", 
  legend = c("No Cluster", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"), 
  col = c("black", "darkgreen", "darkblue", "darkred", "darkorange"), 
  pch = c(20, 13, 17, 18, 23), 
  pt.cex = 2, 
  cex = 1.2, 
  text.col = "black", 
  horiz = F , 
  inset = c(0.1, 0.1), bty="n")

  # # # lty=c(1,1, 1, 1, 1),

  
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



points(xlab.values, prob.estimated.non[1:(repitition+1)], type="p", pch=20)





setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/201-hyper-inf-2-states/")


load(file=paste(".RData", sep=""))

library(rstan)


# # # # fit_simulated <- extract(fit, permuted = TRUE)


# # R<-(length(fit_simulated[[1]])/2)



P_1<-matrix(c(fit_ss$summary[3:6]), nrow=2, ncol=2, byrow=T)

P_2<-matrix(c(fit_ss$summary[7:10]), nrow=2, ncol=2, byrow=T)

prior_pi<-c(fit_ss$summary[1], fit_ss$summary[2])

beta_non_persistent<-fit_ss$summary[14]

alpha_non_persistent<-fit_ss$summary[13]

beta_persistent<-fit_ss$summary[12]

alpha_persistent<-fit_ss$summary[11]


# # # # # # beta_persistent<-fit_ss$summary[16]

# # # alpha_persistent<-fit_ss$summary[15]

# # # beta_non_persistent<-fit_ss$summary[18]

# # # alpha_non_persistent<-fit_ss$summary[17]



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



ratio.states.neighbor<-array(0, dim=c(2, 2, t.record=(min.track.length-1)))



counter<-array(0, dim=c(2, 2, t.record=(min.track.length-1)))


ratio.states<-matrix(0, nrow=(min.track.length-1), ncol=2)


for(t in 2:(min.track.length-1)){
	
	####  the minimum length of trajectories is (min.track.length-1). Note that we remove the last position which means index.i[(min.track.length-1)] 
	
	

for(i in 1:length(cluster.1)){
	
	    index.i<-which(data_cell_motility[, 8]==id[cluster.1[i]])

	    
	    if(classification.matrix[index.i[t], 2]==1 & number.of.neighbors.full.data[index.i[t]]==0) counter[1, 1, t]<-counter[1, 1, t]+1 
	    
	    if(classification.matrix[index.i[t], 2]==1 & number.of.neighbors.full.data[index.i[t]]>0) counter[1, 2, t]<-counter[1, 2, t]+1	
		
			    
	    if(classification.matrix[index.i[t], 2]==2 & number.of.neighbors.full.data[index.i[t]]==0) counter[2, 1, t]<-counter[2, 1, t]+1 
	    
	    if(classification.matrix[index.i[t], 2]==2 & number.of.neighbors.full.data[index.i[t]]>0) counter[2, 2, t]<-counter[2, 2, t]+1	
		

		
		
	}
	
		
}


for (t in 2:(min.track.length-1)){
	
	
	total.time.t<-sum(counter[, , t])
	
	
	for(r in 1:2){ ### states
				
		sum.row<-sum(counter[r, , t])
		
		
		ratio.states[t, r]<-sum.row/total.time.t
		
		for(c in 1:2){  ### whether there exists neighbor or not
			
			ratio.states.neighbor[r, c, t]<-counter[r, c, t]/total.time.t
			
			
		}
	}
	
}


P_independent_of_neighbors<-array(0, dim=c(2, 2, t.record=(min.track.length-1)))


for(t in 1:(min.track.length-1)){
	
	for(r in 1:2){
		
		for(c in 1:2){
			
			P_independent_of_neighbors[r, c, t]<-(P_1[r, c]* ratio.states.neighbor[r, 1, t]+P_2[r, c]* ratio.states.neighbor[r, 2, t])/ratio.states[t, r]
			
		}
		
	}
	
	
	
}



####  note that we do not have angular change for the first position



probability.of.states<-matrix(0, nrow=(min.track.length-1), ncol=2)

probability.of.states[1, ]<-prior_pi


probability.of.states[2, ]<-ratio.states[2, ]


for(t in 3:(min.track.length-1)){
	
	
	probability.of.states[t, ]<-probability.of.states[(t-1), ] %*% P_independent_of_neighbors[, , (t-1)]
	
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



points(xlab.values , empirical.plot, type="l", col="darkgreen")

# axis(1, at=-1:(repitition), labels= xlab.values.axes)



 
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



points(xlab.values, prob.estimated.non[1:(repitition+1)], type="p", pch=13, col="darkgreen")






###################################################################################################################################################################################################################################################




setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/202-hyper-inf-2-states/")


load(file=paste(".RData", sep=""))

library(rstan)


# # # # fit_simulated <- extract(fit, permuted = TRUE)


# # R<-(length(fit_simulated[[1]])/2)



P_1<-matrix(c(fit_ss$summary[3:6]), nrow=2, ncol=2, byrow=T)

P_2<-matrix(c(fit_ss$summary[7:10]), nrow=2, ncol=2, byrow=T)

prior_pi<-c(fit_ss$summary[1], fit_ss$summary[2])

beta_non_persistent<-fit_ss$summary[14]

alpha_non_persistent<-fit_ss$summary[13]

beta_persistent<-fit_ss$summary[12]

alpha_persistent<-fit_ss$summary[11]


# # # # # # beta_persistent<-fit_ss$summary[16]

# # # alpha_persistent<-fit_ss$summary[15]

# # # beta_non_persistent<-fit_ss$summary[18]

# # # alpha_non_persistent<-fit_ss$summary[17]



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



ratio.states.neighbor<-array(0, dim=c(2, 2, t.record=(min.track.length-1)))



counter<-array(0, dim=c(2, 2, t.record=(min.track.length-1)))


ratio.states<-matrix(0, nrow=(min.track.length-1), ncol=2)


for(t in 2:(min.track.length-1)){
	
	####  the minimum length of trajectories is (min.track.length-1). Note that we remove the last position which means index.i[(min.track.length-1)] 
	
	
	

for(i in 1:length(cluster.2)){
	
	    index.i<-which(data_cell_motility[, 8]==id[cluster.2[i]])
	    
	    if(classification.matrix[index.i[t], 2]==1 & number.of.neighbors.full.data[index.i[t]]==0) counter[1, 1, t]<-counter[1, 1, t]+1 
	    
	    if(classification.matrix[index.i[t], 2]==1 & number.of.neighbors.full.data[index.i[t]]>0) counter[1, 2, t]<-counter[1, 2, t]+1	
		
			    
	    if(classification.matrix[index.i[t], 2]==2 & number.of.neighbors.full.data[index.i[t]]==0) counter[2, 1, t]<-counter[2, 1, t]+1 
	    
	    if(classification.matrix[index.i[t], 2]==2 & number.of.neighbors.full.data[index.i[t]]>0) counter[2, 2, t]<-counter[2, 2, t]+1	
		

		
		
	}
	
		
}


for (t in 2:(min.track.length-1)){
	
	
	total.time.t<-sum(counter[, , t])
	
	
	for(r in 1:2){ ### states
				
		sum.row<-sum(counter[r, , t])
		
		
		ratio.states[t, r]<-sum.row/total.time.t
		
		for(c in 1:2){  ### whether there exists neighbor or not
			
			ratio.states.neighbor[r, c, t]<-counter[r, c, t]/total.time.t
			
			
		}
	}
	
}


P_independent_of_neighbors<-array(0, dim=c(2, 2, t.record=(min.track.length-1)))


for(t in 1:(min.track.length-1)){
	
	for(r in 1:2){
		
		for(c in 1:2){
			
			P_independent_of_neighbors[r, c, t]<-(P_1[r, c]* ratio.states.neighbor[r, 1, t]+P_2[r, c]* ratio.states.neighbor[r, 2, t])/ratio.states[t, r]
			
		}
		
	}
	
	
	
}



####  note that we do not have angular change for the first position



probability.of.states<-matrix(0, nrow=(min.track.length-1), ncol=2)

probability.of.states[1, ]<-prior_pi


probability.of.states[2, ]<-ratio.states[2, ]


for(t in 3:(min.track.length-1)){
	
	
	probability.of.states[t, ]<-probability.of.states[(t-1), ] %*% P_independent_of_neighbors[, , (t-1)]
	
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



points(xlab.values , empirical.plot, type="l", col="darkblue")

# axis(1, at=-1:(repitition), labels= xlab.values.axes)



 
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



points(xlab.values, prob.estimated.non[1:(repitition+1)], type="p", pch=17, col="darkblue")






###################################################################################################################################################################################################################################################




setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/203-hyper-inf-2-states/")


load(file=paste(".RData", sep=""))

library(rstan)


# # # # fit_simulated <- extract(fit, permuted = TRUE)


# # R<-(length(fit_simulated[[1]])/2)



P_1<-matrix(c(fit_ss$summary[3:6]), nrow=2, ncol=2, byrow=T)

P_2<-matrix(c(fit_ss$summary[7:10]), nrow=2, ncol=2, byrow=T)

prior_pi<-c(fit_ss$summary[1], fit_ss$summary[2])

beta_non_persistent<-fit_ss$summary[14]

alpha_non_persistent<-fit_ss$summary[13]

beta_persistent<-fit_ss$summary[12]

alpha_persistent<-fit_ss$summary[11]


# # # # # # beta_persistent<-fit_ss$summary[16]

# # # alpha_persistent<-fit_ss$summary[15]

# # # beta_non_persistent<-fit_ss$summary[18]

# # # alpha_non_persistent<-fit_ss$summary[17]



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



ratio.states.neighbor<-array(0, dim=c(2, 2, t.record=(min.track.length-1)))



counter<-array(0, dim=c(2, 2, t.record=(min.track.length-1)))


ratio.states<-matrix(0, nrow=(min.track.length-1), ncol=2)  ###  it is a decoding process. 


for(t in 2:(min.track.length-1)){
	
	####  the minimum length of trajectories is (min.track.length-1). Note that we remove the last position which means index.i[(min.track.length-1)] 
	
	
	for(i in 1: length(cluster.3)){
		
	    index.i<-which(data_cell_motility[, 8]==id[cluster.3[i]])
	    
	    if(classification.matrix[index.i[t], 2]==1 & number.of.neighbors.full.data[index.i[t]]==0) counter[1, 1, t]<-counter[1, 1, t]+1 
	    
	    if(classification.matrix[index.i[t], 2]==1 & number.of.neighbors.full.data[index.i[t]]>0) counter[1, 2, t]<-counter[1, 2, t]+1	
		
			    
	    if(classification.matrix[index.i[t], 2]==2 & number.of.neighbors.full.data[index.i[t]]==0) counter[2, 1, t]<-counter[2, 1, t]+1 
	    
	    if(classification.matrix[index.i[t], 2]==2 & number.of.neighbors.full.data[index.i[t]]>0) counter[2, 2, t]<-counter[2, 2, t]+1	
		

		
		
	}
	
		
}


for (t in 2:(min.track.length-1)){
	
	
	total.time.t<-sum(counter[, , t])
	
	
	for(r in 1:2){ ### states
				
		sum.row<-sum(counter[r, , t])
		
		
		ratio.states[t, r]<-sum.row/total.time.t
		
		for(c in 1:2){  ### whether there exists neighbor or not
			
			ratio.states.neighbor[r, c, t]<-counter[r, c, t]/total.time.t
			
			
		}
	}
	
}




P_independent_of_neighbors<-array(0, dim=c(2, 2, t.record=(min.track.length-1)))


for(t in 1:(min.track.length-1)){
	
	for(r in 1:2){
		
		for(c in 1:2){
			
			P_independent_of_neighbors[r, c, t]<-(P_1[r, c]* ratio.states.neighbor[r, 1, t]+P_2[r, c]* ratio.states.neighbor[r, 2, t])/ratio.states[t, r]
			
		}
		
	}
	
	
	
}



####  note that we do not have angular change for the first position



probability.of.states <- matrix(0, nrow=(min.track.length-1), ncol=2)

probability.of.states[1, ] <- prior_pi


probability.of.states[2, ] <- ratio.states[2, ]


for(t in 3:(min.track.length-1)){
	
	
	probability.of.states[t, ] <- probability.of.states[(t-1), ] %*% P_independent_of_neighbors[, , (t-1)]
	
}




p_non_per <- probability.of.states[2:(min.track.length-1), 1]


p_per <- probability.of.states[2:(min.track.length-1), 2]




### if we consider the same probability for neighbors and states through the time

### What I obtain here is for all cells. The geometric probability should be modified to show the same concept.



repitition<-21

empirical.probabilities.non.persistent <- matrix(0,  nrow= repitition)

waiting.time.matrix.non.persistent <- matrix(0, nrow= repitition)

	
	for(i in 1: length(cluster.3)){
		
	index.i<-which(data_cell_motility[, 8]==id[cluster.3[i]])
	
	for(k in 1: repitition){
		
	  for(j in 1:(length(index.i)-k)){

               if(j>=2){ ####  because we already removed the first position.

               if(classification.matrix[index.i[j-1], 2]!=2 &  (sum(classification.matrix[index.i[j:(j+k-1)], 2]==2)==k) &  classification.matrix[index.i[(j+k)], 2]==1){
              
                  waiting.time.matrix.non.persistent[k] <- waiting.time.matrix.non.persistent[k] + 1
                  
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



points(xlab.values , empirical.plot, type="l", col="darkred")

# axis(1, at=-1:(repitition), labels= xlab.values.axes)



 
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



points(xlab.values, prob.estimated.non[1:(repitition+1)], type="p", pch=18, col="darkred")





###################################################################################################################################################################################################################################################




setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/204-hyper-inf-2-states/")


load(file=paste(".RData", sep=""))

library(rstan)


# # # # fit_simulated <- extract(fit, permuted = TRUE)


# # R<-(length(fit_simulated[[1]])/2)



P_1<-matrix(c(fit_ss$summary[3:6]), nrow=2, ncol=2, byrow=T)

P_2<-matrix(c(fit_ss$summary[7:10]), nrow=2, ncol=2, byrow=T)

prior_pi<-c(fit_ss$summary[1], fit_ss$summary[2])

beta_non_persistent<-fit_ss$summary[14]

alpha_non_persistent<-fit_ss$summary[13]

beta_persistent<-fit_ss$summary[12]

alpha_persistent<-fit_ss$summary[11]


# # # # # # beta_persistent<-fit_ss$summary[16]

# # # alpha_persistent<-fit_ss$summary[15]

# # # beta_non_persistent<-fit_ss$summary[18]

# # # alpha_non_persistent<-fit_ss$summary[17]



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



ratio.states.neighbor<-array(0, dim=c(2, 2, t.record=(min.track.length-1)))



counter<-array(0, dim=c(2, 2, t.record=(min.track.length-1)))


ratio.states<-matrix(0, nrow=(min.track.length-1), ncol=2)


for(t in 2:(min.track.length-1)){
	
	####  the minimum length of trajectories is (min.track.length-1). Note that we remove the last position which means index.i[(min.track.length-1)] 
	
	
	
	for(i in 1: length(cluster.4)){
		
	    index.i<-which(data_cell_motility[, 8]==id[cluster.4[i]])
	    
	    if(classification.matrix[index.i[t], 2]==1 & number.of.neighbors.full.data[index.i[t]]==0) counter[1, 1, t]<-counter[1, 1, t]+1 
	    
	    if(classification.matrix[index.i[t], 2]==1 & number.of.neighbors.full.data[index.i[t]]>0) counter[1, 2, t]<-counter[1, 2, t]+1	
		
			    
	    if(classification.matrix[index.i[t], 2]==2 & number.of.neighbors.full.data[index.i[t]]==0) counter[2, 1, t]<-counter[2, 1, t]+1 
	    
	    if(classification.matrix[index.i[t], 2]==2 & number.of.neighbors.full.data[index.i[t]]>0) counter[2, 2, t]<-counter[2, 2, t]+1	
		

		
		
	}
	
		
}


for (t in 2:(min.track.length-1)){
	
	
	total.time.t<-sum(counter[, , t])
	
	
	for(r in 1:2){ ### states
				
		sum.row<-sum(counter[r, , t])
		
		
		ratio.states[t, r]<-sum.row/total.time.t
		
		for(c in 1:2){  ### whether there exists neighbor or not
			
			ratio.states.neighbor[r, c, t]<-counter[r, c, t]/total.time.t
			
			
		}
	}
	
}


P_independent_of_neighbors<-array(0, dim=c(2, 2, t.record=(min.track.length-1)))


for(t in 1:(min.track.length-1)){
	
	for(r in 1:2){
		
		for(c in 1:2){
			
			P_independent_of_neighbors[r, c, t]<-(P_1[r, c]* ratio.states.neighbor[r, 1, t]+P_2[r, c]* ratio.states.neighbor[r, 2, t])/ratio.states[t, r]
			
		}
		
	}
	
	
	
}



####  note that we do not have angular change for the first position



probability.of.states<-matrix(0, nrow=(min.track.length-1), ncol=2)

probability.of.states[1, ]<-prior_pi


probability.of.states[2, ]<-ratio.states[2, ]


for(t in 3:(min.track.length-1)){
	
	
	probability.of.states[t, ]<-probability.of.states[(t-1), ] %*% P_independent_of_neighbors[, , (t-1)]
	
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



points(xlab.values , empirical.plot, type="l", col="darkorange")

# axis(1, at=-1:(repitition), labels= xlab.values.axes)



 
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



points(xlab.values, prob.estimated.non[1:(repitition+1)], type="p", pch=23, col="darkorange")


###################################################################################################################################################################################################################################################



setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/Position-2-states-NHMM/")


load(file=paste(".RData", sep=""))

library(rstan)


fit_simulated <- extract(fit, permuted = TRUE)


R<-(length(fit_simulated[[1]])/2)



P_1<-matrix(c(fit_ss$summary[3:6]), nrow=2, ncol=2, byrow=T)

P_2<-matrix(c(fit_ss$summary[7:10]), nrow=2, ncol=2, byrow=T)

prior_pi<-c(fit_ss$summary[1], fit_ss$summary[2])

beta_non_persistent<-fit_ss$summary[14]

alpha_non_persistent<-fit_ss$summary[13]

beta_persistent<-fit_ss$summary[12]

alpha_persistent<-fit_ss$summary[11]


# # # # # # beta_persistent<-fit_ss$summary[16]

# # # alpha_persistent<-fit_ss$summary[15]

# # # beta_non_persistent<-fit_ss$summary[18]

# # # alpha_non_persistent<-fit_ss$summary[17]



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



ratio.states.neighbor<-array(0, dim=c(2, 2, t.record=(min.track.length-1)))



counter<-array(0, dim=c(2, 2, t.record=(min.track.length-1)))


ratio.states<-matrix(0, nrow=(min.track.length-1), ncol=2)


for(t in 2:(min.track.length-1)){
	
	####  the minimum length of trajectories is (min.track.length-1). Note that we remove the last position which means index.i[(min.track.length-1)] 
	
	
	for(i in 1: number.of.cells){
		
	    index.i<-which(data_cell_motility[, 8]==id[i])
	    
	    if(classification.matrix[index.i[t], 2]==1 & number.of.neighbors.full.data[index.i[t]]==0) counter[1, 1, t]<-counter[1, 1, t]+1 
	    
	    if(classification.matrix[index.i[t], 2]==1 & number.of.neighbors.full.data[index.i[t]]>0) counter[1, 2, t]<-counter[1, 2, t]+1	
		
			    
	    if(classification.matrix[index.i[t], 2]==2 & number.of.neighbors.full.data[index.i[t]]==0) counter[2, 1, t]<-counter[2, 1, t]+1 
	    
	    if(classification.matrix[index.i[t], 2]==2 & number.of.neighbors.full.data[index.i[t]]>0) counter[2, 2, t]<-counter[2, 2, t]+1	
		

		
		
	}
	
		
}


for (t in 2:(min.track.length-1)){
	
	
	total.time.t<-sum(counter[, , t])
	
	
	for(r in 1:2){ ### states
				
		sum.row<-sum(counter[r, , t])
		
		
		ratio.states[t, r]<-sum.row/total.time.t
		
		for(c in 1:2){  ### whether there exists neighbor or not
			
			ratio.states.neighbor[r, c, t]<-counter[r, c, t]/total.time.t
			
			
		}
	}
	
}


P_independent_of_neighbors<-array(0, dim=c(2, 2, t.record=(min.track.length-1)))


for(t in 1:(min.track.length-1)){
	
	for(r in 1:2){
		
		for(c in 1:2){
			
			P_independent_of_neighbors[r, c, t]<-(P_1[r, c]* ratio.states.neighbor[r, 1, t]+P_2[r, c]* ratio.states.neighbor[r, 2, t])/ratio.states[t, r]
			
		}
		
	}
	
	
	
}



####  note that we do not have angular change for the first position



probability.of.states<-matrix(0, nrow=(min.track.length-1), ncol=2)

probability.of.states[1, ]<-prior_pi


probability.of.states[2, ]<-ratio.states[2, ]


for(t in 3:(min.track.length-1)){
	
	
	probability.of.states[t, ]<-probability.of.states[(t-1), ] %*% P_independent_of_neighbors[, , (t-1)]
	
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

      
setwd("/Users/e2torkas/Desktop/Clustering-independent-of-direction/Revision-Aug-7/")


       


plot(xlab.values[3:(repitition+2)] , log(freq.non.persistent[1: repitition]), xlab="Trials",  main="(b)",  type="p", xlim=c(0, repitition), lwd=0.7, xaxt = "n", pch=20, ylab="log(Frequency)", cex.lab=1.5, cex.main=1.5, cex.axis=1.5)


axis(1, at=0:(repitition), labels= xlab.values.axes[2:(repitition+2)], cex.lab=1.5)

legend("bottomleft", 
  legend = c("No Cluster"), 
  col = c("black"), 
  pch = c(20), 
  pt.cex = 2, 
  cex = 1.2, 
  text.col = "black", 
  horiz = F , 
  inset = c(0.1, 0.1))



regression.line<-lm(log(freq.non.persistent[1: repitition])~xlab.values[3:(repitition+2)])

abline(b=regression.line[[1]][2], a=regression.line[[1]][1])

setwd("/Users/e2torkas/Desktop/Clustering-independent-of-direction/Revision-Aug-7/")


dev.print(pdf, 'Rplot-non-persistent-2-states-NHMM-Frequency.pdf', width=15, height=6, paper='special')

# setwd("/Users/e2torkas/Desktop/Clustering-independent-of-direction/Revision-May-28/")

dev.copy(jpeg, filename="Rplot-non-persistent-2-states-NHMM-Frequency.jpeg");
          
      dev.off()
      






