R=201








###  setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Motility-Data")

###  setwd("/Users/elahehtorkashvand/Desktop/Cell-Motility-Project-Analysis/Motility-Data")


# Simulation on constraint bayes of the small area covariate subject to the functional measurement error:

set.seed(5)

date()


# # # #####################################################################
# # # #####################################################################
# # # ###  clusters


cluster.1<-c(1,   2,   6,   8,   9,  12,  17,  18,  19,  20,  21,  22,  46,  52,  54,  57,  60,  62,  68,  69,  77,  83,  84,  85,  86,  90,  91, 104, 105, 113, 124, 127, 128, 132, 135, 150, 151, 152, 153, 157, 160, 162, 163, 166, 169, 171, 173, 178, 183, 187, 191, 192, 193, 199, 206, 207, 216, 218, 234, 238, 241, 244, 247, 248, 252, 253, 256, 257, 258, 260, 266, 267, 273, 276, 277, 283, 289, 293, 296, 298, 302, 303, 304, 310, 313, 321, 322, 324, 325, 332, 341, 345, 346,  3,   7,  10,  13,  15,  16,  23,  24,  28,  32,  35,  36,  37,  39,  45,  53,  55,  71,  73,  75,  78,  92,  93, 103, 108, 109, 110, 111, 130, 133, 134, 137, 139, 143, 144, 146, 148, 161, 165, 172, 174, 197, 201, 202, 204, 210, 211, 217, 221, 228, 230, 231, 240, 242, 243, 251, 254, 255, 263, 264, 265, 269, 271, 281, 287, 294, 295, 299, 300, 315, 318, 323, 327, 329, 330, 331, 335, 336, 338, 339, 342, 343, 349, 352)



cluster.2<-c( 4,   5,  11,  14,  25,  26,  29,  30,  31,  33,  34,  38,  40,  41,  47 , 48 , 50 , 51 , 58,  59,  61,  63,  64,  65,  66,  67,  70,  72,  74,  79,  80,  81,  82,  87,  88,  94,  95,  96 , 97,  98, 100, 101, 102, 106, 107, 114, 115, 116, 118, 120 ,121, 122, 123, 125, 126, 129, 136, 140, 141, 142, 145, 147, 149, 154, 155, 156, 158, 164, 167, 168, 170, 175, 179, 180 ,181, 182, 184, 185, 186, 188 ,189, 190 ,194, 200, 203, 205, 208, 209, 213, 214, 219, 220, 222, 223, 224, 225 ,227, 233, 236, 237, 239 ,250, 259, 261, 262, 268, 270, 272, 274, 275, 278, 279, 280, 282, 284, 286, 288, 290 ,292, 297, 305, 307 ,308, 311, 312, 316, 320, 326, 328, 333, 334, 337, 340, 344, 348, 350, 351, 27,  42,  43,  44,  49,  56 , 76,  89,  99, 112, 117, 119, 131, 138, 159 ,176 ,177 ,195 ,196 ,198 ,212, 215, 226, 229, 232, 235, 245, 246, 249, 285, 291, 301, 306, 309, 314, 317, 319, 347)


cluster.3<-c(1, 2, 3)

cluster.4<-c(5, 6, 7)
 
 

cluster<-list("1"=cluster.1, "2"=cluster.2, "3"=cluster.3, "4"=cluster.4)


#################################################################################

    
                              
                 number.of.neighbors.cluster<-NULL

                 angular.change.cluster<-NULL
                 

                
				data_cell_motility_1<-read.csv("Position_1.csv", header=TRUE)

				dim(data_cell_motility_1)

				#########   In order to have an understanding of the number of T-cells observed in the study, we use the following comment

				number.of.cells_1<-length(unique(data_cell_motility_1[, 8]))

    
                

                # data_cell_motility_2<-read.csv("Position_2.csv", header=TRUE)

				# dim(data_cell_motility_2)

				# data_cell_motility_2[, 8]<-data_cell_motility_2[, 8]+(number.of.cells_1)

				#####   In order to have an understanding of the number of T-cells observed in the study, we use the following comment

				# number.of.cells_2<-length(unique(data_cell_motility_2[, 8]))



                
                data_cell_motility_3<-read.csv("Position_3.csv", header=TRUE)

				dim(data_cell_motility_3)

				data_cell_motility_3[, 8]<-data_cell_motility_3[, 8]+(number.of.cells_1)

				#####   In order to have an understanding of the number of T-cells observed in the study, we use the following comment

                
                number.of.cells_3<-length(unique(data_cell_motility_3[, 8]))
                



                
                data_cell_motility_6<-read.csv("Position_6.csv", header=TRUE)

				dim(data_cell_motility_6)

				data_cell_motility_6[, 8]<-data_cell_motility_6[, 8]+(number.of.cells_1+number.of.cells_3)

				#####   In order to have an understanding of the number of T-cells observed in the study, we use the following comment

				number.of.cells_6<-length(unique(data_cell_motility_6[, 8]))



    
    
                
                data_cell_motility<-rbind(data_cell_motility_1, data_cell_motility_3,  data_cell_motility_6)

				number.of.cells<-number.of.cells_1+number.of.cells_3+number.of.cells_6

	            #########   In order to have an understanding of the number of T-cells observed in the study, we use the following comment

             
				number.of.cells<-length(unique(data_cell_motility[, 8]))

				id<-seq(1000000000, (1000000000+number.of.cells-1), 1)
				
				track.length<-NULL
				
                for(i in 1: number.of.cells){
                	
                	track.length[i]<-length(which(data_cell_motility[, 8]==id[i]))
                	
                }

                number.of.positions.total<-sum(track.length)

				
				angular.change.complete.full.data<-matrix(0, nrow= number.of.positions.total, ncol=1)
				
				
				number.of.neighbors.full.data<-matrix(0, nrow= number.of.positions.total, ncol=1)

                
                radius.complete.full.data<-matrix(0, nrow= number.of.positions.total, ncol=1)






######  We are going to consider two separate states for each trajectory using angular change. These states are persistent and non-persistent.
##########   Feb 13, 2016


####  Note that as we are interested to see how movements on different directions interact with each other, X, Y, and Z are modeled in the same model.
 #####  for the office Mac I should change it to e2torkas
 
############   Cell Motility Project-Sep 12-2016   ###############
####  WARNING: for the purpose of dependeny, experiments should be kept separated.

 
date()


# # install.packages("mgcv", dependencies=T)


# install.packages("boot", dependencies=T)


# install.packages("MotilityLab", dependencies=T)


# install.packages("moments", dependencies=T)


# install.packages("cluster", dependencies=T)


# install.packages("scatterplot3d", dependencies=T)


# install.packages("rgl", dependencies=T)


# install.packages("coda", dependencies=T)


# install.packages("rjags", dependencies=T)


# install.packages("R2WinBUGS", dependencies=T)


# install.packages("dclone", dependencies=T)


# install.packages("stats", dependencies=T)


# install.packages("nlme", dependencies=T)


# install.packages("MASS", dependencies=T)


# install.packages("splines", dependencies=T)


# install.packages("lattice", dependencies=T)


# install.packages("xtable", dependencies=T)


# install.packages("rstan", dependencies=T)


# install.packages("Rcpp", dependencies=T)


# install.packages("parallel", dependencies=T)


# install.packages("mvtnorm", dependencies=T)



require("mgcv")


library(boot)


require("MotilityLab")


library(moments)


library(cluster)


library(scatterplot3d)


library(rgl)


library(coda)


library(stats)


library(nlme)


library(MASS)


library(splines)


library(lattice)


library(xtable)


library(rstan)


library(Rcpp)


library(parallel)


library(mvtnorm)
		

library(CircStats)	

	
	   tracks<-list() 
		
###  This gives the trajectories in space. I suggest you to run this code. By looking at the trajectories, it make sense to cluster them in different groups.

       track.length<-NULL
              
       psi.trajectory<-NULL       
       
       radius.trajectory<-NULL
       
       theta.trajectory<-NULL
       
       phi.trajectory<-NULL
       
       #### We are looking for a model that gives us the best displacement fit!

       #### displacement<-NULL
       
       
       
       
       increment.x<-increment.y<-increment.z<-NULL

       
       

                      
  #####  stan file
  
                      
                  
  stanmodelcode="data {
	
  int<lower=1> N;        // Number of points from 2 to track length-1 of the i'th for a cluster. 
          
  int<lower=1> K;          // number of states

  int<lower=1> J;          // number of cells in the cluster.
      
  vector[N] angular_change_cluster;    // an array containning the angular change of all trajectories.

  # # # vector[N] radius_cluster;     // an array containning the radius (second-to-the-end) of all trajectories.

  vector[N+J] number_of_neighbor_cluster;  // we will use it as a covariate to determine two transition matrices.
       
  vector<lower=0>[K] alpha_prior;   // dirichlet distribution parameters

  vector[K] alpha_P_1[K];   // dirichlet distribution parameters for P_1

  vector[K] alpha_P_2[K];   // dirichlet distribution parameters for P_2
  
  int<lower=1> counter_n[J];  // where to start the markovian process
  
  int<lower=1> b_final_f[J];  // where to end the markovian process
  
  int<lower=1>  a_initial_f[J];  // first position of each trajectory
  

      
}


parameters {
	
	

    simplex[K] prior_pi;    // k is the number of observed models we consider.
    
    simplex[K] P[2, K];    //  transition matrix 
   
    real<lower=0, upper=4> alpha_per_angle;  
      
    real<lower=0, upper=10> beta_per_angle;
        
    real<lower=0, upper=10> alpha_non_per_angle;
    
    real<lower=0, upper=4> beta_non_per_angle;
   
    
   

}



model {
	  

  matrix[K, K] P_transformed[K];   //saving P as a matrix format.

  vector[N+J] prior_pi_alternative[K];
  
  vector[N] prior_pi_modeling[K];
   
  row_vector[N] models[K];    // coordinates 
  
  row_vector[K] a;

  row_vector[K] b;
  
  int t;
  
  vector[K] transfer;
  
  
  
  
 	
  for(k in 1:K){

           target+= dirichlet_lpdf(P[1, k]| alpha_P_1[k]);
  
  }

  
  
  for(k in 1:K){
  
           target+= dirichlet_lpdf(P[2, k]| alpha_P_2[k]);

   }

  
     
  
     

  for(l in 1:2){
  
   for(m in 1:K){
  	
  	 for(n in 1:K){
  		
  		 P_transformed[l, m, n] = P[l, m, n];
  		
  	 }
  	
    }    	

   }
  
  
     target+= dirichlet_lpdf(prior_pi| alpha_prior);   
  
	
   for(j in 1:J){	

     # target+= normal_lpdf(alpha_r[j]| b_0, sigma_alpha);
  	
     # target+= chi_square_lpdf(radius[a_initial_f[j]:b_final_f[j]]| alpha_r[j]); 

     
     for(k in 1:K){

           prior_pi_alternative[k, a_initial_f[j]] = prior_pi[k];


     }
  	 

     for(l in counter_n[j]: b_final_f[j]){


           for(k in 1:K){

                  b[k] = prior_pi_alternative[k, (l-1)];

           }



     if(number_of_neighbor_cluster[l-1]<=0){
       
           a = b * P_transformed[1];
          
     }


     if(number_of_neighbor_cluster[l-1]>0){
       
           a = b * P_transformed[2];
          
     }



     for(k in 1:K){
           
           prior_pi_alternative[k, l] = a[k];

     }

  	
  }

 } 
 
 t = 1;
 
 for(j in 1:J){
 	
 	for(l in counter_n[j]: b_final_f[j]){
 		
 		for(k in 1:K){
 		
 		      prior_pi_modeling[k, t] = prior_pi_alternative[k, l];
 		
 				
} 	

    t = t+1;

}
 	
 	
 }
 

  


##   Introducing no priors means we consider flat prior for the parameter. Note that the domain of the flat prior comes from the bound we have already defined for parameters.
   


   for(n in 1:N){
   	
   	           	transfer[1] = log(prior_pi_modeling[1, n]) +  beta_lpdf( (angular_change_cluster[n]/pi())| alpha_non_per_angle, beta_non_per_angle); 
   	          	 	          	
   	          	transfer[2] = log(prior_pi_modeling[2, n]) +  beta_lpdf((angular_change_cluster[n]/pi())| alpha_per_angle, beta_per_angle);  	   
   	          

                target += log_sum_exp( transfer );

   }

                       
  
 }"
                                
                                
                                
                                
###################################################################################################################
###############################################################  Spherical representation and angular change 


cluster<-list("1"=cluster.1, "2"=cluster.2, "3"=cluster.3, "4"=cluster.4)


 Mode<-function(x) {
 	
  ux <- unique(x)
  
  ux[which.max(tabulate(match(x, ux)))]
  
}





          
    radius.phi.theta.psi<-function(trajectory, trajectory.num){
					
					length.trajectory<-length(trajectory[, 1])
					
					length.sub.trajectories<-NULL
					
					theta<-NULL
					
					phi<-NULL
					
					for(i in 1: (length.trajectory-1)){
						
					   	length.sub.trajectories[i]<-sqrt((trajectory[(i+1), 2]-trajectory[i, 2])^2+(trajectory[(i+1), 3]-trajectory[i, 3])^2+(trajectory[(i+1), 4]-trajectory[i, 4])^2)
					   	
					   	# theta[i]<-acos((trajectory[(i+1), 4]-trajectory[i, 4])/length.sub.trajectories[i])
					   	
					   	#####   Considering the revised theta under my assumption
					   	
					    theta[i]<-acos(abs((trajectory[(i+1), 4]-trajectory[i, 4])/length.sub.trajectories[i]))
				
					   	
					   	# if(theta[i]>(pi/2)) theta[i]<-pi-theta[i]
					   	
					   	phi[i]<-atan((trajectory[(i+1), 3]-trajectory[i, 3])/(trajectory[(i+1), 2]-trajectory[i, 2]))
					   	
					   	if((trajectory[(i+1), 3]-trajectory[i, 3])<0 & (trajectory[(i+1), 2]-trajectory[i, 2])<0) phi[i]<-pi+phi[i]
					   	
					   	if((trajectory[(i+1), 3]-trajectory[i, 3])>0 & (trajectory[(i+1), 2]-trajectory[i, 2])<0) phi[i]<-pi+phi[i]
					   		   	
				        if(((trajectory[(i+1), 3]-trajectory[i, 3])==0) & ((trajectory[(i+1), 2]-trajectory[i, 2])==0)){
				        	
				        	if(trajectory[i, 3]>0 & trajectory[i, 2]>0) phi[i]<-pi/4
				        	
				        	if(trajectory[i, 3]>0 & trajectory[i, 2]<0) phi[i]<-3*pi/4
				        	
				        	if(trajectory[i, 3]<0 & trajectory[i, 2]<0) phi[i]<--3*pi/4
				        	
				        	if(trajectory[i, 3]<0 & trajectory[i, 2]>0) phi[i]<--pi/4
				        	}
					   			
						}
						
						angle.between.sub.trajectories<-NULL
					
					for(i in 1: (length.trajectory-2)){
						
						u<-c((trajectory[(i+1), 2]-trajectory[i, 2]), (trajectory[(i+1), 3]-trajectory[i, 3]), (trajectory[(i+1), 4]-trajectory[i, 4]))
						
						# print("u"); print(u)
						
						v<-c((trajectory[(i+2), 2]-trajectory[i+1, 2]), (trajectory[(i+2), 3]-trajectory[i+1, 3]), (trajectory[(i+2), 4]-trajectory[i+1, 4]))
						# print("v"); print(v)
						
						u.cross.product.v<-c((u[2]*v[3]-u[3]*v[2]), (u[3]*v[1]-u[1]*v[3]), (u[1]*v[2]-u[2]*v[1]))
						
						u.dot.product.v<-u[1]*v[1]+u[2]*v[2]+u[3]*v[3]
						
						length.u.cross.product.v<-sqrt((u.cross.product.v[1])^2+(u.cross.product.v[2])^2+(u.cross.product.v[3])^2)
						
						length.u<-sqrt((u[1])^2+(u[2])^2+(u[3])^2)
						
						length.v<-sqrt((v[1])^2+(v[2])^2+(v[3])^2)
						
						# angle.between.sub.trajectories[i]<-asin(length.u.cross.product.v/(length.u* length.v))
						
						angle.between.sub.trajectories[i]<-acos(min(c(1, max(c(-1, u.dot.product.v/(length.u* length.v))))))	
												
					}
				
						
					increments.radius.theta.phi.psi<-list(radius= length.sub.trajectories, theta.component= theta, phi.component= phi, psi.component= angle.between.sub.trajectories)
					    
					    
					    if(sum(is.na(length.sub.trajectories))>0) {print("length.trajectory"); print(length.sub.trajectories)}
					    	    
					    
					    if(sum(is.na(theta)>0)){ print("length(theta)"); print(theta)}
					    
					    
					    if(sum(is.na(phi)>0)) {print("length(phi.component)"); print(phi)}
					    
					    
					    if(sum(is.na(angle.between.sub.trajectories)>0)) {
					    	
					    	print("length(psi.component)"); 
					    	
					    	print(angle.between.sub.trajectories);
					    	
					    	print("trajectory.num");
					    	
					    	print(trajectory.num)
					    	
					    	}
				
				
						
						return(increments.radius.theta.phi.psi)
						
	
}
				
				
								
			###  This gives the trajectories in space. I suggest you to run this code. By looking at the trajectories, it make sense to cluster them in different groups.
				
                         
		       		        
		      total.phi<-NULL
			    
			  total.psi<-NULL
			    
			  total.theta<-NULL
			      
			  total.radius<-NULL  
				
              increment.x.cell<-NULL
		       
		      increment.y.cell<-NULL
		        
		      increment.z.cell<-NULL 
		        
		    
		    
		    
           				                              
                
                track.length<-NULL
                           
                              
                id_1<-seq(1000000000, (1000000000+number.of.cells_1-1), 1)

				
				for(i in 1: number.of.cells_1){
				
				    
				    index.i<-which(data_cell_motility[, 8]==id_1[i])

				    
				    track.length[i]<-length(index.i)
				    
				    				    				    
				    data_i<-matrix(c(data_cell_motility[index.i, 7], data_cell_motility[index.i, 1], data_cell_motility[index.i, 2], data_cell_motility[index.i, 3]), nrow=length(index.i), ncol=4)
				    				
				    x<-tracks(data_i)
				        				        
                                  u<-radius.phi.theta.psi(trajectory=x$'1', trajectory.num=i) 

                                  # r.array<-c(0.25, 0.5, 1, 2, 3, 4, 5, 10, 15, 20) 
                                  
                                  r.array<-c(20, 22, 25, 30, 35, 40, 45, 50, 60, 65)

                                  count.r<-2

                                  r.neighbor.sphere<-r.array[count.r]

                                  number.of.neighbors<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j.x<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j.y<-matrix(0, nrow=(track.length[i]), ncol=1)	
				    
				    min.distance.i.from.j.z<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    ###  Here we deal with the positions. Note that neighbors can be defined for each position.
				    
				  
                   	
                   	           for(k in 1: (track.length[i])){

				    		    distance.i.from.j<-NULL
				    		    
				    		    distance.i.from.j.x<-NULL
				    		    
				    		    distance.i.from.j.y<-NULL
				    		    
				    		    distance.i.from.j.z<-NULL
				    		    
				    		    l<-1
				    		    
				    		    t.i<-data_cell_motility[index.i[k], 7]
				    		    
				    		    center<-c(data_cell_motility[index.i[k], 1], data_cell_motility[index.i[k], 2], data_cell_motility[index.i[k], 3])
				    		    
						    	for( j in 1: number.of.cells_1){
						    		
						    		if(j!=i){
						    			
						    	        index.j<-which(data_cell_motility[, 8]==id_1[j])
						    	        
						    	        t.j<-data_cell_motility[index.j, 7]
						    	        
						    			t.inter<-intersect(t.i, t.j)
						    			
						    									    			
						    			if(length(t.inter)>0){
						    				
						    				t.initial<-t.inter
						    				
								    		t.initial.index<-index.j[which(data_cell_motility[index.j, 7]==t.initial)]
								    		
								    		distance.i.from.j[l]<-sqrt(sum((center-data_cell_motility[t.initial.index, 1:3])^2))
								    		
								    		distance.i.from.j.x[l]<-center[1]-data_cell_motility[t.initial.index, 1]
								    		
								    		distance.i.from.j.y[l]<-center[2]-data_cell_motility[t.initial.index, 2]
								    		
								    		distance.i.from.j.z[l]<-center[3]-data_cell_motility[t.initial.index, 3]	
								    										    	
								    		t.count<-0
								    		
					    			    	if(distance.i.from.j[l]<= r.neighbor.sphere){
					    			    			
					    			    		     number.of.neighbors[k]<-number.of.neighbors[k]+1
					    			    		
					    			    	}
					    			    		
					    			    		l<-l+1
					    			    		
					    			    		   	}
					    			    		   	
					    			    if(l==1){
					    			    	
					    			    	distance.i.from.j[l]<-0   ####  This zero means i and j are not in the same time span!
					    			    	
					    			    	distance.i.from.j.x[l]<-0
					    			    	
					    			    	distance.i.from.j.y[l]<-0
					    			    	
					    			    	distance.i.from.j.z[l]<-0
					    			    	
					    			    }		 
					    			      	
					                 }
					
				    			}
				    			
				    			distance.i.from.j<-distance.i.from.j[abs(distance.i.from.j)>0]
				    			
				    			min.distance.i.from.j[k]<-min(distance.i.from.j)
				    			
				    			distance.i.from.j.x<-distance.i.from.j.x[abs(distance.i.from.j.x)>0]
				    			
				    			min.distance.i.from.j.x[k]<-distance.i.from.j.x[which(abs(distance.i.from.j.x)==min(abs(distance.i.from.j.x)))[1]]
				    			
				    			distance.i.from.j.y<-distance.i.from.j.y[abs(distance.i.from.j.y)>0]
				    			
				    			min.distance.i.from.j.y[k]<-distance.i.from.j.y[which(abs(distance.i.from.j.y)==min(abs(distance.i.from.j.y)))[1]]
				    			
				    			distance.i.from.j.z<-distance.i.from.j.z[abs(distance.i.from.j.z)>0]
				    			
				    			min.distance.i.from.j.z[k]<-distance.i.from.j.z[which(abs(distance.i.from.j.z)==min(abs(distance.i.from.j.z)))[1]]
				    			
}


                       
				    
				    increment.x<-increment.y<-increment.z<-NULL

                                  for(j in 2: length(index.i)){
	        	
	        	                increment.x[j-1]<-data_cell_motility[index.i[j], 1]-data_cell_motility[index.i[j-1], 1]
	        	
	        	                increment.y[j-1]<-data_cell_motility[index.i[j], 2]-data_cell_motility[index.i[j-1], 2]
	        	
	        	                increment.z[j-1]<-data_cell_motility[index.i[j], 3]-data_cell_motility[index.i[j-1], 3]
	        	
	                           }

                   increment.trajectory<-matrix(c(increment.x, increment.y, increment.z), ncol=3, nrow=(track.length[i]-1))
                   
                   
                   number.of.neighbors.full.data[index.i]<-number.of.neighbors

                   angular.change.complete.full.data[index.i[2:(track.length[i]-1)]]<-u$psi.component

                   radius.complete.full.data[index.i[2:(track.length[i]-1)]]<-(u$radius[2:(track.length[i]-1)])

               
                   }
                   
                   
                
                
                data_cell_motility<-rbind(data_cell_motility_1, data_cell_motility_3,  data_cell_motility_6)

				number.of.cells<-number.of.cells_1+number.of.cells_3+number.of.cells_6



                id_3<-seq(1000000000, (1000000000+number.of.cells_3-1), 1)+number.of.cells_1

				
                track.length<-NULL
                
                                
				
				for(i in 1: number.of.cells_3){


                    # i<-cluster[[1]][cluster.counter]
				
				    
				    index.i<-which(data_cell_motility[, 8]==id_3[i])

				    
				    track.length[i]<-length(index.i)
				    
				    				    				    
				    data_i<-matrix(c(data_cell_motility[index.i, 7], data_cell_motility[index.i, 1], data_cell_motility[index.i, 2], data_cell_motility[index.i, 3]), nrow=length(index.i), ncol=4)
				    				
				    x<-tracks(data_i)
				        				        
                                  u<-radius.phi.theta.psi(trajectory=x$'1', trajectory.num=i) 

                                  # r.array<-c(0.25, 0.5, 1, 2, 3, 4, 5, 10, 15, 20) 
                                  
                                  r.array<-c(20, 22, 25, 30, 35, 40, 45, 50, 60, 65)

                                  count.r<-2

                                  r.neighbor.sphere<-r.array[count.r]

                                  number.of.neighbors<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j.x<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j.y<-matrix(0, nrow=(track.length[i]), ncol=1)	
				    
				    min.distance.i.from.j.z<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    ###  Here we deal with the positions. Note that neighbors can be defined for each position.
				    
				  
                   	
                   	           for(k in 1: (track.length[i])){

				    		    distance.i.from.j<-NULL
				    		    
				    		    distance.i.from.j.x<-NULL
				    		    
				    		    distance.i.from.j.y<-NULL
				    		    
				    		    distance.i.from.j.z<-NULL
				    		    
				    		    l<-1
				    		    
				    		    t.i<-data_cell_motility[index.i[k], 7]
				    		    
				    		    center<-c(data_cell_motility[index.i[k], 1], data_cell_motility[index.i[k], 2], data_cell_motility[index.i[k], 3])
				    		    
						    	for( j in 1: number.of.cells_3){
						    		
						    		if(j!=i){
						    			
						    	        index.j<-which(data_cell_motility[, 8]==id_3[j])
						    	        
						    	        t.j<-data_cell_motility[index.j, 7]
						    	        
						    			t.inter<-intersect(t.i, t.j)
						    			
						    									    			
						    			if(length(t.inter)>0){
						    				
						    				t.initial<-t.inter
						    				
								    		t.initial.index<-index.j[which(data_cell_motility[index.j, 7]==t.initial)]
								    		
								    		distance.i.from.j[l]<-sqrt(sum((center-data_cell_motility[t.initial.index, 1:3])^2))
								    		
								    		distance.i.from.j.x[l]<-center[1]-data_cell_motility[t.initial.index, 1]
								    		
								    		distance.i.from.j.y[l]<-center[2]-data_cell_motility[t.initial.index, 2]
								    		
								    		distance.i.from.j.z[l]<-center[3]-data_cell_motility[t.initial.index, 3]	
								    										    	
								    		t.count<-0
								    		
					    			    	if(distance.i.from.j[l]<= r.neighbor.sphere){
					    			    			
					    			    		     number.of.neighbors[k]<-number.of.neighbors[k]+1
					    			    		
					    			    	}
					    			    		
					    			    		l<-l+1
					    			    		
					    			    		   	}
					    			    		   	
					    			    if(l==1){
					    			    	
					    			    	distance.i.from.j[l]<-0   ####  This zero means i and j are not in the same time span!
					    			    	
					    			    	distance.i.from.j.x[l]<-0
					    			    	
					    			    	distance.i.from.j.y[l]<-0
					    			    	
					    			    	distance.i.from.j.z[l]<-0
					    			    	
					    			    }		 
					    			      	
					                 }
					
				    			}
				    			
				    			distance.i.from.j<-distance.i.from.j[abs(distance.i.from.j)>0]
				    			
				    			min.distance.i.from.j[k]<-min(distance.i.from.j)
				    			
				    			distance.i.from.j.x<-distance.i.from.j.x[abs(distance.i.from.j.x)>0]
				    			
				    			min.distance.i.from.j.x[k]<-distance.i.from.j.x[which(abs(distance.i.from.j.x)==min(abs(distance.i.from.j.x)))[1]]
				    			
				    			distance.i.from.j.y<-distance.i.from.j.y[abs(distance.i.from.j.y)>0]
				    			
				    			min.distance.i.from.j.y[k]<-distance.i.from.j.y[which(abs(distance.i.from.j.y)==min(abs(distance.i.from.j.y)))[1]]
				    			
				    			distance.i.from.j.z<-distance.i.from.j.z[abs(distance.i.from.j.z)>0]
				    			
				    			min.distance.i.from.j.z[k]<-distance.i.from.j.z[which(abs(distance.i.from.j.z)==min(abs(distance.i.from.j.z)))[1]]
				    			
}


                       
				    
				    increment.x<-increment.y<-increment.z<-NULL

                                  for(j in 2: length(index.i)){
	        	
	        	                increment.x[j-1]<-data_cell_motility[index.i[j], 1]-data_cell_motility[index.i[j-1], 1]
	        	
	        	                increment.y[j-1]<-data_cell_motility[index.i[j], 2]-data_cell_motility[index.i[j-1], 2]
	        	
	        	                increment.z[j-1]<-data_cell_motility[index.i[j], 3]-data_cell_motility[index.i[j-1], 3]
	        	
	                           }

                   increment.trajectory<-matrix(c(increment.x, increment.y, increment.z), ncol=3, nrow=(track.length[i]-1))
                   
                   
  
                   number.of.neighbors.full.data[index.i]<-number.of.neighbors

                   angular.change.complete.full.data[index.i[2:(track.length[i]-1)]]<-u$psi.component

                   radius.complete.full.data[index.i[2:(track.length[i]-1)]]<-(u$radius[2:(track.length[i]-1)])

               
                   }


                
                data_cell_motility<-rbind(data_cell_motility_1, data_cell_motility_3,  data_cell_motility_6)

				number.of.cells<-number.of.cells_1+number.of.cells_3+number.of.cells_6


				
                track.length<-NULL
                

                id_6<-seq(1000000000, (1000000000+number.of.cells_6-1), 1)+number.of.cells_3+number.of.cells_1

				
				for(i in 1: number.of.cells_6){
				
				    
				    index.i<-which(data_cell_motility[, 8]==id_6[i])

				    
				    track.length[i]<-length(index.i)
				    
				    				    				    
				    data_i<-matrix(c(data_cell_motility[index.i, 7], data_cell_motility[index.i, 1], data_cell_motility[index.i, 2], data_cell_motility[index.i, 3]), nrow=length(index.i), ncol=4)
				    				
				    x<-tracks(data_i)
				        				        
                                  u<-radius.phi.theta.psi(trajectory=x$'1', trajectory.num=i) 

                                  # r.array<-c(0.25, 0.5, 1, 2, 3, 4, 5, 10, 15, 20) 
                                  
                                  r.array<-c(20, 22, 25, 30, 35, 40, 45, 50, 60, 65)

                                  count.r<-2

                                  r.neighbor.sphere<-r.array[count.r]

                                  number.of.neighbors<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j.x<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j.y<-matrix(0, nrow=(track.length[i]), ncol=1)	
				    
				    min.distance.i.from.j.z<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    ###  Here we deal with the positions. Note that neighbors can be defined for each position.
				    
				  
                   	
                   	           for(k in 1: (track.length[i])){

				    		    distance.i.from.j<-NULL
				    		    
				    		    distance.i.from.j.x<-NULL
				    		    
				    		    distance.i.from.j.y<-NULL
				    		    
				    		    distance.i.from.j.z<-NULL
				    		    
				    		    l<-1
				    		    
				    		    t.i<-data_cell_motility[index.i[k], 7]
				    		    
				    		    center<-c(data_cell_motility[index.i[k], 1], data_cell_motility[index.i[k], 2], data_cell_motility[index.i[k], 3])
				    		    
						    	for( j in 1: number.of.cells_6){
						    		
						    		if(j!=i){
						    			
						    	        index.j<-which(data_cell_motility[, 8]==id_6[j])
						    	        
						    	        t.j<-data_cell_motility[index.j, 7]
						    	        
						    			t.inter<-intersect(t.i, t.j)
						    			
						    									    			
						    			if(length(t.inter)>0){
						    				
						    				t.initial<-t.inter
						    				
								    		t.initial.index<-index.j[which(data_cell_motility[index.j, 7]==t.initial)]
								    		
								    		distance.i.from.j[l]<-sqrt(sum((center-data_cell_motility[t.initial.index, 1:3])^2))
								    		
								    		distance.i.from.j.x[l]<-center[1]-data_cell_motility[t.initial.index, 1]
								    		
								    		distance.i.from.j.y[l]<-center[2]-data_cell_motility[t.initial.index, 2]
								    		
								    		distance.i.from.j.z[l]<-center[3]-data_cell_motility[t.initial.index, 3]	
								    										    	
								    		t.count<-0
								    		
					    			    	if(distance.i.from.j[l]<= r.neighbor.sphere){
					    			    			
					    			    		     number.of.neighbors[k]<-number.of.neighbors[k]+1
					    			    		
					    			    	}
					    			    		
					    			    		l<-l+1
					    			    		
					    			    		   	}
					    			    		   	
					    			    if(l==1){
					    			    	
					    			    	distance.i.from.j[l]<-0   ####  This zero means i and j are not in the same time span!
					    			    	
					    			    	distance.i.from.j.x[l]<-0
					    			    	
					    			    	distance.i.from.j.y[l]<-0
					    			    	
					    			    	distance.i.from.j.z[l]<-0
					    			    	
					    			    }		 
					    			      	
					                 }
					
				    			}
				    			
				    			distance.i.from.j<-distance.i.from.j[abs(distance.i.from.j)>0]
				    			
				    			min.distance.i.from.j[k]<-min(distance.i.from.j)
				    			
				    			distance.i.from.j.x<-distance.i.from.j.x[abs(distance.i.from.j.x)>0]
				    			
				    			min.distance.i.from.j.x[k]<-distance.i.from.j.x[which(abs(distance.i.from.j.x)==min(abs(distance.i.from.j.x)))[1]]
				    			
				    			distance.i.from.j.y<-distance.i.from.j.y[abs(distance.i.from.j.y)>0]
				    			
				    			min.distance.i.from.j.y[k]<-distance.i.from.j.y[which(abs(distance.i.from.j.y)==min(abs(distance.i.from.j.y)))[1]]
				    			
				    			distance.i.from.j.z<-distance.i.from.j.z[abs(distance.i.from.j.z)>0]
				    			
				    			min.distance.i.from.j.z[k]<-distance.i.from.j.z[which(abs(distance.i.from.j.z)==min(abs(distance.i.from.j.z)))[1]]
				    			
}


                       
				    
				    increment.x<-increment.y<-increment.z<-NULL

                                  for(j in 2: length(index.i)){
	        	
	        	                increment.x[j-1]<-data_cell_motility[index.i[j], 1]-data_cell_motility[index.i[j-1], 1]
	        	
	        	                increment.y[j-1]<-data_cell_motility[index.i[j], 2]-data_cell_motility[index.i[j-1], 2]
	        	
	        	                increment.z[j-1]<-data_cell_motility[index.i[j], 3]-data_cell_motility[index.i[j-1], 3]
	        	
	                           }

                   increment.trajectory<-matrix(c(increment.x, increment.y, increment.z), ncol=3, nrow=(track.length[i]-1))
                   
                   
  
                   number.of.neighbors.full.data[index.i]<-number.of.neighbors

                   angular.change.complete.full.data[index.i[2:(track.length[i]-1)]]<-u$psi.component

                   radius.complete.full.data[index.i[2:(track.length[i]-1)]]<-(u$radius[2:(track.length[i]-1)])


               
                   }

				
				
				
				
				
				
				
				index.neighbor.cluster.1<-NULL
				
				index.angular.change.cluster.1<-NULL
				
				for(i in 1:length(cluster.1)){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.1[i]])
					
					
					 index.neighbor.cluster.1<-c(index.neighbor.cluster.1, index.i[1:(length(index.i)-1)])
					 
					 index.angular.change.cluster.1<-c(index.angular.change.cluster.1, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				index.neighbor.cluster.2<-NULL
				
				index.angular.change.cluster.2<-NULL
				
				for(i in 1:length(cluster.2)){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.2[i]])
					
					
					 index.neighbor.cluster.2<-c(index.neighbor.cluster.2, index.i[1:(length(index.i)-1)])
					 
					 index.angular.change.cluster.2<-c(index.angular.change.cluster.2, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				index.neighbor.cluster.3<-NULL
				
				index.angular.change.cluster.3<-NULL
				
				for(i in 1:length(cluster.3)){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.3[i]])
					
					
					 index.neighbor.cluster.3<-c(index.neighbor.cluster.3, index.i[1:(length(index.i)-1)])
					 
					 index.angular.change.cluster.3<-c(index.angular.change.cluster.3, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				index.neighbor.cluster.4<-NULL
				
				index.angular.change.cluster.4<-NULL
				
				for(i in 1:length(cluster.4)){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.4[i]])
					
					
					 index.neighbor.cluster.4<-c(index.neighbor.cluster.4, index.i[1:(length(index.i)-1)])
					 
					 index.angular.change.cluster.4<-c(index.angular.change.cluster.4, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				
				data_cell_cluster_1<-data_cell_motility[index.neighbor.cluster.1, ]
				
				data_cell_cluster_2<-data_cell_motility[index.neighbor.cluster.2, ]
				
				data_cell_cluster_3<-data_cell_motility[index.neighbor.cluster.3, ]
				
				data_cell_cluster_4<-data_cell_motility[index.neighbor.cluster.4, ]

########  Note that you have already removed the end points. Also, the length of trajectories is shorter.
          
                
		        track.length.cluster.1<-NULL
	             
		        a.cluster.1<-NULL
		        
		        b.cluster.1<-NULL
		        
		        for(i in 1: length(cluster.1)){
	                    
	                     index.i<-which(data_cell_cluster_1[, 8]==id[cluster.1[i]])
		                 
		                 track.length.cluster.1[i]<-length(index.i)
		                 
		                 if(i==1) {
		                 	
		                 	a.cluster.1[i]<-1
		                 	
		                 	b.cluster.1[i]<-track.length.cluster.1[i]
		                 	
		                 	}
		                 if(i>1){
		                 	
		                 	a.cluster.1[i]<-sum(track.length.cluster.1[1:(i-1)])+1
		                 	
		                 	b.cluster.1[i]<-sum(track.length.cluster.1[1:i])
		                 			                 	
		                 	}
		                  
                        
                     }
		                
                        total.n.cluster.1<-sum(track.length.cluster.1)-length(track.length.cluster.1)
                       
                        counter.cluster.1<-a.cluster.1+1

######################################################                          
                
		        track.length.cluster.2<-NULL
	             
		        a.cluster.2<-NULL
		        
		        b.cluster.2<-NULL
		        
		        for(i in 1: length(cluster.2)){
	                    
	                     index.i<-which(data_cell_cluster_2[, 8]==id[cluster.2[i]])
		                 
		                 track.length.cluster.2[i]<-length(index.i)
		                 
		                 if(i==1) {
		                 	
		                 	a.cluster.2[i]<-1
		                 	
		                 	b.cluster.2[i]<-track.length.cluster.2[i]
		                 	
		                 	}
		                 if(i>1){
		                 	
		                 	a.cluster.2[i]<-sum(track.length.cluster.2[1:(i-1)])+1
		                 	
		                 	b.cluster.2[i]<-sum(track.length.cluster.2[1:i])
		                 			                 	
		                 	}
		                  
                        
                     }
		                
                        total.n.cluster.2<-sum(track.length.cluster.2)-length(track.length.cluster.2)
                       
                        counter.cluster.2<-a.cluster.2+1


################################################                                 
                
		        track.length.cluster.3<-NULL
	             
		        a.cluster.3<-NULL
		        
		        b.cluster.3<-NULL
		        
		        for(i in 1: length(cluster.3)){
	                    
	                     index.i<-which(data_cell_cluster_3[, 8]==id[cluster.3[i]])
		                 
		                 track.length.cluster.3[i]<-length(index.i)
		                 
		                 if(i==1) {
		                 	
		                 	a.cluster.3[i]<-1
		                 	
		                 	b.cluster.3[i]<-track.length.cluster.3[i]
		                 	
		                 	}
		                 if(i>1){
		                 	
		                 	a.cluster.3[i]<-sum(track.length.cluster.3[1:(i-1)])+1
		                 	
		                 	b.cluster.3[i]<-sum(track.length.cluster.3[1:i])
		                 			                 	
		                 	}
		                  
                        
                     }
		                
                        total.n.cluster.3<-sum(track.length.cluster.3)-length(track.length.cluster.3)
                       
                        counter.cluster.3<-a.cluster.3+1



###########################################################
                                   
                
		        track.length.cluster.4<-NULL
	             
		        a.cluster.4<-NULL
		        
		        b.cluster.4<-NULL
		        
		        for(i in 1: length(cluster.4)){
	                    
	                     index.i<-which(data_cell_cluster_4[, 8]==id[cluster.4[i]])
		                 
		                 track.length.cluster.4[i]<-length(index.i)
		                 
		                 if(i==1) {
		                 	
		                 	a.cluster.4[i]<-1
		                 	
		                 	b.cluster.4[i]<-track.length.cluster.4[i]
		                 	
		                 	}
		                 if(i>1){
		                 	
		                 	a.cluster.4[i]<-sum(track.length.cluster.4[1:(i-1)])+1
		                 	
		                 	b.cluster.4[i]<-sum(track.length.cluster.4[1:i])
		                 			                 	
		                 	}
		                  
                        
                     }
		                
                        total.n.cluster.4<-sum(track.length.cluster.4)-length(track.length.cluster.4)
                       
                        counter.cluster.4<-a.cluster.4+1


################################################


                    length.angular.change.up.to.cluster.1<-sum(track.length.cluster.1)-length(track.length.cluster.1)
                    
                    length.angular.change.up.to.cluster.2<-sum(track.length.cluster.2)-length(track.length.cluster.2)+length.angular.change.up.to.cluster.1
                    
                    length.angular.change.up.to.cluster.3<-sum(track.length.cluster.3)-length(track.length.cluster.3)+length.angular.change.up.to.cluster.2
                    
                    length.angular.change.up.to.cluster.4<-sum(track.length.cluster.4)-length(track.length.cluster.4)+length.angular.change.up.to.cluster.3



                    t<-which(angular.change.complete.full.data==0)

                    angular.change.complete.full.data[t]<-0.00001


                    psi.cluster<-angular.change.complete.full.data[index.angular.change.cluster.1]
                    
                    alpha_1_prior<-100*(mean(psi.cluster>=(pi/2))+0.01)   

                    alpha_2_prior<-100*(mean((psi.cluster<(pi/2)))+0.01)
                    
 



                    number.of.neighbors.full.data.matrix<-matrix(c(number.of.neighbors.full.data, rep(0, length(number.of.neighbors.full.data))), ncol=2, nrow=length(number.of.neighbors.full.data), byrow=F)

                    number.of.neighbors.full.data.matrix[index.angular.change.cluster.1, 2]<-1

                    number.of.neighbors.full.data.matrix[index.angular.change.cluster.2, 2]<-2

                    number.of.neighbors.full.data.matrix[index.angular.change.cluster.3, 2]<-3

                    number.of.neighbors.full.data.matrix[index.angular.change.cluster.4, 2]<-4







                    

                    t<-which(number.of.neighbors.full.data.matrix[, 2]==1)                   

                    psi.cluster<-angular.change.complete.full.data[index.angular.change.cluster.1]
                    
                    alpha_1_prior<-100*(mean(psi.cluster>=(pi/2))+0.01)   

                    alpha_2_prior<-100*(mean((psi.cluster<(pi/2)))+0.01)
                    
 



                    number.of.neighbors.full.data.matrix<-matrix(c(number.of.neighbors.full.data, rep(0, length(number.of.neighbors.full.data))), ncol=2, nrow=length(number.of.neighbors.full.data), byrow=F)

                    number.of.neighbors.full.data.matrix[index.angular.change.cluster.1, 2]<-1

                    number.of.neighbors.full.data.matrix[index.angular.change.cluster.2, 2]<-2

                    number.of.neighbors.full.data.matrix[index.angular.change.cluster.3, 2]<-3

                    number.of.neighbors.full.data.matrix[index.angular.change.cluster.4, 2]<-4







                    

                    t<-which(number.of.neighbors.full.data.matrix[, 1]==0 & number.of.neighbors.full.data.matrix[, 2]==1)

                    psi.cluster.no.neighbor<-angular.change.complete.full.data[t]


                    n.transition.1.1<-0

                    n.transition.1.2<-0
   
                    n.transition.2.1<-0

                    n.transition.2.2<-0
                   


				
				for(i in 1:length(cluster.1)){
					
					
					    index.i<-which(data_cell_motility[, 8]==id[cluster.1[i]])
							                
		                k<-track.length[i]<-length(index.i)

                    for(l in 3:(k-1)){



                        if( (angular.change.complete.full.data[index.i[l-1]]>=(pi/2)) & (angular.change.complete.full.data[index.i[l]]>=(pi/2)) & number.of.neighbors.full.data.matrix[index.i[l-1], 1]==0 ) n.transition.1.1<-n.transition.1.1+1 

                        if( (angular.change.complete.full.data[index.i[l-1]]>=(pi/2)) & (angular.change.complete.full.data[index.i[l]]<(pi/2)) & number.of.neighbors.full.data.matrix[index.i[l-1], 1]==0 ) n.transition.1.2<-n.transition.1.2+1 

 


                        if( (angular.change.complete.full.data[index.i[l-1]]<(pi/2)) & (angular.change.complete.full.data[index.i[l]]>=(pi/2)) & number.of.neighbors.full.data.matrix[index.i[l-1], 1]==0 ) n.transition.2.1<-n.transition.2.1+1 

                        if( (angular.change.complete.full.data[index.i[l-1]]<(pi/2)) & (angular.change.complete.full.data[index.i[l]]<(pi/2)) & number.of.neighbors.full.data.matrix[index.i[l-1], 1]==0 ) n.transition.2.2<-n.transition.2.2+1 

 
 

                     }
                     
                     }
                     
                     



                    alpha_P_1_1<-100*(n.transition.1.1/(n.transition.1.1+n.transition.1.2))

                    alpha_P_1_2<-100*(n.transition.1.2/(n.transition.1.1+n.transition.1.2))
                    
 
 
                    alpha_P_2_1<-100*(n.transition.2.1/(n.transition.2.1+n.transition.2.2))

                    alpha_P_2_2<-100*(n.transition.2.2/(n.transition.2.1+n.transition.2.2))


                 
                    alpha_P_1_informative<-matrix(c(alpha_P_1_1, alpha_P_1_2, alpha_P_2_1, alpha_P_2_2), nrow=2, ncol=2, byrow=T)




                    t<-which(number.of.neighbors.full.data.matrix[, 1]!=0)

                    psi.cluster.with.neighbor<-angular.change.complete.full.data[t]

                    


                    n.transition.1.1<-0

                    n.transition.1.2<-0
                    
                    n.transition.2.1<-0

                    n.transition.2.2<-0
                    


				
				for(i in 1:length(cluster.1)){
					
					
					    index.i<-which(data_cell_motility[, 8]==id[cluster.1[i]])
							                
		                k<-track.length[i]<-length(index.i)

                    for(l in 3:(k-1)){
              	
                    	
                        if( (angular.change.complete.full.data[index.i[l-1]]>=(pi/2)) & (angular.change.complete.full.data[index.i[l]]>=(pi/2)) & number.of.neighbors.full.data.matrix[index.i[l-1], 1]!=0 ) n.transition.1.1<-n.transition.1.1+1 

                        if( (angular.change.complete.full.data[index.i[l-1]]>=(pi/2)) & (angular.change.complete.full.data[index.i[l]]<(pi/2)) & number.of.neighbors.full.data.matrix[index.i[l-1], 1]!=0 ) n.transition.1.2<-n.transition.1.2+1 

 


                        if( (angular.change.complete.full.data[index.i[l-1]]<(pi/2)) & (angular.change.complete.full.data[index.i[l]]>=(pi/2)) & number.of.neighbors.full.data.matrix[index.i[l-1], 1]!=0 ) n.transition.2.1<-n.transition.2.1+1 

                        if( (angular.change.complete.full.data[index.i[l-1]]<(pi/2)) & (angular.change.complete.full.data[index.i[l]]<(pi/2)) & number.of.neighbors.full.data.matrix[index.i[l-1], 1]!=0 ) n.transition.2.2<-n.transition.2.2+1 


                     }
                     
                     }


                    alpha_P_1_1<-100*(n.transition.1.1/(n.transition.1.1+n.transition.1.2))

                    alpha_P_1_2<-100*(n.transition.1.2/(n.transition.1.1+n.transition.1.2))
                    
 
 
                    alpha_P_2_1<-100*(n.transition.2.1/(n.transition.2.1+n.transition.2.2))

                    alpha_P_2_2<-100*(n.transition.2.2/(n.transition.2.1+n.transition.2.2))

                 
                    alpha_P_2_informative<-matrix(c(alpha_P_1_1, alpha_P_1_2, alpha_P_2_1, alpha_P_2_2), nrow=2, ncol=2, byrow=T)












   
                       dat.cluster.1 <- list(N = sum(track.length.cluster.1)-length(track.length.cluster.1), angular_change_cluster= angular.change.complete.full.data[index.angular.change.cluster.1], number_of_neighbor_cluster= number.of.neighbors.full.data[index.neighbor.cluster.1], alpha_prior=c(alpha_1_prior, alpha_2_prior), alpha_P_1=alpha_P_1_informative, alpha_P_2=alpha_P_2_informative, K=2, counter_n=counter.cluster.1, a_initial_f=a.cluster.1, b_final_f=b.cluster.1, J=length(cluster[[1]]))
                    
                    
                     dat.cluster.2 <- list(N = sum(track.length.cluster.2)-length(track.length.cluster.2), angular_change_cluster= angular.change.complete.full.data[index.angular.change.cluster.2], number_of_neighbor_cluster= number.of.neighbors.full.data[index.neighbor.cluster.2], alpha_prior=c(alpha_1_prior, alpha_2_prior), alpha_P_1=alpha_P_1_informative, alpha_P_2=alpha_P_2_informative, K=2, counter_n=counter.cluster.2, a_initial_f=a.cluster.2, b_final_f=b.cluster.2, J=length(cluster[[2]]))



                     dat.cluster.3 <- list(N = sum(track.length.cluster.3)-length(track.length.cluster.3), angular_change_cluster= angular.change.complete.full.data[index.angular.change.cluster.3], number_of_neighbor_cluster= number.of.neighbors.full.data[index.neighbor.cluster.3], alpha_prior=c(alpha_1_prior, alpha_2_prior), alpha_P_1=alpha_P_1_informative, alpha_P_2=alpha_P_2_informative, K=2, counter_n=counter.cluster.3, a_initial_f=a.cluster.3, b_final_f=b.cluster.3, J=length(cluster[[3]]))

                  
                     dat.cluster.4 <- list(N = sum(track.length.cluster.4)-length(track.length.cluster.4), angular_change_cluster= angular.change.complete.full.data[index.angular.change.cluster.4], number_of_neighbor_cluster= number.of.neighbors.full.data[index.neighbor.cluster.4], alpha_prior=c(alpha_1_prior, alpha_2_prior), alpha_P_1=alpha_P_1_informative, alpha_P_2=alpha_P_2_informative, K=2, counter_n=counter.cluster.4, a_initial_f=a.cluster.4, b_final_f=b.cluster.4, J=length(cluster[[4]]))



                    save(list = ls(all=TRUE), file =".RData-22")


                    date() 
   

    
                    rstan_options(auto_write = TRUE)
                    
                    
                    options(mc.cores = parallel::detectCores())


                    
                    fit <- stan(model_code = stanmodelcode, model_name = "example",
                data = dat.cluster.1, iter = 3000, chains = 8, verbose = TRUE, control = list(adapt_delta = 0.9999, stepsize = 0.001, max_treedepth =15), cores=8)
                
                
                    fit_ss <- summary(fit) # fit_ss is a list


                    print(fit_ss$summary)
                    
                    
                    save(list = ls(all=TRUE), file = ".RData-22")


                    date() 

   

                                          
                
                   





                    
                   



                    
                   




                    
                   