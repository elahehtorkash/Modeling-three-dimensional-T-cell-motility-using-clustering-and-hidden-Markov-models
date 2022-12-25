library(rstan)

remove(list=ls())


setwd("/global/scratch/torkashe/Clustering-Results/cluster-1-3-states-HMM/")

load(".RData")




cluster.c<-1


fit.ext<-extract(fit)


Rep.angle<-(length(fit.ext[[1]])/3)


counter.c<-as.matrix(rep(cluster.c, Rep.angle))


prior_pi<-as.matrix(fit.ext[[1]])


P_1<-fit.ext[[2]]


# P_2<-fit.ext[[2]][, 2, , ]


alpha_per<-as.matrix(fit.ext[[3]])


beta_per<-as.matrix(fit.ext[[4]])


alpha_non_per<-as.matrix(fit.ext[[5]])


beta_non_per<-as.matrix(fit.ext[[6]])


alpha_transition<-beta_transition <-as.matrix(fit.ext[[7]])

parameter.fit.angle<-cbind(prior_pi, P_1[, 1, 1], P_1[, 1, 2], P_1[, 1, 3], P_1[, 2, 1], P_1[, 2, 2], P_1[, 2, 3], P_1[, 3, 1], P_1[, 3, 2], P_1[, 3, 3], alpha_per, beta_per, alpha_non_per, beta_non_per, alpha_transition , beta_transition, counter.c)

###################################################

######################################################################################################################################################################################################################################################################################  Held-Out Data  #########################################################

   
                

                data_cell_motility_2<-read.csv("Position_2.csv", header=TRUE)

				dim(data_cell_motility_2)


				  In order to have an understanding of the number of T-cells observed in the study, we use the following comment


				number.of.cells_2<-length(unique(data_cell_motility_2[, 8]))



                
                data_cell_motility_4<-read.csv("Position_Longest.csv", header=TRUE)

				dim(data_cell_motility_4)

				data_cell_motility_4[, 8]<-data_cell_motility_4[, 8]+(number.of.cells_2)
				

				  In order to have an understanding of the number of T-cells observed in the study, we use the following comment

                
                number.of.cells_4<-length(unique(data_cell_motility_4[, 8]))
                


    
    
                
                data_cell_motility<-rbind(data_cell_motility_2, data_cell_motility_4)

				number.of.cells<-number.of.cells_2+number.of.cells_4

	              In order to have an understanding of the number of T-cells observed in the study, we use the following comment

             
				number.of.cells<-length(unique(data_cell_motility[, 8]))

				id<-seq(1000000000, (1000000000+number.of.cells-1), 1)
				
				track.length<-NULL
				
                for(i in 1: number.of.cells){
                	
                	track.length[i]<-length(which(data_cell_motility[, 8]==id[i]))
                	
                }

                number.of.positions.total<-sum(track.length)
                
                
                
##########################################
##########################################
				
				
				angular.change.complete.full.data<-matrix(0, nrow= number.of.positions.total, ncol=1)
				
				
				number.of.neighbors.full.data<-matrix(0, nrow= number.of.positions.total, ncol=1)

                
                radius.complete.full.data<-matrix(0, nrow= number.of.positions.total, ncol=1)


						
		tracks<-list() 
				    
				    
				    												
        This gives the trajectories in space. I suggest you to run this code. By looking at the trajectories, it make sense to cluster them in different groups.

       track.length<-NULL
              
       hierarchical.clustering.criteria<-matrix(0, nrow= number.of.cells, ncol=length(points.for.estimation.angle)+1)
      
       radius.trajectory<-NULL
       
       theta.trajectory<-NULL
       
       phi.trajectory<-NULL
              
       r.domain.ratio<-rep(0, number.of.cells)
       
       p.persistent.nonpersistent<-matrix(0, nrow= number.of.cells, ncol=2)
                     
        edited theta and phi
       
       theta.breaks<-seq(0, 3.5, 0.5)
       
       theta.breaks<-seq(0, 1.6, 0.4)
       
       phi.breaks<-seq(-3.5, 3.5, 0.5)
       
       x.cell<-y.cell<-z.cell<-time.cell<-NULL 
       
       increment.x<-increment.y<-increment.z<-NULL

       B.boot<-10000
                                          
       persistent.probability<-matrix(0, nrow= number.of.cells, ncol=length(r.array))
       
       conditional.persistent.probability<-matrix(0, nrow= number.of.cells, ncol=length(r.array))
                             
                      
#######################################################################################################
###################################################  Spherical representation and angular change 

               
    radius.phi.theta.psi<-function(trajectory, trajectory.num){
					
					length.trajectory<-length(trajectory[, 1])
					
					length.sub.trajectories<-NULL
					
					theta<-NULL
					
					phi<-NULL
					
					for(i in 1: (length.trajectory-1)){
						
					   	length.sub.trajectories[i]<-sqrt((trajectory[(i+1), 2]-trajectory[i, 2])^2+(trajectory[(i+1), 3]-trajectory[i, 3])^2+(trajectory[(i+1), 4]-trajectory[i, 4])^2)
					   	
					   	theta[i]<-acos((trajectory[(i+1), 4]-trajectory[i, 4])/length.sub.trajectories[i])
					   	
					   	  Considering the revised theta under my assumption
					   	
					    theta[i]<-acos(abs((trajectory[(i+1), 4]-trajectory[i, 4])/length.sub.trajectories[i]))
				
					   	
					   	if(theta[i]>(pi/2)) theta[i]<-pi-theta[i]
					   	
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
						
						print("u"); print(u)
						
						v<-c((trajectory[(i+2), 2]-trajectory[i+1, 2]), (trajectory[(i+2), 3]-trajectory[i+1, 3]), (trajectory[(i+2), 4]-trajectory[i+1, 4]))
						print("v"); print(v)
						
						u.cross.product.v<-c((u[2]*v[3]-u[3]*v[2]), (u[3]*v[1]-u[1]*v[3]), (u[1]*v[2]-u[2]*v[1]))
						
						u.dot.product.v<-u[1]*v[1]+u[2]*v[2]+u[3]*v[3]
						
						length.u.cross.product.v<-sqrt((u.cross.product.v[1])^2+(u.cross.product.v[2])^2+(u.cross.product.v[3])^2)
						
						length.u<-sqrt((u[1])^2+(u[2])^2+(u[3])^2)
						
						length.v<-sqrt((v[1])^2+(v[2])^2+(v[3])^2)
						
						angle.between.sub.trajectories[i]<-asin(length.u.cross.product.v/(length.u* length.v))
						
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


 
								
				 This gives the trajectories in space. I suggest you to run this code. By looking at the trajectories, it make sense to cluster them in different groups.
				
				track.length<-NULL
				       
		        x.cell<-NULL
		       
		        y.cell<-NULL
		        
		        z.cell<-NULL 
		       
		        time.cell<-NULL
		        
		        total.phi<-NULL
			    
			    total.psi<-NULL
			    
			    total.theta<-NULL
			      
			    total.radius<-NULL  
				
                increment.x.cell<-NULL
		       
		        increment.y.cell<-NULL
		        
		        increment.z.cell<-NULL 

                times.present.in.study<-NULL
                        
                    	
                count.r<-3
                   	
                r.neighbor.sphere<-r.array[count.r]

               
				
				for(i in 1: number.of.cells){
				
				    
				    index.i<-which(data_cell_motility[, 8]==id[i])
				    
				    track.length[i]<-length(index.i)
				    
				    times.present.in.study<-c(times.present.in.study, length(index.i))
				    
				    x.cell<-c(x.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 1])
				    
				    y.cell<-c(y.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 2])
				    
				    z.cell<-c(z.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 3])
				    
				    time.cell<-c(time.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 7])
				    
				    data_i<-matrix(c(data_cell_motility[index.i, 7], data_cell_motility[index.i, 1], data_cell_motility[index.i, 2], data_cell_motility[index.i, 3]), nrow=length(index.i), ncol=4)
				    
				    				    
        
        #############################  Finding the corresponding clustering criteria
        
        
                    increment.code<-matrix(0, nrow=3, ncol=1)
        
                    increment.x<-increment.x.total.simulated[index.i[1:(l-1)]]
                    
                    increment.y<-increment.y.total.simulated[index.i[1:(l-1)]]
                    
                    increment.z<-increment.z.total.simulated[index.i[1:(l-1)]]
                    
                    
                    				    				    
				    increment.x<-increment.y<-increment.z<-NULL
				    
				    increment.code<-matrix(0, nrow=3, ncol=1)

                    for(j in 2: length(index.i)){
	        	
	        	             increment.x[j-1]<-data_cell_motility[index.i[j], 1]-data_cell_motility[index.i[j-1], 1]
	        	
	        	             increment.y[j-1]<-data_cell_motility[index.i[j], 2]-data_cell_motility[index.i[j-1], 2]
	        	
	        	             increment.z[j-1]<-data_cell_motility[index.i[j], 3]-data_cell_motility[index.i[j-1], 3]
	        	
	                }

        
        	                
	                	                                
	                B=10000; 
	                
	                n = length(increment.x)
	                
	                
	                normal.samples = matrix( rnorm(n*B, mean(increment.x), sd(increment.x)), B, n)
                   
                    normal.statistics = apply(normal.samples, 1, quantile)
                    
                    first.quantile<-normal.statistics[2, ]
                    
                    third.quantile<-normal.statistics[4, ]
                    
                    min.boot<-normal.statistics[1, ]
                    
                    max.boot<-normal.statistics[5, ]
                    
                    cm_normal<-(third.quantile-first.quantile)/(max.boot-min.boot)
                    
                    interval = quantile(cm_normal, c(0.025, 0.975))
                    

	                cm_x_increments<-(summary(increment.x)[5]-summary(increment.x)[2])/(summary(increment.x)[6]-summary(increment.x)[1])

                     |  cm_b > interval[2]	                
	                
 	                if(cm_x_increments < interval[1] ){
	                	
	                	increment.code[1]<-1
	                	
	                }
	                
	                
	                
	                
	                normal.samples = matrix( rnorm(n*B, mean(increment.y), sd(increment.y)), B, n)
                   
                    normal.statistics = apply(normal.samples, 1, quantile)
                    
                    first.quantile<-normal.statistics[2, ]
                    
                    third.quantile<-normal.statistics[4, ]
                    
                    min.boot<-normal.statistics[1, ]
                    
                    max.boot<-normal.statistics[5, ]
                    
                    cm_normal<-(third.quantile-first.quantile)/(max.boot-min.boot)
                    
                    interval = quantile(cm_normal, c(0.025, 0.975))	                
	                
	                cm_y_increments<-(summary(increment.y)[5]-summary(increment.y)[2])/(summary(increment.y)[6]-summary(increment.y)[1])
	                
	                
	                if(cm_y_increments < interval[1]  ){
	                	
	                	increment.code[2]<-1
	                	
	                }

	                	                
	                
	               	                                
	                B=10000; 
	                
	                n = length(increment.z)
	                
	                
	                normal.samples = matrix( rnorm(n*B, mean(increment.z), sd(increment.z)), B, n)
                   
                    normal.statistics = apply(normal.samples, 1, quantile)
                    
                    first.quantile<-normal.statistics[2, ]
                    
                    third.quantile<-normal.statistics[4, ]
                    
                    min.boot<-normal.statistics[1, ]
                    
                    max.boot<-normal.statistics[5, ]
                    
                    cm_normal<-(third.quantile-first.quantile)/(max.boot-min.boot)
	                
	                interval = quantile(cm_normal, c(0.025, 0.975))
	                
	                cm_z_increments<-(summary(increment.z)[5]-summary(increment.z)[2])/(summary(increment.z)[6]-summary(increment.z)[1])
	                
	                
	                if(cm_z_increments < interval[1]  ){
	                	
	                	increment.code[3]<-1
	                	
	                }

	                
	                
                    increment.trajectory<-matrix(c(increment.x, increment.y, increment.z), ncol=3, nrow=(track.length[i]-1))
       
                    if(sum(increment.code)>0) r.domain.ratio[i]<-1
                                      
                 	x<-tracks(data_i)
                 	
                 	u<-radius.phi.theta.psi(trajectory=x$'1', trajectory.num=i)
                 	
                 	
                 	
                   angular.change.complete.full.data[index.i[2:(track.length[i]-1)]]<-u$psi.component

                   radius.complete.full.data[index.i[2:(track.length[i]-1)]]<-(u$radius[2:(track.length[i]-1)])

        
                    
                    p.persistent.nonpersistent[i, ]<-c(sum(u$psi.component<=pi/2)/length(u$psi.component), (1-(sum(u$psi.component<=pi/2)/length(u$psi.component))))

        
                   
 
 }            
 
 
 
                     The other way around is to use the median of the new dataset in order to make clusters.                 
                    
                    n_c<-NULL
							
												
					cluster.1<-which(p.persistent.nonpersistent[, 1] >= median(p.persistent.nonpersistent[, 1]) & r.domain.ratio==1 &  track.length>=29)   ### persistent and long
					
					n_c[1]<-length(cluster.1)
					
					cluster.2<-which(p.persistent.nonpersistent[, 1] >= median(p.persistent.nonpersistent[, 1]) & r.domain.ratio==0 &  track.length>=29)   ### persistent and short 
					
					n_c[2]<-length(cluster.2)
					
					cluster.3<-which(p.persistent.nonpersistent[, 1] < median(p.persistent.nonpersistent[, 1]) & r.domain.ratio==1 &  track.length>=29)   ### non-persistent and long
					
					n_c[3]<-length(cluster.3)
					
					cluster.4<-which(p.persistent.nonpersistent[, 1] < median(p.persistent.nonpersistent[, 1]) & r.domain.ratio==0 &  track.length>=29)   ### non-persistent and short
					
					n_c[4]<-length(cluster.4)
					

                    print("n_c")
                    
                    print(n_c)
                    
                    
                    
                    
                    
###########################################################################################################################################################################################################################################  neighbors

		       		        
		      total.phi<-NULL
			    
			  total.psi<-NULL
			    
			  total.theta<-NULL
			      
			  total.radius<-NULL  
				
              increment.x.cell<-NULL
		       
		      increment.y.cell<-NULL
		        
		      increment.z.cell<-NULL 
		        
		    
		    
           				                              
                
                track.length<-NULL
                           
                              
                id_2<-seq(1000000000, (1000000000+number.of.cells_2-1), 1)

				
				for(i in 1: number.of.cells_2){
				
				    
				    index.i<-which(data_cell_motility[, 8]==id_2[i])

				    
				    track.length[i]<-length(index.i)
				    
				    				    				    
				    data_i<-matrix(c(data_cell_motility[index.i, 7], data_cell_motility[index.i, 1], data_cell_motility[index.i, 2], data_cell_motility[index.i, 3]), nrow=length(index.i), ncol=4)
				    				
				    x<-tracks(data_i)
				        				        
                                  u<-radius.phi.theta.psi(trajectory=x$'1', trajectory.num=i) 

                                  r.array<-c(20, 22, 25, 30, 35, 40, 45, 50, 60, 65)

                                  count.r<-2

                                  r.neighbor.sphere<-r.array[count.r]

                                  number.of.neighbors<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j.x<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j.y<-matrix(0, nrow=(track.length[i]), ncol=1)	
				    
				    min.distance.i.from.j.z<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				     Here we deal with the positions. Note that neighbors can be defined for each position.
				    
				  
                   	
                   	           for(k in 1: (track.length[i])){

				    		    distance.i.from.j<-NULL
				    		    
				    		    distance.i.from.j.x<-NULL
				    		    
				    		    distance.i.from.j.y<-NULL
				    		    
				    		    distance.i.from.j.z<-NULL
				    		    
				    		    l<-1
				    		    
				    		    t.i<-data_cell_motility[index.i[k], 7]
				    		    
				    		    center<-c(data_cell_motility[index.i[k], 1], data_cell_motility[index.i[k], 2], data_cell_motility[index.i[k], 3])
				    		    
						    	for( j in 1: number.of.cells_2){
						    		
						    		if(j!=i){
						    			
						    	        index.j<-which(data_cell_motility[, 8]==id_2[j])
						    	        
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
				    			
				    			distance.i.from.j.x<-distance.i.from.j[abs(distance.i.from.j.x)>0]
				    			
				    			min.distance.i.from.j.x[k]<-distance.i.from.j.x[which(abs(distance.i.from.j.x)==min(abs(distance.i.from.j.x)))[1]]
				    			
				    			distance.i.from.j.y<-distance.i.from.j[abs(distance.i.from.j.y)>0]
				    			
				    			min.distance.i.from.j.y[k]<-distance.i.from.j.y[which(abs(distance.i.from.j.y)==min(abs(distance.i.from.j.y)))[1]]
				    			
				    			distance.i.from.j.z<-distance.i.from.j[abs(distance.i.from.j.z)>0]
				    			
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
                   
                   
                
                
                data_cell_motility<-rbind(data_cell_motility_2, data_cell_motility_4)

				number.of.cells<-number.of.cells_2+number.of.cells_4



                id_4<-seq(1000000000, (1000000000+number.of.cells_4-1), 1)+number.of.cells_2

				
                track.length<-NULL
                
                                
				
				for(i in 1: number.of.cells_4){


                    i<-cluster[[1]][cluster.counter]
				
				    
				    index.i<-which(data_cell_motility[, 8]==id_4[i])

				    
				    track.length[i]<-length(index.i)
				    
				    				    				    
				    data_i<-matrix(c(data_cell_motility[index.i, 7], data_cell_motility[index.i, 1], data_cell_motility[index.i, 2], data_cell_motility[index.i, 3]), nrow=length(index.i), ncol=4)
				    				
				    x<-tracks(data_i)
				        				        
                                  u<-radius.phi.theta.psi(trajectory=x$'1', trajectory.num=i) 

                                  r.array<-c(20, 22, 25, 30, 35, 40, 45, 50, 60, 65)

                                  count.r<-2

                                  r.neighbor.sphere<-r.array[count.r]

                                  number.of.neighbors<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j.x<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j.y<-matrix(0, nrow=(track.length[i]), ncol=1)	
				    
				    min.distance.i.from.j.z<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				     Here we deal with the positions. Note that neighbors can be defined for each position.
				    
				  
                   	
                   	           for(k in 1: (track.length[i])){

				    		    distance.i.from.j<-NULL
				    		    
				    		    distance.i.from.j.x<-NULL
				    		    
				    		    distance.i.from.j.y<-NULL
				    		    
				    		    distance.i.from.j.z<-NULL
				    		    
				    		    l<-1
				    		    
				    		    t.i<-data_cell_motility[index.i[k], 7]
				    		    
				    		    center<-c(data_cell_motility[index.i[k], 1], data_cell_motility[index.i[k], 2], data_cell_motility[index.i[k], 3])
				    		    
						    	for( j in 1: number.of.cells_4){
						    		
						    		if(j!=i){
						    			
						    	        index.j<-which(data_cell_motility[, 8]==id_4[j])
						    	        
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
				    			
				    			distance.i.from.j.x<-distance.i.from.j[abs(distance.i.from.j.x)>0]
				    			
				    			min.distance.i.from.j.x[k]<-distance.i.from.j.x[which(abs(distance.i.from.j.x)==min(abs(distance.i.from.j.x)))[1]]
				    			
				    			distance.i.from.j.y<-distance.i.from.j[abs(distance.i.from.j.y)>0]
				    			
				    			min.distance.i.from.j.y[k]<-distance.i.from.j.y[which(abs(distance.i.from.j.y)==min(abs(distance.i.from.j.y)))[1]]
				    			
				    			distance.i.from.j.z<-distance.i.from.j[abs(distance.i.from.j.z)>0]
				    			
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


                    
                  number.of.neighbors.full.data <- number.of.neighbors.full.data  
                  
                  
                  angular.change.complete.full.data <- angular.change.complete.full.data
                  
                  
                  radius.complete.full.data <- radius.complete.full.data



                  data_cell_motility <- data_cell_motility


                  cluster<-list("1"=cluster.1, "2"=cluster.2, "3"=cluster.3, "4"=cluster.4)

 

#######################################################################################################################################################################################################################################

                				
				
				index.angular.change<-NULL
				
				for(i in 1: number.of.cells){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[i])
					
					
					 
					 index.angular.change<-c(index.angular.change, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				
				index.angular.change.cluster.1<-NULL
				
				for(i in 1:length(cluster.1){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.1[i]])
					
					
					 
					 index.angular.change.cluster.1<-c(index.angular.change.cluster.1, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				
				index.angular.change.cluster.2<-NULL
				
				for(i in 1:length(cluster.2)){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.2[i]])
					
					
					 
					 index.angular.change.cluster.2<-c(index.angular.change.cluster.2, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				
				index.angular.change.cluster.3<-NULL
				
				for(i in 1:length(cluster.3)){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.3[i]])
					
					
					 
					 index.angular.change.cluster.3<-c(index.angular.change.cluster.3, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				
				index.angular.change.cluster.4<-NULL
				
				for(i in 1:length(cluster.4)){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.4[i]])
					
					
					 
					 index.angular.change.cluster.4<-c(index.angular.change.cluster.4, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				index.angular.change.cluster<-c(index.angular.change.cluster.1, index.angular.change.cluster.2, index.angular.change.cluster.3, index.angular.change.cluster.4)


###########################################################################################################################################################################################################################


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


devience.infor.criterion<-function(parameter, parameter.bayes.f){
	
	t<-apply(parameter.fit.angle, 1, log.likelihood.f)
	
	Bayes<-log.likelihood.f(parameter.bayes.f)


     u<-2*(Bayes-mean(t))
    
    
    DIC.cluster<- (-2*(Bayes-u))



    print("DIC.cluster")


    print(DIC.cluster)

	
	
}


##########################################################################

# # # # # # DIC.cluster<-devience.infor.criterion(parameter.fit.angle, parameter.bayes)

# # # print(DIC.cluster)

###########################################################################################################################################################################################



##########################################################################


WAIC.f<-function(parameter){
	
	t<-apply(parameter, 1, log.likelihood.f)
	
	likelihood<-exp(t)
	
	mean.likelihood<-mean(likelihood)	
	
	if(round(mean.likelihood, 4)==0) mean.likelihood<-10^-4

    u <-2*(log(mean.likelihood)-mean(t))
    
    WAIC.cluster<-(-2*(log(mean.likelihood)-u))


    print("WAIC.cluster")


    print(WAIC.cluster)

	
	
}


########################################################################

WAIC.cluster<-WAIC.f(parameter.fit.angle)

print(WAIC.cluster)



###########################################################################################################################################################################################

library(rstan)

remove(list=ls())


setwd("/global/scratch/torkashe/Clustering-Results/cluster-2-3-states-HMM/")

load(".RData")




cluster.c<-2



fit.ext<-extract(fit)


Rep.angle<-(length(fit.ext[[1]])/3)


counter.c<-as.matrix(rep(cluster.c, Rep.angle))


prior_pi<-as.matrix(fit.ext[[1]])


P_1<-fit.ext[[2]]


# P_2<-fit.ext[[2]][, 2, , ]


alpha_per<-as.matrix(fit.ext[[3]])


beta_per<-as.matrix(fit.ext[[4]])


alpha_non_per<-as.matrix(fit.ext[[5]])


beta_non_per<-as.matrix(fit.ext[[6]])


alpha_transition<-beta_transition <-as.matrix(fit.ext[[7]])

parameter.fit.angle<-cbind(prior_pi, P_1[, 1, 1], P_1[, 1, 2], P_1[, 1, 3], P_1[, 2, 1], P_1[, 2, 2], P_1[, 2, 3], P_1[, 3, 1], P_1[, 3, 2], P_1[, 3, 3], alpha_per, beta_per, alpha_non_per, beta_non_per, alpha_transition , beta_transition, counter.c)



###################################################



######################################################################################################################################################################################################################################################################################  Held-Out Data  #########################################################

   
                

                data_cell_motility_2<-read.csv("Position_2.csv", header=TRUE)

				dim(data_cell_motility_2)


				  In order to have an understanding of the number of T-cells observed in the study, we use the following comment


				number.of.cells_2<-length(unique(data_cell_motility_2[, 8]))



                
                data_cell_motility_4<-read.csv("Position_Longest.csv", header=TRUE)

				dim(data_cell_motility_4)

				data_cell_motility_4[, 8]<-data_cell_motility_4[, 8]+(number.of.cells_2)
				

				  In order to have an understanding of the number of T-cells observed in the study, we use the following comment

                
                number.of.cells_4<-length(unique(data_cell_motility_4[, 8]))
                


    
    
                
                data_cell_motility<-rbind(data_cell_motility_2, data_cell_motility_4)

				number.of.cells<-number.of.cells_2+number.of.cells_4

	              In order to have an understanding of the number of T-cells observed in the study, we use the following comment

             
				number.of.cells<-length(unique(data_cell_motility[, 8]))

				id<-seq(1000000000, (1000000000+number.of.cells-1), 1)
				
				track.length<-NULL
				
                for(i in 1: number.of.cells){
                	
                	track.length[i]<-length(which(data_cell_motility[, 8]==id[i]))
                	
                }

                number.of.positions.total<-sum(track.length)
                
                
                
##########################################
##########################################
				
				
				angular.change.complete.full.data<-matrix(0, nrow= number.of.positions.total, ncol=1)
				
				
				number.of.neighbors.full.data<-matrix(0, nrow= number.of.positions.total, ncol=1)

                
                radius.complete.full.data<-matrix(0, nrow= number.of.positions.total, ncol=1)


						
		tracks<-list() 
				    
				    
				    												
        This gives the trajectories in space. I suggest you to run this code. By looking at the trajectories, it make sense to cluster them in different groups.

       track.length<-NULL
              
       hierarchical.clustering.criteria<-matrix(0, nrow= number.of.cells, ncol=length(points.for.estimation.angle)+1)
      
       radius.trajectory<-NULL
       
       theta.trajectory<-NULL
       
       phi.trajectory<-NULL
              
       r.domain.ratio<-rep(0, number.of.cells)
       
       p.persistent.nonpersistent<-matrix(0, nrow= number.of.cells, ncol=2)
                     
        edited theta and phi
       
       theta.breaks<-seq(0, 3.5, 0.5)
       
       theta.breaks<-seq(0, 1.6, 0.4)
       
       phi.breaks<-seq(-3.5, 3.5, 0.5)
       
       x.cell<-y.cell<-z.cell<-time.cell<-NULL 
       
       increment.x<-increment.y<-increment.z<-NULL

       B.boot<-10000
                                          
       persistent.probability<-matrix(0, nrow= number.of.cells, ncol=length(r.array))
       
       conditional.persistent.probability<-matrix(0, nrow= number.of.cells, ncol=length(r.array))
                             
                      
#######################################################################################################
###################################################  Spherical representation and angular change 

               
    radius.phi.theta.psi<-function(trajectory, trajectory.num){
					
					length.trajectory<-length(trajectory[, 1])
					
					length.sub.trajectories<-NULL
					
					theta<-NULL
					
					phi<-NULL
					
					for(i in 1: (length.trajectory-1)){
						
					   	length.sub.trajectories[i]<-sqrt((trajectory[(i+1), 2]-trajectory[i, 2])^2+(trajectory[(i+1), 3]-trajectory[i, 3])^2+(trajectory[(i+1), 4]-trajectory[i, 4])^2)
					   	
					   	theta[i]<-acos((trajectory[(i+1), 4]-trajectory[i, 4])/length.sub.trajectories[i])
					   	
					   	  Considering the revised theta under my assumption
					   	
					    theta[i]<-acos(abs((trajectory[(i+1), 4]-trajectory[i, 4])/length.sub.trajectories[i]))
				
					   	
					   	if(theta[i]>(pi/2)) theta[i]<-pi-theta[i]
					   	
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
						
						print("u"); print(u)
						
						v<-c((trajectory[(i+2), 2]-trajectory[i+1, 2]), (trajectory[(i+2), 3]-trajectory[i+1, 3]), (trajectory[(i+2), 4]-trajectory[i+1, 4]))
						print("v"); print(v)
						
						u.cross.product.v<-c((u[2]*v[3]-u[3]*v[2]), (u[3]*v[1]-u[1]*v[3]), (u[1]*v[2]-u[2]*v[1]))
						
						u.dot.product.v<-u[1]*v[1]+u[2]*v[2]+u[3]*v[3]
						
						length.u.cross.product.v<-sqrt((u.cross.product.v[1])^2+(u.cross.product.v[2])^2+(u.cross.product.v[3])^2)
						
						length.u<-sqrt((u[1])^2+(u[2])^2+(u[3])^2)
						
						length.v<-sqrt((v[1])^2+(v[2])^2+(v[3])^2)
						
						angle.between.sub.trajectories[i]<-asin(length.u.cross.product.v/(length.u* length.v))
						
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


 
								
				 This gives the trajectories in space. I suggest you to run this code. By looking at the trajectories, it make sense to cluster them in different groups.
				
				track.length<-NULL
				       
		        x.cell<-NULL
		       
		        y.cell<-NULL
		        
		        z.cell<-NULL 
		       
		        time.cell<-NULL
		        
		        total.phi<-NULL
			    
			    total.psi<-NULL
			    
			    total.theta<-NULL
			      
			    total.radius<-NULL  
				
                increment.x.cell<-NULL
		       
		        increment.y.cell<-NULL
		        
		        increment.z.cell<-NULL 

                times.present.in.study<-NULL
                        
                    	
                count.r<-3
                   	
                r.neighbor.sphere<-r.array[count.r]

               
				
				for(i in 1: number.of.cells){
				
				    
				    index.i<-which(data_cell_motility[, 8]==id[i])
				    
				    track.length[i]<-length(index.i)
				    
				    times.present.in.study<-c(times.present.in.study, length(index.i))
				    
				    x.cell<-c(x.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 1])
				    
				    y.cell<-c(y.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 2])
				    
				    z.cell<-c(z.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 3])
				    
				    time.cell<-c(time.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 7])
				    
				    data_i<-matrix(c(data_cell_motility[index.i, 7], data_cell_motility[index.i, 1], data_cell_motility[index.i, 2], data_cell_motility[index.i, 3]), nrow=length(index.i), ncol=4)
				    
				    				    
        
        #############################  Finding the corresponding clustering criteria
        
        
                    increment.code<-matrix(0, nrow=3, ncol=1)
        
                    increment.x<-increment.x.total.simulated[index.i[1:(l-1)]]
                    
                    increment.y<-increment.y.total.simulated[index.i[1:(l-1)]]
                    
                    increment.z<-increment.z.total.simulated[index.i[1:(l-1)]]
                    
                    
                    				    				    
				    increment.x<-increment.y<-increment.z<-NULL
				    
				    increment.code<-matrix(0, nrow=3, ncol=1)

                    for(j in 2: length(index.i)){
	        	
	        	             increment.x[j-1]<-data_cell_motility[index.i[j], 1]-data_cell_motility[index.i[j-1], 1]
	        	
	        	             increment.y[j-1]<-data_cell_motility[index.i[j], 2]-data_cell_motility[index.i[j-1], 2]
	        	
	        	             increment.z[j-1]<-data_cell_motility[index.i[j], 3]-data_cell_motility[index.i[j-1], 3]
	        	
	                }

        
        	                
	                	                                
	                B=10000; 
	                
	                n = length(increment.x)
	                
	                
	                normal.samples = matrix( rnorm(n*B, mean(increment.x), sd(increment.x)), B, n)
                   
                    normal.statistics = apply(normal.samples, 1, quantile)
                    
                    first.quantile<-normal.statistics[2, ]
                    
                    third.quantile<-normal.statistics[4, ]
                    
                    min.boot<-normal.statistics[1, ]
                    
                    max.boot<-normal.statistics[5, ]
                    
                    cm_normal<-(third.quantile-first.quantile)/(max.boot-min.boot)
                    
                    interval = quantile(cm_normal, c(0.025, 0.975))
                    

	                cm_x_increments<-(summary(increment.x)[5]-summary(increment.x)[2])/(summary(increment.x)[6]-summary(increment.x)[1])

                     |  cm_b > interval[2]	                
	                
 	                if(cm_x_increments < interval[1] ){
	                	
	                	increment.code[1]<-1
	                	
	                }
	                
	                
	                
	                
	                normal.samples = matrix( rnorm(n*B, mean(increment.y), sd(increment.y)), B, n)
                   
                    normal.statistics = apply(normal.samples, 1, quantile)
                    
                    first.quantile<-normal.statistics[2, ]
                    
                    third.quantile<-normal.statistics[4, ]
                    
                    min.boot<-normal.statistics[1, ]
                    
                    max.boot<-normal.statistics[5, ]
                    
                    cm_normal<-(third.quantile-first.quantile)/(max.boot-min.boot)
                    
                    interval = quantile(cm_normal, c(0.025, 0.975))	                
	                
	                cm_y_increments<-(summary(increment.y)[5]-summary(increment.y)[2])/(summary(increment.y)[6]-summary(increment.y)[1])
	                
	                
	                if(cm_y_increments < interval[1]  ){
	                	
	                	increment.code[2]<-1
	                	
	                }

	                	                
	                
	               	                                
	                B=10000; 
	                
	                n = length(increment.z)
	                
	                
	                normal.samples = matrix( rnorm(n*B, mean(increment.z), sd(increment.z)), B, n)
                   
                    normal.statistics = apply(normal.samples, 1, quantile)
                    
                    first.quantile<-normal.statistics[2, ]
                    
                    third.quantile<-normal.statistics[4, ]
                    
                    min.boot<-normal.statistics[1, ]
                    
                    max.boot<-normal.statistics[5, ]
                    
                    cm_normal<-(third.quantile-first.quantile)/(max.boot-min.boot)
	                
	                interval = quantile(cm_normal, c(0.025, 0.975))
	                
	                cm_z_increments<-(summary(increment.z)[5]-summary(increment.z)[2])/(summary(increment.z)[6]-summary(increment.z)[1])
	                
	                
	                if(cm_z_increments < interval[1]  ){
	                	
	                	increment.code[3]<-1
	                	
	                }

	                
	                
                    increment.trajectory<-matrix(c(increment.x, increment.y, increment.z), ncol=3, nrow=(track.length[i]-1))
       
                    if(sum(increment.code)>0) r.domain.ratio[i]<-1
                                      
                 	x<-tracks(data_i)
                 	
                 	u<-radius.phi.theta.psi(trajectory=x$'1', trajectory.num=i)
                 	
                 	
                 	
                   angular.change.complete.full.data[index.i[2:(track.length[i]-1)]]<-u$psi.component

                   radius.complete.full.data[index.i[2:(track.length[i]-1)]]<-(u$radius[2:(track.length[i]-1)])

        
                    
                    p.persistent.nonpersistent[i, ]<-c(sum(u$psi.component<=pi/2)/length(u$psi.component), (1-(sum(u$psi.component<=pi/2)/length(u$psi.component))))

        
                   
 
 }            
 
 
 
                     The other way around is to use the median of the new dataset in order to make clusters.                 
                    
                    n_c<-NULL
							
												
					cluster.1<-which(p.persistent.nonpersistent[, 1] >= median(p.persistent.nonpersistent[, 1]) & r.domain.ratio==1 &  track.length>=29)   ### persistent and long
					
					n_c[1]<-length(cluster.1)
					
					cluster.2<-which(p.persistent.nonpersistent[, 1] >= median(p.persistent.nonpersistent[, 1]) & r.domain.ratio==0 &  track.length>=29)   ### persistent and short 
					
					n_c[2]<-length(cluster.2)
					
					cluster.3<-which(p.persistent.nonpersistent[, 1] < median(p.persistent.nonpersistent[, 1]) & r.domain.ratio==1 &  track.length>=29)   ### non-persistent and long
					
					n_c[3]<-length(cluster.3)
					
					cluster.4<-which(p.persistent.nonpersistent[, 1] < median(p.persistent.nonpersistent[, 1]) & r.domain.ratio==0 &  track.length>=29)   ### non-persistent and short
					
					n_c[4]<-length(cluster.4)
					

                    print("n_c")
                    
                    print(n_c)
                    
                    
                    
                    
                    
###########################################################################################################################################################################################################################################  neighbors

		       		        
		      total.phi<-NULL
			    
			  total.psi<-NULL
			    
			  total.theta<-NULL
			      
			  total.radius<-NULL  
				
              increment.x.cell<-NULL
		       
		      increment.y.cell<-NULL
		        
		      increment.z.cell<-NULL 
		        
		    
		    
           				                              
                
                track.length<-NULL
                           
                              
                id_2<-seq(1000000000, (1000000000+number.of.cells_2-1), 1)

				
				for(i in 1: number.of.cells_2){
				
				    
				    index.i<-which(data_cell_motility[, 8]==id_2[i])

				    
				    track.length[i]<-length(index.i)
				    
				    				    				    
				    data_i<-matrix(c(data_cell_motility[index.i, 7], data_cell_motility[index.i, 1], data_cell_motility[index.i, 2], data_cell_motility[index.i, 3]), nrow=length(index.i), ncol=4)
				    				
				    x<-tracks(data_i)
				        				        
                                  u<-radius.phi.theta.psi(trajectory=x$'1', trajectory.num=i) 

                                  r.array<-c(20, 22, 25, 30, 35, 40, 45, 50, 60, 65)

                                  count.r<-2

                                  r.neighbor.sphere<-r.array[count.r]

                                  number.of.neighbors<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j.x<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j.y<-matrix(0, nrow=(track.length[i]), ncol=1)	
				    
				    min.distance.i.from.j.z<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				     Here we deal with the positions. Note that neighbors can be defined for each position.
				    
				  
                   	
                   	           for(k in 1: (track.length[i])){

				    		    distance.i.from.j<-NULL
				    		    
				    		    distance.i.from.j.x<-NULL
				    		    
				    		    distance.i.from.j.y<-NULL
				    		    
				    		    distance.i.from.j.z<-NULL
				    		    
				    		    l<-1
				    		    
				    		    t.i<-data_cell_motility[index.i[k], 7]
				    		    
				    		    center<-c(data_cell_motility[index.i[k], 1], data_cell_motility[index.i[k], 2], data_cell_motility[index.i[k], 3])
				    		    
						    	for( j in 1: number.of.cells_2){
						    		
						    		if(j!=i){
						    			
						    	        index.j<-which(data_cell_motility[, 8]==id_2[j])
						    	        
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
				    			
				    			distance.i.from.j.x<-distance.i.from.j[abs(distance.i.from.j.x)>0]
				    			
				    			min.distance.i.from.j.x[k]<-distance.i.from.j.x[which(abs(distance.i.from.j.x)==min(abs(distance.i.from.j.x)))[1]]
				    			
				    			distance.i.from.j.y<-distance.i.from.j[abs(distance.i.from.j.y)>0]
				    			
				    			min.distance.i.from.j.y[k]<-distance.i.from.j.y[which(abs(distance.i.from.j.y)==min(abs(distance.i.from.j.y)))[1]]
				    			
				    			distance.i.from.j.z<-distance.i.from.j[abs(distance.i.from.j.z)>0]
				    			
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
                   
                   
                
                
                data_cell_motility<-rbind(data_cell_motility_2, data_cell_motility_4)

				number.of.cells<-number.of.cells_2+number.of.cells_4



                id_4<-seq(1000000000, (1000000000+number.of.cells_4-1), 1)+number.of.cells_2

				
                track.length<-NULL
                
                                
				
				for(i in 1: number.of.cells_4){


                    i<-cluster[[1]][cluster.counter]
				
				    
				    index.i<-which(data_cell_motility[, 8]==id_4[i])

				    
				    track.length[i]<-length(index.i)
				    
				    				    				    
				    data_i<-matrix(c(data_cell_motility[index.i, 7], data_cell_motility[index.i, 1], data_cell_motility[index.i, 2], data_cell_motility[index.i, 3]), nrow=length(index.i), ncol=4)
				    				
				    x<-tracks(data_i)
				        				        
                                  u<-radius.phi.theta.psi(trajectory=x$'1', trajectory.num=i) 

                                  r.array<-c(20, 22, 25, 30, 35, 40, 45, 50, 60, 65)

                                  count.r<-2

                                  r.neighbor.sphere<-r.array[count.r]

                                  number.of.neighbors<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j.x<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j.y<-matrix(0, nrow=(track.length[i]), ncol=1)	
				    
				    min.distance.i.from.j.z<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				     Here we deal with the positions. Note that neighbors can be defined for each position.
				    
				  
                   	
                   	           for(k in 1: (track.length[i])){

				    		    distance.i.from.j<-NULL
				    		    
				    		    distance.i.from.j.x<-NULL
				    		    
				    		    distance.i.from.j.y<-NULL
				    		    
				    		    distance.i.from.j.z<-NULL
				    		    
				    		    l<-1
				    		    
				    		    t.i<-data_cell_motility[index.i[k], 7]
				    		    
				    		    center<-c(data_cell_motility[index.i[k], 1], data_cell_motility[index.i[k], 2], data_cell_motility[index.i[k], 3])
				    		    
						    	for( j in 1: number.of.cells_4){
						    		
						    		if(j!=i){
						    			
						    	        index.j<-which(data_cell_motility[, 8]==id_4[j])
						    	        
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
				    			
				    			distance.i.from.j.x<-distance.i.from.j[abs(distance.i.from.j.x)>0]
				    			
				    			min.distance.i.from.j.x[k]<-distance.i.from.j.x[which(abs(distance.i.from.j.x)==min(abs(distance.i.from.j.x)))[1]]
				    			
				    			distance.i.from.j.y<-distance.i.from.j[abs(distance.i.from.j.y)>0]
				    			
				    			min.distance.i.from.j.y[k]<-distance.i.from.j.y[which(abs(distance.i.from.j.y)==min(abs(distance.i.from.j.y)))[1]]
				    			
				    			distance.i.from.j.z<-distance.i.from.j[abs(distance.i.from.j.z)>0]
				    			
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


                    
                  number.of.neighbors.full.data <- number.of.neighbors.full.data  
                  
                  
                  angular.change.complete.full.data <- angular.change.complete.full.data
                  
                  
                  radius.complete.full.data <- radius.complete.full.data



                  data_cell_motility <- data_cell_motility


                  cluster<-list("1"=cluster.1, "2"=cluster.2, "3"=cluster.3, "4"=cluster.4)


#######################################################################################################################################################################################################################################

                				
				
				index.angular.change<-NULL
				
				for(i in 1: number.of.cells){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[i])
					
					
					 
					 index.angular.change<-c(index.angular.change, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				
				index.angular.change.cluster.1<-NULL
				
				for(i in 1:length(cluster.1){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.1[i]])
					
					
					 
					 index.angular.change.cluster.1<-c(index.angular.change.cluster.1, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				
				index.angular.change.cluster.2<-NULL
				
				for(i in 1:length(cluster.2)){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.2[i]])
					
					
					 
					 index.angular.change.cluster.2<-c(index.angular.change.cluster.2, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				
				index.angular.change.cluster.3<-NULL
				
				for(i in 1:length(cluster.3)){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.3[i]])
					
					
					 
					 index.angular.change.cluster.3<-c(index.angular.change.cluster.3, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				
				index.angular.change.cluster.4<-NULL
				
				for(i in 1:length(cluster.4)){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.4[i]])
					
					
					 
					 index.angular.change.cluster.4<-c(index.angular.change.cluster.4, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				index.angular.change.cluster<-c(index.angular.change.cluster.1, index.angular.change.cluster.2, index.angular.change.cluster.3, index.angular.change.cluster.4)


###########################################################################################################################################################################################################################




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


devience.infor.criterion<-function(parameter, parameter.bayes.f){
	
	t<-apply(parameter.fit.angle, 1, log.likelihood.f)
	
	Bayes<-log.likelihood.f(parameter.bayes.f)


     u<-2*(Bayes-mean(t))
    
    
    DIC.cluster<- (-2*(Bayes-u))



    print("DIC.cluster")


    print(DIC.cluster)

	
	
}


##########################################################################

# # # # DIC.cluster<-devience.infor.criterion(parameter.fit.angle, parameter.bayes)

# # print(DIC.cluster)

###########################################################################################################################################################################################



##########################################################################


WAIC.f<-function(parameter){
	
	t<-apply(parameter, 1, log.likelihood.f)
	
	likelihood<-exp(t)
	
	mean.likelihood<-mean(likelihood)	
	
	if(round(mean.likelihood, 4)==0) mean.likelihood<-10^-4

    u <-2*(log(mean.likelihood)-mean(t))
    
    WAIC.cluster<-(-2*(log(mean.likelihood)-u))


    print("WAIC.cluster")


    print(WAIC.cluster)

	
	
}


########################################################################

WAIC.cluster<-WAIC.f(parameter.fit.angle)

print(WAIC.cluster)



###########################################################################################################################################################################################



library(rstan)

remove(list=ls())



setwd("/global/scratch/torkashe/Clustering-Results/cluster-3-3-states-HMM/")

load(".RData")




cluster.c<-3



fit.ext<-extract(fit)


Rep.angle<-(length(fit.ext[[1]])/3)


counter.c<-as.matrix(rep(cluster.c, Rep.angle))


prior_pi<-as.matrix(fit.ext[[1]])


P_1<-fit.ext[[2]]


# P_2<-fit.ext[[2]][, 2, , ]


alpha_per<-as.matrix(fit.ext[[3]])


beta_per<-as.matrix(fit.ext[[4]])


alpha_non_per<-as.matrix(fit.ext[[5]])


beta_non_per<-as.matrix(fit.ext[[6]])


alpha_transition<-beta_transition <-as.matrix(fit.ext[[7]])

parameter.fit.angle<-cbind(prior_pi, P_1[, 1, 1], P_1[, 1, 2], P_1[, 1, 3], P_1[, 2, 1], P_1[, 2, 2], P_1[, 2, 3], P_1[, 3, 1], P_1[, 3, 2], P_1[, 3, 3], alpha_per, beta_per, alpha_non_per, beta_non_per, alpha_transition , beta_transition, counter.c)



###################################################


######################################################################################################################################################################################################################################################################################  Held-Out Data  #########################################################

   
                

                data_cell_motility_2<-read.csv("Position_2.csv", header=TRUE)

				dim(data_cell_motility_2)


				  In order to have an understanding of the number of T-cells observed in the study, we use the following comment


				number.of.cells_2<-length(unique(data_cell_motility_2[, 8]))



                
                data_cell_motility_4<-read.csv("Position_Longest.csv", header=TRUE)

				dim(data_cell_motility_4)

				data_cell_motility_4[, 8]<-data_cell_motility_4[, 8]+(number.of.cells_2)
				

				  In order to have an understanding of the number of T-cells observed in the study, we use the following comment

                
                number.of.cells_4<-length(unique(data_cell_motility_4[, 8]))
                


    
    
                
                data_cell_motility<-rbind(data_cell_motility_2, data_cell_motility_4)

				number.of.cells<-number.of.cells_2+number.of.cells_4

	              In order to have an understanding of the number of T-cells observed in the study, we use the following comment

             
				number.of.cells<-length(unique(data_cell_motility[, 8]))

				id<-seq(1000000000, (1000000000+number.of.cells-1), 1)
				
				track.length<-NULL
				
                for(i in 1: number.of.cells){
                	
                	track.length[i]<-length(which(data_cell_motility[, 8]==id[i]))
                	
                }

                number.of.positions.total<-sum(track.length)
                
                
                
##########################################
##########################################
				
				
				angular.change.complete.full.data<-matrix(0, nrow= number.of.positions.total, ncol=1)
				
				
				number.of.neighbors.full.data<-matrix(0, nrow= number.of.positions.total, ncol=1)

                
                radius.complete.full.data<-matrix(0, nrow= number.of.positions.total, ncol=1)


						
		tracks<-list() 
				    
				    
				    												
        This gives the trajectories in space. I suggest you to run this code. By looking at the trajectories, it make sense to cluster them in different groups.

       track.length<-NULL
              
       hierarchical.clustering.criteria<-matrix(0, nrow= number.of.cells, ncol=length(points.for.estimation.angle)+1)
      
       radius.trajectory<-NULL
       
       theta.trajectory<-NULL
       
       phi.trajectory<-NULL
              
       r.domain.ratio<-rep(0, number.of.cells)
       
       p.persistent.nonpersistent<-matrix(0, nrow= number.of.cells, ncol=2)
                     
        edited theta and phi
       
       theta.breaks<-seq(0, 3.5, 0.5)
       
       theta.breaks<-seq(0, 1.6, 0.4)
       
       phi.breaks<-seq(-3.5, 3.5, 0.5)
       
       x.cell<-y.cell<-z.cell<-time.cell<-NULL 
       
       increment.x<-increment.y<-increment.z<-NULL

       B.boot<-10000
                                          
       persistent.probability<-matrix(0, nrow= number.of.cells, ncol=length(r.array))
       
       conditional.persistent.probability<-matrix(0, nrow= number.of.cells, ncol=length(r.array))
                             
                      
#######################################################################################################
###################################################  Spherical representation and angular change 

               
    radius.phi.theta.psi<-function(trajectory, trajectory.num){
					
					length.trajectory<-length(trajectory[, 1])
					
					length.sub.trajectories<-NULL
					
					theta<-NULL
					
					phi<-NULL
					
					for(i in 1: (length.trajectory-1)){
						
					   	length.sub.trajectories[i]<-sqrt((trajectory[(i+1), 2]-trajectory[i, 2])^2+(trajectory[(i+1), 3]-trajectory[i, 3])^2+(trajectory[(i+1), 4]-trajectory[i, 4])^2)
					   	
					   	theta[i]<-acos((trajectory[(i+1), 4]-trajectory[i, 4])/length.sub.trajectories[i])
					   	
					   	  Considering the revised theta under my assumption
					   	
					    theta[i]<-acos(abs((trajectory[(i+1), 4]-trajectory[i, 4])/length.sub.trajectories[i]))
				
					   	
					   	if(theta[i]>(pi/2)) theta[i]<-pi-theta[i]
					   	
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
						
						print("u"); print(u)
						
						v<-c((trajectory[(i+2), 2]-trajectory[i+1, 2]), (trajectory[(i+2), 3]-trajectory[i+1, 3]), (trajectory[(i+2), 4]-trajectory[i+1, 4]))
						print("v"); print(v)
						
						u.cross.product.v<-c((u[2]*v[3]-u[3]*v[2]), (u[3]*v[1]-u[1]*v[3]), (u[1]*v[2]-u[2]*v[1]))
						
						u.dot.product.v<-u[1]*v[1]+u[2]*v[2]+u[3]*v[3]
						
						length.u.cross.product.v<-sqrt((u.cross.product.v[1])^2+(u.cross.product.v[2])^2+(u.cross.product.v[3])^2)
						
						length.u<-sqrt((u[1])^2+(u[2])^2+(u[3])^2)
						
						length.v<-sqrt((v[1])^2+(v[2])^2+(v[3])^2)
						
						angle.between.sub.trajectories[i]<-asin(length.u.cross.product.v/(length.u* length.v))
						
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


 
								
				 This gives the trajectories in space. I suggest you to run this code. By looking at the trajectories, it make sense to cluster them in different groups.
				
				track.length<-NULL
				       
		        x.cell<-NULL
		       
		        y.cell<-NULL
		        
		        z.cell<-NULL 
		       
		        time.cell<-NULL
		        
		        total.phi<-NULL
			    
			    total.psi<-NULL
			    
			    total.theta<-NULL
			      
			    total.radius<-NULL  
				
                increment.x.cell<-NULL
		       
		        increment.y.cell<-NULL
		        
		        increment.z.cell<-NULL 

                times.present.in.study<-NULL
                        
                    	
                count.r<-3
                   	
                r.neighbor.sphere<-r.array[count.r]

               
				
				for(i in 1: number.of.cells){
				
				    
				    index.i<-which(data_cell_motility[, 8]==id[i])
				    
				    track.length[i]<-length(index.i)
				    
				    times.present.in.study<-c(times.present.in.study, length(index.i))
				    
				    x.cell<-c(x.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 1])
				    
				    y.cell<-c(y.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 2])
				    
				    z.cell<-c(z.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 3])
				    
				    time.cell<-c(time.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 7])
				    
				    data_i<-matrix(c(data_cell_motility[index.i, 7], data_cell_motility[index.i, 1], data_cell_motility[index.i, 2], data_cell_motility[index.i, 3]), nrow=length(index.i), ncol=4)
				    
				    				    
        
        #############################  Finding the corresponding clustering criteria
        
        
                    increment.code<-matrix(0, nrow=3, ncol=1)
        
                    increment.x<-increment.x.total.simulated[index.i[1:(l-1)]]
                    
                    increment.y<-increment.y.total.simulated[index.i[1:(l-1)]]
                    
                    increment.z<-increment.z.total.simulated[index.i[1:(l-1)]]
                    
                    
                    				    				    
				    increment.x<-increment.y<-increment.z<-NULL
				    
				    increment.code<-matrix(0, nrow=3, ncol=1)

                    for(j in 2: length(index.i)){
	        	
	        	             increment.x[j-1]<-data_cell_motility[index.i[j], 1]-data_cell_motility[index.i[j-1], 1]
	        	
	        	             increment.y[j-1]<-data_cell_motility[index.i[j], 2]-data_cell_motility[index.i[j-1], 2]
	        	
	        	             increment.z[j-1]<-data_cell_motility[index.i[j], 3]-data_cell_motility[index.i[j-1], 3]
	        	
	                }

        
        	                
	                	                                
	                B=10000; 
	                
	                n = length(increment.x)
	                
	                
	                normal.samples = matrix( rnorm(n*B, mean(increment.x), sd(increment.x)), B, n)
                   
                    normal.statistics = apply(normal.samples, 1, quantile)
                    
                    first.quantile<-normal.statistics[2, ]
                    
                    third.quantile<-normal.statistics[4, ]
                    
                    min.boot<-normal.statistics[1, ]
                    
                    max.boot<-normal.statistics[5, ]
                    
                    cm_normal<-(third.quantile-first.quantile)/(max.boot-min.boot)
                    
                    interval = quantile(cm_normal, c(0.025, 0.975))
                    

	                cm_x_increments<-(summary(increment.x)[5]-summary(increment.x)[2])/(summary(increment.x)[6]-summary(increment.x)[1])

                     |  cm_b > interval[2]	                
	                
 	                if(cm_x_increments < interval[1] ){
	                	
	                	increment.code[1]<-1
	                	
	                }
	                
	                
	                
	                
	                normal.samples = matrix( rnorm(n*B, mean(increment.y), sd(increment.y)), B, n)
                   
                    normal.statistics = apply(normal.samples, 1, quantile)
                    
                    first.quantile<-normal.statistics[2, ]
                    
                    third.quantile<-normal.statistics[4, ]
                    
                    min.boot<-normal.statistics[1, ]
                    
                    max.boot<-normal.statistics[5, ]
                    
                    cm_normal<-(third.quantile-first.quantile)/(max.boot-min.boot)
                    
                    interval = quantile(cm_normal, c(0.025, 0.975))	                
	                
	                cm_y_increments<-(summary(increment.y)[5]-summary(increment.y)[2])/(summary(increment.y)[6]-summary(increment.y)[1])
	                
	                
	                if(cm_y_increments < interval[1]  ){
	                	
	                	increment.code[2]<-1
	                	
	                }

	                	                
	                
	               	                                
	                B=10000; 
	                
	                n = length(increment.z)
	                
	                
	                normal.samples = matrix( rnorm(n*B, mean(increment.z), sd(increment.z)), B, n)
                   
                    normal.statistics = apply(normal.samples, 1, quantile)
                    
                    first.quantile<-normal.statistics[2, ]
                    
                    third.quantile<-normal.statistics[4, ]
                    
                    min.boot<-normal.statistics[1, ]
                    
                    max.boot<-normal.statistics[5, ]
                    
                    cm_normal<-(third.quantile-first.quantile)/(max.boot-min.boot)
	                
	                interval = quantile(cm_normal, c(0.025, 0.975))
	                
	                cm_z_increments<-(summary(increment.z)[5]-summary(increment.z)[2])/(summary(increment.z)[6]-summary(increment.z)[1])
	                
	                
	                if(cm_z_increments < interval[1]  ){
	                	
	                	increment.code[3]<-1
	                	
	                }

	                
	                
                    increment.trajectory<-matrix(c(increment.x, increment.y, increment.z), ncol=3, nrow=(track.length[i]-1))
       
                    if(sum(increment.code)>0) r.domain.ratio[i]<-1
                                      
                 	x<-tracks(data_i)
                 	
                 	u<-radius.phi.theta.psi(trajectory=x$'1', trajectory.num=i)
                 	
                 	
                 	
                   angular.change.complete.full.data[index.i[2:(track.length[i]-1)]]<-u$psi.component

                   radius.complete.full.data[index.i[2:(track.length[i]-1)]]<-(u$radius[2:(track.length[i]-1)])

        
                    
                    p.persistent.nonpersistent[i, ]<-c(sum(u$psi.component<=pi/2)/length(u$psi.component), (1-(sum(u$psi.component<=pi/2)/length(u$psi.component))))

        
                   
 
 }            
 
 
 
                     The other way around is to use the median of the new dataset in order to make clusters.                 
                    
                    n_c<-NULL
							
												
					cluster.1<-which(p.persistent.nonpersistent[, 1] >= median(p.persistent.nonpersistent[, 1]) & r.domain.ratio==1 &  track.length>=29)   ### persistent and long
					
					n_c[1]<-length(cluster.1)
					
					cluster.2<-which(p.persistent.nonpersistent[, 1] >= median(p.persistent.nonpersistent[, 1]) & r.domain.ratio==0 &  track.length>=29)   ### persistent and short 
					
					n_c[2]<-length(cluster.2)
					
					cluster.3<-which(p.persistent.nonpersistent[, 1] < median(p.persistent.nonpersistent[, 1]) & r.domain.ratio==1 &  track.length>=29)   ### non-persistent and long
					
					n_c[3]<-length(cluster.3)
					
					cluster.4<-which(p.persistent.nonpersistent[, 1] < median(p.persistent.nonpersistent[, 1]) & r.domain.ratio==0 &  track.length>=29)   ### non-persistent and short
					
					n_c[4]<-length(cluster.4)
					

                    print("n_c")
                    
                    print(n_c)
                    
                    
                    
                    
                    
###########################################################################################################################################################################################################################################  neighbors

		       		        
		      total.phi<-NULL
			    
			  total.psi<-NULL
			    
			  total.theta<-NULL
			      
			  total.radius<-NULL  
				
              increment.x.cell<-NULL
		       
		      increment.y.cell<-NULL
		        
		      increment.z.cell<-NULL 
		        
		    
		    
           				                              
                
                track.length<-NULL
                           
                              
                id_2<-seq(1000000000, (1000000000+number.of.cells_2-1), 1)

				
				for(i in 1: number.of.cells_2){
				
				    
				    index.i<-which(data_cell_motility[, 8]==id_2[i])

				    
				    track.length[i]<-length(index.i)
				    
				    				    				    
				    data_i<-matrix(c(data_cell_motility[index.i, 7], data_cell_motility[index.i, 1], data_cell_motility[index.i, 2], data_cell_motility[index.i, 3]), nrow=length(index.i), ncol=4)
				    				
				    x<-tracks(data_i)
				        				        
                                  u<-radius.phi.theta.psi(trajectory=x$'1', trajectory.num=i) 

                                  r.array<-c(20, 22, 25, 30, 35, 40, 45, 50, 60, 65)

                                  count.r<-2

                                  r.neighbor.sphere<-r.array[count.r]

                                  number.of.neighbors<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j.x<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j.y<-matrix(0, nrow=(track.length[i]), ncol=1)	
				    
				    min.distance.i.from.j.z<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				     Here we deal with the positions. Note that neighbors can be defined for each position.
				    
				  
                   	
                   	           for(k in 1: (track.length[i])){

				    		    distance.i.from.j<-NULL
				    		    
				    		    distance.i.from.j.x<-NULL
				    		    
				    		    distance.i.from.j.y<-NULL
				    		    
				    		    distance.i.from.j.z<-NULL
				    		    
				    		    l<-1
				    		    
				    		    t.i<-data_cell_motility[index.i[k], 7]
				    		    
				    		    center<-c(data_cell_motility[index.i[k], 1], data_cell_motility[index.i[k], 2], data_cell_motility[index.i[k], 3])
				    		    
						    	for( j in 1: number.of.cells_2){
						    		
						    		if(j!=i){
						    			
						    	        index.j<-which(data_cell_motility[, 8]==id_2[j])
						    	        
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
				    			
				    			distance.i.from.j.x<-distance.i.from.j[abs(distance.i.from.j.x)>0]
				    			
				    			min.distance.i.from.j.x[k]<-distance.i.from.j.x[which(abs(distance.i.from.j.x)==min(abs(distance.i.from.j.x)))[1]]
				    			
				    			distance.i.from.j.y<-distance.i.from.j[abs(distance.i.from.j.y)>0]
				    			
				    			min.distance.i.from.j.y[k]<-distance.i.from.j.y[which(abs(distance.i.from.j.y)==min(abs(distance.i.from.j.y)))[1]]
				    			
				    			distance.i.from.j.z<-distance.i.from.j[abs(distance.i.from.j.z)>0]
				    			
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
                   
                   
                
                
                data_cell_motility<-rbind(data_cell_motility_2, data_cell_motility_4)

				number.of.cells<-number.of.cells_2+number.of.cells_4



                id_4<-seq(1000000000, (1000000000+number.of.cells_4-1), 1)+number.of.cells_2

				
                track.length<-NULL
                
                                
				
				for(i in 1: number.of.cells_4){


                    i<-cluster[[1]][cluster.counter]
				
				    
				    index.i<-which(data_cell_motility[, 8]==id_4[i])

				    
				    track.length[i]<-length(index.i)
				    
				    				    				    
				    data_i<-matrix(c(data_cell_motility[index.i, 7], data_cell_motility[index.i, 1], data_cell_motility[index.i, 2], data_cell_motility[index.i, 3]), nrow=length(index.i), ncol=4)
				    				
				    x<-tracks(data_i)
				        				        
                                  u<-radius.phi.theta.psi(trajectory=x$'1', trajectory.num=i) 

                                  r.array<-c(20, 22, 25, 30, 35, 40, 45, 50, 60, 65)

                                  count.r<-2

                                  r.neighbor.sphere<-r.array[count.r]

                                  number.of.neighbors<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j.x<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j.y<-matrix(0, nrow=(track.length[i]), ncol=1)	
				    
				    min.distance.i.from.j.z<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				     Here we deal with the positions. Note that neighbors can be defined for each position.
				    
				  
                   	
                   	           for(k in 1: (track.length[i])){

				    		    distance.i.from.j<-NULL
				    		    
				    		    distance.i.from.j.x<-NULL
				    		    
				    		    distance.i.from.j.y<-NULL
				    		    
				    		    distance.i.from.j.z<-NULL
				    		    
				    		    l<-1
				    		    
				    		    t.i<-data_cell_motility[index.i[k], 7]
				    		    
				    		    center<-c(data_cell_motility[index.i[k], 1], data_cell_motility[index.i[k], 2], data_cell_motility[index.i[k], 3])
				    		    
						    	for( j in 1: number.of.cells_4){
						    		
						    		if(j!=i){
						    			
						    	        index.j<-which(data_cell_motility[, 8]==id_4[j])
						    	        
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
				    			
				    			distance.i.from.j.x<-distance.i.from.j[abs(distance.i.from.j.x)>0]
				    			
				    			min.distance.i.from.j.x[k]<-distance.i.from.j.x[which(abs(distance.i.from.j.x)==min(abs(distance.i.from.j.x)))[1]]
				    			
				    			distance.i.from.j.y<-distance.i.from.j[abs(distance.i.from.j.y)>0]
				    			
				    			min.distance.i.from.j.y[k]<-distance.i.from.j.y[which(abs(distance.i.from.j.y)==min(abs(distance.i.from.j.y)))[1]]
				    			
				    			distance.i.from.j.z<-distance.i.from.j[abs(distance.i.from.j.z)>0]
				    			
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


                    
                  number.of.neighbors.full.data <- number.of.neighbors.full.data  
                  
                  
                  angular.change.complete.full.data <- angular.change.complete.full.data
                  
                  
                  radius.complete.full.data <- radius.complete.full.data



                  data_cell_motility <- data_cell_motility
                  
                  
                  cluster<-list("1"=cluster.1, "2"=cluster.2, "3"=cluster.3, "4"=cluster.4)



#######################################################################################################################################################################################################################################

                				
				
				index.angular.change<-NULL
				
				for(i in 1: number.of.cells){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[i])
					
					
					 
					 index.angular.change<-c(index.angular.change, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				
				index.angular.change.cluster.1<-NULL
				
				for(i in 1:length(cluster.1){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.1[i]])
					
					
					 
					 index.angular.change.cluster.1<-c(index.angular.change.cluster.1, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				
				index.angular.change.cluster.2<-NULL
				
				for(i in 1:length(cluster.2)){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.2[i]])
					
					
					 
					 index.angular.change.cluster.2<-c(index.angular.change.cluster.2, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				
				index.angular.change.cluster.3<-NULL
				
				for(i in 1:length(cluster.3)){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.3[i]])
					
					
					 
					 index.angular.change.cluster.3<-c(index.angular.change.cluster.3, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				
				index.angular.change.cluster.4<-NULL
				
				for(i in 1:length(cluster.4)){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.4[i]])
					
					
					 
					 index.angular.change.cluster.4<-c(index.angular.change.cluster.4, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				index.angular.change.cluster<-c(index.angular.change.cluster.1, index.angular.change.cluster.2, index.angular.change.cluster.3, index.angular.change.cluster.4)


###########################################################################################################################################################################################################################



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


devience.infor.criterion<-function(parameter, parameter.bayes.f){
	
	t<-apply(parameter.fit.angle, 1, log.likelihood.f)
	
	Bayes<-log.likelihood.f(parameter.bayes.f)


     u<-2*(Bayes-mean(t))
    
    
    DIC.cluster<- (-2*(Bayes-u))



    print("DIC.cluster")


    print(DIC.cluster)

	
	
}


##########################################################################

# # # # DIC.cluster<-devience.infor.criterion(parameter.fit.angle, parameter.bayes)

# # print(DIC.cluster)

###########################################################################################################################################################################################



##########################################################################


WAIC.f<-function(parameter){
	
	t<-apply(parameter, 1, log.likelihood.f)
	
	likelihood<-exp(t)
	
	mean.likelihood<-mean(likelihood)	
	
	if(round(mean.likelihood, 4)==0) mean.likelihood<-10^-4

    u <-2*(log(mean.likelihood)-mean(t))
    
    WAIC.cluster<-(-2*(log(mean.likelihood)-u))


    print("WAIC.cluster")


    print(WAIC.cluster)

	
	
}


########################################################################

WAIC.cluster<-WAIC.f(parameter.fit.angle)

print(WAIC.cluster)



###########################################################################################################################################################################################


library(rstan)

remove(list=ls())



setwd("/global/scratch/torkashe/Clustering-Results/cluster-4-3-states-HMM/")

load(".RData")




cluster.c<-4


fit.ext<-extract(fit)


Rep.angle<-(length(fit.ext[[1]])/3)


counter.c<-as.matrix(rep(cluster.c, Rep.angle))


prior_pi<-as.matrix(fit.ext[[1]])


P_1<-fit.ext[[2]]


# P_2<-fit.ext[[2]][, 2, , ]


alpha_per<-as.matrix(fit.ext[[3]])


beta_per<-as.matrix(fit.ext[[4]])


alpha_non_per<-as.matrix(fit.ext[[5]])


beta_non_per<-as.matrix(fit.ext[[6]])


alpha_transition<-beta_transition <-as.matrix(fit.ext[[7]])

parameter.fit.angle<-cbind(prior_pi, P_1[, 1, 1], P_1[, 1, 2], P_1[, 1, 3], P_1[, 2, 1], P_1[, 2, 2], P_1[, 2, 3], P_1[, 3, 1], P_1[, 3, 2], P_1[, 3, 3], alpha_per, beta_per, alpha_non_per, beta_non_per, alpha_transition , beta_transition, counter.c)



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


devience.infor.criterion<-function(parameter, parameter.bayes.f){
	
	t<-apply(parameter.fit.angle, 1, log.likelihood.f)
	
	Bayes<-log.likelihood.f(parameter.bayes.f)


     u<-2*(Bayes-mean(t))
    
    
    DIC.cluster<- (-2*(Bayes-u))



    print("DIC.cluster")


    print(DIC.cluster)

	
	
}


##########################################################################

# # # # DIC.cluster<-devience.infor.criterion(parameter.fit.angle, parameter.bayes)

# # print(DIC.cluster)

###########################################################################################################################################################################################



##########################################################################


WAIC.f<-function(parameter){
	
	t<-apply(parameter, 1, log.likelihood.f)
	
	likelihood<-exp(t)
	
	mean.likelihood<-mean(likelihood)	
	
	if(round(mean.likelihood, 4)==0) mean.likelihood<-10^-4

    u <-2*(log(mean.likelihood)-mean(t))
    
    WAIC.cluster<-(-2*(log(mean.likelihood)-u))


    print("WAIC.cluster")


    print(WAIC.cluster)

	
	
}


########################################################################

WAIC.cluster<-WAIC.f(parameter.fit.angle)

print(WAIC.cluster)



###########################################################################################################################################################################################



library(rstan)

remove(list=ls())


setwd("/global/scratch/torkashe/Clustering-Results/2-clusters-3-states-HMM-1/")

load(".RData")




cluster.c<-1


fit.ext<-extract(fit)


Rep.angle<-(length(fit.ext[[1]])/3)


counter.c<-as.matrix(rep(cluster.c, Rep.angle))


prior_pi<-as.matrix(fit.ext[[1]])


P_1<-fit.ext[[2]]


# P_2<-fit.ext[[2]][, 2, , ]


alpha_per<-as.matrix(fit.ext[[3]])


beta_per<-as.matrix(fit.ext[[4]])


alpha_non_per<-as.matrix(fit.ext[[5]])


beta_non_per<-as.matrix(fit.ext[[6]])


 alpha_transition<-beta_transition <-as.matrix(fit.ext[[7]])

parameter.fit.angle<-cbind(prior_pi, P_1[, 1, 1], P_1[, 1, 2], P_1[, 1, 3], P_1[, 2, 1], P_1[, 2, 2], P_1[, 2, 3], P_1[, 3, 1], P_1[, 3, 2], P_1[, 3, 3], alpha_per, beta_per, alpha_non_per, beta_non_per, alpha_transition , beta_transition, counter.c)



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


devience.infor.criterion<-function(parameter, parameter.bayes.f){
	
	t<-apply(parameter.fit.angle, 1, log.likelihood.f)
	
	Bayes<-log.likelihood.f(parameter.bayes.f)


     u<-2*(Bayes-mean(t))
    
    
    DIC.cluster<- (-2*(Bayes-u))



    print("DIC.cluster")


    print(DIC.cluster)

	
	
}


##########################################################################

# # # # # # DIC.cluster<-devience.infor.criterion(parameter.fit.angle, parameter.bayes)

# # # print(DIC.cluster)

###########################################################################################################################################################################################



##########################################################################


WAIC.f<-function(parameter){
	
	t<-apply(parameter, 1, log.likelihood.f)
	
	likelihood<-exp(t)
	
	mean.likelihood<-mean(likelihood)	
	
	if(round(mean.likelihood, 4)==0) mean.likelihood<-10^-4

    u <-2*(log(mean.likelihood)-mean(t))
    
    WAIC.cluster<-(-2*(log(mean.likelihood)-u))


    print("WAIC.cluster")


    print(WAIC.cluster)

	
	
}


########################################################################

WAIC.cluster<-WAIC.f(parameter.fit.angle)

print(WAIC.cluster)



###########################################################################################################################################################################################

library(rstan)

remove(list=ls())


setwd("/global/scratch/torkashe/Clustering-Results/2-clusters-3-states-HMM-2/")

load(".RData")




cluster.c<-2



fit.ext<-extract(fit)


Rep.angle<-(length(fit.ext[[1]])/3)


counter.c<-as.matrix(rep(cluster.c, Rep.angle))


prior_pi<-as.matrix(fit.ext[[1]])


P_1<-fit.ext[[2]]


# P_2<-fit.ext[[2]][, 2, , ]


alpha_per<-as.matrix(fit.ext[[3]])


beta_per<-as.matrix(fit.ext[[4]])


alpha_non_per<-as.matrix(fit.ext[[5]])


beta_non_per<-as.matrix(fit.ext[[6]])


alpha_transition<-beta_transition <-as.matrix(fit.ext[[7]])

parameter.fit.angle<-cbind(prior_pi, P_1[, 1, 1], P_1[, 1, 2], P_1[, 1, 3], P_1[, 2, 1], P_1[, 2, 2], P_1[, 2, 3], P_1[, 3, 1], P_1[, 3, 2], P_1[, 3, 3], alpha_per, beta_per, alpha_non_per, beta_non_per, alpha_transition , beta_transition, counter.c)


######################################################################################################################################################################################################################################################################################  Held-Out Data  #########################################################

   
                

                data_cell_motility_2<-read.csv("Position_2.csv", header=TRUE)

				dim(data_cell_motility_2)


				  In order to have an understanding of the number of T-cells observed in the study, we use the following comment


				number.of.cells_2<-length(unique(data_cell_motility_2[, 8]))



                
                data_cell_motility_4<-read.csv("Position_Longest.csv", header=TRUE)

				dim(data_cell_motility_4)

				data_cell_motility_4[, 8]<-data_cell_motility_4[, 8]+(number.of.cells_2)
				

				  In order to have an understanding of the number of T-cells observed in the study, we use the following comment

                
                number.of.cells_4<-length(unique(data_cell_motility_4[, 8]))
                


    
    
                
                data_cell_motility<-rbind(data_cell_motility_2, data_cell_motility_4)

				number.of.cells<-number.of.cells_2+number.of.cells_4

	              In order to have an understanding of the number of T-cells observed in the study, we use the following comment

             
				number.of.cells<-length(unique(data_cell_motility[, 8]))

				id<-seq(1000000000, (1000000000+number.of.cells-1), 1)
				
				track.length<-NULL
				
                for(i in 1: number.of.cells){
                	
                	track.length[i]<-length(which(data_cell_motility[, 8]==id[i]))
                	
                }

                number.of.positions.total<-sum(track.length)
                
                
                
##########################################
##########################################
				
				
				angular.change.complete.full.data<-matrix(0, nrow= number.of.positions.total, ncol=1)
				
				
				number.of.neighbors.full.data<-matrix(0, nrow= number.of.positions.total, ncol=1)

                
                radius.complete.full.data<-matrix(0, nrow= number.of.positions.total, ncol=1)


						
		tracks<-list() 
				    
				    
				    												
        This gives the trajectories in space. I suggest you to run this code. By looking at the trajectories, it make sense to cluster them in different groups.

       track.length<-NULL
              
       hierarchical.clustering.criteria<-matrix(0, nrow= number.of.cells, ncol=length(points.for.estimation.angle)+1)
      
       radius.trajectory<-NULL
       
       theta.trajectory<-NULL
       
       phi.trajectory<-NULL
              
       r.domain.ratio<-rep(0, number.of.cells)
       
       p.persistent.nonpersistent<-matrix(0, nrow= number.of.cells, ncol=2)
                     
        edited theta and phi
       
       theta.breaks<-seq(0, 3.5, 0.5)
       
       theta.breaks<-seq(0, 1.6, 0.4)
       
       phi.breaks<-seq(-3.5, 3.5, 0.5)
       
       x.cell<-y.cell<-z.cell<-time.cell<-NULL 
       
       increment.x<-increment.y<-increment.z<-NULL

       B.boot<-10000
                                          
       persistent.probability<-matrix(0, nrow= number.of.cells, ncol=length(r.array))
       
       conditional.persistent.probability<-matrix(0, nrow= number.of.cells, ncol=length(r.array))
                             
                      
#######################################################################################################
###################################################  Spherical representation and angular change 

               
    radius.phi.theta.psi<-function(trajectory, trajectory.num){
					
					length.trajectory<-length(trajectory[, 1])
					
					length.sub.trajectories<-NULL
					
					theta<-NULL
					
					phi<-NULL
					
					for(i in 1: (length.trajectory-1)){
						
					   	length.sub.trajectories[i]<-sqrt((trajectory[(i+1), 2]-trajectory[i, 2])^2+(trajectory[(i+1), 3]-trajectory[i, 3])^2+(trajectory[(i+1), 4]-trajectory[i, 4])^2)
					   	
					   	theta[i]<-acos((trajectory[(i+1), 4]-trajectory[i, 4])/length.sub.trajectories[i])
					   	
					   	  Considering the revised theta under my assumption
					   	
					    theta[i]<-acos(abs((trajectory[(i+1), 4]-trajectory[i, 4])/length.sub.trajectories[i]))
				
					   	
					   	if(theta[i]>(pi/2)) theta[i]<-pi-theta[i]
					   	
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
						
						print("u"); print(u)
						
						v<-c((trajectory[(i+2), 2]-trajectory[i+1, 2]), (trajectory[(i+2), 3]-trajectory[i+1, 3]), (trajectory[(i+2), 4]-trajectory[i+1, 4]))
						print("v"); print(v)
						
						u.cross.product.v<-c((u[2]*v[3]-u[3]*v[2]), (u[3]*v[1]-u[1]*v[3]), (u[1]*v[2]-u[2]*v[1]))
						
						u.dot.product.v<-u[1]*v[1]+u[2]*v[2]+u[3]*v[3]
						
						length.u.cross.product.v<-sqrt((u.cross.product.v[1])^2+(u.cross.product.v[2])^2+(u.cross.product.v[3])^2)
						
						length.u<-sqrt((u[1])^2+(u[2])^2+(u[3])^2)
						
						length.v<-sqrt((v[1])^2+(v[2])^2+(v[3])^2)
						
						angle.between.sub.trajectories[i]<-asin(length.u.cross.product.v/(length.u* length.v))
						
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


 
								
				 This gives the trajectories in space. I suggest you to run this code. By looking at the trajectories, it make sense to cluster them in different groups.
				
				track.length<-NULL
				       
		        x.cell<-NULL
		       
		        y.cell<-NULL
		        
		        z.cell<-NULL 
		       
		        time.cell<-NULL
		        
		        total.phi<-NULL
			    
			    total.psi<-NULL
			    
			    total.theta<-NULL
			      
			    total.radius<-NULL  
				
                increment.x.cell<-NULL
		       
		        increment.y.cell<-NULL
		        
		        increment.z.cell<-NULL 

                times.present.in.study<-NULL
                        
                    	
                count.r<-3
                   	
                r.neighbor.sphere<-r.array[count.r]

               
				
				for(i in 1: number.of.cells){
				
				    
				    index.i<-which(data_cell_motility[, 8]==id[i])
				    
				    track.length[i]<-length(index.i)
				    
				    times.present.in.study<-c(times.present.in.study, length(index.i))
				    
				    x.cell<-c(x.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 1])
				    
				    y.cell<-c(y.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 2])
				    
				    z.cell<-c(z.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 3])
				    
				    time.cell<-c(time.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 7])
				    
				    data_i<-matrix(c(data_cell_motility[index.i, 7], data_cell_motility[index.i, 1], data_cell_motility[index.i, 2], data_cell_motility[index.i, 3]), nrow=length(index.i), ncol=4)
				    
				    				    
        
        #############################  Finding the corresponding clustering criteria
        
        
                    increment.code<-matrix(0, nrow=3, ncol=1)
        
                    increment.x<-increment.x.total.simulated[index.i[1:(l-1)]]
                    
                    increment.y<-increment.y.total.simulated[index.i[1:(l-1)]]
                    
                    increment.z<-increment.z.total.simulated[index.i[1:(l-1)]]
                    
                    
                    				    				    
				    increment.x<-increment.y<-increment.z<-NULL
				    
				    increment.code<-matrix(0, nrow=3, ncol=1)

                    for(j in 2: length(index.i)){
	        	
	        	             increment.x[j-1]<-data_cell_motility[index.i[j], 1]-data_cell_motility[index.i[j-1], 1]
	        	
	        	             increment.y[j-1]<-data_cell_motility[index.i[j], 2]-data_cell_motility[index.i[j-1], 2]
	        	
	        	             increment.z[j-1]<-data_cell_motility[index.i[j], 3]-data_cell_motility[index.i[j-1], 3]
	        	
	                }

        
        	                
	                	                                
	                B=10000; 
	                
	                n = length(increment.x)
	                
	                
	                normal.samples = matrix( rnorm(n*B, mean(increment.x), sd(increment.x)), B, n)
                   
                    normal.statistics = apply(normal.samples, 1, quantile)
                    
                    first.quantile<-normal.statistics[2, ]
                    
                    third.quantile<-normal.statistics[4, ]
                    
                    min.boot<-normal.statistics[1, ]
                    
                    max.boot<-normal.statistics[5, ]
                    
                    cm_normal<-(third.quantile-first.quantile)/(max.boot-min.boot)
                    
                    interval = quantile(cm_normal, c(0.025, 0.975))
                    

	                cm_x_increments<-(summary(increment.x)[5]-summary(increment.x)[2])/(summary(increment.x)[6]-summary(increment.x)[1])

                     |  cm_b > interval[2]	                
	                
 	                if(cm_x_increments < interval[1] ){
	                	
	                	increment.code[1]<-1
	                	
	                }
	                
	                
	                
	                
	                normal.samples = matrix( rnorm(n*B, mean(increment.y), sd(increment.y)), B, n)
                   
                    normal.statistics = apply(normal.samples, 1, quantile)
                    
                    first.quantile<-normal.statistics[2, ]
                    
                    third.quantile<-normal.statistics[4, ]
                    
                    min.boot<-normal.statistics[1, ]
                    
                    max.boot<-normal.statistics[5, ]
                    
                    cm_normal<-(third.quantile-first.quantile)/(max.boot-min.boot)
                    
                    interval = quantile(cm_normal, c(0.025, 0.975))	                
	                
	                cm_y_increments<-(summary(increment.y)[5]-summary(increment.y)[2])/(summary(increment.y)[6]-summary(increment.y)[1])
	                
	                
	                if(cm_y_increments < interval[1]  ){
	                	
	                	increment.code[2]<-1
	                	
	                }

	                	                
	                
	               	                                
	                B=10000; 
	                
	                n = length(increment.z)
	                
	                
	                normal.samples = matrix( rnorm(n*B, mean(increment.z), sd(increment.z)), B, n)
                   
                    normal.statistics = apply(normal.samples, 1, quantile)
                    
                    first.quantile<-normal.statistics[2, ]
                    
                    third.quantile<-normal.statistics[4, ]
                    
                    min.boot<-normal.statistics[1, ]
                    
                    max.boot<-normal.statistics[5, ]
                    
                    cm_normal<-(third.quantile-first.quantile)/(max.boot-min.boot)
	                
	                interval = quantile(cm_normal, c(0.025, 0.975))
	                
	                cm_z_increments<-(summary(increment.z)[5]-summary(increment.z)[2])/(summary(increment.z)[6]-summary(increment.z)[1])
	                
	                
	                if(cm_z_increments < interval[1]  ){
	                	
	                	increment.code[3]<-1
	                	
	                }

	                
	                
                    increment.trajectory<-matrix(c(increment.x, increment.y, increment.z), ncol=3, nrow=(track.length[i]-1))
       
                    if(sum(increment.code)>0) r.domain.ratio[i]<-1
                                      
                 	x<-tracks(data_i)
                 	
                 	u<-radius.phi.theta.psi(trajectory=x$'1', trajectory.num=i)
                 	
                 	
                 	
                   angular.change.complete.full.data[index.i[2:(track.length[i]-1)]]<-u$psi.component

                   radius.complete.full.data[index.i[2:(track.length[i]-1)]]<-(u$radius[2:(track.length[i]-1)])

        
                    
                    p.persistent.nonpersistent[i, ]<-c(sum(u$psi.component<=pi/2)/length(u$psi.component), (1-(sum(u$psi.component<=pi/2)/length(u$psi.component))))

        
                   
 
 }            
 
 
 
                     The other way around is to use the median of the new dataset in order to make clusters.                 
                    
                    n_c<-NULL
							
												
					cluster.1<-which(p.persistent.nonpersistent[, 1] >= median(p.persistent.nonpersistent[, 1]) & r.domain.ratio==1 &  track.length>=29)   ### persistent and long
					
					n_c[1]<-length(cluster.1)
					
					cluster.2<-which(p.persistent.nonpersistent[, 1] >= median(p.persistent.nonpersistent[, 1]) & r.domain.ratio==0 &  track.length>=29)   ### persistent and short 
					
					n_c[2]<-length(cluster.2)
					
					cluster.3<-which(p.persistent.nonpersistent[, 1] < median(p.persistent.nonpersistent[, 1]) & r.domain.ratio==1 &  track.length>=29)   ### non-persistent and long
					
					n_c[3]<-length(cluster.3)
					
					cluster.4<-which(p.persistent.nonpersistent[, 1] < median(p.persistent.nonpersistent[, 1]) & r.domain.ratio==0 &  track.length>=29)   ### non-persistent and short
					
					n_c[4]<-length(cluster.4)
					

                    print("n_c")
                    
                    print(n_c)
                    
                    
                    
                    
                    
###########################################################################################################################################################################################################################################  neighbors

		       		        
		      total.phi<-NULL
			    
			  total.psi<-NULL
			    
			  total.theta<-NULL
			      
			  total.radius<-NULL  
				
              increment.x.cell<-NULL
		       
		      increment.y.cell<-NULL
		        
		      increment.z.cell<-NULL 
		        
		    
		    
           				                              
                
                track.length<-NULL
                           
                              
                id_2<-seq(1000000000, (1000000000+number.of.cells_2-1), 1)

				
				for(i in 1: number.of.cells_2){
				
				    
				    index.i<-which(data_cell_motility[, 8]==id_2[i])

				    
				    track.length[i]<-length(index.i)
				    
				    				    				    
				    data_i<-matrix(c(data_cell_motility[index.i, 7], data_cell_motility[index.i, 1], data_cell_motility[index.i, 2], data_cell_motility[index.i, 3]), nrow=length(index.i), ncol=4)
				    				
				    x<-tracks(data_i)
				        				        
                                  u<-radius.phi.theta.psi(trajectory=x$'1', trajectory.num=i) 

                                  r.array<-c(20, 22, 25, 30, 35, 40, 45, 50, 60, 65)

                                  count.r<-2

                                  r.neighbor.sphere<-r.array[count.r]

                                  number.of.neighbors<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j.x<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j.y<-matrix(0, nrow=(track.length[i]), ncol=1)	
				    
				    min.distance.i.from.j.z<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				     Here we deal with the positions. Note that neighbors can be defined for each position.
				    
				  
                   	
                   	           for(k in 1: (track.length[i])){

				    		    distance.i.from.j<-NULL
				    		    
				    		    distance.i.from.j.x<-NULL
				    		    
				    		    distance.i.from.j.y<-NULL
				    		    
				    		    distance.i.from.j.z<-NULL
				    		    
				    		    l<-1
				    		    
				    		    t.i<-data_cell_motility[index.i[k], 7]
				    		    
				    		    center<-c(data_cell_motility[index.i[k], 1], data_cell_motility[index.i[k], 2], data_cell_motility[index.i[k], 3])
				    		    
						    	for( j in 1: number.of.cells_2){
						    		
						    		if(j!=i){
						    			
						    	        index.j<-which(data_cell_motility[, 8]==id_2[j])
						    	        
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
				    			
				    			distance.i.from.j.x<-distance.i.from.j[abs(distance.i.from.j.x)>0]
				    			
				    			min.distance.i.from.j.x[k]<-distance.i.from.j.x[which(abs(distance.i.from.j.x)==min(abs(distance.i.from.j.x)))[1]]
				    			
				    			distance.i.from.j.y<-distance.i.from.j[abs(distance.i.from.j.y)>0]
				    			
				    			min.distance.i.from.j.y[k]<-distance.i.from.j.y[which(abs(distance.i.from.j.y)==min(abs(distance.i.from.j.y)))[1]]
				    			
				    			distance.i.from.j.z<-distance.i.from.j[abs(distance.i.from.j.z)>0]
				    			
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
                   
                   
                
                
                data_cell_motility<-rbind(data_cell_motility_2, data_cell_motility_4)

				number.of.cells<-number.of.cells_2+number.of.cells_4



                id_4<-seq(1000000000, (1000000000+number.of.cells_4-1), 1)+number.of.cells_2

				
                track.length<-NULL
                
                                
				
				for(i in 1: number.of.cells_4){


                    i<-cluster[[1]][cluster.counter]
				
				    
				    index.i<-which(data_cell_motility[, 8]==id_4[i])

				    
				    track.length[i]<-length(index.i)
				    
				    				    				    
				    data_i<-matrix(c(data_cell_motility[index.i, 7], data_cell_motility[index.i, 1], data_cell_motility[index.i, 2], data_cell_motility[index.i, 3]), nrow=length(index.i), ncol=4)
				    				
				    x<-tracks(data_i)
				        				        
                                  u<-radius.phi.theta.psi(trajectory=x$'1', trajectory.num=i) 

                                  r.array<-c(20, 22, 25, 30, 35, 40, 45, 50, 60, 65)

                                  count.r<-2

                                  r.neighbor.sphere<-r.array[count.r]

                                  number.of.neighbors<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j.x<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    min.distance.i.from.j.y<-matrix(0, nrow=(track.length[i]), ncol=1)	
				    
				    min.distance.i.from.j.z<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				     Here we deal with the positions. Note that neighbors can be defined for each position.
				    
				  
                   	
                   	           for(k in 1: (track.length[i])){

				    		    distance.i.from.j<-NULL
				    		    
				    		    distance.i.from.j.x<-NULL
				    		    
				    		    distance.i.from.j.y<-NULL
				    		    
				    		    distance.i.from.j.z<-NULL
				    		    
				    		    l<-1
				    		    
				    		    t.i<-data_cell_motility[index.i[k], 7]
				    		    
				    		    center<-c(data_cell_motility[index.i[k], 1], data_cell_motility[index.i[k], 2], data_cell_motility[index.i[k], 3])
				    		    
						    	for( j in 1: number.of.cells_4){
						    		
						    		if(j!=i){
						    			
						    	        index.j<-which(data_cell_motility[, 8]==id_4[j])
						    	        
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
				    			
				    			distance.i.from.j.x<-distance.i.from.j[abs(distance.i.from.j.x)>0]
				    			
				    			min.distance.i.from.j.x[k]<-distance.i.from.j.x[which(abs(distance.i.from.j.x)==min(abs(distance.i.from.j.x)))[1]]
				    			
				    			distance.i.from.j.y<-distance.i.from.j[abs(distance.i.from.j.y)>0]
				    			
				    			min.distance.i.from.j.y[k]<-distance.i.from.j.y[which(abs(distance.i.from.j.y)==min(abs(distance.i.from.j.y)))[1]]
				    			
				    			distance.i.from.j.z<-distance.i.from.j[abs(distance.i.from.j.z)>0]
				    			
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


                    
                  number.of.neighbors.full.data <- number.of.neighbors.full.data  
                  
                  
                  angular.change.complete.full.data <- angular.change.complete.full.data
                  
                  
                  radius.complete.full.data <- radius.complete.full.data



                  data_cell_motility <- data_cell_motility
                  
                  
                  cluster<-list("1"=cluster.1, "2"=cluster.2, "3"=cluster.3, "4"=cluster.4)



#######################################################################################################################################################################################################################################

                				
				
				index.angular.change<-NULL
				
				for(i in 1: number.of.cells){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[i])
					
					
					 
					 index.angular.change<-c(index.angular.change, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				
				index.angular.change.cluster.1<-NULL
				
				for(i in 1:length(cluster.1){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.1[i]])
					
					
					 
					 index.angular.change.cluster.1<-c(index.angular.change.cluster.1, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				
				index.angular.change.cluster.2<-NULL
				
				for(i in 1:length(cluster.2)){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.2[i]])
					
					
					 
					 index.angular.change.cluster.2<-c(index.angular.change.cluster.2, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				
				index.angular.change.cluster.3<-NULL
				
				for(i in 1:length(cluster.3)){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.3[i]])
					
					
					 
					 index.angular.change.cluster.3<-c(index.angular.change.cluster.3, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				
				index.angular.change.cluster.4<-NULL
				
				for(i in 1:length(cluster.4)){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.4[i]])
					
					
					 
					 index.angular.change.cluster.4<-c(index.angular.change.cluster.4, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				index.angular.change.cluster<-c(index.angular.change.cluster.1, index.angular.change.cluster.2, index.angular.change.cluster.3, index.angular.change.cluster.4)


###########################################################################################################################################################################################################################



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


devience.infor.criterion<-function(parameter, parameter.bayes.f){
	
	t<-apply(parameter.fit.angle, 1, log.likelihood.f)
	
	Bayes<-log.likelihood.f(parameter.bayes.f)


     u<-2*(Bayes-mean(t))
    
    
    DIC.cluster<- (-2*(Bayes-u))



    print("DIC.cluster")


    print(DIC.cluster)

	
	
}


##########################################################################

# # # # DIC.cluster<-devience.infor.criterion(parameter.fit.angle, parameter.bayes)

# # print(DIC.cluster)

###########################################################################################################################################################################################



##########################################################################


WAIC.f<-function(parameter){
	
	t<-apply(parameter, 1, log.likelihood.f)
	
	likelihood<-exp(t)
	
	mean.likelihood<-mean(likelihood)	
	
	if(round(mean.likelihood, 4)==0) mean.likelihood<-10^-4

    u <-2*(log(mean.likelihood)-mean(t))
    
    WAIC.cluster<-(-2*(log(mean.likelihood)-u))


    print("WAIC.cluster")


    print(WAIC.cluster)

	
	
}


########################################################################

WAIC.cluster<-WAIC.f(parameter.fit.angle)

print(WAIC.cluster)



###########################################################################################################################################################################################
