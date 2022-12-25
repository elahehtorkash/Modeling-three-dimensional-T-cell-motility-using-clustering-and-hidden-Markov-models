#############################################    Diagram of Angular Change ########################################
################   2 states NHMM

remove(list=ls())




###  setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Motility-Data")

setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/cluster-evaluation-2-states-HHMM/")

##########   FEB 5, 2016

date()

set.seed(5)

####  Note that as we are interested to see how movements on different directions interact with each other, X, Y, and Z are modeled in the same model.
 
############   Cell Motility Project-Sep 12-2016   ###############
####  WARNING: for the purpose of dependeny, experiments should be kept separated.


require("mgcv")


library(boot)


require("MotilityLab")


library(moments)


library(cluster)


library(scatterplot3d)


library(rgl)


library(coda)


library(rjags)


library(R2WinBUGS)


library(dclone)


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


library(matlib)

library(sde)

library("rstan")

library("MotilityLab")

library(ggplot2)



# # # # # # 

# # # k_0<-4


# # # mean.log.mean.predictive.cluster<-matrix(0, nrow=k_0, ncol=5)


		
# # # #################################################################################################
# # # ################################ the radius we search for the number of cells in neighborhood 				

                
		
# # # ##################################################################################################################################################################################################################################################################################################  Held-Out Data  #########################################################

   
                

                # # # data_cell_motility_2<-read.csv("Position_2.csv", header=TRUE)

				# # # dim(data_cell_motility_2)


				# # # ##   In order to have an understanding of the number of T-cells observed in the study, we use the following comment


				# # # number.of.cells_2<-length(unique(data_cell_motility_2[, 8]))



                
                # # # data_cell_motility_4<-read.csv("Position_Longest.csv", header=TRUE)

				# # # dim(data_cell_motility_4)

				# # # data_cell_motility_4[, 8]<-data_cell_motility_4[, 8]+(number.of.cells_2)
				

				# # # #####   In order to have an understanding of the number of T-cells observed in the study, we use the following comment

                
                # # # number.of.cells_4<-length(unique(data_cell_motility_4[, 8]))
                


    
    
                
                # # # data_cell_motility<-rbind(data_cell_motility_2, data_cell_motility_4)

				# # # number.of.cells<-number.of.cells_2+number.of.cells_4

	            # # # #########   In order to have an understanding of the number of T-cells observed in the study, we use the following comment

             
				# # # number.of.cells<-length(unique(data_cell_motility[, 8]))

				# # # id<-seq(1000000000, (1000000000+number.of.cells-1), 1)
				
				# # # track.length<-NULL
				
                # # # for(i in 1: number.of.cells){
                	
                	# # # track.length[i]<-length(which(data_cell_motility[, 8]==id[i]))
                	
                # # # }

                # # # number.of.positions.total<-sum(track.length)
                
                
                
# # # ######################################################
# # # ######################################################
				
				
				# # # angular.change.complete.full.data<-matrix(0, nrow= number.of.positions.total, ncol=1)
				
				
				# # # number.of.neighbors.full.data<-matrix(0, nrow= number.of.positions.total, ncol=1)

                
                # # # radius.complete.full.data<-matrix(0, nrow= number.of.positions.total, ncol=1)


						
		# # # tracks<-list() 
				    
				    
				    												
       # # # ###  This gives the trajectories in space. I suggest you to run this code. By looking at the trajectories, it make sense to cluster them in different groups.

       # # # track.length<-NULL
              
       # # # # hierarchical.clustering.criteria<-matrix(0, nrow= number.of.cells, ncol=length(points.for.estimation.angle)+1)
      
       # # # radius.trajectory<-NULL
       
       # # # theta.trajectory<-NULL
       
       # # # phi.trajectory<-NULL
              
       # # # r.domain.ratio<-rep(0, number.of.cells)
       
       # # # p.persistent.nonpersistent<-matrix(0, nrow= number.of.cells, ncol=2)
                     
       # # # ####  edited theta and phi
       
       # # # theta.breaks<-seq(0, 3.5, 0.5)
       
       # # # # theta.breaks<-seq(0, 1.6, 0.4)
       
       # # # phi.breaks<-seq(-3.5, 3.5, 0.5)
       
       # # # x.cell<-y.cell<-z.cell<-time.cell<-NULL 
       
       # # # increment.x<-increment.y<-increment.z<-NULL

       # # # B.boot<-10000
                                          
       # # # persistent.probability<-matrix(0, nrow= number.of.cells, ncol=length(r.array))
       
       # # # conditional.persistent.probability<-matrix(0, nrow= number.of.cells, ncol=length(r.array))
                             
                      
# # # ###################################################################################################################
# # # ###############################################################  Spherical representation and angular change 

               
    # # # radius.phi.theta.psi<-function(trajectory, trajectory.num){
					
					# # # length.trajectory<-length(trajectory[, 1])
					
					# # # length.sub.trajectories<-NULL
					
					# # # theta<-NULL
					
					# # # phi<-NULL
					
					# # # for(i in 1: (length.trajectory-1)){
						
					   	# # # length.sub.trajectories[i]<-sqrt((trajectory[(i+1), 2]-trajectory[i, 2])^2+(trajectory[(i+1), 3]-trajectory[i, 3])^2+(trajectory[(i+1), 4]-trajectory[i, 4])^2)
					   	
					   	# # # # theta[i]<-acos((trajectory[(i+1), 4]-trajectory[i, 4])/length.sub.trajectories[i])
					   	
					   	# # # #####   Considering the revised theta under my assumption
					   	
					    # # # theta[i]<-acos(abs((trajectory[(i+1), 4]-trajectory[i, 4])/length.sub.trajectories[i]))
				
					   	
					   	# # # # if(theta[i]>(pi/2)) theta[i]<-pi-theta[i]
					   	
					   	# # # phi[i]<-atan((trajectory[(i+1), 3]-trajectory[i, 3])/(trajectory[(i+1), 2]-trajectory[i, 2]))
					   	
					   	# # # if((trajectory[(i+1), 3]-trajectory[i, 3])<0 & (trajectory[(i+1), 2]-trajectory[i, 2])<0) phi[i]<-pi+phi[i]
					   	
					   	# # # if((trajectory[(i+1), 3]-trajectory[i, 3])>0 & (trajectory[(i+1), 2]-trajectory[i, 2])<0) phi[i]<-pi+phi[i]
					   		   	
				        # # # if(((trajectory[(i+1), 3]-trajectory[i, 3])==0) & ((trajectory[(i+1), 2]-trajectory[i, 2])==0)){
				        	
				        	# # # if(trajectory[i, 3]>0 & trajectory[i, 2]>0) phi[i]<-pi/4
				        	
				        	# # # if(trajectory[i, 3]>0 & trajectory[i, 2]<0) phi[i]<-3*pi/4
				        	
				        	# # # if(trajectory[i, 3]<0 & trajectory[i, 2]<0) phi[i]<--3*pi/4
				        	
				        	# # # if(trajectory[i, 3]<0 & trajectory[i, 2]>0) phi[i]<--pi/4
				        	# # # }
					   			
						# # # }
						
						# # # angle.between.sub.trajectories<-NULL
					
					# # # for(i in 1: (length.trajectory-2)){
						
						# # # u<-c((trajectory[(i+1), 2]-trajectory[i, 2]), (trajectory[(i+1), 3]-trajectory[i, 3]), (trajectory[(i+1), 4]-trajectory[i, 4]))
						
						# # # # print("u"); print(u)
						
						# # # v<-c((trajectory[(i+2), 2]-trajectory[i+1, 2]), (trajectory[(i+2), 3]-trajectory[i+1, 3]), (trajectory[(i+2), 4]-trajectory[i+1, 4]))
						# # # # print("v"); print(v)
						
						# # # u.cross.product.v<-c((u[2]*v[3]-u[3]*v[2]), (u[3]*v[1]-u[1]*v[3]), (u[1]*v[2]-u[2]*v[1]))
						
						# # # u.dot.product.v<-u[1]*v[1]+u[2]*v[2]+u[3]*v[3]
						
						# # # length.u.cross.product.v<-sqrt((u.cross.product.v[1])^2+(u.cross.product.v[2])^2+(u.cross.product.v[3])^2)
						
						# # # length.u<-sqrt((u[1])^2+(u[2])^2+(u[3])^2)
						
						# # # length.v<-sqrt((v[1])^2+(v[2])^2+(v[3])^2)
						
						# # # # angle.between.sub.trajectories[i]<-asin(length.u.cross.product.v/(length.u* length.v))
						
						# # # angle.between.sub.trajectories[i]<-acos(min(c(1, max(c(-1, u.dot.product.v/(length.u* length.v))))))	
												
					# # # }
				
						
					# # # increments.radius.theta.phi.psi<-list(radius= length.sub.trajectories, theta.component= theta, phi.component= phi, psi.component= angle.between.sub.trajectories)
					    
					    
					    # # # if(sum(is.na(length.sub.trajectories))>0) {print("length.trajectory"); print(length.sub.trajectories)}
					    	    
					    
					    # # # if(sum(is.na(theta)>0)){ print("length(theta)"); print(theta)}
					    
					    
					    # # # if(sum(is.na(phi)>0)) {print("length(phi.component)"); print(phi)}
					    
					    
					    # # # if(sum(is.na(angle.between.sub.trajectories)>0)) {
					    	
					    	# # # print("length(psi.component)"); 
					    	
					    	# # # print(angle.between.sub.trajectories);
					    	
					    	# # # print("trajectory.num");
					    	
					    	# # # print(trajectory.num)
					    	
					    	# # # }
				
				
						
						# # # return(increments.radius.theta.phi.psi)
						
	
# # # }


 
								
				# # # ###  This gives the trajectories in space. I suggest you to run this code. By looking at the trajectories, it make sense to cluster them in different groups.
				
				# # # track.length<-NULL
				       
		        # # # x.cell<-NULL
		       
		        # # # y.cell<-NULL
		        
		        # # # z.cell<-NULL 
		       
		        # # # time.cell<-NULL
		        
		        # # # total.phi<-NULL
			    
			    # # # total.psi<-NULL
			    
			    # # # total.theta<-NULL
			      
			    # # # total.radius<-NULL  
				
                # # # increment.x.cell<-NULL
		       
		        # # # increment.y.cell<-NULL
		        
		        # # # increment.z.cell<-NULL 

                # # # times.present.in.study<-NULL
                    
                    
                # # # r.array<-c(20, 22, 25, 30, 35, 40, 45, 50, 60, 65)
    
                    	
                # # # count.r<-2
                   	
                # # # r.neighbor.sphere<-r.array[count.r]

               
				
				# # # for(i in 1: number.of.cells){
				
				    
				    # # # index.i<-which(data_cell_motility[, 8]==id[i])
				    
				    # # # track.length[i]<-length(index.i)
				    
				    # # # times.present.in.study<-c(times.present.in.study, length(index.i))
				    
				    # # # x.cell<-c(x.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 1])
				    
				    # # # y.cell<-c(y.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 2])
				    
				    # # # z.cell<-c(z.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 3])
				    
				    # # # time.cell<-c(time.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 7])
				    
				    # # # data_i<-matrix(c(data_cell_motility[index.i, 7], data_cell_motility[index.i, 1], data_cell_motility[index.i, 2], data_cell_motility[index.i, 3]), nrow=length(index.i), ncol=4)
				    
				    				    
        
        # # # #########################################  Finding the corresponding clustering criteria
        
        
                    # # # increment.code<-matrix(0, nrow=3, ncol=1)
        
# # # # # # #                     increment.x<-increment.x.total.simulated[index.i[1:(l-1)]]
                    
                    # # # # # increment.y<-increment.y.total.simulated[index.i[1:(l-1)]]
                    
                    # # # # # increment.z<-increment.z.total.simulated[index.i[1:(l-1)]]
                    
                    
                    				    				    
				    # # # increment.x<-increment.y<-increment.z<-NULL
				    
				    # # # increment.code<-matrix(0, nrow=3, ncol=1)

                    # # # for(j in 2: length(index.i)){
	        	
	        	             # # # increment.x[j-1]<-data_cell_motility[index.i[j], 1]-data_cell_motility[index.i[j-1], 1]
	        	
	        	             # # # increment.y[j-1]<-data_cell_motility[index.i[j], 2]-data_cell_motility[index.i[j-1], 2]
	        	
	        	             # # # increment.z[j-1]<-data_cell_motility[index.i[j], 3]-data_cell_motility[index.i[j-1], 3]
	        	
	                # # # }

        
        	                
	                	                                
	                # # # B=10000; 
	                
	                # # # n = length(increment.x)
	                
	                
	                # # # normal.samples = matrix( rnorm(n*B, mean(increment.x), sd(increment.x)), B, n)
                   
                    # # # normal.statistics = apply(normal.samples, 1, quantile)
                    
                    # # # first.quantile<-normal.statistics[2, ]
                    
                    # # # third.quantile<-normal.statistics[4, ]
                    
                    # # # min.boot<-normal.statistics[1, ]
                    
                    # # # max.boot<-normal.statistics[5, ]
                    
                    # # # cm_normal<-(third.quantile-first.quantile)/(max.boot-min.boot)
                    
                    # # # interval = quantile(cm_normal, c(0.025, 0.975))
                    

	                # # # cm_x_increments<-(summary(increment.x)[5]-summary(increment.x)[2])/(summary(increment.x)[6]-summary(increment.x)[1])

                    # # # ###  |  cm_b > interval[2]	                
	                
 	                # # # if(cm_x_increments < interval[1] ){
	                	
	                	# # # increment.code[1]<-1
	                	
	                # # # }
	                
	                
	                
	                
	                # # # normal.samples = matrix( rnorm(n*B, mean(increment.y), sd(increment.y)), B, n)
                   
                    # # # normal.statistics = apply(normal.samples, 1, quantile)
                    
                    # # # first.quantile<-normal.statistics[2, ]
                    
                    # # # third.quantile<-normal.statistics[4, ]
                    
                    # # # min.boot<-normal.statistics[1, ]
                    
                    # # # max.boot<-normal.statistics[5, ]
                    
                    # # # cm_normal<-(third.quantile-first.quantile)/(max.boot-min.boot)
                    
                    # # # interval = quantile(cm_normal, c(0.025, 0.975))	                
	                
	                # # # cm_y_increments<-(summary(increment.y)[5]-summary(increment.y)[2])/(summary(increment.y)[6]-summary(increment.y)[1])
	                
	                
	                # # # if(cm_y_increments < interval[1]  ){
	                	
	                	# # # increment.code[2]<-1
	                	
	                # # # }

	                	                
	                
	               	                                
	                # # # B=10000; 
	                
	                # # # n = length(increment.z)
	                
	                
	                # # # normal.samples = matrix( rnorm(n*B, mean(increment.z), sd(increment.z)), B, n)
                   
                    # # # normal.statistics = apply(normal.samples, 1, quantile)
                    
                    # # # first.quantile<-normal.statistics[2, ]
                    
                    # # # third.quantile<-normal.statistics[4, ]
                    
                    # # # min.boot<-normal.statistics[1, ]
                    
                    # # # max.boot<-normal.statistics[5, ]
                    
                    # # # cm_normal<-(third.quantile-first.quantile)/(max.boot-min.boot)
	                
	                # # # interval = quantile(cm_normal, c(0.025, 0.975))
	                
	                # # # cm_z_increments<-(summary(increment.z)[5]-summary(increment.z)[2])/(summary(increment.z)[6]-summary(increment.z)[1])
	                
	                
	                # # # if(cm_z_increments < interval[1]  ){
	                	
	                	# # # increment.code[3]<-1
	                	
	                # # # }

	                
	                
                    # # # increment.trajectory<-matrix(c(increment.x, increment.y, increment.z), ncol=3, nrow=(track.length[i]-1))
       
                    # # # if(sum(increment.code)>0) r.domain.ratio[i]<-1
                                      
                 	# # # x<-tracks(data_i)
                 	
                 	# # # u<-radius.phi.theta.psi(trajectory=x$'1', trajectory.num=i)
                 	
                 	
                 	
                   # # # angular.change.complete.full.data[index.i[2:(track.length[i]-1)]]<-u$psi.component

                   # # # radius.complete.full.data[index.i[2:(track.length[i]-1)]]<-(u$radius[2:(track.length[i]-1)])

        
                    
                    # # # p.persistent.nonpersistent[i, ]<-c(sum(u$psi.component<=pi/2)/length(u$psi.component), (1-(sum(u$psi.component<=pi/2)/length(u$psi.component))))

        
                   
 
 # # # }            
 
 
 
                    # # # ####  The other way around is to use the median of the new dataset in order to make clusters.                 
                    
                    # # # n_c<-NULL
							
												
					# # # cluster.1.held.out<-which(p.persistent.nonpersistent[, 1] >= median(p.persistent.nonpersistent[, 1]) & r.domain.ratio==1 &  track.length>=29)   ### persistent and long
					
					# # # n_c[1]<-length(cluster.1.held.out)
					
					# # # cluster.2.held.out<-which(p.persistent.nonpersistent[, 1] >= median(p.persistent.nonpersistent[, 1]) & r.domain.ratio==0 &  track.length>=29)   ### persistent and short 
					
					# # # n_c[2]<-length(cluster.2.held.out)
					
					# # # cluster.3.held.out<-which(p.persistent.nonpersistent[, 1] < median(p.persistent.nonpersistent[, 1]) & r.domain.ratio==1 &  track.length>=29)   ### non-persistent and long
					
					# # # n_c[3]<-length(cluster.3.held.out)
					
					# # # cluster.4.held.out<-which(p.persistent.nonpersistent[, 1] < median(p.persistent.nonpersistent[, 1]) & r.domain.ratio==0 &  track.length>=29)   ### non-persistent and short
					
					# # # n_c[4]<-length(cluster.4.held.out)
					

                    # # # print("n_c")
                    
                    # # # print(n_c)
                    
                    
                    
                    
                    
# # # #######################################################################################################################################################################################################################################################  neighbors

		       		        
		      # # # total.phi<-NULL
			    
			  # # # total.psi<-NULL
			    
			  # # # total.theta<-NULL
			      
			  # # # total.radius<-NULL  
				
              # # # increment.x.cell<-NULL
		       
		      # # # increment.y.cell<-NULL
		        
		      # # # increment.z.cell<-NULL 
		        
		    
		    
           				                              
                
                # # # track.length<-NULL
                           
                              
                # # # id_2<-seq(1000000000, (1000000000+number.of.cells_2-1), 1)

				
				# # # for(i in 1: number.of.cells_2){
				
				    
				    # # # index.i<-which(data_cell_motility[, 8]==id_2[i])

				    
				    # # # track.length[i]<-length(index.i)
				    
				    				    				    
				    # # # data_i<-matrix(c(data_cell_motility[index.i, 7], data_cell_motility[index.i, 1], data_cell_motility[index.i, 2], data_cell_motility[index.i, 3]), nrow=length(index.i), ncol=4)
				    				
				    # # # x<-tracks(data_i)
				        				        
                                  # # # u<-radius.phi.theta.psi(trajectory=x$'1', trajectory.num=i) 

                                  # # # r.array<-c(20, 22, 25, 30, 35, 40, 45, 50, 60, 65)


                                  # # # count.r<-2

                                  # # # r.neighbor.sphere<-r.array[count.r]

                                  # # # number.of.neighbors<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    # # # min.distance.i.from.j<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    # # # min.distance.i.from.j.x<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    # # # min.distance.i.from.j.y<-matrix(0, nrow=(track.length[i]), ncol=1)	
				    
				    # # # min.distance.i.from.j.z<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    # # # ###  Here we deal with the positions. Note that neighbors can be defined for each position.
				    
				  
                   	
                   	           # # # for(k in 1: (track.length[i])){

				    		    # # # distance.i.from.j<-NULL
				    		    
				    		    # # # distance.i.from.j.x<-NULL
				    		    
				    		    # # # distance.i.from.j.y<-NULL
				    		    
				    		    # # # distance.i.from.j.z<-NULL
				    		    
				    		    # # # l<-1
				    		    
				    		    # # # t.i<-data_cell_motility[index.i[k], 7]
				    		    
				    		    # # # center<-c(data_cell_motility[index.i[k], 1], data_cell_motility[index.i[k], 2], data_cell_motility[index.i[k], 3])
				    		    
						    	# # # for( j in 1: number.of.cells_2){
						    		
						    		# # # if(j!=i){
						    			
						    	        # # # index.j<-which(data_cell_motility[, 8]==id_2[j])
						    	        
						    	        # # # t.j<-data_cell_motility[index.j, 7]
						    	        
						    			# # # t.inter<-intersect(t.i, t.j)
						    			
						    									    	
						    									    			
						    			# # # if(length(t.inter)>0){
						    				
						    				# # # t.initial<-t.inter
						    				
								    		# # # t.initial.index<-index.j[which(data_cell_motility[index.j, 7]==t.initial)]
								    		
								    		# # # distance.i.from.j[l]<-sqrt(sum((center-data_cell_motility[t.initial.index, 1:3])^2))
								    		
								    		# # # distance.i.from.j.x[l]<-center[1]-data_cell_motility[t.initial.index, 1]
								    		
								    		# # # distance.i.from.j.y[l]<-center[2]-data_cell_motility[t.initial.index, 2]
								    		
								    		# # # distance.i.from.j.z[l]<-center[3]-data_cell_motility[t.initial.index, 3]	
								    										    	
								    		# # # t.count<-0
								    		
					    			    	# # # if(distance.i.from.j[l]<= r.neighbor.sphere){
					    			    			
					    			    		     # # # number.of.neighbors[k]<-number.of.neighbors[k]+1
					    			    		
					    			    	# # # }
					    			    		
					    			    		# # # l<-l+1
					    			    		
					    			    		   	# # # }
					    			    		   	
					    			    # # # if(l==1){
					    			    	
					    			    	# # # distance.i.from.j[l]<-0   ####  This zero means i and j are not in the same time span!
					    			    	
					    			    	# # # distance.i.from.j.x[l]<-0
					    			    	
					    			    	# # # distance.i.from.j.y[l]<-0
					    			    	
					    			    	# # # distance.i.from.j.z[l]<-0
					    			    	
					    			    # # # }		 
					    			      	
					                 # # # }
					
				    			# # # }
				    			
				    			# # # distance.i.from.j<-distance.i.from.j[abs(distance.i.from.j)>0]
				    			
				    			# # # min.distance.i.from.j[k]<-min(distance.i.from.j)
				    			
				    			# # # distance.i.from.j.x<-distance.i.from.j.x[abs(distance.i.from.j.x)>0]
				    			
				    			# # # min.distance.i.from.j.x[k]<-distance.i.from.j.x[which(abs(distance.i.from.j.x)==min(abs(distance.i.from.j.x)))[1]]
				    			
				    			# # # distance.i.from.j.y<-distance.i.from.j.y[abs(distance.i.from.j.y)>0]
				    			
				    			# # # min.distance.i.from.j.y[k]<-distance.i.from.j.y[which(abs(distance.i.from.j.y)==min(abs(distance.i.from.j.y)))[1]]
				    			
				    			# # # distance.i.from.j.z<-distance.i.from.j.z[abs(distance.i.from.j.z)>0]
				    			
				    			# # # min.distance.i.from.j.z[k]<-distance.i.from.j.z[which(abs(distance.i.from.j.z)==min(abs(distance.i.from.j.z)))[1]]
				    			
# # # }



                       
				    
				    # # # increment.x<-increment.y<-increment.z<-NULL

                                  # # # for(j in 2: length(index.i)){
	        	
	        	                # # # increment.x[j-1]<-data_cell_motility[index.i[j], 1]-data_cell_motility[index.i[j-1], 1]
	        	
	        	                # # # increment.y[j-1]<-data_cell_motility[index.i[j], 2]-data_cell_motility[index.i[j-1], 2]
	        	
	        	                # # # increment.z[j-1]<-data_cell_motility[index.i[j], 3]-data_cell_motility[index.i[j-1], 3]
	        	
	                           # # # }

                   # # # increment.trajectory<-matrix(c(increment.x, increment.y, increment.z), ncol=3, nrow=(track.length[i]-1))
                   
                   
                   # # # number.of.neighbors.full.data[index.i]<-number.of.neighbors

                   # # # angular.change.complete.full.data[index.i[2:(track.length[i]-1)]]<-u$psi.component

                   # # # radius.complete.full.data[index.i[2:(track.length[i]-1)]]<-(u$radius[2:(track.length[i]-1)])

               
                   # # # }
                   
                   
                
                
                # # # data_cell_motility<-rbind(data_cell_motility_2, data_cell_motility_4)

				# # # number.of.cells<-number.of.cells_2+number.of.cells_4



                # # # id_4<-seq(1000000000, (1000000000+number.of.cells_4-1), 1)+number.of.cells_2

				
                # # # track.length<-NULL
                
                                
				
				# # # for(i in 1: number.of.cells_4){


                    # # # # i<-cluster[[1]][cluster.counter]
				
				    
				    # # # index.i<-which(data_cell_motility[, 8]==id_4[i])

				    
				    # # # track.length[i]<-length(index.i)
				    
				    				    				    
				    # # # data_i<-matrix(c(data_cell_motility[index.i, 7], data_cell_motility[index.i, 1], data_cell_motility[index.i, 2], data_cell_motility[index.i, 3]), nrow=length(index.i), ncol=4)
				    				
				    # # # x<-tracks(data_i)
				        				        
                                  # # # u<-radius.phi.theta.psi(trajectory=x$'1', trajectory.num=i) 

                                  # # # r.array<-c(20, 22, 25, 30, 35, 40, 45, 50, 60, 65)

                                  # # # count.r<-2

                                  # # # r.neighbor.sphere<-r.array[count.r]

                                  # # # number.of.neighbors<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    # # # min.distance.i.from.j<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    # # # min.distance.i.from.j.x<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    # # # min.distance.i.from.j.y<-matrix(0, nrow=(track.length[i]), ncol=1)	
				    
				    # # # min.distance.i.from.j.z<-matrix(0, nrow=(track.length[i]), ncol=1)
				    
				    # # # ###  Here we deal with the positions. Note that neighbors can be defined for each position.
				    
				  
                   	
                   	           # # # for(k in 1: (track.length[i])){

				    		    # # # distance.i.from.j<-NULL
				    		    
				    		    # # # distance.i.from.j.x<-NULL
				    		    
				    		    # # # distance.i.from.j.y<-NULL
				    		    
				    		    # # # distance.i.from.j.z<-NULL
				    		    
				    		    # # # l<-1
				    		    
				    		    # # # t.i<-data_cell_motility[index.i[k], 7]
				    		    
				    		    # # # center<-c(data_cell_motility[index.i[k], 1], data_cell_motility[index.i[k], 2], data_cell_motility[index.i[k], 3])
				    		    
						    	# # # for( j in 1: number.of.cells_4){
						    		
						    		# # # if(j!=i){
						    			
						    	        # # # index.j<-which(data_cell_motility[, 8]==id_4[j])
						    	        
						    	        # # # t.j<-data_cell_motility[index.j, 7]
						    	        
						    			# # # t.inter<-intersect(t.i, t.j)
						    			
						    									    	
						    									    			
						    			# # # if(length(t.inter)>0){
						    				
						    				# # # t.initial<-t.inter
						    				
								    		# # # t.initial.index<-index.j[which(data_cell_motility[index.j, 7]==t.initial)]
								    		
								    		# # # distance.i.from.j[l]<-sqrt(sum((center-data_cell_motility[t.initial.index, 1:3])^2))
								    		
								    		# # # distance.i.from.j.x[l]<-center[1]-data_cell_motility[t.initial.index, 1]
								    		
								    		# # # distance.i.from.j.y[l]<-center[2]-data_cell_motility[t.initial.index, 2]
								    		
								    		# # # distance.i.from.j.z[l]<-center[3]-data_cell_motility[t.initial.index, 3]	
								    										    	
								    		# # # t.count<-0
								    		
					    			    	# # # if(distance.i.from.j[l]<= r.neighbor.sphere){
					    			    			
					    			    		     # # # number.of.neighbors[k]<-number.of.neighbors[k]+1
					    			    		
					    			    	# # # }
					    			    		
					    			    		# # # l<-l+1
					    			    		
					    			    		   	# # # }
					    			    		   	
					    			    # # # if(l==1){
					    			    	
					    			    	# # # distance.i.from.j[l]<-0   ####  This zero means i and j are not in the same time span!
					    			    	
					    			    	# # # distance.i.from.j.x[l]<-0
					    			    	
					    			    	# # # distance.i.from.j.y[l]<-0
					    			    	
					    			    	# # # distance.i.from.j.z[l]<-0
					    			    	
					    			    # # # }		 
					    			      	
					                 # # # }
					
				    			# # # }
				    			
				    			# # # distance.i.from.j<-distance.i.from.j[abs(distance.i.from.j)>0]
				    			
				    			# # # min.distance.i.from.j[k]<-min(distance.i.from.j)
				    			
				    			# # # distance.i.from.j.x<-distance.i.from.j.x[abs(distance.i.from.j.x)>0]
				    			
				    			# # # min.distance.i.from.j.x[k]<-distance.i.from.j.x[which(abs(distance.i.from.j.x)==min(abs(distance.i.from.j.x)))[1]]
				    			
				    			# # # distance.i.from.j.y<-distance.i.from.j.y[abs(distance.i.from.j.y)>0]
				    			
				    			# # # min.distance.i.from.j.y[k]<-distance.i.from.j.y[which(abs(distance.i.from.j.y)==min(abs(distance.i.from.j.y)))[1]]
				    			
				    			# # # distance.i.from.j.z<-distance.i.from.j.z[abs(distance.i.from.j.z)>0]
				    			
				    			# # # min.distance.i.from.j.z[k]<-distance.i.from.j.z[which(abs(distance.i.from.j.z)==min(abs(distance.i.from.j.z)))[1]]
				    			
# # # }


                       
				    
				    # # # increment.x<-increment.y<-increment.z<-NULL

                                  # # # for(j in 2: length(index.i)){
	        	
	        	                # # # increment.x[j-1]<-data_cell_motility[index.i[j], 1]-data_cell_motility[index.i[j-1], 1]
	        	
	        	                # # # increment.y[j-1]<-data_cell_motility[index.i[j], 2]-data_cell_motility[index.i[j-1], 2]
	        	
	        	                # # # increment.z[j-1]<-data_cell_motility[index.i[j], 3]-data_cell_motility[index.i[j-1], 3]
	        	
	                           # # # }

                   # # # increment.trajectory<-matrix(c(increment.x, increment.y, increment.z), ncol=3, nrow=(track.length[i]-1))
                   
                   
  
                   # # # number.of.neighbors.full.data[index.i]<-number.of.neighbors

                   # # # angular.change.complete.full.data[index.i[2:(track.length[i]-1)]]<-u$psi.component

                   # # # radius.complete.full.data[index.i[2:(track.length[i]-1)]]<-(u$radius[2:(track.length[i]-1)])

               
                   # # # }
                   
                   
                                     
                    setwd("/Users/e2torkas/Desktop/AOAS-R-Programs/")
                    
                    
                    # load(".RHeldOutDataFourClusters")
                   
                    load(".RHeldOutData")
 


                    
                  number.of.neighbors.full.data.held.out<-number.of.neighbors.full.data  
                  
                  
                  angular.change.complete.full.data.held.out<-angular.change.complete.full.data
                  
                  
                  radius.complete.full.data.held.out<-radius.complete.full.data



                  data_cell_motility_held_out<-data_cell_motility


###################################################################################################################################################################################################################################################

                				
				
				index.angular.change.held.out<-NULL
				
				for(i in 1: number.of.cells){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[i])
					
					
					 
					 index.angular.change.held.out<-c(index.angular.change.held.out, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				
				index.angular.change.cluster.held.out.1<-NULL
				
				for(i in 1:length(cluster.1.held.out)){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.1.held.out[i]])
					
					
					 
					 index.angular.change.cluster.held.out.1<-c(index.angular.change.cluster.held.out.1, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				
				index.angular.change.cluster.held.out.2<-NULL
				
				for(i in 1:length(cluster.2.held.out)){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.2.held.out[i]])
					
					
					 
					 index.angular.change.cluster.held.out.2<-c(index.angular.change.cluster.held.out.2, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				
				index.angular.change.cluster.held.out.3<-NULL
				
				for(i in 1:length(cluster.3.held.out)){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.3.held.out[i]])
					
					
					 
					 index.angular.change.cluster.held.out.3<-c(index.angular.change.cluster.held.out.3, index.i[2:(length(index.i)-1)])
					
					
				}
				
				
				
				
				index.angular.change.cluster.held.out.4<-NULL
				
				for(i in 1:length(cluster.4.held.out)){
					
					
					 index.i<-which(data_cell_motility[, 8]==id[cluster.4.held.out[i]])
					
					
					 
					 index.angular.change.cluster.held.out.4<-c(index.angular.change.cluster.held.out.4, index.i[2:(length(index.i)-1)])
					
					
				}
				
				


#######################################################################################################################################################################################################################################





setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/2-clusters-NHMM-1/")


load(file=paste(".RData", sep=""))





index.angular.change <- index.angular.change.cluster.1





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
###  Defining the ratio for the purpose of classification

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

1:black <- non-persistent

2:darkred <- persistent


index.persistent<-which(classification.matrix[index.angular.change, 2]==2)

index.non.persistent<-which(classification.matrix[index.angular.change, 2]==1)


For now, as the program includes the radius devided by max for others.

# setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/cluster-2-hyper-informative-second-prior-2-states/")
# load(file=paste(".Rdata", sep=""))



persist<-radius.complete.full.data[index.angular.change[index.persistent]]

non.persist<-radius.complete.full.data[index.angular.change[index.non.persistent]]


		
##################################################################################################################################################################################################################################################################################################  Held-Out Data  #########################################################

                data_cell_motility<-data_cell_motility_held_out
             
				number.of.cells<-length(unique(data_cell_motility[, 8]))

				id<-seq(1000000000, (1000000000+number.of.cells-1), 1)
				
				track.length<-NULL
				
                for(i in 1: number.of.cells){
                	
                	track.length[i]<-length(which(data_cell_motility[, 8]==id[i]))
                	
                }

                number.of.positions.total<-sum(track.length)
                
                
				number.of.neighbors.full.data<-matrix(number.of.neighbors.full.data.held.out, nrow= number.of.positions.total)
				
				
				angular.change.complete.full.data<-matrix(angular.change.complete.full.data.held.out, nrow= number.of.positions.total)
							
				
                index.angular.change<- c(index.angular.change.cluster.held.out.1, index.angular.change.cluster.held.out.2)
				

                held.out.points <- angular.change.complete.full.data[index.angular.change]/pi

# #################################################
# #################################################


load(".Rdata")



alpha_non_persistent_radius<-fit_radius_ss$summary[1]

beta_non_persistent_radius<-fit_radius_ss$summary[2]


alpha_persistent_radius<-fit_radius_ss$summary[3]

beta_persistent_radius<-fit_radius_ss$summary[4]



index.angular.change <- index.angular.change.cluster.1

             
       
      
     setwd("/Users/e2torkas/Desktop/Clustering-independent-of-direction/Revision-Sep-11/")
       
     
    
      pdf("Rplot-Angular-Change-2-states-NHMM-2-clusters.pdf", width=10, height=8, paper='special') 
 
 
      
     
     # i<-2

     # m <- rbind(cbind(c(1, 1, 1), c(2, 3, 4)), c(5))

     m <- cbind(rbind(c(1, 2), c(3, 4)))
     
     # m <- cbind(rbind(c(1, 2), c(3, 4)))
     
     # par(mar=c(4,7,2,1)) 
     
     # m <- cbind(c(1, 1, 1), c(2, 3, 4))
     
     # m <- rbind(c(5), cbind(c(1, 1, 1), c(2, 3, 4)))

     layout(m)
     

     par(lwd=2.5, mar=c(5, 5, 4, 2))

# par(mfrow=c(2, 2), lwd=2, mar=c(5, 5, 4, 2))




                 
                  hist(held.out.points, freq=F, xlim=c(0, 1), breaks=10, xlab=expression(paste(eta, sep="")), main="HP", ylim=c(0, 4), cex.lab=1.5, cex.axis=2.25, cex.main=1.7, col="#999999")

                  curve(weight[1]*dbeta(x, alpha_non_persistent, beta_non_persistent)+weight[2]*dbeta(x, alpha_persistent, beta_persistent), add=T)

   
                  curve(weight[1]*dbeta(x, alpha_non_persistent, beta_non_persistent), lty=3, add=T, col="darkred")
                  
                  
                  curve(weight[2]*dbeta(x, alpha_persistent, beta_persistent), lty=3, add=T, col="yellow")
                  
                  
                  
                           
                probDist <- (weight[1]*pbeta(angular.change.complete.full.data[index.angular.change]/pi, alpha_non_persistent, beta_non_persistent)+weight[2]*pbeta(angular.change.complete.full.data[index.angular.change]/pi, alpha_persistent, beta_persistent))
                     
                probTraining <- ppoints(length(held.out.points)) 
               
                                 
                length.min<-min(c(length(probDist), length(probTraining)))
                
                
                x <- sample(probDist, length.min, replace=F)
                
                
                y <- sample(probTraining, length.min, replace=F)
                
                
                x <- sort(x)
                
                
                y <- sort(y)
                
                
                
 
                plot(x, y, main="HP P-P Plot", xlab=" ", ylab=" ", cex.lab=1.5, cex.axis=1, cex.main=1.7, cex=0.05)



                abline(0,1, col="darkblue") 
                 

# ###################################################################################################################################################################################################################################################





setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/2-clusters-NHMM-2/")



load(file=paste(".RData", sep=""))





index.angular.change <- index.angular.change.cluster.2





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
###  Defining the ratio for the purpose of classification

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

1:black <- non-persistent

2:darkred <- persistent


index.persistent<-which(classification.matrix[index.angular.change, 2]==2)

index.non.persistent<-which(classification.matrix[index.angular.change, 2]==1)


For now, as the program includes the radius devided by max for others.

# setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Clustering-Results/cluster-2-hyper-informative-second-prior-2-states/")
# load(file=paste(".Rdata", sep=""))



persist<-radius.complete.full.data[index.angular.change[index.persistent]]

non.persist<-radius.complete.full.data[index.angular.change[index.non.persistent]]


		
##################################################################################################################################################################################################################################################################################################  Held-Out Data  #########################################################

                data_cell_motility<-data_cell_motility_held_out
             
				number.of.cells<-length(unique(data_cell_motility[, 8]))

				id<-seq(1000000000, (1000000000+number.of.cells-1), 1)
				
				track.length<-NULL
				
                for(i in 1: number.of.cells){
                	
                	track.length[i]<-length(which(data_cell_motility[, 8]==id[i]))
                	
                }

                number.of.positions.total<-sum(track.length)
                
                
				number.of.neighbors.full.data<-matrix(number.of.neighbors.full.data.held.out, nrow= number.of.positions.total)
				
				
				angular.change.complete.full.data<-matrix(angular.change.complete.full.data.held.out, nrow= number.of.positions.total)
							
				
                index.angular.change<-c(index.angular.change.cluster.held.out.3, index.angular.change.cluster.held.out.4)
				

                held.out.points <- angular.change.complete.full.data[index.angular.change]/pi

# #################################################
# #################################################


load(".Rdata")



alpha_non_persistent_radius<-fit_radius_ss$summary[1]

beta_non_persistent_radius<-fit_radius_ss$summary[2]


alpha_persistent_radius<-fit_radius_ss$summary[3]

beta_persistent_radius<-fit_radius_ss$summary[4]



index.angular.change <- index.angular.change.cluster.2




                 
                  hist(held.out.points, freq=F, xlim=c(0, 1), breaks=10, xlab=expression(paste(eta, sep="")), main="NP", ylim=c(0, 4), cex.lab=1.5, cex.axis=2.25, cex.main=1.7, col="#999999")

                  curve(weight[1]*dbeta(x, alpha_non_persistent, beta_non_persistent)+weight[2]*dbeta(x, alpha_persistent, beta_persistent), add=T)

   
                  curve(weight[1]*dbeta(x, alpha_non_persistent, beta_non_persistent), lty=3, add=T, col="darkred")
                  
                  
                  curve(weight[2]*dbeta(x, alpha_persistent, beta_persistent), lty=3, add=T, col="yellow")
                  
                  
                  
                           
                probDist <- (weight[1]*pbeta(angular.change.complete.full.data[index.angular.change]/pi, alpha_non_persistent, beta_non_persistent)+weight[2]*pbeta(angular.change.complete.full.data[index.angular.change]/pi, alpha_persistent, beta_persistent))
                     
                probTraining <- ppoints(length(held.out.points)) 
               
                                 
                length.min<-min(c(length(probDist), length(probTraining)))
                
                
                x <- sample(probDist, length.min, replace=F)
                
                
                y <- sample(probTraining, length.min, replace=F)
                
                
                x <- sort(x)
                
                
                y <- sort(y)
                
                
                
 
                plot(x, y, main="NP P-P Plot", xlab=" ", ylab=" ", cex.lab=1.5, cex.axis=1, cex.main=1.7, cex=0.05)



                abline(0,1, col="darkblue") 
                 
                 

                
                dev.off()
      
 
 

