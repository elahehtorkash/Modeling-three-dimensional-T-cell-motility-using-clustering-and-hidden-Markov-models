######  Graphs of AoAS

remove(list=ls())

setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Motility-Data")


require("mgcv")

library("ggmap")

require("lattice")

require("MotilityLab")

library(moments)

# library(VineCopula)

library(cluster)

library(scatterplot3d)


cluster.1<-c(1,   2,   6,   8,   9,  11,  12,  17,  18,  19,  20,  21,  22,  46,  50,  52,  54,  57,  58,  60,  62,  68,  69,  72,  74,  77,  80,  83,  84,  85,  86,  90,  91, 104, 105, 113, 120, 121, 128, 136, 137, 142, 156, 157, 168, 169, 171, 174, 177, 182, 193, 195, 199, 203, 211, 217, 227, 231, 233, 238,
  240, 241, 242, 243, 248, 252, 254, 266, 271, 279, 282, 283, 287, 288, 296, 298, 302, 306, 308, 309, 310, 314, 317, 332, 333, 334, 335, 339, 342, 344,
 345, 348, 351, 353, 355, 357, 360, 364, 365, 369, 373, 374, 375, 381, 386, 388, 389, 390, 398, 400, 401, 404, 407, 416, 420, 423, 426, 429, 430, 434,
 435, 438, 439, 440, 442, 448, 449, 455, 459, 462, 465, 468, 471, 475, 478, 480, 484, 485, 486, 492, 494, 495, 498, 503, 504, 507, 514, 522, 523, 526,
  527, 528, 532)



cluster.2<-c(3,   7,  10,  13,  15,  16,  23,  24,  28,  32,  35,  36,  37,  39,  45,  53,  55,  71,  73,  75,  78,  92,  93, 103, 108, 109, 110, 111, 112, 123,
  129, 147, 149, 151, 172, 178, 181, 187, 188, 194, 197, 201, 204, 210, 214, 219, 244, 267, 268, 273, 285, 297, 312, 315, 316, 319, 320, 321, 325, 326,
  328, 330, 343, 347, 354, 356 , 379, 383, 384, 392, 393, 399, 403, 410, 412, 413, 422, 424, 425, 433, 436, 437, 445, 446, 447, 451, 453, 458, 463, 467,
  469, 476, 477, 481, 482, 497, 500, 505, 506, 509, 511, 512, 513, 517, 518, 520, 521, 524, 525, 531, 534)


cluster.3<-c( 4,   5,  14,  25,  26,  29,  30,  31,  33,  34,  38,  40,  41,  47,  48  ,51,  59,  61,  63,  64,  65,  66,  67,  70,  79,  81,  82,  87,  88,  94,
  95,  96,  97,  98, 100, 101, 102, 106, 107, 114, 115, 116, 118, 122, 124, 125, 126, 127, 130, 131, 132, 134, 135, 138, 139, 140, 141, 143, 144, 146,
   150, 152, 153, 154, 158, 160, 161, 162 ,163, 164, 165, 167, 170, 175, 180, 184, 186, 190, 191, 196, 200, 202, 205, 206, 207, 208, 212, 215, 216, 218,
   220, 221, 222, 223, 224, 225, 226, 230, 232, 234, 235, 236, 239, 247, 250, 251, 253, 255, 256, 258, 259, 260, 261, 262, 263, 264, 265, 270, 272, 274,
 276, 277, 278, 284, 286, 289, 290, 291, 292, 295, 300, 303, 304, 305, 307, 311, 318, 322, 323, 324, 327, 329, 331, 336, 337, 338, 340, 346, 349, 350,
  352, 361, 362, 363, 366, 367, 368, 370 ,371 ,372, 376, 382, 385, 387, 391, 395, 396, 402, 405 ,406, 409, 415, 418 ,419, 421, 432, 441, 443, 444, 450,
  452, 454, 456, 457, 460, 461, 464, 466, 470, 472, 474, 479, 487, 489, 490, 493, 502 ,508, 510, 515, 516 ,519, 530, 533)


cluster.4<-c(27,  42,  43,  44,  49,  56 , 76,  89,  99, 117, 119, 133, 145 ,148, 155 ,159, 166 ,173, 176, 179, 183, 185, 189, 192, 198, 209, 213, 228, 229, 237, 245,
  246, 249, 257, 269, 275, 280, 281, 293, 294, 299, 301, 313, 341, 358, 359, 377, 378, 380, 394, 397, 408, 411, 414, 417, 427, 428, 431, 473, 483, 488, 491,
 496, 499, 501, 529)
 
 
 
 
##########################################################################################################################
##########################################################################################################################

 
#######   clsuters, Joel comments


cluster.1<-c(1,   2,   6,   8,   9,  12,  17,  18,  19,  20,  21,  22,  46,  52,  54,  57,  60,  62,  68,  69,  77,  83,  84,  85,  86,  90,  91, 104, 105, 113, 124, 127, 128, 132, 135, 150, 151, 152, 153, 157, 160, 162, 163, 166, 169, 171, 173, 178, 183, 187, 191, 192, 193, 199, 206, 207, 216, 218, 234, 238, 241, 244, 247, 248, 252, 253, 256, 257, 258, 260, 266, 267, 273, 276, 277, 283, 289, 293, 296, 298, 302, 303, 304, 310, 313, 321, 322, 324, 325, 332, 341, 345, 346)			                                     


cluster.2<-c( 3,   7,  10,  13,  15,  16,  23,  24,  28,  32,  35,  36,  37,  39,  45,  53,  55,  71,  73,  75,  78,  92,  93, 103, 108, 109, 110, 111, 130, 133, 134, 137, 139, 143, 144, 146, 148, 161, 165, 172, 174, 197, 201, 202, 204, 210, 211, 217, 221, 228, 230, 231, 240, 242, 243, 251, 254, 255, 263, 264, 265, 269, 271, 281, 287, 294, 295, 299, 300, 315, 318, 323, 327, 329, 330, 331, 335, 336, 338, 339, 342, 343, 349, 352)



cluster.3<-c( 4,   5,  11,  14,  25,  26,  29,  30,  31,  33,  34,  38,  40,  41,  47 , 48 , 50 , 51 , 58,  59,  61,  63,  64,  65,  66,  67,  70,  72,  74,  79,  80,  81,  82,  87,  88,  94,  95,  96 , 97,  98, 100, 101, 102, 106, 107, 114, 115, 116, 118, 120 ,121, 122, 123, 125, 126, 129, 136, 140, 141, 142, 145, 147, 149, 154, 155, 156, 158, 164, 167, 168, 170, 175, 179, 180 ,181, 182, 184, 185, 186, 188 ,189, 190 ,194, 200, 203, 205, 208, 209, 213, 214, 219, 220, 222, 223, 224, 225 ,227, 233, 236, 237, 239 ,250, 259, 261, 262, 268, 270, 272, 274, 275, 278, 279, 280, 282, 284, 286, 288, 290 ,292, 297, 305, 307 ,308, 311, 312, 316, 320, 326, 328, 333, 334, 337, 340, 344, 348, 350, 351)


cluster.4<-c(27,  42,  43,  44,  49,  56 , 76,  89,  99, 112, 117, 119, 131, 138, 159 ,176 ,177 ,195 ,196 ,198 ,212, 215, 226, 229, 232, 235, 245, 246, 249, 285, 291, 301, 306, 309, 314, 317, 319, 347)


 
 




        data_cell_motility<-read.csv("Position_1.csv", header=TRUE)

        dim(data_cell_motility)


#########   In order to have an understanding of the number of T-cells observed in the study, we use the following comment

        number.of.cells<-length(unique(data_cell_motility[, 8]))
###  [1] 117

		id<-seq(1000000000, (1000000000+ number.of.cells-1), 1)
		
		times.present.in.study<-NULL
		
		tracks<-list()      


            
    radius.phi.theta.psi<-function(trajectory, trajectory.num=1){
					
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



setwd("/Users/e2torkas/Desktop/Clustering-independent-of-direction/Revision-Sep-11/")





par(mar=c(4,4.5,2,1))

par(oma=c(0,0,0,0))


###  Cluster 1: 

# # # #       quartz()

      # # pWidth = 10

      # # pHeight = 8

      # # plot.window(c(0,pWidth),
             # # c(0,pHeight))  
             
     
      pdf("cluster_1.pdf", width=10, height=10, paper='special') 
 
       
      par(lwd=2.5) 
       
     
     i<-3

     m <- rbind(cbind(c(1, 1, 1), c(2, 3, 4)), c(5))
     
     # par(mar=c(4,7,2,1)) 
     
     # m <- cbind(c(1, 1, 1), c(2, 3, 4))
     
     # m <- rbind(c(5), cbind(c(1, 1, 1), c(2, 3, 4)))

     layout(m)
     
     ### layout.show(5)  ### to see how it works
     
     index<-which(data_cell_motility[, 8]==id[i])
	  
	  data_i<-matrix(c(data_cell_motility[index, 7], data_cell_motility[index, 1], data_cell_motility[index, 2], data_cell_motility[index, 3]), nrow=length(index), ncol=4)
                   	
                   	x<-tracks(data_i)
                 	
                   u<-radius.phi.theta.psi(trajectory=x$'1')
                 
                   increment.x<-increment.y<-increment.z<-NULL
				    
				    increment.code<-matrix(0, nrow=3, ncol=1)

                    for(j in 2: length(index)){
	        	
	        	             increment.x[j-1]<-data_cell_motility[index[j], 1]-data_cell_motility[index[j-1], 1]
	        	
	        	             increment.y[j-1]<-data_cell_motility[index[j], 2]-data_cell_motility[index[j-1], 2]
	        	
	        	             increment.z[j-1]<-data_cell_motility[index[j], 3]-data_cell_motility[index[j-1], 3]
	        	
	                }
   	
	    					                                                                    
	 x.min<-min(data_cell_motility[index, 1])-5
		
	 x.max<-max(data_cell_motility[index, 1])+5
			
	 t.x<-0.5*(x.min+x.max)				
			
	 y.min<-min(data_cell_motility[index, 2])-5
						
	 y.max<-max(data_cell_motility[index, 2])+5
			
	 t.y<-0.5*(y.min+y.max)
			
	 z.min<-0
			
	 z.max<-max(data_cell_motility[index, 3])+5
	 
	 
	 # pdf( "cluster_1.pdf", width = 11, height = 8 )
	 
	 	 
	 s3d <- scatterplot3d(data_cell_motility[index, 1:3], pch = 16, 
              xlab = expression(paste("X (", mu, "m)", sep="" )),
              ylab = expression(paste("Y (", mu, "m)", sep="" )),
              zlab = expression(paste("Z (", mu, "m)", sep="" )), angle=70, type="l", cex.lab=1, cex.axis=1, cex.main=2,  xlim=c(x.min, x.max), ylim=c(y.min, y.max), zlim=c(z.min, z.max), color="darkred", main="(a)")	 
	 			

     text(s3d$xyz.convert(x=data_cell_motility[index[1], 1], y=data_cell_motility[index[1], 2], z=data_cell_motility[index[1], 3]), labels="Initial Point", pos=1, col="black") 
     
     
     text(s3d$xyz.convert(x=data_cell_motility[index[length(index)], 1], y=data_cell_motility[index[length(index)], 2], z=data_cell_motility[index[length(index)],3]), labels="Final Point", pos=1, col="black")
     
     
       hist(increment.x, ylab=" ", xlab=expression(paste("X (", mu, "m)", sep="" )), breaks=15, border="darkgreen", 
     col="green", cex.lab=2, cex.axis=1.5, cex.main=2,  freq=F, main="(b)", ylim=c(0, 0.55))
     
      hist(increment.y, ylab=" ", xlab=expression(paste("Y (", mu, "m)", sep="" )), breaks=15, border="darkblue", 
     col="blue", cex.lab=2, cex.axis=1.5, cex.main=2,  freq=F, main="(c)", ylim=c(0,0.55))
     
      hist(increment.z, ylab=" ", xlab=expression(paste("Z (", mu, "m)", sep="" )), breaks=15, border="darkred", 
     col="red", cex.lab=2, cex.axis=1.5, cex.main=2,  freq=F, main="(d)", ylim=c(0,0.55)) 
         
          
       plot( u$psi.component, xlab="Time step" , ylab=expression(paste(psi)), type="l", ylim=c(0, pi), pch=20, col="darkred", cex.lab=2, cex.axis=2, cex.main=2,  main="(e)", yaxt="n")
                   
                   
      # x.breaks<-length(u$psi.component) 
                   
      # seq.x.breaks<-seq(1, x.breaks, 1)
                   
      # grid(x.breaks)
                   
      abline(h=pi/2, lty=3)

      axis(side=2, at=c(0, pi/2, pi), labels=c(0, expression(paste(pi/2)), expression(paste(pi))), cex.axis=1.5,  las=2, ylab=expression(paste(psi)))


# #      dev.copy(jpeg,
          # filename="cluster_1.jpeg");
          
          dev.off()


par(mar=c(4,4.5,2,1))

par(oma=c(0,0,0,0))

###  Cluster 2:  

      # # # quartz()

      # # # pWidth = 10

      # # # pHeight = 8

      # # # plot.window(c(0,pWidth),
             # # # c(0,pHeight))         
        
     
      pdf("cluster_2.pdf", width=10, height=10, paper='special')    
        
     par(lwd=2.5)
     
     
     i<-6

     m <- rbind(cbind(c(1, 1, 1), c(2, 3, 4)), c(5))
     
     # par(mar=c(4,7,2,1)) 
     
     # m <- cbind(c(1, 1, 1), c(2, 3, 4))
     
     # m <- rbind(c(5), cbind(c(1, 1, 1), c(2, 3, 4)))

     layout(m)
     
     ### layout.show(5)  ### to see how it works
     
     index<-which(data_cell_motility[, 8]==id[i])
	  
	  data_i<-matrix(c(data_cell_motility[index, 7], data_cell_motility[index, 1], data_cell_motility[index, 2], data_cell_motility[index, 3]), nrow=length(index), ncol=4)
                   	
                   	x<-tracks(data_i)
                 	
                   u<-radius.phi.theta.psi(trajectory=x$'1')
                 
                   increment.x<-increment.y<-increment.z<-NULL
				    
				    increment.code<-matrix(0, nrow=3, ncol=1)

                    for(j in 2: length(index)){
	        	
	        	             increment.x[j-1]<-data_cell_motility[index[j], 1]-data_cell_motility[index[j-1], 1]
	        	
	        	             increment.y[j-1]<-data_cell_motility[index[j], 2]-data_cell_motility[index[j-1], 2]
	        	
	        	             increment.z[j-1]<-data_cell_motility[index[j], 3]-data_cell_motility[index[j-1], 3]
	        	
	                }
   	
	    					                                                                    
	 x.min<-min(data_cell_motility[index, 1])-5
		
	 x.max<-max(data_cell_motility[index, 1])+5
			
	 t.x<-0.5*(x.min+x.max)				
			
	 y.min<-min(data_cell_motility[index, 2])-5
						
	 y.max<-max(data_cell_motility[index, 2])+5
			
	 t.y<-0.5*(y.min+y.max)
			
	 z.min<-0
			
	 z.max<-max(data_cell_motility[index, 3])+5
	 
 	 
	 # pdf( "cluster_1.pdf", width = 11, height = 8 )
	 
	 	 
	 s3d <- scatterplot3d(data_cell_motility[index, 1:3], pch = 16, 
              xlab = expression(paste("X (", mu, "m)", sep="" )),
              ylab = expression(paste("Y (", mu, "m)", sep="" )),
              zlab = expression(paste("Z (", mu, "m)", sep="" )), angle=70, type="l", cex.lab=1, cex.axis=1, cex.main=2,  xlim=c(x.min, x.max), ylim=c(y.min, y.max), zlim=c(z.min, z.max), color="darkred", main="(a)")	 
	 			

     text(s3d$xyz.convert(x=data_cell_motility[index[1], 1], y=data_cell_motility[index[1], 2], z=data_cell_motility[index[1], 3]), labels="Initial Point", pos=1, col="black") 
     
     
     text(s3d$xyz.convert(x=data_cell_motility[index[length(index)], 1], y=data_cell_motility[index[length(index)], 2], z=data_cell_motility[index[length(index)],3]), labels="Final Point", pos=1, col="black")
     
     
       hist(increment.x, ylab=" ", xlab=expression(paste("X (", mu, "m)", sep="" )), breaks=15, border="darkgreen", 
     col="green", cex.lab=2, cex.axis=1.5, cex.main=2,  freq=F, main="(b)", ylim=c(0, 0.55))
     
      hist(increment.y, ylab=" ", xlab=expression(paste("Y (", mu, "m)", sep="" )), breaks=15, border="darkblue", 
     col="blue", cex.lab=2, cex.axis=1.5, cex.main=2,  freq=F, main="(c)", ylim=c(0, 0.55))
     
      hist(increment.z, ylab=" ", xlab=expression(paste("Z (", mu, "m)", sep="" )), breaks=15, border="darkred", 
     col="red", cex.lab=2, cex.axis=1.5, cex.main=2,  freq=F, main="(d)", ylim=c(0, 0.55)) 
         
         
       plot( u$psi.component, xlab="Time step" , ylab=expression(paste(psi)), type="l", ylim=c(0, pi), pch=20, col="darkred", cex.lab=2, cex.axis=2, cex.main=2,  main="(e)", yaxt="n")
                   
                   
      # x.breaks<-length(u$psi.component) 
                   
      # seq.x.breaks<-seq(1, x.breaks, 1)
                   
      # grid(x.breaks)
                   
      abline(h=pi/2, lty=3)

      axis(side=2, at=c(0, pi/2, pi), labels=c(0, expression(paste(pi/2)), expression(paste(pi))), cex.axis=1.5,  las=2, ylab=expression(paste(psi)))




      
# # # #       dev.copy(jpeg,
          # # filename="cluster_2.jpeg");

      dev.off()
      
    

par(mar=c(4,4.5,2,1))

par(oma=c(0,0,0,0))
      
     ###  Cluster 3:  
     
     
      pdf("cluster_3.pdf", width=10, height=10, paper='special') 
 

      par(lwd=2.5)
     
      # # # quartz()

      # # # pWidth = 10

      # # # pHeight = 8

      # # # plot.window(c(0,pWidth),
             # # # c(0,pHeight))                
     
     i<-1

     m <- rbind(cbind(c(1, 1, 1), c(2, 3, 4)), c(5))
     
     # par(mar=c(4,7,2,1)) 
     
     # m <- cbind(c(1, 1, 1), c(2, 3, 4))
     
     # m <- rbind(c(5), cbind(c(1, 1, 1), c(2, 3, 4)))

     layout(m)
     
     ### layout.show(5)  ### to see how it works
     
     index<-which(data_cell_motility[, 8]==id[i])
	  
	  data_i<-matrix(c(data_cell_motility[index, 7], data_cell_motility[index, 1], data_cell_motility[index, 2], data_cell_motility[index, 3]), nrow=length(index), ncol=4)
                   	
                   	x<-tracks(data_i)
                 	
                   u<-radius.phi.theta.psi(trajectory=x$'1')
                 
                   increment.x<-increment.y<-increment.z<-NULL
				    
				    increment.code<-matrix(0, nrow=3, ncol=1)

                    for(j in 2: length(index)){
	        	
	        	             increment.x[j-1]<-data_cell_motility[index[j], 1]-data_cell_motility[index[j-1], 1]
	        	
	        	             increment.y[j-1]<-data_cell_motility[index[j], 2]-data_cell_motility[index[j-1], 2]
	        	
	        	             increment.z[j-1]<-data_cell_motility[index[j], 3]-data_cell_motility[index[j-1], 3]
	        	
	                }
   	
	    					                                                                    
	 x.min<-min(data_cell_motility[index, 1])-5
		
	 x.max<-max(data_cell_motility[index, 1])+5
			
	 t.x<-0.5*(x.min+x.max)				
			
	 y.min<-min(data_cell_motility[index, 2])-5
						
	 y.max<-max(data_cell_motility[index, 2])+5
			
	 t.y<-0.5*(y.min+y.max)
			
	 z.min<-0
			
	 z.max<-max(data_cell_motility[index, 3])+5
	 
	 
	 
	 # pdf( "cluster_1.pdf", width = 11, height = 8 )
	 	 
	 s3d <- scatterplot3d(data_cell_motility[index, 1:3], pch = 16, 
              xlab = expression(paste("X (", mu, "m)", sep="" )),
              ylab = expression(paste("Y (", mu, "m)", sep="" )),
              zlab = expression(paste("Z (", mu, "m)", sep="" )), angle=70, type="l", cex.lab=1, cex.axis=1, cex.main=2,  xlim=c(x.min, x.max), ylim=c(y.min, y.max), zlim=c(z.min, z.max), color="darkred", main="(a)")	 
	 			

     text(s3d$xyz.convert(x=data_cell_motility[index[1], 1], y=data_cell_motility[index[1], 2], z=data_cell_motility[index[1], 3]), labels="Initial Point", pos=1, col="black") 
     
     
     text(s3d$xyz.convert(x=data_cell_motility[index[length(index)], 1], y=data_cell_motility[index[length(index)], 2], z=data_cell_motility[index[length(index)],3]), labels="Final Point", pos=1, col="black")
     
     
       hist(increment.x, ylab=" ", xlab=expression(paste("X (", mu, "m)", sep="" )), breaks=15, border="darkgreen", 
     col="green", cex.lab=2, cex.axis=1.5, cex.main=2,  freq=F, main="(b)", ylim=c(0, 0.55))
     
      hist(increment.y, ylab=" ", xlab=expression(paste("Y (", mu, "m)", sep="" )), breaks=15, border="darkblue", 
     col="blue", cex.lab=2, cex.axis=1.5, cex.main=2,  freq=F, main="(c)", ylim=c(0, 0.55))
     
      hist(increment.z, ylab=" ", xlab=expression(paste("Z (", mu, "m)", sep="" )), breaks=15, border="darkred", 
     col="red", cex.lab=2, cex.axis=1.5, cex.main=2,  freq=F, main="(d)", ylim=c(0, 0.8)) 

         
       plot( u$psi.component, xlab="Time step" , ylab=expression(paste(psi)), type="l", ylim=c(0, pi), pch=20, col="darkred", cex.lab=2, cex.axis=2, cex.main=2,  main="(e)", yaxt="n")
                   
                   
      # x.breaks<-length(u$psi.component) 
                   
      # seq.x.breaks<-seq(1, x.breaks, 1)
                   
      # grid(x.breaks)
                   
      abline(h=pi/2, lty=3)

      axis(side=2, at=c(0, pi/2, pi), labels=c(0, expression(paste(pi/2)), expression(paste(pi))), cex.axis=1.5,  las=2, ylab=expression(paste(psi)))


      
# # # #       dev.copy(jpeg,
          # # filename="cluster_3.jpeg");

      dev.off()
      
      
      
	par(mar=c(4,4.5,2,1))
	
	par(oma=c(0,0,0,0))


      
     ###  Cluster 4:          
     
     
      # # # quartz()

# # # #       pWidth = 10

      # # pHeight = 8

      # # plot.window(c(0,pWidth),
             # # c(0,pHeight))  
             
             
     

      pdf("cluster_4.pdf", width=10, height=10, paper='special') 
 
 
      par(lwd=2.5)
     
     i<-2

     m <- rbind(cbind(c(1, 1, 1), c(2, 3, 4)), c(5))
     
     # par(mar=c(4,7,2,1)) 
     
     # m <- cbind(c(1, 1, 1), c(2, 3, 4))
     
     # m <- rbind(c(5), cbind(c(1, 1, 1), c(2, 3, 4)))

     layout(m)
     
     ### layout.show(5)  ### to see how it works
     
     index<-which(data_cell_motility[, 8]==id[i])
	  
	  data_i<-matrix(c(data_cell_motility[index, 7], data_cell_motility[index, 1], data_cell_motility[index, 2], data_cell_motility[index, 3]), nrow=length(index), ncol=4)
                   	
                   	x<-tracks(data_i)
                 	
                   u<-radius.phi.theta.psi(trajectory=x$'1')
                 
                   increment.x<-increment.y<-increment.z<-NULL
				    
				    increment.code<-matrix(0, nrow=3, ncol=1)

                    for(j in 2: length(index)){
	        	
	        	             increment.x[j-1]<-data_cell_motility[index[j], 1]-data_cell_motility[index[j-1], 1]
	        	
	        	             increment.y[j-1]<-data_cell_motility[index[j], 2]-data_cell_motility[index[j-1], 2]
	        	
	        	             increment.z[j-1]<-data_cell_motility[index[j], 3]-data_cell_motility[index[j-1], 3]
	        	
	                }
   	
	    					                                                                    
	 x.min<-min(data_cell_motility[index, 1])-5
		
	 x.max<-max(data_cell_motility[index, 1])+5
			
	 t.x<-0.5*(x.min+x.max)				
			
	 y.min<-min(data_cell_motility[index, 2])-5
						
	 y.max<-max(data_cell_motility[index, 2])+5
			
	 t.y<-0.5*(y.min+y.max)
			
	 z.min<-0
			
	 z.max<-max(data_cell_motility[index, 3])+5
	 
	 
	 # pdf( "cluster_1.pdf", width = 11, height = 8 )
	 
	 	 
	 s3d <- scatterplot3d(data_cell_motility[index, 1:3], pch = 16, 
              xlab = expression(paste("X (", mu, "m)", sep="" )),
              ylab = expression(paste("Y (", mu, "m)", sep="" )),
              zlab = expression(paste("Z (", mu, "m)", sep="" )), angle=70, type="l", cex.lab=1, cex.axis=1, cex.main=2,  xlim=c(x.min, x.max), ylim=c(y.min, y.max), zlim=c(z.min, z.max), color="darkred", main="(a)")	 
	 			

     text(s3d$xyz.convert(x=data_cell_motility[index[1], 1], y=data_cell_motility[index[1], 2], z=data_cell_motility[index[1], 3]), labels="Initial Point", pos=1, col="black") 
     
     
     text(s3d$xyz.convert(x=data_cell_motility[index[length(index)], 1], y=data_cell_motility[index[length(index)], 2], z=data_cell_motility[index[length(index)],3]), labels="Final Point", pos=1, col="black")
     
     
       hist(increment.x, ylab=" ", xlab=expression(paste("X (", mu, "m)", sep="" )), breaks=15, border="darkgreen", 
     col="green", cex.lab=2, cex.axis=1.5, cex.main=2,  freq=F, main="(b)", ylim=c(0, 0.55))
     
      hist(increment.y, ylab=" ", xlab=expression(paste("Y (", mu, "m)", sep="" )), breaks=15, border="darkblue", 
     col="blue", cex.lab=2, cex.axis=1.5, cex.main=2,  freq=F, main="(c)", ylim=c(0, 0.55))
     
      hist(increment.z, ylab=" ", xlab=expression(paste("Z (", mu, "m)", sep="" )), breaks=15, border="darkred", 
     col="red", cex.lab=2, cex.axis=1.5, cex.main=2,  freq=F, main="(d)", ylim=c(0, 0.55)) 
         
       plot( u$psi.component, xlab="Time step" , ylab=expression(paste(psi)), type="l", ylim=c(0, pi), pch=20, col="darkred", cex.lab=2, cex.axis=2, cex.main=2,  main="(e)", yaxt="n")
                   
                   
      # x.breaks<-length(u$psi.component) 
                   
      # seq.x.breaks<-seq(1, x.breaks, 1)
                   
      # grid(x.breaks)
                   
      abline(h=pi/2, lty=3)

      axis(side=2, at=c(0, pi/2, pi), labels=c(0, expression(paste(pi/2)), expression(paste(pi))), cex.axis=1.5,  las=2, ylab=expression(paste(psi)))

      
# #       dev.copy(jpeg,
          # filename="cluster_4.jpeg");
          
      dev.off()
      
      
      # install.packages("plotrix", dependencies=T)
      
      library(plotrix)
      
      quartz()

      pWidth = 8

      pHeight = 8

      plot.window(c(0,pWidth),
             c(0,pHeight))   
             
             
      plot(NA, xlim=c(0, 8), ylim=c(0, 5), xlab="X", ylab="Y")

      # arrows(0, 0, 4, 0, col=4, lwd=2)
      
      
      # arrows(0, 0, 3, 3, col=3, lwd=2)
      
      
      # arrows(0, 0, 3, 0, col=2, lwd=2, lty=3)
      
      # lines(c(3, 3), c(0, 3),  lty=2)
      
      
      vecs <- data.frame(vname=c("a","b","transb", ""),
                   x0=c(0, 0, 0, 4),y0=c(0, 0, 0, 0), x1=c(4, 3, 3, 7), y1=c(0, 3, 0, 3),  
                   col=c(4, 3, 1, 2))
      with( vecs, mapply("arrows", x0, y0, x1, y1, col=c(4, 3, 1, 2), lwd=3, lty=c(1, 1, 3, 1)) )
      
      
      with(vecs, mapply('text', x=x1[1:4]-.1, y=y1[1:4]+.1, 
  labels=expression(a, b, b*cos(theta), b ) ))
      
  
      lines(c(3, 3), c(0, 3),  lty=2)
             
             
      alpha_angle= 3.5*pi/12

      draw.arc(0, 0, 0.8, angle2=alpha_angle, col="darkred")

      text(1, 0.3, expression(alpha),cex=1.5)       
             
             
      dev.copy(jpeg,
          filename="persistency_concept.jpeg");
          
      dev.off()       
      
      
      
      
      
      
      #########################################   radius determination   #########################################################
      
      
      clusters<-list("1"=cluster.1, "2"=cluster.2, "3"=cluster.3, "4"=cluster.4)
      
      Q.3.list<-list()
      
      median.list<-list()
      
      for(l in 1:4){
      
	      name <- paste("", l, "",sep="")
	      
	      tmp <- list(third.quartile.radius[clusters[[l]]])
	      
	      Q.3.list[name]<-tmp
	      
	      name <- paste("", l, "",sep="")
	      
	      tmp.m <- list(median.radius[clusters[[l]]])
	      
	      median.list[name]<-tmp.m
 
               
}    


k_0<-4  

par(mfrow=c(2, 2))

for(l in 1: k_0){
	
	hist(Q.3.list[[l]], xlim=c(0, 1), freq=F, xlab= paste("Cluster ", l, sep=""), cex.lab=1.5, cex.axis=2.25, cex.main=2.25, lwd=5, main=" ",  breaks=seq(0, 1, 0.1), ylim=c(0, 4), border="black", col="#999999")
	
	# axis(side=1, at=seq(0, 1, 0.1), labels=seq(0, 1, 0.1), las=1)
	
}



k_0<-4  

par(mfrow=c(2, 2))

for(l in 1: k_0){
	
	hist(median.list[[l]], xlim=c(0, 1), freq=F, xlab= paste("Cluster ", l, sep=""), cex.lab=1.5, cex.axis=2, cex.main=2.5, lwd=3, main=" ",  breaks=seq(0, 1, 0.1), xaxt="n", col="#999999")
	
	axis(side=1, at=seq(0, 1, 0.1), labels=seq(0, 1, 0.1), las=1)
	
}


par(mfrow=c(2, 2))



