remove(list=ls())


setwd("/Users/e2torkas/Desktop/Cell-Motility-Project-Analysis/Motility-Data/")



library(matlib)

library(sde)

s.angle=function(x1,x2){ #this calculates the second angle going from increment x1 to x2

if( t(x1)%*%x2==sqrt(t(x1)%*%x1)*sqrt(t(x2)%*%x2) | t(x1)%*%x2==-sqrt(t(x1)%*%x1)*sqrt(t(x2)%*%x2)){return(0)}

	if( x1[3] != 0 ){
	A =cbind(x1,as.matrix(c(1,0,0)),as.matrix(c(0,1,0)))
	A.qr=qr(A)
      Qm=qr.Q(A.qr)
	t2=t(x2)%*%Qm[,2]
	t3=t(x2)%*%Qm[,3]
	if(t2<0){ theta=(atan(t3/t2)+pi)}
		if(t2>=0){  
			if(t3>0){theta=atan(t3/t2)} 
			if(t3<=0){theta= atan(t3/t2)+2*pi}
		}
	}	
	
	if( (x1[3]==0)& (x1[1]!=0)){
	A =cbind(x1,as.matrix(c(0,0,1)),as.matrix(c(0,1,0)))
	A.qr=qr(A)
      Qm=qr.Q(A.qr)
	t2=t(x2)%*%Qm[,2]
	t3=t(x2)%*%Qm[,3]
	if(t2<0){ theta=(atan(t3/t2)+pi)}
		if(t2>=0){  
			if(t3>0){theta=atan(t3/t2)} 
			if(t3<=0){theta= atan(t3/t2)+2*pi}
		}
	}
	if( (x1[3]==0)& (x1[1]==0)){
	A =cbind(x1,as.matrix(c(0,0,1)),as.matrix(c(1,0,0)))
	A.qr=qr(A)
      Qm=qr.Q(A.qr)
	t2=t(x2)%*%Qm[,2]
	t3=t(x2)%*%Qm[,3]
	if(t2<0){ theta=(atan(t3/t2)+pi)}
		if(t2>=0){  
			if(t3>0){theta=atan(t3/t2)} 
			if(t3<=0){theta= atan(t3/t2)+2*pi}
		}
	}

theta
}








############################################################################################################################
############################################################################################################################


                r.array<-c(0.25, 0.5, 1, 2, 3, 4, 5, 10, 15, 20)



                
				data_cell_motility_1<-read.csv("Position_1.csv", header=TRUE)

				dim(data_cell_motility_1)

				#########   In order to have an understanding of the number of T-cells observed in the study, we use the following comment

				number.of.cells_1<-length(unique(data_cell_motility_1[, 8]))


                
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
		        
		        k_0=4
		        
		        times.present.in.study<-NULL
		        
		        
# # # # cluster.1<-c(1,   2,   6,   8,   9,  11,  12,  17,  18,  19,  20,  21,  22,  46,  50,  52,  54,  57,  58,  60,  62,  68,  69,  72,  74,  77,  80,  83,  84,  85,  86,  90,  91, 104, 105, 113, 120, 121, 128, 136, 137, 142, 156, 157, 168, 169, 171, 174, 177, 182, 193, 195, 199, 203, 211, 217, 227, 231, 233, 238,
  # # # # 240, 241, 242, 243, 248, 252, 254, 266, 271, 279, 282, 283, 287, 288, 296, 298, 302, 306, 308, 309, 310, 314, 317, 332, 333, 334, 335, 339, 342, 344,
 # # # # 345, 348, 351, 353, 355, 357, 360, 364, 365, 369, 373, 374, 375, 381, 386, 388, 389, 390, 398, 400, 401, 404, 407, 416, 420, 423, 426, 429, 430, 434,
 # # # # 435, 438, 439, 440, 442, 448, 449, 455, 459, 462, 465, 468, 471, 475, 478, 480, 484, 485, 486, 492, 494, 495, 498, 503, 504, 507, 514, 522, 523, 526,
  # # # # 527, 528, 532)



# # # # cluster.2<-c(3,   7,  10,  13,  15,  16,  23,  24,  28,  32,  35,  36,  37,  39,  45,  53,  55,  71,  73,  75,  78,  92,  93, 103, 108, 109, 110, 111, 112, 123,
  # # # # 129, 147, 149, 151, 172, 178, 181, 187, 188, 194, 197, 201, 204, 210, 214, 219, 244, 267, 268, 273, 285, 297, 312, 315, 316, 319, 320, 321, 325, 326,
  # # # # 328, 330, 343, 347, 354, 356 , 379, 383, 384, 392, 393, 399, 403, 410, 412, 413, 422, 424, 425, 433, 436, 437, 445, 446, 447, 451, 453, 458, 463, 467,
  # # # # 469, 476, 477, 481, 482, 497, 500, 505, 506, 509, 511, 512, 513, 517, 518, 520, 521, 524, 525, 531, 534)


# # # # cluster.3<-c( 4,   5,  14,  25,  26,  29,  30,  31,  33,  34,  38,  40,  41,  47,  48  ,51,  59,  61,  63,  64,  65,  66,  67,  70,  79,  81,  82,  87,  88,  94,
  # # # # 95,  96,  97,  98, 100, 101, 102, 106, 107, 114, 115, 116, 118, 122, 124, 125, 126, 127, 130, 131, 132, 134, 135, 138, 139, 140, 141, 143, 144, 146,
   # # # # 150, 152, 153, 154, 158, 160, 161, 162 ,163, 164, 165, 167, 170, 175, 180, 184, 186, 190, 191, 196, 200, 202, 205, 206, 207, 208, 212, 215, 216, 218,
   # # # # 220, 221, 222, 223, 224, 225, 226, 230, 232, 234, 235, 236, 239, 247, 250, 251, 253, 255, 256, 258, 259, 260, 261, 262, 263, 264, 265, 270, 272, 274,
 # # # # 276, 277, 278, 284, 286, 289, 290, 291, 292, 295, 300, 303, 304, 305, 307, 311, 318, 322, 323, 324, 327, 329, 331, 336, 337, 338, 340, 346, 349, 350,
  # # # # 352, 361, 362, 363, 366, 367, 368, 370 ,371 ,372, 376, 382, 385, 387, 391, 395, 396, 402, 405 ,406, 409, 415, 418 ,419, 421, 432, 441, 443, 444, 450,
  # # # # 452, 454, 456, 457, 460, 461, 464, 466, 470, 472, 474, 479, 487, 489, 490, 493, 502 ,508, 510, 515, 516 ,519, 530, 533)


# # # # cluster.4<-c(27,  42,  43,  44,  49,  56 , 76,  89,  99, 117, 119, 133, 145 ,148, 155 ,159, 166 ,173, 176, 179, 183, 185, 189, 192, 198, 209, 213, 228, 229, 237, 245,
  # # # # 246, 249, 257, 269, 275, 280, 281, 293, 294, 299, 301, 313, 341, 358, 359, 377, 378, 380, 394, 397, 408, 411, 414, 417, 427, 428, 431, 473, 483, 488, 491,
 # # # # 496, 499, 501, 529)
 
 
 ##########################################################################################################################
############################# Clusters are obtained from the real data analysis


cluster.1<-c(1,   2,   6,   8,   9,  12,  17,  18,  19,  20,  21,  22,  46,  52,  54,  57,  60,  62,  68,  69,  77,  83,  84,  85,  86,  90,  91, 104, 105, 113, 124, 127, 128, 132, 135, 150, 151, 152, 153, 157, 160, 162, 163, 166, 169, 171, 173, 178, 183, 187, 191, 192, 193, 199, 206, 207, 216, 218, 234, 238, 241, 244, 247, 248, 252, 253, 256, 257, 258, 260, 266, 267, 273, 277, 283, 289, 293, 296, 298, 302, 303, 304, 310, 313, 321, 322, 324, 325, 332, 341, 345, 346, 35, 133, 217, 221, 271)			                                     


cluster.2<-c( 3,   7,  10,  13,  15,  16,  23,  24,  28,  32,  36,  37,  39,  45,  53,  55,  71,  73,  75,  78,  92,  93, 103, 108, 109, 110, 111, 130, 134, 137, 139, 143, 144, 146, 148, 161, 165, 172, 174, 197, 201, 202, 204, 210, 211, 228, 230, 231, 240, 242, 243, 251, 254, 255, 263, 264, 265, 269, 276, 281, 287, 294, 295, 299, 300, 315, 318, 323, 327, 329, 330, 331, 335, 336, 338, 339, 342, 343, 349, 352)



cluster.3<-c( 4,   5,  11,  14,  25,  26,  29,  30,  31,  33,  34,  38,  40,  41,  47 , 48 , 50 , 51 , 58,  59,  61,  63,  64,  65,  66,  67,  70,  72,  74,  79,  80,  81,  82,  87,  88,  94,  95,  96 , 97,  98, 100, 101, 102, 106, 107, 114, 115, 116, 118, 120 ,121, 122, 123, 125, 126, 129, 136, 140, 141, 142, 145, 147, 149, 154, 155, 156, 158, 164, 167, 168, 170, 175, 179, 180 ,181, 182, 184, 185, 186, 188 ,189, 190 ,194, 200, 203, 205, 208, 209, 213, 214, 219, 220, 222, 223, 224, 225 ,227, 233, 236, 237, 239 ,250, 259, 261, 262, 268, 270, 272, 274, 275, 278, 279, 280, 282, 284, 286, 288, 290 ,292, 297, 305, 307 ,308, 311, 312, 316, 320, 326, 328, 333, 334, 337, 340, 344, 348, 350, 351)


cluster.4<-c(27,  42,  43,  44,  49,  56 , 76,  89,  99, 112, 117, 119, 131, 138, 159 ,176 ,177 ,195 ,196 ,198 ,212, 215, 226, 229, 232, 235, 245, 246, 249, 285, 291, 301, 306, 309, 314, 317, 319, 347)


 
 

cluster<-list("1"=cluster.1, "2"=cluster.2, "3"=cluster.3, "4"=cluster.4)

 
 
 
                        
                    	
                count.r<-3
                   	
                r.neighbor.sphere<-r.array[count.r]

                theta.cluster<-list()
                
                
                for(l in 1:k_0){
                	
                	name <- paste("", l, "",sep="")
                	
                	theta.cluster.array<-NULL


				for(i in 1: length(cluster[[l]])){
				
				     
				    theta.cell<-NULL 
				    
				    index.i<-which(data_cell_motility[, 8]==id[i])
				    
				    track.length[i]<-length(index.i)
				    
				    times.present.in.study<-c(times.present.in.study, length(index.i))
				    
				    x.cell<-c(x.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 1])
				    
				    y.cell<-c(y.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 2])
				    
				    z.cell<-c(z.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 3])
				    
				    time.cell<-c(time.cell, data_cell_motility[index.i[1:(track.length[i]-1)], 7])
				    
				    data_i<-matrix(c(data_cell_motility[index.i, 7], data_cell_motility[index.i, 1], data_cell_motility[index.i, 2], data_cell_motility[index.i, 3]), nrow=length(index.i), ncol=4)
				    
				    				    
				    increment.x<-increment.y<-increment.z<-NULL
				    
				    increment.code<-matrix(0, nrow=3, ncol=1)

                    for(j in 2: length(index.i)){
	        	
	        	             increment.x[j-1]<-data_cell_motility[index.i[j], 1]-data_cell_motility[index.i[j-1], 1]
	        	
	        	             increment.y[j-1]<-data_cell_motility[index.i[j], 2]-data_cell_motility[index.i[j-1], 2]
	        	
	        	             increment.z[j-1]<-data_cell_motility[index.i[j], 3]-data_cell_motility[index.i[j-1], 3]
	        	
	                }
	                
	                
	                for(j in 1:(length(index.i)-2)){
	                	
	                	
	                	t.1<-c(increment.x[j], increment.y[j], increment.z[j])
	                	
	                	t.2<-c(increment.x[j+1], increment.y[j+1], increment.z[j+1])
	                	
	                	theta<-s.angle(t.1, t.2)
	                	
	                	
	                	theta.cell<-c(theta.cell, theta)
	                	
	                	
	                	
	                }
	                
	                
	                theta.cluster.array<-c(theta.cluster.array, theta.cell)
	                
	                
	                
	                
}

                    tmp <- list(theta.cluster.array)

                    theta.cluster[name]<-tmp


}

# xlab= expression(paste(theta)),

             
    
    setwd("/Users/e2torkas/Desktop/Clustering-independent-of-direction/Revision-Jun-11/")


    pdf("Rplot-no-cluster-second-angle.pdf", width=10, height=5, paper='special') 
 
         



     
      # # # quartz()

      # # # pWidth = 10

      # # # pHeight = 8

      # # # plot.window(c(0,pWidth),
             # # # c(0,pHeight))        
     
     # i<-2

     # m <- rbind(cbind(c(1, 1, 1), c(2, 3, 4)), c(5))

     # m <- cbind(rbind(c(1, 2), c(3, 4), c(5, 6), c(7, 8)))
     
     # m <- cbind(rbind(c(1, 2), c(3, 4)))
     
     # par(mar=c(4,7,2,1)) 
     
     # m <- cbind(c(1, 1, 1), c(2, 3, 4))
     
     # m <- rbind(c(5), cbind(c(1, 1, 1), c(2, 3, 4)))

     layout(m)
     

     par(lwd=2, mar=c(5, 5, 4, 2))

# par(mfrow=c(2, 2), lwd=2, mar=c(5, 5, 4, 2))


	
	hist(theta.cluster.array, xlim=c(0, 2*pi), freq=F, xlab= expression(paste(theta)), cex.lab=1.5, cex.axis=1.5, cex.main=2,  xaxt="n", breaks=seq(0, 2*pi, pi/5), main="", col="#999999")
	
	
	axis(side=1, at=c(0, pi, 2*pi), labels=c(0, expression(paste(pi)), expression(paste(2*pi))), las=1, cex.axis=1.5)
	
	

# # # # # # setwd("/Users/e2torkas/Desktop/Clustering-independent-of-direction/Revision-May-28/")

# # # dev.copy(jpeg,
          # # # filename="Rplot-no-cluster-second-angle.jpeg");
          
      dev.off()
      
 
	
	

   
# # # # # #       quartz()

      # # # pWidth = 10

      # # # pHeight = 8

      # # # plot.window(c(0,pWidth),
             # # # c(0,pHeight))        
             
             
    
    setwd("/Users/e2torkas/Desktop/Clustering-independent-of-direction/Revision-Jun-11/")


    pdf("Rplot-second-angle-clusters.pdf", width=10, height=10, paper='special') 
 
         
             
     
     # i<-2

     # m <- rbind(cbind(c(1, 1, 1), c(2, 3, 4)), c(5))

     # m <- cbind(rbind(c(1, 2), c(3, 4), c(5, 6), c(7, 8)))
     
     m <- cbind(rbind(c(1, 2), c(3, 4)))
     
     # par(mar=c(4,7,2,1)) 
     
     # m <- cbind(c(1, 1, 1), c(2, 3, 4))
     
     # m <- rbind(c(5), cbind(c(1, 1, 1), c(2, 3, 4)))

     layout(m)
     

     par(lwd=2, mar=c(5, 5, 4, 2))


# for(l in 1: k_0){
	
	hist(theta.cluster[[1]], xlim=c(0, 2*pi), freq=F, xlab= paste("HP/L", sep=""), cex.lab=1.5, cex.axis=1.5, cex.main=2,  xaxt="n", breaks=seq(0, 2*pi, pi/5), main="", col="#999999")
	
	
	axis(side=1, at=c(0, pi, 2*pi), labels=c(0, expression(paste(pi)), expression(paste(2*pi))), las=1, cex.axis=1.5)
	
	
	
	hist(theta.cluster[[2]], xlim=c(0, 2*pi), freq=F, xlab= paste("HP/S", sep=""), cex.lab=1.5, cex.axis=1.5, cex.main=2,  xaxt="n", breaks=seq(0, 2*pi, pi/5), main="", col="#999999")
	
	
	axis(side=1, at=c(0, pi, 2*pi), labels=c(0, expression(paste(pi)), expression(paste(2*pi))), las=1, cex.axis=1.5)



hist(theta.cluster[[3]], xlim=c(0, 2*pi), freq=F, xlab= paste("NP/L", sep=""), cex.lab=1.5, cex.axis=1.5, cex.main=2,  xaxt="n", breaks=seq(0, 2*pi, pi/5), main="", col="#999999")
	
	
	axis(side=1, at=c(0, pi, 2*pi), labels=c(0, expression(paste(pi)), expression(paste(2*pi))), las=1, cex.axis=1.5)



hist(theta.cluster[[4]], xlim=c(0, 2*pi), freq=F, xlab= paste("NP/S", sep=""), cex.lab=1.5, cex.axis=1.5, cex.main=2,  xaxt="n", breaks=seq(0, 2*pi, pi/5), main="", col="#999999")
	
	
	axis(side=1, at=c(0, pi, 2*pi), labels=c(0, expression(paste(pi)), expression(paste(2*pi))), las=1, cex.axis=1.5)

	
# # # setwd("/Users/e2torkas/Desktop/Clustering-independent-of-direction/Revision-May-28/")

# # # dev.copy(jpeg,
          # # # filename="Rplot-second-angle-clusters.jpeg");
          
      dev.off()
      

# }