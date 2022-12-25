########  to work with held-out data: in the first step we call the simulated values from the fitted models. In the next steps, we plot them against the values obtained from held-out data.
remove(list=ls())


setwd("/sbb/scratch/torkashe/Clustering-Results/measures-of-motility-2-clusters-ET-2-states-HHMM/")


############################################################################################################################################################################################################################################################################################################################################################################



library(matlib)

library(sde)

library("rstan")

library("MotilityLab")

library(ggplot2)


####  Functions to generate trajectories.

##########   Generating trajectories


new.inc=function(r,phi,theta,x1){ 
  
  fake<-sum(x1^2)
  
  if(fake ==0){return("Error, first increment has zero length")}
  if( x1[3] != 0 ){
    A =cbind(x1,as.matrix(c(1,0,0)),as.matrix(c(0,1,0)))
    A.qr=qr(A)
    Qm=qr.Q(A.qr)
    u=Qm[,2]
    v=Qm[,3]
    c1=r*cos(phi)
    r2=sqrt(r^2-c1^2)
    c2=r2*cos(theta)
    c3=r2*sin(theta)
    x2=c1*x1/sqrt(fake)+c2*u+c3*v	
  }	
  if( (x1[3]==0)& (x1[1]!=0)){
    A =cbind(x1,as.matrix(c(0,0,1)),as.matrix(c(0,1,0)))
    A.qr=qr(A)
    Qm=qr.Q(A.qr)
    u=Qm[,2]
    v=Qm[,3]
    c1=r*cos(phi)
    r2=sqrt(r^2-c1^2)
    c2=r2*cos(theta)
    c3=r2*sin(theta)
    x2=c1*x1/sqrt(fake)+c2*u+c3*v
  }
  if( (x1[3]==0)& (x1[1]==0)){
    A =cbind(x1,as.matrix(c(0,0,1)),as.matrix(c(1,0,0)))
    A.qr=qr(A)
    Qm=qr.Q(A.qr)
    u=Qm[,2]
    v=Qm[,3]
    c1=r*cos(phi)
    r2=sqrt(r^2-c1^2)
    c2=r2*cos(theta)
    c3=r2*sin(theta)
    x2=c1*x1/sqrt(fake)+c2*u+c3*v
  }
  
  x2
}


s.angle=function(x1,x2){ #this calculates the second angle going from increment x1 to x2
  
  fake1<-sum(x1^2)
  
  fake2<-sum(x2^2)
  
  if( t(x1)%*%x2==sqrt(fake1)*sqrt(fake2) | t(x1)%*%x2==-sqrt(fake1)*sqrt(fake2)){return(0)}
  
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




###################################################################################################################
###############################################################  Spherical representation and angular change 


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






##################################################################################################################################################################################################################################################################################################  Held-Out Data  #########################################################




data_cell_motility_2<-read.csv("Position_2.csv", header=TRUE)

dim(data_cell_motility_2)


##   In order to have an understanding of the number of T-cells observed in the study, we use the following comment


number.of.cells_2<-length(unique(data_cell_motility_2[, 8]))




data_cell_motility_4<-read.csv("Position_Longest.csv", header=TRUE)

dim(data_cell_motility_4)

data_cell_motility_4[, 8]<-data_cell_motility_4[, 8]+(number.of.cells_2)


#####   In order to have an understanding of the number of T-cells observed in the study, we use the following comment


number.of.cells_4<-length(unique(data_cell_motility_4[, 8]))






data_cell_motility<-rbind(data_cell_motility_2, data_cell_motility_4)

number.of.cells<-number.of.cells_2+number.of.cells_4

#########   In order to have an understanding of the number of T-cells observed in the study, we use the following comment


number.of.cells<-length(unique(data_cell_motility[, 8]))

id<-seq(1000000000, (1000000000+number.of.cells-1), 1)

track.length<-NULL

for(i in 1: number.of.cells){
  
  track.length[i]<-length(which(data_cell_motility[, 8]==id[i]))
  
}

number.of.positions.total<-sum(track.length)



######################################################
######################################################
######################################################
########   measures obtained from the real data



mean.displacement.track.real<-matrix(0, nrow= number.of.positions.total)

mean.hurst.track.real<-matrix(0, nrow= number.of.positions.total)

mean.squareDisplacement.track.real<-matrix(0, nrow= number.of.positions.total)

mean.outreachRatio.track.real<-matrix(0, nrow= number.of.positions.total)

mean.straightness.track.real<-matrix(0, nrow= number.of.positions.total)

mean.max.displacement.ratio.track.real<-matrix(0, nrow= number.of.positions.total)

mean.meanTurningAngle.track.real<-matrix(0, nrow= number.of.positions.total)

mean.asphericity.track.real<-matrix(0, nrow= number.of.positions.total)

increment.x.total<-increment.y.total<-increment.z.total<-NULL

index.measures<-index.meanTurningAngle <-NULL

theta.cell.total<-matrix(0, nrow= number.of.positions.total)



for(i in 1: number.of.cells){
  
  
  hurst.track<-NULL
  
  displacement.track<-NULL
  
  squareDisplacement.track<-NULL
  
  outreachRatio.track<-NULL
  
  straightness.track<-NULL
  
  max.displacement.ratio.track<-NULL
  
  meanTurningAngle.track<-NULL
  
  asphericity.track<-NULL
  
  
  
  
  index.i<-which(data_cell_motility[, 8]==id[i])
  
  l<-length(index.i)
  
  index.measures<-c(index.measures, index.i[2:l])
  
  index.meanTurningAngle<-c(index.meanTurningAngle, index.i[2:(l-1)])
  
  
  data_i<-matrix(c(data_cell_motility[index.i, 7], data_cell_motility[index.i, 1], data_cell_motility[index.i, 2], data_cell_motility[index.i, 3]), nrow=length(index.i), ncol=4)
  
  
  
  increment.x<-increment.y<-increment.z<-NULL
  
  increment.x[1]<-data_cell_motility[index.i[1], 1]
  
  increment.y[1]<-data_cell_motility[index.i[1], 2]
  
  increment.z[1]<-data_cell_motility[index.i[1], 3]
  
  for(j in 2: length(index.i)){
    
    increment.x[j]<-data_cell_motility[index.i[j], 1]-data_cell_motility[index.i[j-1], 1]
    
    increment.y[j]<-data_cell_motility[index.i[j], 2]-data_cell_motility[index.i[j-1], 2]
    
    increment.z[j]<-data_cell_motility[index.i[j], 3]-data_cell_motility[index.i[j-1], 3]
    
  }
  
  theta.cell<-NULL
  
  
  for(j in 1:(length(index.i)-1)){
    
    
    t.1<-c(increment.x[j], increment.y[j], increment.z[j])
    
    t.2<-c(increment.x[j+1], increment.y[j+1], increment.z[j+1])
    
    theta<-s.angle(t.1, t.2)
    
    
    theta.cell<-c(theta.cell, theta)
    
    
    
  }
  
  theta.cell.total[index.i[1:(l-1)]]<-theta.cell
  
  increment.x.total<-c(increment.x[2:(l-1)], increment.x.total)
  
  increment.y.total<-c(increment.y[2:(l-1)], increment.y.total)
  
  increment.z.total<-c(increment.z[2:(l-1)], increment.z.total)
  
  
  
  for(k in 2:l){
    
    
    if(k>2) hurst.track[k-1]<-hurstExponent(data_i[1:k, ])
    
    displacement.track[k-1]<-displacement(data_i[1:k, ])
    
    squareDisplacement.track[k-1]<-squareDisplacement(data_i[1:k, ])
    
    outreachRatio.track[k-1]<-outreachRatio(data_i[1:k, ])
    
    straightness.track[k-1]<-straightness(data_i[1:k, ])
    
    max.displacement.ratio.track[k-1]<- displacementRatio(data_i[1:k, ])
    
    if(k>2) meanTurningAngle.track[k-1]<-meanTurningAngle(data_i[1:k, ])
    
    asphericity.track[k-1]<-asphericity(data_i[1:k, ])
    
    
  }
  
  mean.hurst.track.real[index.i[3:l]]<-hurst.track[2:(l-1)]
  
  mean.displacement.track.real[index.i[2:l]]<-displacement.track
  
  mean.squareDisplacement.track.real[index.i[2:l]]<-squareDisplacement.track
  
  mean.outreachRatio.track.real[index.i[2:l]]<-outreachRatio.track
  
  mean.straightness.track.real[index.i[2:l]]<-straightness.track
  
  mean.max.displacement.ratio.track.real[index.i[2:l]]<-max.displacement.ratio.track
  
  mean.meanTurningAngle.track.real[index.i[2:(l-1)]]<-meanTurningAngle.track[2:(l-1)]
  
  mean.asphericity.track.real[index.i[2:l]]<-asphericity.track
  
  
}


time.msd<-21




separated.data.based.on.time.asphericity<-matrix(0, nrow= time.msd, ncol=number.of.cells)

measure.greg.asphericity.data <-matrix(0, nrow=time.msd)



separated.data.based.on.time.max.displacement.ratio<-matrix(0, nrow= time.msd, ncol=number.of.cells)

measure.greg.max.displacement.ratio.data <-matrix(0, nrow=time.msd)




separated.data.based.on.time.hurst<-matrix(0, nrow= time.msd, ncol=number.of.cells)

measure.greg.hurst.data <-matrix(0, nrow=time.msd)




separated.data.based.on.time.displacement<-matrix(0, nrow= time.msd, ncol=number.of.cells)

measure.greg.displacement.data <-matrix(0, nrow=time.msd)





separated.data.based.on.time.squared.displacement<-matrix(0, nrow= time.msd, ncol=number.of.cells)

measure.greg.squared.displacement.data <-matrix(0, nrow=time.msd)




separated.data.based.on.time.straightness<-matrix(0, nrow= time.msd, ncol=number.of.cells)

measure.greg.straightness.data <-matrix(0, nrow=time.msd)



separated.data.based.on.time.outreachratio<-matrix(0, nrow= time.msd, ncol=number.of.cells)

measure.greg.outreachratio.data <-matrix(0, nrow=time.msd)



true.radius<-matrix(0, nrow= time.msd, ncol=number.of.cells)


data.sd.fence.displacement.track<-matrix(0, nrow= time.msd)

data.sd.fence.hurst.track<-matrix(0, nrow= time.msd)

data.sd.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

data.sd.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

data.sd.fence.straightness.track<-matrix(0, nrow= time.msd)

data.sd.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

data.sd.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

data.sd.fence.asphericity.track<-matrix(0, nrow= time.msd)







data.lower.fence.displacement.track<-matrix(0, nrow= time.msd)

data.lower.fence.hurst.track<-matrix(0, nrow= time.msd)

data.lower.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

data.lower.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

data.lower.fence.straightness.track<-matrix(0, nrow= time.msd)

data.lower.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

data.lower.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

data.lower.fence.asphericity.track<-matrix(0, nrow= time.msd)





data.upper.fence.displacement.track<-matrix(0, nrow= time.msd)

data.upper.fence.hurst.track<-matrix(0, nrow= time.msd)

data.upper.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

data.upper.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

data.upper.fence.straightness.track<-matrix(0, nrow= time.msd)

data.upper.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

data.upper.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

data.upper.fence.asphericity.track<-matrix(0, nrow= time.msd)







time.msd<-21


counter.time.data<-matrix(0, nrow=time.msd)




for(i in 1:number.of.cells){
  
  
  
  index.i<-which(data_cell_motility[, 8]==id[i])
  
  
  
  # # # # for(i in 1:length(cluster.1[cluster.1>249])){
  
  # # t<-cluster.1[cluster.1>249]-249; 
  
  # # index.i<-which(data_cell_motility[, 8]==id[t[i]])
  
  
  
  
  l<-length(index.i)
  
  for(k in 2: 21){
    
    # t<-0
    
    
    
    
    separated.data.based.on.time.displacement[k, i]<-mean.displacement.track.real[index.i[k]]
    
    measure.greg.displacement.data[k]<-measure.greg.displacement.data[k]+ mean.displacement.track.real[index.i[k]]
    
    
    
    separated.data.based.on.time.hurst[k, i]<-mean.hurst.track.real[index.i[k]]
    
    measure.greg.hurst.data[k]<-measure.greg.hurst.data[k]+ mean.hurst.track.real[index.i[k]]
    
    
    
    
    separated.data.based.on.time.squared.displacement[k, i]<-mean.squareDisplacement.track.real[index.i[k]]
    
    measure.greg.squared.displacement.data[k] <-measure.greg.squared.displacement.data[k]+ mean.squareDisplacement.track.real[index.i[k]]
    
    
    
    
    
    
    
    
    
    
    separated.data.based.on.time.straightness[k, i]<-mean.straightness.track.real[index.i[k]]
    
    measure.greg.straightness.data[k] <-measure.greg.straightness.data[k]+mean.straightness.track.real[index.i[k]]
    
    
    
    
    
    separated.data.based.on.time.asphericity[k, i]<-mean.asphericity.track.real[index.i[k]]
    
    measure.greg.asphericity.data[k] <-measure.greg.asphericity.data[k]+mean.asphericity.track.real[index.i[k]]
    
    
    
    
    
    
    
    separated.data.based.on.time.max.displacement.ratio[k, i]<-mean.max.displacement.ratio.track.real[index.i[k]]
    
    measure.greg.max.displacement.ratio.data[k] <-measure.greg.max.displacement.ratio.data[k]+ mean.max.displacement.ratio.track.real[index.i[k]]
    
    
    
    
    
    separated.data.based.on.time.outreachratio[k, i]<-mean.outreachRatio.track.real[index.i[k]]
    
    measure.greg.outreachratio.data[k] <-measure.greg.outreachratio.data[k]+ mean.outreachRatio.track.real[index.i[k]]
    
    
    counter.time.data[k]<-counter.time.data[k]+1
    
    
  }	
  
  
  
  # # # # # # 		             for(k in 30: time.msd){
  
  # # # if(l>=k){
  
  # # # # t<-0
  
  
  
  
  # # # separated.data.based.on.time.displacement[k, i]<-mean.displacement.track.real[index.i[k]]
  
  # # # measure.greg.displacement.data[k]<-measure.greg.displacement.data[k]+ mean.displacement.track.real[index.i[k]]
  
  
  
  
  
  
  
  
  # # # separated.data.based.on.time.squared.displacement[k, i]<-mean.squareDisplacement.track.real[index.i[k]]
  
  # # # measure.greg.squared.displacement.data[k] <-measure.greg.squared.displacement.data[k]+ mean.squareDisplacement.track.real[index.i[k]]
  
  
  
  
  
  
  
  # # # separated.data.based.on.time.straightness[k, i]<-mean.straightness.track.real[index.i[k]]
  
  # # # measure.greg.straightness.data[k] <-measure.greg.straightness.data[k]+mean.straightness.track.real[index.i[k]]
  
  
  
  
  
  # # # separated.data.based.on.time.outreachratio[k, i]<-mean.outreachRatio.track.real[index.i[k]]
  
  # # # measure.greg.outreachratio.data[k] <-measure.greg.outreachratio.data[k]+ mean.outreachRatio.track.real[index.i[k]]
  
  
  # # # counter.time.data[k]<-counter.time.data[k]+1
  
  
  # # # }             
  
  # # # }		
  
  
}  



for(k in 2:time.msd){
  
  
  
  t<-which(separated.data.based.on.time.displacement[k, ]==0)
  
  if(length(t)>0) values<-separated.data.based.on.time.displacement[k, -t]
  
  if(length(t)==0) values<-separated.data.based.on.time.displacement[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  data.lower.fence.displacement.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  data.upper.fence.displacement.track[k]<-summary(values)[5]+1.5*iqr
  
  data.sd.fence.displacement.track[k]<-sd(values)
  
  
  ########################################################################################################################## 
  
  
  
  
  t<-which(separated.data.based.on.time.hurst[k, ]==0)
  
  if(length(t)>0) values<-separated.data.based.on.time.hurst[k, -t]
  
  if(length(t)==0) values<-separated.data.based.on.time.hurst[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  data.lower.fence.hurst.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  data.upper.fence.hurst.track[k]<-summary(values)[5]+1.5*iqr
  
  data.sd.fence.hurst.track[k]<-sd(values)
  
  
  ########################################################################################################################## 
  
  
  
  
  
  
  t<-which(separated.data.based.on.time.max.displacement.ratio[k, ]==0)
  
  if(length(t)>0) values<-separated.data.based.on.time.max.displacement.ratio[k, -t]
  
  if(length(t)==0) values<-separated.data.based.on.time.max.displacement.ratio[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  data.lower.fence.max.displacement.ratio.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  data.upper.fence.max.displacement.ratio.track[k]<-summary(values)[5]+1.5*iqr
  
  data.sd.fence.max.displacement.ratio.track[k]<-sd(values)
  
  
  ##########################################################################################################################   
  
  
  
  
  
  t<-which(separated.data.based.on.time.asphericity[k, ]==0)
  
  if(length(t)>0) values<-separated.data.based.on.time.asphericity[k, -t]
  
  if(length(t)==0) values<-separated.data.based.on.time.asphericity[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  data.lower.fence.asphericity.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  data.upper.fence.asphericity.track[k]<-summary(values)[5]+1.5*iqr
  
  data.sd.fence.asphericity.track[k]<-sd(values)
  
  
  ##########################################################################################################################   
  
  
  
  
  
  
  
  
  t<-which(separated.data.based.on.time.squared.displacement[k, ]==0)
  
  if(length(t)>0) values<-separated.data.based.on.time.squared.displacement[k, -t]
  
  if(length(t)==0) values<-separated.data.based.on.time.squared.displacement[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  data.lower.fence.squareDisplacement.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  data.upper.fence.squareDisplacement.track[k]<-summary(values)[5]+1.5*iqr
  
  data.sd.fence.squareDisplacement.track[k]<-sd(values)
  
  
  
  ##########################################################################################################################
  
  
  
  
  
  t<-which(separated.data.based.on.time.straightness[k, ]==0)
  
  if(length(t)>0) values<-separated.data.based.on.time.straightness[k, -t]
  
  if(length(t)==0) values<-separated.data.based.on.time.straightness[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  data.lower.fence.straightness.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  data.upper.fence.straightness.track[k]<-summary(values)[5]+1.5*iqr
  
  data.sd.fence.straightness.track[k]<-sd(values)
  
  
  
  ##########################################################################################################################
  
  
  
  
  
  
  t<-which(separated.data.based.on.time.outreachratio[k, ]==0)
  
  if(length(t)>0) values<-separated.data.based.on.time.outreachratio[k, -t]
  
  if(length(t)==0) values<-separated.data.based.on.time.outreachratio[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  data.lower.fence.outreachRatio.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  data.upper.fence.outreachRatio.track[k]<-summary(values)[5]+1.5*iqr
  
  data.sd.fence.outreachRatio.track[k]<-sd(values)
  
  
}	 



time.msd<-21




M.G.D.asphericity<-measure.greg.asphericity.data[2: time.msd]/counter.time.data[2: time.msd]



M.G.D.max.displacement.ratio<-measure.greg.max.displacement.ratio.data[2: time.msd]/counter.time.data[2: time.msd]


M.G.D.hurst<-measure.greg.hurst.data[3: time.msd]/counter.time.data[3: time.msd]



M.G.D.displacement<-measure.greg.displacement.data[2: time.msd]/counter.time.data[2: time.msd]



M.G.D.squared.displacement<-measure.greg.squared.displacement.data[2: time.msd]/counter.time.data[2: time.msd]



M.G.D.outreachratio<-measure.greg.outreachratio.data[2: time.msd]/counter.time.data[2: time.msd]



M.G.D.straightness<-measure.greg.straightness.data[2: time.msd]/counter.time.data[2: time.msd]



##############################################

setwd("/sbb/scratch/torkashe/Spherical-Coordinates-Motility-Metrics/Asphericity/")



xlab.values<-seq(1, (time.msd-1), 1)



quartz()

pWidth = 10

pHeight = 8

plot.window(c(0,pWidth),
            c(0,pHeight))        

# i<-2

# m <- rbind(cbind(c(1, 1, 1), c(2, 3, 4)), c(5))

# m <- cbind(rbind(c(1, 2), c(3, 4), c(5, 6), c(7, 8)))

# # # m <- cbind(c(1, 2))

# # # # par(mar=c(4,7,2,1)) 

# # # # m <- cbind(c(1, 1, 1), c(2, 3, 4))

# # # # m <- rbind(c(5), cbind(c(1, 1, 1), c(2, 3, 4)))

# # # layout(m)


# par(lwd=2, mar=c(5, 5, 4, 2))




par(lwd=2.5, cex.axis=2, mar=c(5, 5, 4, 2) + 0.1)






ylab<-c(10, 100, 1000)


# ylab.values<-seq(0.6, 0.8, .05)


ylab.values<-c(0.60, 0.65, 0.70, 0.75)


x.lab<-seq(2, (time.msd-1), 1)


xlab.values<-log(x.lab)/log(10)


plot(x.lab, M.G.D.asphericity, col="black", xaxt = "n", xlab="Step", ylab="Asphericity", main="",   ylim=c(min(M.G.D.asphericity+ data.sd.fence.asphericity.track[3: time.msd]/sqrt(counter.time.data[3: time.msd]))-.1, max(M.G.D.asphericity+ data.sd.fence.asphericity.track[3: time.msd]/sqrt(counter.time.data[3: time.msd]))+.1),  pch=24, cex.lab=2, cex.main=2.5, type="l", cex=1.5, lwd=1.5)

axis(1, at= x.lab, labels= x.lab)


arrows(x.lab, M.G.D.asphericity-data.sd.fence.asphericity.track[3: time.msd]/sqrt(counter.time.data[3: time.msd]), x.lab, M.G.D.asphericity+ data.sd.fence.asphericity.track[3: time.msd]/sqrt(counter.time.data[3: time.msd]), length=0.05, angle=90, code=3, col="black")

axis(4, at= ylab.values, labels= ylab.values)




###################################################################################################################################################################################################################################################


remove(list=ls())		

##########################################################################################################################    

setwd("/sbb/scratch/torkashe/Clustering-Results/measures-of-motility-4-clusters-ET-2-states-NHMM")

load(".RMData")




simulated.mean.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)



simulated.lower.fence.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)




simulated.upper.fence.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)





simulated.mean.hurst.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.asphericity.track.observation<-matrix(0, nrow= number.of.positions.total)




simulated.lower.fence.hurst.track.observation<-matrix(0, nrow= number.of.positions.total)  	   

simulated.lower.fence.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.asphericity.observation<-matrix(0, nrow= number.of.positions.total)




simulated.upper.fence.hurst.track.observation<-matrix(0, nrow= number.of.positions.total) 	   

simulated.upper.fence.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.asphericity.observation<-matrix(0, nrow= number.of.positions.total)









for(i in 1: number.of.positions.total){
  
  
  
  
  simulated.mean.hurst.track.inside[i]<-mean(mean.hurst.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.hurst.track.inside[1:(r-1), i])[5]-summary(mean.hurst.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.hurst.track.inside[i]<-summary(mean.hurst.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.hurst.track.inside[i]<-summary(mean.hurst.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.displacement.track.inside[i]<-mean(mean.displacement.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.displacement.track.inside[1:(r-1), i])[5]-summary(mean.displacement.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.displacement.track.inside[i]<-summary(mean.displacement.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.displacement.track.inside[i]<-summary(mean.displacement.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.squareDisplacement.track.inside[i]<-mean(mean.squareDisplacement.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[5]-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.squareDisplacement.track.inside[i]<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.squareDisplacement.track.inside[i]<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  simulated.mean.outreachRatio.track.inside[i]<-mean(mean.outreachRatio.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[5]-summary(mean.outreachRatio.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.outreachRatio.track.inside[i]<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.outreachRatio.track.inside[i]<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.straightness.track.inside[i]<-mean(mean.straightness.track.inside[1:(r-1), i])
  
  
  
  
  
  iqr<-summary(mean.straightness.track.inside[1:(r-1), i])[5]-summary(mean.straightness.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.straightness.track.inside[i]<-summary(mean.straightness.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.straightness.track.inside[i]<-summary(mean.straightness.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  simulated.mean.max.displacement.ratio.track.inside[i]<-mean(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])
  
  
  
  
  iqr<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[5]-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.max.displacement.ratio.track.inside[i]<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.max.displacement.ratio.track.inside[i]<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  
  simulated.mean.meanTurningAngle.track.inside[i]<-mean(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[5]-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.meanTurningAngle.track.inside[i]<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.meanTurningAngle.track.inside[i]<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  simulated.mean.asphericity.track.inside[i]<-mean(mean.simulated.asphericity.inside[1:(r-1), i])
  
  
  
  
  iqr<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[5]-summary(mean.simulated.asphericity.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.asphericity.track.inside[i]<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.asphericity.track.inside[i]<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  ##############################################
  ##############################################    
  
  
  simulated.mean.hurst.track.observation[i]<-mean(mean.hurst.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.hurst.track.observation[1:(r-1), i])[5]-summary(mean.hurst.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.hurst.track.observation[i]<-summary(mean.hurst.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.hurst.track.observation[i]<-summary(mean.hurst.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  simulated.mean.displacement.track.observation[i]<-mean(mean.displacement.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.displacement.track.observation[1:(r-1), i])[5]-summary(mean.displacement.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.displacement.track.observation[i]<-summary(mean.displacement.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.displacement.track.observation[i]<-summary(mean.displacement.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.squareDisplacement.track.observation[i]<-mean(mean.squareDisplacement.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[5]-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.squareDisplacement.track.observation[i]<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.squareDisplacement.track.observation[i]<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  simulated.mean.outreachRatio.track.observation[i]<-mean(mean.outreachRatio.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[5]-summary(mean.outreachRatio.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.outreachRatio.track.observation[i]<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.outreachRatio.track.observation[i]<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.straightness.track.observation[i]<-mean(mean.straightness.track.observation[1:(r-1), i])
  
  
  
  
  
  iqr<-summary(mean.straightness.track.observation[1:(r-1), i])[5]-summary(mean.straightness.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.straightness.track.observation[i]<-summary(mean.straightness.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.straightness.track.observation[i]<-summary(mean.straightness.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  
  
  
  
  
}	   

###################################################################################################################################################################################################################################################




time.msd<-21


counter.time.simulated<-matrix(0, nrow=time.msd)


measure.greg.hurst.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.hurst<-matrix(0, nrow= time.msd, ncol=number.of.cells)



measure.greg.displacement.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.displacement<-matrix(0, nrow= time.msd, ncol=number.of.cells)


measure.greg.max.displacement.ratio.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.max.displacement.ratio<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.squared.displacement.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.squared.displacement<-matrix(0, nrow= time.msd, ncol=number.of.cells)



measure.greg.asphericity.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.asphericity<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.straightness.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.straightness<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.outreachratio.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.outreachratio<-matrix(0, nrow= time.msd, ncol=number.of.cells)



true.radius<-matrix(0, nrow= time.msd, ncol=number.of.cells)






simulated.lower.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.asphericity.track<-matrix(0, nrow= time.msd)




simulated.upper.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.asphericity.track<-matrix(0, nrow= time.msd)



simulated.sd.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.asphericity.track<-matrix(0, nrow= time.msd)





for(i in 1:number.of.cells){
  
  
  
  index.i<-which(data_cell_motility[, 8]==id[i])
  
  
  
  # # # # for(i in 1:length(cluster.1[cluster.1>249])){
  
  # # t<-cluster.1[cluster.1>249]-249; 
  
  # # index.i<-which(data_cell_motility[, 8]==id[t[i]])
  
  
  
  
  l<-length(index.i)
  
  for(k in 2: time.msd){
    
    
    
    measure.greg.asphericity.simulated[k]<-measure.greg.asphericity.simulated[k]+ simulated.mean.asphericity.track.inside[index.i[k]]
    
    separated.sim.based.on.time.asphericity[k, i]<-simulated.mean.asphericity.track.inside[index.i[k]]
    
    
    
    
    # t<-0
    
    
    measure.greg.displacement.simulated[k]<-measure.greg.displacement.simulated[k]+ simulated.mean.displacement.track.inside[index.i[k]]
    
    separated.sim.based.on.time.displacement[k, i]<-simulated.mean.displacement.track.inside[index.i[k]]
    
    
    
    measure.greg.hurst.simulated[k]<-measure.greg.hurst.simulated[k]+ simulated.mean.hurst.track.inside[index.i[k]]
    
    separated.sim.based.on.time.hurst[k, i]<-simulated.mean.hurst.track.inside[index.i[k]]
    
    
    
    
    measure.greg.max.displacement.ratio.simulated[k]<-measure.greg.max.displacement.ratio.simulated[k]+ simulated.mean.max.displacement.ratio.track.inside[index.i[k]]
    
    separated.sim.based.on.time.max.displacement.ratio[k, i]<-simulated.mean.max.displacement.ratio.track.inside[index.i[k]]
    
    
    
    
    
    measure.greg.squared.displacement.simulated[k] <-measure.greg.squared.displacement.simulated[k]+ simulated.mean.squareDisplacement.track.inside[index.i[k]]
    
    separated.sim.based.on.time.squared.displacement[k, i]<-simulated.mean.squareDisplacement.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    measure.greg.straightness.simulated[k] <-measure.greg.straightness.simulated[k] +simulated.mean.straightness.track.inside[index.i[k]]
    
    separated.sim.based.on.time.straightness[k, i]<-simulated.mean.straightness.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    
    measure.greg.outreachratio.simulated[k] <-measure.greg.outreachratio.simulated[k]+simulated.mean.outreachRatio.track.inside[index.i[k]]
    
    separated.sim.based.on.time.outreachratio[k, i]<-simulated.mean.outreachRatio.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    counter.time.simulated[k]<-counter.time.simulated[k]+1
    
    
  }	
  
  
  
  # # # # # # 		            for(k in 30: time.msd){
  
  # # # if(l>=k){
  
  # # # # t<-0
  
  
  # # # measure.greg.displacement.simulated[k]<-measure.greg.displacement.simulated[k]+ simulated.mean.displacement.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.displacement[k, i]<-simulated.mean.displacement.track.inside[index.i[k]]
  
  
  
  
  # # # measure.greg.hurst.simulated[k]<-measure.greg.hurst.simulated[k]+ simulated.mean.hurst.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.hurst[k, i]<-simulated.mean.hurst.track.inside[index.i[k]]
  
  
  
  
  
  # # # measure.greg.squared.displacement.simulated[k] <-measure.greg.squared.displacement.simulated[k]+ simulated.mean.squareDisplacement.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.squared.displacement[k, i]<-simulated.mean.squareDisplacement.track.inside[index.i[k]]
  
  
  
  
  
  
  # # # measure.greg.straightness.simulated[k] <-measure.greg.straightness.simulated[k] +simulated.mean.straightness.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.straightness[k, i]<-simulated.mean.straightness.track.inside[index.i[k]]
  
  
  
  
  
  
  
  # # # measure.greg.outreachratio.simulated[k] <-measure.greg.outreachratio.simulated[k]+simulated.mean.outreachRatio.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.outreachratio[k, i]<-simulated.mean.outreachRatio.track.inside[index.i[k]]
  
  
  # # # counter.time.simulated[k]<-counter.time.simulated[k]+1
  
  
  # # # }             
  
  # # # }		
  
  
}  




for(k in 2:time.msd){
  
  
  t<-which(separated.sim.based.on.time.asphericity[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.asphericity[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.asphericity[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.asphericity.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.asphericity.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.asphericity.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  t<-which(separated.sim.based.on.time.max.displacement.ratio[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.max.displacement.ratio[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.max.displacement.ratio[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.max.displacement.ratio.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.max.displacement.ratio.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.max.displacement.ratio.track[k]<-sd(values)
  
  
  
  ##########################################################################################################################
  
  
  
  
  
  
  t<-which(separated.sim.based.on.time.hurst[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.hurst[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.hurst[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.hurst.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.hurst.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.hurst.track[k]<-sd(values)
  
  
  
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  
  
  t<-which(separated.sim.based.on.time.displacement[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.displacement[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.displacement[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.displacement.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.displacement.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.displacement.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  t<-which(separated.sim.based.on.time.squared.displacement[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.squared.displacement[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.squared.displacement[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.squareDisplacement.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.squareDisplacement.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.squareDisplacement.track[k]<-sd(values)
  
  
  
  
  ##########################################################################################################################
  
  
  
  t<-which(separated.sim.based.on.time.straightness[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.straightness[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.straightness[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.straightness.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.straightness.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.straightness.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  t<-which(separated.sim.based.on.time.outreachratio[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.outreachratio[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.outreachratio[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.outreachRatio.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.outreachRatio.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.outreachRatio.track[k]<-sd(values)
  
  
  
  
  
  
  
}	  	   







time.msd<-21



M.G.S.hurst<-measure.greg.hurst.simulated[3: time.msd]/counter.time.simulated[3: time.msd]




M.G.S.diplacement<-measure.greg.displacement.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.max.diplacement.ratio<-measure.greg.max.displacement.ratio.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.asphericity<-measure.greg.asphericity.simulated[2: time.msd]/counter.time.simulated[2: time.msd]




M.G.S.squared.diplacement<-measure.greg.squared.displacement.simulated[2: time.msd]/counter.time.simulated[2: time.msd]




M.G.S.outreachratio<-measure.greg.outreachratio.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.straightness<-measure.greg.straightness.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



##############################################

setwd("/sbb/scratch/torkashe/Spherical-Coordinates-Motility-Metrics/Asphericity/")



ylab<-c(10, 100, 1000)

ylab.values<-log(ylab)/log(10)

x.lab<-seq(2, (time.msd-1), 1)

xlab.values<-log(x.lab)/log(10)

# points(M.G.S.max.diplacement.ratio, col="violetred", pch=20, type="p")

points(x.lab, M.G.S.asphericity, col="darkred", pch=20, type="l", cex=1.5, lwd=1.5)

# points(xlab , M.G.S.asphericity, col="violetred", pch=20, type="p")



# arrows(seq(1, (time.msd-1), 1), M.G.S.diplacement-simulated.sd.fence.displacement.track[2: time.msd], seq(1, (time.msd-1), 1), M.G.S.diplacement+simulated.sd.fence.displacement.track[2: time.msd], length=0.05, angle=90, code=3, col="darkgreen")



# # # # u<-ecdf(M.G.S.diplacement)

# # plot(u, add=T, col="darkred") 




###################################################################################################################################################################################################################################################


remove(list=ls())		

##########################################################################################################################    

setwd("/sbb/scratch/torkashe/Clustering-Results/measures-of-motility-2-clusters-ET-2-states-NHMM")

load(".RMData")




simulated.mean.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)



simulated.lower.fence.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)




simulated.upper.fence.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)





simulated.mean.hurst.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.asphericity.track.observation<-matrix(0, nrow= number.of.positions.total)




simulated.lower.fence.hurst.track.observation<-matrix(0, nrow= number.of.positions.total)  	   

simulated.lower.fence.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.asphericity.observation<-matrix(0, nrow= number.of.positions.total)




simulated.upper.fence.hurst.track.observation<-matrix(0, nrow= number.of.positions.total) 	   

simulated.upper.fence.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.asphericity.observation<-matrix(0, nrow= number.of.positions.total)









for(i in 1: number.of.positions.total){
  
  
  
  
  simulated.mean.hurst.track.inside[i]<-mean(mean.hurst.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.hurst.track.inside[1:(r-1), i])[5]-summary(mean.hurst.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.hurst.track.inside[i]<-summary(mean.hurst.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.hurst.track.inside[i]<-summary(mean.hurst.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.displacement.track.inside[i]<-mean(mean.displacement.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.displacement.track.inside[1:(r-1), i])[5]-summary(mean.displacement.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.displacement.track.inside[i]<-summary(mean.displacement.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.displacement.track.inside[i]<-summary(mean.displacement.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.squareDisplacement.track.inside[i]<-mean(mean.squareDisplacement.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[5]-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.squareDisplacement.track.inside[i]<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.squareDisplacement.track.inside[i]<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  simulated.mean.outreachRatio.track.inside[i]<-mean(mean.outreachRatio.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[5]-summary(mean.outreachRatio.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.outreachRatio.track.inside[i]<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.outreachRatio.track.inside[i]<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.straightness.track.inside[i]<-mean(mean.straightness.track.inside[1:(r-1), i])
  
  
  
  
  
  iqr<-summary(mean.straightness.track.inside[1:(r-1), i])[5]-summary(mean.straightness.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.straightness.track.inside[i]<-summary(mean.straightness.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.straightness.track.inside[i]<-summary(mean.straightness.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  simulated.mean.max.displacement.ratio.track.inside[i]<-mean(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])
  
  
  
  
  iqr<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[5]-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.max.displacement.ratio.track.inside[i]<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.max.displacement.ratio.track.inside[i]<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  
  simulated.mean.meanTurningAngle.track.inside[i]<-mean(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[5]-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.meanTurningAngle.track.inside[i]<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.meanTurningAngle.track.inside[i]<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  simulated.mean.asphericity.track.inside[i]<-mean(mean.simulated.asphericity.inside[1:(r-1), i])
  
  
  
  
  iqr<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[5]-summary(mean.simulated.asphericity.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.asphericity.track.inside[i]<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.asphericity.track.inside[i]<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  ##############################################
  ##############################################    
  
  
  simulated.mean.hurst.track.observation[i]<-mean(mean.hurst.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.hurst.track.observation[1:(r-1), i])[5]-summary(mean.hurst.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.hurst.track.observation[i]<-summary(mean.hurst.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.hurst.track.observation[i]<-summary(mean.hurst.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  simulated.mean.displacement.track.observation[i]<-mean(mean.displacement.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.displacement.track.observation[1:(r-1), i])[5]-summary(mean.displacement.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.displacement.track.observation[i]<-summary(mean.displacement.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.displacement.track.observation[i]<-summary(mean.displacement.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.squareDisplacement.track.observation[i]<-mean(mean.squareDisplacement.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[5]-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.squareDisplacement.track.observation[i]<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.squareDisplacement.track.observation[i]<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  simulated.mean.outreachRatio.track.observation[i]<-mean(mean.outreachRatio.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[5]-summary(mean.outreachRatio.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.outreachRatio.track.observation[i]<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.outreachRatio.track.observation[i]<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.straightness.track.observation[i]<-mean(mean.straightness.track.observation[1:(r-1), i])
  
  
  
  
  
  iqr<-summary(mean.straightness.track.observation[1:(r-1), i])[5]-summary(mean.straightness.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.straightness.track.observation[i]<-summary(mean.straightness.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.straightness.track.observation[i]<-summary(mean.straightness.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  
  
  
  
  
}	   

###################################################################################################################################################################################################################################################




time.msd<-21


counter.time.simulated<-matrix(0, nrow=time.msd)


measure.greg.hurst.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.hurst<-matrix(0, nrow= time.msd, ncol=number.of.cells)



measure.greg.displacement.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.displacement<-matrix(0, nrow= time.msd, ncol=number.of.cells)


measure.greg.max.displacement.ratio.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.max.displacement.ratio<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.squared.displacement.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.squared.displacement<-matrix(0, nrow= time.msd, ncol=number.of.cells)



measure.greg.asphericity.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.asphericity<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.straightness.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.straightness<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.outreachratio.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.outreachratio<-matrix(0, nrow= time.msd, ncol=number.of.cells)



true.radius<-matrix(0, nrow= time.msd, ncol=number.of.cells)






simulated.lower.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.asphericity.track<-matrix(0, nrow= time.msd)




simulated.upper.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.asphericity.track<-matrix(0, nrow= time.msd)



simulated.sd.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.asphericity.track<-matrix(0, nrow= time.msd)





for(i in 1:number.of.cells){
  
  
  
  index.i<-which(data_cell_motility[, 8]==id[i])
  
  
  
  # # # # for(i in 1:length(cluster.1[cluster.1>249])){
  
  # # t<-cluster.1[cluster.1>249]-249; 
  
  # # index.i<-which(data_cell_motility[, 8]==id[t[i]])
  
  
  
  
  l<-length(index.i)
  
  for(k in 2: time.msd){
    
    
    
    measure.greg.asphericity.simulated[k]<-measure.greg.asphericity.simulated[k]+ simulated.mean.asphericity.track.inside[index.i[k]]
    
    separated.sim.based.on.time.asphericity[k, i]<-simulated.mean.asphericity.track.inside[index.i[k]]
    
    
    
    
    # t<-0
    
    
    measure.greg.displacement.simulated[k]<-measure.greg.displacement.simulated[k]+ simulated.mean.displacement.track.inside[index.i[k]]
    
    separated.sim.based.on.time.displacement[k, i]<-simulated.mean.displacement.track.inside[index.i[k]]
    
    
    
    measure.greg.hurst.simulated[k]<-measure.greg.hurst.simulated[k]+ simulated.mean.hurst.track.inside[index.i[k]]
    
    separated.sim.based.on.time.hurst[k, i]<-simulated.mean.hurst.track.inside[index.i[k]]
    
    
    
    
    measure.greg.max.displacement.ratio.simulated[k]<-measure.greg.max.displacement.ratio.simulated[k]+ simulated.mean.max.displacement.ratio.track.inside[index.i[k]]
    
    separated.sim.based.on.time.max.displacement.ratio[k, i]<-simulated.mean.max.displacement.ratio.track.inside[index.i[k]]
    
    
    
    
    
    measure.greg.squared.displacement.simulated[k] <-measure.greg.squared.displacement.simulated[k]+ simulated.mean.squareDisplacement.track.inside[index.i[k]]
    
    separated.sim.based.on.time.squared.displacement[k, i]<-simulated.mean.squareDisplacement.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    measure.greg.straightness.simulated[k] <-measure.greg.straightness.simulated[k] +simulated.mean.straightness.track.inside[index.i[k]]
    
    separated.sim.based.on.time.straightness[k, i]<-simulated.mean.straightness.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    
    measure.greg.outreachratio.simulated[k] <-measure.greg.outreachratio.simulated[k]+simulated.mean.outreachRatio.track.inside[index.i[k]]
    
    separated.sim.based.on.time.outreachratio[k, i]<-simulated.mean.outreachRatio.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    counter.time.simulated[k]<-counter.time.simulated[k]+1
    
    
  }	
  
  
  
  # # # # # # 		            for(k in 30: time.msd){
  
  # # # if(l>=k){
  
  # # # # t<-0
  
  
  # # # measure.greg.displacement.simulated[k]<-measure.greg.displacement.simulated[k]+ simulated.mean.displacement.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.displacement[k, i]<-simulated.mean.displacement.track.inside[index.i[k]]
  
  
  
  
  # # # measure.greg.hurst.simulated[k]<-measure.greg.hurst.simulated[k]+ simulated.mean.hurst.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.hurst[k, i]<-simulated.mean.hurst.track.inside[index.i[k]]
  
  
  
  
  
  # # # measure.greg.squared.displacement.simulated[k] <-measure.greg.squared.displacement.simulated[k]+ simulated.mean.squareDisplacement.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.squared.displacement[k, i]<-simulated.mean.squareDisplacement.track.inside[index.i[k]]
  
  
  
  
  
  
  # # # measure.greg.straightness.simulated[k] <-measure.greg.straightness.simulated[k] +simulated.mean.straightness.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.straightness[k, i]<-simulated.mean.straightness.track.inside[index.i[k]]
  
  
  
  
  
  
  
  # # # measure.greg.outreachratio.simulated[k] <-measure.greg.outreachratio.simulated[k]+simulated.mean.outreachRatio.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.outreachratio[k, i]<-simulated.mean.outreachRatio.track.inside[index.i[k]]
  
  
  # # # counter.time.simulated[k]<-counter.time.simulated[k]+1
  
  
  # # # }             
  
  # # # }		
  
  
}  




for(k in 2:time.msd){
  
  
  t<-which(separated.sim.based.on.time.asphericity[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.asphericity[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.asphericity[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.asphericity.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.asphericity.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.asphericity.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  t<-which(separated.sim.based.on.time.max.displacement.ratio[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.max.displacement.ratio[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.max.displacement.ratio[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.max.displacement.ratio.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.max.displacement.ratio.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.max.displacement.ratio.track[k]<-sd(values)
  
  
  
  ##########################################################################################################################
  
  
  
  
  
  
  t<-which(separated.sim.based.on.time.hurst[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.hurst[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.hurst[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.hurst.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.hurst.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.hurst.track[k]<-sd(values)
  
  
  
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  
  
  t<-which(separated.sim.based.on.time.displacement[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.displacement[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.displacement[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.displacement.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.displacement.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.displacement.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  t<-which(separated.sim.based.on.time.squared.displacement[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.squared.displacement[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.squared.displacement[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.squareDisplacement.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.squareDisplacement.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.squareDisplacement.track[k]<-sd(values)
  
  
  
  
  ##########################################################################################################################
  
  
  
  t<-which(separated.sim.based.on.time.straightness[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.straightness[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.straightness[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.straightness.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.straightness.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.straightness.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  t<-which(separated.sim.based.on.time.outreachratio[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.outreachratio[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.outreachratio[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.outreachRatio.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.outreachRatio.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.outreachRatio.track[k]<-sd(values)
  
  
  
  
  
  
  
}	  	   








time.msd<-21



M.G.S.hurst<-measure.greg.hurst.simulated[3: time.msd]/counter.time.simulated[3: time.msd]




M.G.S.diplacement<-measure.greg.displacement.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.max.diplacement.ratio<-measure.greg.max.displacement.ratio.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.asphericity<-measure.greg.asphericity.simulated[2: time.msd]/counter.time.simulated[2: time.msd]




M.G.S.squared.diplacement<-measure.greg.squared.displacement.simulated[2: time.msd]/counter.time.simulated[2: time.msd]




M.G.S.outreachratio<-measure.greg.outreachratio.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.straightness<-measure.greg.straightness.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



##############################################


setwd("/sbb/scratch/torkashe/Spherical-Coordinates-Motility-Metrics/Asphericity/")


ylab<-c(10, 100, 1000)

ylab.values<-log(ylab)/log(10)

x.lab<-seq(2, (time.msd-1), 1)

xlab.values<-log(x.lab)/log(10)

# points(M.G.S.max.diplacement.ratio, col="violetred", pch=20, type="p")

points(x.lab, M.G.S.asphericity, col="darkred", pch=11, type="l", cex=1.5, lwd=1.5, lty=3)

# points(xlab , M.G.S.asphericity, col="violetred", pch=20, type="p")



# arrows(seq(1, (time.msd-1), 1), M.G.S.diplacement-simulated.sd.fence.displacement.track[2: time.msd], seq(1, (time.msd-1), 1), M.G.S.diplacement+simulated.sd.fence.displacement.track[2: time.msd], length=0.05, angle=90, code=3, col="darkgreen")



# # # # u<-ecdf(M.G.S.diplacement)

# # plot(u, add=T, col="darkred") 


###################################################################################################################################################################################################################################################


remove(list=ls())		

##########################################################################################################################    

setwd("/sbb/scratch/torkashe/Clustering-Results/measures-of-motility-3-clusters-mokh-2-states-NHMM/")

load(".RMData")



simulated.mean.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)



simulated.lower.fence.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)




simulated.upper.fence.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)





simulated.mean.hurst.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.asphericity.track.observation<-matrix(0, nrow= number.of.positions.total)




simulated.lower.fence.hurst.track.observation<-matrix(0, nrow= number.of.positions.total)  	   

simulated.lower.fence.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.asphericity.observation<-matrix(0, nrow= number.of.positions.total)




simulated.upper.fence.hurst.track.observation<-matrix(0, nrow= number.of.positions.total) 	   

simulated.upper.fence.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.asphericity.observation<-matrix(0, nrow= number.of.positions.total)









for(i in 1: number.of.positions.total){
  
  
  
  
  simulated.mean.hurst.track.inside[i]<-mean(mean.hurst.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.hurst.track.inside[1:(r-1), i])[5]-summary(mean.hurst.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.hurst.track.inside[i]<-summary(mean.hurst.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.hurst.track.inside[i]<-summary(mean.hurst.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.displacement.track.inside[i]<-mean(mean.displacement.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.displacement.track.inside[1:(r-1), i])[5]-summary(mean.displacement.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.displacement.track.inside[i]<-summary(mean.displacement.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.displacement.track.inside[i]<-summary(mean.displacement.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.squareDisplacement.track.inside[i]<-mean(mean.squareDisplacement.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[5]-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.squareDisplacement.track.inside[i]<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.squareDisplacement.track.inside[i]<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  simulated.mean.outreachRatio.track.inside[i]<-mean(mean.outreachRatio.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[5]-summary(mean.outreachRatio.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.outreachRatio.track.inside[i]<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.outreachRatio.track.inside[i]<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.straightness.track.inside[i]<-mean(mean.straightness.track.inside[1:(r-1), i])
  
  
  
  
  
  iqr<-summary(mean.straightness.track.inside[1:(r-1), i])[5]-summary(mean.straightness.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.straightness.track.inside[i]<-summary(mean.straightness.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.straightness.track.inside[i]<-summary(mean.straightness.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  simulated.mean.max.displacement.ratio.track.inside[i]<-mean(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])
  
  
  
  
  iqr<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[5]-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.max.displacement.ratio.track.inside[i]<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.max.displacement.ratio.track.inside[i]<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  
  simulated.mean.meanTurningAngle.track.inside[i]<-mean(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[5]-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.meanTurningAngle.track.inside[i]<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.meanTurningAngle.track.inside[i]<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  simulated.mean.asphericity.track.inside[i]<-mean(mean.simulated.asphericity.inside[1:(r-1), i])
  
  
  
  
  iqr<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[5]-summary(mean.simulated.asphericity.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.asphericity.track.inside[i]<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.asphericity.track.inside[i]<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  ##############################################
  ##############################################    
  
  
  simulated.mean.hurst.track.observation[i]<-mean(mean.hurst.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.hurst.track.observation[1:(r-1), i])[5]-summary(mean.hurst.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.hurst.track.observation[i]<-summary(mean.hurst.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.hurst.track.observation[i]<-summary(mean.hurst.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  simulated.mean.displacement.track.observation[i]<-mean(mean.displacement.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.displacement.track.observation[1:(r-1), i])[5]-summary(mean.displacement.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.displacement.track.observation[i]<-summary(mean.displacement.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.displacement.track.observation[i]<-summary(mean.displacement.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.squareDisplacement.track.observation[i]<-mean(mean.squareDisplacement.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[5]-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.squareDisplacement.track.observation[i]<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.squareDisplacement.track.observation[i]<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  simulated.mean.outreachRatio.track.observation[i]<-mean(mean.outreachRatio.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[5]-summary(mean.outreachRatio.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.outreachRatio.track.observation[i]<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.outreachRatio.track.observation[i]<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.straightness.track.observation[i]<-mean(mean.straightness.track.observation[1:(r-1), i])
  
  
  
  
  
  iqr<-summary(mean.straightness.track.observation[1:(r-1), i])[5]-summary(mean.straightness.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.straightness.track.observation[i]<-summary(mean.straightness.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.straightness.track.observation[i]<-summary(mean.straightness.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  
  
  
  
  
}	   

###################################################################################################################################################################################################################################################




time.msd<-21


counter.time.simulated<-matrix(0, nrow=time.msd)


measure.greg.hurst.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.hurst<-matrix(0, nrow= time.msd, ncol=number.of.cells)



measure.greg.displacement.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.displacement<-matrix(0, nrow= time.msd, ncol=number.of.cells)


measure.greg.max.displacement.ratio.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.max.displacement.ratio<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.squared.displacement.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.squared.displacement<-matrix(0, nrow= time.msd, ncol=number.of.cells)



measure.greg.asphericity.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.asphericity<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.straightness.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.straightness<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.outreachratio.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.outreachratio<-matrix(0, nrow= time.msd, ncol=number.of.cells)



true.radius<-matrix(0, nrow= time.msd, ncol=number.of.cells)






simulated.lower.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.asphericity.track<-matrix(0, nrow= time.msd)




simulated.upper.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.asphericity.track<-matrix(0, nrow= time.msd)



simulated.sd.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.asphericity.track<-matrix(0, nrow= time.msd)





for(i in 1:number.of.cells){
  
  
  
  index.i<-which(data_cell_motility[, 8]==id[i])
  
  
  
  # # # # for(i in 1:length(cluster.1[cluster.1>249])){
  
  # # t<-cluster.1[cluster.1>249]-249; 
  
  # # index.i<-which(data_cell_motility[, 8]==id[t[i]])
  
  
  
  
  l<-length(index.i)
  
  for(k in 2: time.msd){
    
    
    
    measure.greg.asphericity.simulated[k]<-measure.greg.asphericity.simulated[k]+ simulated.mean.asphericity.track.inside[index.i[k]]
    
    separated.sim.based.on.time.asphericity[k, i]<-simulated.mean.asphericity.track.inside[index.i[k]]
    
    
    
    
    # t<-0
    
    
    measure.greg.displacement.simulated[k]<-measure.greg.displacement.simulated[k]+ simulated.mean.displacement.track.inside[index.i[k]]
    
    separated.sim.based.on.time.displacement[k, i]<-simulated.mean.displacement.track.inside[index.i[k]]
    
    
    
    measure.greg.hurst.simulated[k]<-measure.greg.hurst.simulated[k]+ simulated.mean.hurst.track.inside[index.i[k]]
    
    separated.sim.based.on.time.hurst[k, i]<-simulated.mean.hurst.track.inside[index.i[k]]
    
    
    
    
    measure.greg.max.displacement.ratio.simulated[k]<-measure.greg.max.displacement.ratio.simulated[k]+ simulated.mean.max.displacement.ratio.track.inside[index.i[k]]
    
    separated.sim.based.on.time.max.displacement.ratio[k, i]<-simulated.mean.max.displacement.ratio.track.inside[index.i[k]]
    
    
    
    
    
    measure.greg.squared.displacement.simulated[k] <-measure.greg.squared.displacement.simulated[k]+ simulated.mean.squareDisplacement.track.inside[index.i[k]]
    
    separated.sim.based.on.time.squared.displacement[k, i]<-simulated.mean.squareDisplacement.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    measure.greg.straightness.simulated[k] <-measure.greg.straightness.simulated[k] +simulated.mean.straightness.track.inside[index.i[k]]
    
    separated.sim.based.on.time.straightness[k, i]<-simulated.mean.straightness.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    
    measure.greg.outreachratio.simulated[k] <-measure.greg.outreachratio.simulated[k]+simulated.mean.outreachRatio.track.inside[index.i[k]]
    
    separated.sim.based.on.time.outreachratio[k, i]<-simulated.mean.outreachRatio.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    counter.time.simulated[k]<-counter.time.simulated[k]+1
    
    
  }	
  
  
  
  # # # # # # 		            for(k in 30: time.msd){
  
  # # # if(l>=k){
  
  # # # # t<-0
  
  
  # # # measure.greg.displacement.simulated[k]<-measure.greg.displacement.simulated[k]+ simulated.mean.displacement.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.displacement[k, i]<-simulated.mean.displacement.track.inside[index.i[k]]
  
  
  
  
  # # # measure.greg.hurst.simulated[k]<-measure.greg.hurst.simulated[k]+ simulated.mean.hurst.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.hurst[k, i]<-simulated.mean.hurst.track.inside[index.i[k]]
  
  
  
  
  
  # # # measure.greg.squared.displacement.simulated[k] <-measure.greg.squared.displacement.simulated[k]+ simulated.mean.squareDisplacement.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.squared.displacement[k, i]<-simulated.mean.squareDisplacement.track.inside[index.i[k]]
  
  
  
  
  
  
  # # # measure.greg.straightness.simulated[k] <-measure.greg.straightness.simulated[k] +simulated.mean.straightness.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.straightness[k, i]<-simulated.mean.straightness.track.inside[index.i[k]]
  
  
  
  
  
  
  
  # # # measure.greg.outreachratio.simulated[k] <-measure.greg.outreachratio.simulated[k]+simulated.mean.outreachRatio.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.outreachratio[k, i]<-simulated.mean.outreachRatio.track.inside[index.i[k]]
  
  
  # # # counter.time.simulated[k]<-counter.time.simulated[k]+1
  
  
  # # # }             
  
  # # # }		
  
  
}  




for(k in 2:time.msd){
  
  
  t<-which(separated.sim.based.on.time.asphericity[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.asphericity[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.asphericity[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.asphericity.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.asphericity.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.asphericity.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  t<-which(separated.sim.based.on.time.max.displacement.ratio[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.max.displacement.ratio[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.max.displacement.ratio[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.max.displacement.ratio.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.max.displacement.ratio.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.max.displacement.ratio.track[k]<-sd(values)
  
  
  
  ##########################################################################################################################
  
  
  
  
  
  
  t<-which(separated.sim.based.on.time.hurst[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.hurst[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.hurst[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.hurst.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.hurst.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.hurst.track[k]<-sd(values)
  
  
  
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  
  
  t<-which(separated.sim.based.on.time.displacement[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.displacement[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.displacement[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.displacement.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.displacement.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.displacement.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  t<-which(separated.sim.based.on.time.squared.displacement[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.squared.displacement[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.squared.displacement[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.squareDisplacement.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.squareDisplacement.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.squareDisplacement.track[k]<-sd(values)
  
  
  
  
  ##########################################################################################################################
  
  
  
  t<-which(separated.sim.based.on.time.straightness[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.straightness[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.straightness[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.straightness.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.straightness.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.straightness.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  t<-which(separated.sim.based.on.time.outreachratio[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.outreachratio[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.outreachratio[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.outreachRatio.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.outreachRatio.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.outreachRatio.track[k]<-sd(values)
  
  
  
  
  
  
  
}	  	   








time.msd<-21



M.G.S.hurst<-measure.greg.hurst.simulated[3: time.msd]/counter.time.simulated[3: time.msd]




M.G.S.diplacement<-measure.greg.displacement.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.max.diplacement.ratio<-measure.greg.max.displacement.ratio.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.asphericity<-measure.greg.asphericity.simulated[2: time.msd]/counter.time.simulated[2: time.msd]




M.G.S.squared.diplacement<-measure.greg.squared.displacement.simulated[2: time.msd]/counter.time.simulated[2: time.msd]




M.G.S.outreachratio<-measure.greg.outreachratio.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.straightness<-measure.greg.straightness.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



##############################################


setwd("/sbb/scratch/torkashe/Spherical-Coordinates-Motility-Metrics/Asphericity/")



ylab<-c(10, 100, 1000)

ylab.values<-log(ylab)/log(10)

x.lab<-seq(2, (time.msd-1), 1)

xlab.values<-log(x.lab)/log(10)




points(x.lab, M.G.S.asphericity, col="yellow", pch=20, type="l", cex=1.5, lwd=1.5)


# arrows(seq(1, (time.msd-1), 1), M.G.S.diplacement-simulated.sd.fence.displacement.track[2: time.msd], seq(1, (time.msd-1), 1), M.G.S.diplacement+simulated.sd.fence.displacement.track[2: time.msd], length=0.05, angle=90, code=3, col="darkblue")



# # # # u<-ecdf(M.G.S.diplacement)

# # plot(u, add=T, col="darkred") 

###################################################################################################################################################################################################################################################


remove(list=ls())		

##########################################################################################################################    

setwd("/sbb/scratch/torkashe/Clustering-Results/measures-of-motility-3-clusters-mokh-2-states-HHMM/")

load(".RMData")



simulated.mean.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)



simulated.lower.fence.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)




simulated.upper.fence.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)





simulated.mean.hurst.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.asphericity.track.observation<-matrix(0, nrow= number.of.positions.total)




simulated.lower.fence.hurst.track.observation<-matrix(0, nrow= number.of.positions.total)  	   

simulated.lower.fence.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.asphericity.observation<-matrix(0, nrow= number.of.positions.total)




simulated.upper.fence.hurst.track.observation<-matrix(0, nrow= number.of.positions.total) 	   

simulated.upper.fence.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.asphericity.observation<-matrix(0, nrow= number.of.positions.total)









for(i in 1: number.of.positions.total){
  
  
  
  
  simulated.mean.hurst.track.inside[i]<-mean(mean.hurst.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.hurst.track.inside[1:(r-1), i])[5]-summary(mean.hurst.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.hurst.track.inside[i]<-summary(mean.hurst.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.hurst.track.inside[i]<-summary(mean.hurst.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.displacement.track.inside[i]<-mean(mean.displacement.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.displacement.track.inside[1:(r-1), i])[5]-summary(mean.displacement.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.displacement.track.inside[i]<-summary(mean.displacement.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.displacement.track.inside[i]<-summary(mean.displacement.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.squareDisplacement.track.inside[i]<-mean(mean.squareDisplacement.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[5]-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.squareDisplacement.track.inside[i]<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.squareDisplacement.track.inside[i]<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  simulated.mean.outreachRatio.track.inside[i]<-mean(mean.outreachRatio.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[5]-summary(mean.outreachRatio.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.outreachRatio.track.inside[i]<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.outreachRatio.track.inside[i]<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.straightness.track.inside[i]<-mean(mean.straightness.track.inside[1:(r-1), i])
  
  
  
  
  
  iqr<-summary(mean.straightness.track.inside[1:(r-1), i])[5]-summary(mean.straightness.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.straightness.track.inside[i]<-summary(mean.straightness.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.straightness.track.inside[i]<-summary(mean.straightness.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  simulated.mean.max.displacement.ratio.track.inside[i]<-mean(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])
  
  
  
  
  iqr<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[5]-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.max.displacement.ratio.track.inside[i]<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.max.displacement.ratio.track.inside[i]<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  
  simulated.mean.meanTurningAngle.track.inside[i]<-mean(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[5]-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.meanTurningAngle.track.inside[i]<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.meanTurningAngle.track.inside[i]<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  simulated.mean.asphericity.track.inside[i]<-mean(mean.simulated.asphericity.inside[1:(r-1), i])
  
  
  
  
  iqr<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[5]-summary(mean.simulated.asphericity.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.asphericity.track.inside[i]<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.asphericity.track.inside[i]<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  ##############################################
  ##############################################    
  
  
  simulated.mean.hurst.track.observation[i]<-mean(mean.hurst.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.hurst.track.observation[1:(r-1), i])[5]-summary(mean.hurst.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.hurst.track.observation[i]<-summary(mean.hurst.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.hurst.track.observation[i]<-summary(mean.hurst.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  simulated.mean.displacement.track.observation[i]<-mean(mean.displacement.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.displacement.track.observation[1:(r-1), i])[5]-summary(mean.displacement.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.displacement.track.observation[i]<-summary(mean.displacement.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.displacement.track.observation[i]<-summary(mean.displacement.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.squareDisplacement.track.observation[i]<-mean(mean.squareDisplacement.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[5]-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.squareDisplacement.track.observation[i]<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.squareDisplacement.track.observation[i]<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  simulated.mean.outreachRatio.track.observation[i]<-mean(mean.outreachRatio.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[5]-summary(mean.outreachRatio.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.outreachRatio.track.observation[i]<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.outreachRatio.track.observation[i]<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.straightness.track.observation[i]<-mean(mean.straightness.track.observation[1:(r-1), i])
  
  
  
  
  
  iqr<-summary(mean.straightness.track.observation[1:(r-1), i])[5]-summary(mean.straightness.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.straightness.track.observation[i]<-summary(mean.straightness.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.straightness.track.observation[i]<-summary(mean.straightness.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  
  
  
  
  
}	   

###################################################################################################################################################################################################################################################




time.msd<-21


counter.time.simulated<-matrix(0, nrow=time.msd)


measure.greg.hurst.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.hurst<-matrix(0, nrow= time.msd, ncol=number.of.cells)



measure.greg.displacement.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.displacement<-matrix(0, nrow= time.msd, ncol=number.of.cells)


measure.greg.max.displacement.ratio.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.max.displacement.ratio<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.squared.displacement.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.squared.displacement<-matrix(0, nrow= time.msd, ncol=number.of.cells)



measure.greg.asphericity.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.asphericity<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.straightness.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.straightness<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.outreachratio.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.outreachratio<-matrix(0, nrow= time.msd, ncol=number.of.cells)



true.radius<-matrix(0, nrow= time.msd, ncol=number.of.cells)






simulated.lower.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.asphericity.track<-matrix(0, nrow= time.msd)




simulated.upper.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.asphericity.track<-matrix(0, nrow= time.msd)



simulated.sd.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.asphericity.track<-matrix(0, nrow= time.msd)





for(i in 1:number.of.cells){
  
  
  
  index.i<-which(data_cell_motility[, 8]==id[i])
  
  
  
  # # # # for(i in 1:length(cluster.1[cluster.1>249])){
  
  # # t<-cluster.1[cluster.1>249]-249; 
  
  # # index.i<-which(data_cell_motility[, 8]==id[t[i]])
  
  
  
  
  l<-length(index.i)
  
  for(k in 2: time.msd){
    
    
    
    measure.greg.asphericity.simulated[k]<-measure.greg.asphericity.simulated[k]+ simulated.mean.asphericity.track.inside[index.i[k]]
    
    separated.sim.based.on.time.asphericity[k, i]<-simulated.mean.asphericity.track.inside[index.i[k]]
    
    
    
    
    # t<-0
    
    
    measure.greg.displacement.simulated[k]<-measure.greg.displacement.simulated[k]+ simulated.mean.displacement.track.inside[index.i[k]]
    
    separated.sim.based.on.time.displacement[k, i]<-simulated.mean.displacement.track.inside[index.i[k]]
    
    
    
    measure.greg.hurst.simulated[k]<-measure.greg.hurst.simulated[k]+ simulated.mean.hurst.track.inside[index.i[k]]
    
    separated.sim.based.on.time.hurst[k, i]<-simulated.mean.hurst.track.inside[index.i[k]]
    
    
    
    
    measure.greg.max.displacement.ratio.simulated[k]<-measure.greg.max.displacement.ratio.simulated[k]+ simulated.mean.max.displacement.ratio.track.inside[index.i[k]]
    
    separated.sim.based.on.time.max.displacement.ratio[k, i]<-simulated.mean.max.displacement.ratio.track.inside[index.i[k]]
    
    
    
    
    
    measure.greg.squared.displacement.simulated[k] <-measure.greg.squared.displacement.simulated[k]+ simulated.mean.squareDisplacement.track.inside[index.i[k]]
    
    separated.sim.based.on.time.squared.displacement[k, i]<-simulated.mean.squareDisplacement.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    measure.greg.straightness.simulated[k] <-measure.greg.straightness.simulated[k] +simulated.mean.straightness.track.inside[index.i[k]]
    
    separated.sim.based.on.time.straightness[k, i]<-simulated.mean.straightness.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    
    measure.greg.outreachratio.simulated[k] <-measure.greg.outreachratio.simulated[k]+simulated.mean.outreachRatio.track.inside[index.i[k]]
    
    separated.sim.based.on.time.outreachratio[k, i]<-simulated.mean.outreachRatio.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    counter.time.simulated[k]<-counter.time.simulated[k]+1
    
    
  }	
  
  
  
  # # # # # # 		            for(k in 30: time.msd){
  
  # # # if(l>=k){
  
  # # # # t<-0
  
  
  # # # measure.greg.displacement.simulated[k]<-measure.greg.displacement.simulated[k]+ simulated.mean.displacement.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.displacement[k, i]<-simulated.mean.displacement.track.inside[index.i[k]]
  
  
  
  
  # # # measure.greg.hurst.simulated[k]<-measure.greg.hurst.simulated[k]+ simulated.mean.hurst.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.hurst[k, i]<-simulated.mean.hurst.track.inside[index.i[k]]
  
  
  
  
  
  # # # measure.greg.squared.displacement.simulated[k] <-measure.greg.squared.displacement.simulated[k]+ simulated.mean.squareDisplacement.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.squared.displacement[k, i]<-simulated.mean.squareDisplacement.track.inside[index.i[k]]
  
  
  
  
  
  
  # # # measure.greg.straightness.simulated[k] <-measure.greg.straightness.simulated[k] +simulated.mean.straightness.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.straightness[k, i]<-simulated.mean.straightness.track.inside[index.i[k]]
  
  
  
  
  
  
  
  # # # measure.greg.outreachratio.simulated[k] <-measure.greg.outreachratio.simulated[k]+simulated.mean.outreachRatio.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.outreachratio[k, i]<-simulated.mean.outreachRatio.track.inside[index.i[k]]
  
  
  # # # counter.time.simulated[k]<-counter.time.simulated[k]+1
  
  
  # # # }             
  
  # # # }		
  
  
}  




for(k in 2:time.msd){
  
  
  t<-which(separated.sim.based.on.time.asphericity[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.asphericity[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.asphericity[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.asphericity.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.asphericity.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.asphericity.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  t<-which(separated.sim.based.on.time.max.displacement.ratio[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.max.displacement.ratio[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.max.displacement.ratio[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.max.displacement.ratio.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.max.displacement.ratio.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.max.displacement.ratio.track[k]<-sd(values)
  
  
  
  ##########################################################################################################################
  
  
  
  
  
  
  t<-which(separated.sim.based.on.time.hurst[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.hurst[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.hurst[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.hurst.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.hurst.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.hurst.track[k]<-sd(values)
  
  
  
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  
  
  t<-which(separated.sim.based.on.time.displacement[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.displacement[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.displacement[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.displacement.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.displacement.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.displacement.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  t<-which(separated.sim.based.on.time.squared.displacement[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.squared.displacement[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.squared.displacement[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.squareDisplacement.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.squareDisplacement.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.squareDisplacement.track[k]<-sd(values)
  
  
  
  
  ##########################################################################################################################
  
  
  
  t<-which(separated.sim.based.on.time.straightness[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.straightness[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.straightness[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.straightness.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.straightness.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.straightness.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  t<-which(separated.sim.based.on.time.outreachratio[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.outreachratio[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.outreachratio[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.outreachRatio.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.outreachRatio.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.outreachRatio.track[k]<-sd(values)
  
  
  
  
  
  
  
}	  	   








time.msd<-21



M.G.S.hurst<-measure.greg.hurst.simulated[3: time.msd]/counter.time.simulated[3: time.msd]




M.G.S.diplacement<-measure.greg.displacement.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.max.diplacement.ratio<-measure.greg.max.displacement.ratio.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.asphericity<-measure.greg.asphericity.simulated[2: time.msd]/counter.time.simulated[2: time.msd]




M.G.S.squared.diplacement<-measure.greg.squared.displacement.simulated[2: time.msd]/counter.time.simulated[2: time.msd]




M.G.S.outreachratio<-measure.greg.outreachratio.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.straightness<-measure.greg.straightness.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



##############################################


setwd("/sbb/scratch/torkashe/Spherical-Coordinates-Motility-Metrics/Asphericity/")



ylab<-c(10, 100, 1000)

ylab.values<-log(ylab)/log(10)

x.lab<-seq(2, (time.msd-1), 1)

xlab.values<-log(x.lab)/log(10)




points(x.lab, M.G.S.asphericity, col="yellow", pch=11, type="l", cex=1.5, lwd=1.5)


# arrows(seq(1, (time.msd-1), 1), M.G.S.diplacement-simulated.sd.fence.displacement.track[2: time.msd], seq(1, (time.msd-1), 1), M.G.S.diplacement+simulated.sd.fence.displacement.track[2: time.msd], length=0.05, angle=90, code=3, col="darkblue")



# # # # u<-ecdf(M.G.S.diplacement)

# # plot(u, add=T, col="darkred") 



##################################################################################################################################################################################################################################################

remove(list=ls())


##########################################################################################################################    

setwd("/sbb/scratch/torkashe/Clustering-Results/measures-of-motility-4-clusters-spectral-2-states-NHMM/")

load(".RMData")




simulated.mean.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)



simulated.lower.fence.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)




simulated.upper.fence.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)





simulated.mean.hurst.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.asphericity.track.observation<-matrix(0, nrow= number.of.positions.total)




simulated.lower.fence.hurst.track.observation<-matrix(0, nrow= number.of.positions.total)  	   

simulated.lower.fence.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.asphericity.observation<-matrix(0, nrow= number.of.positions.total)




simulated.upper.fence.hurst.track.observation<-matrix(0, nrow= number.of.positions.total) 	   

simulated.upper.fence.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.asphericity.observation<-matrix(0, nrow= number.of.positions.total)









for(i in 1: number.of.positions.total){
  
  
  
  
  simulated.mean.hurst.track.inside[i]<-mean(mean.hurst.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.hurst.track.inside[1:(r-1), i])[5]-summary(mean.hurst.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.hurst.track.inside[i]<-summary(mean.hurst.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.hurst.track.inside[i]<-summary(mean.hurst.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.displacement.track.inside[i]<-mean(mean.displacement.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.displacement.track.inside[1:(r-1), i])[5]-summary(mean.displacement.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.displacement.track.inside[i]<-summary(mean.displacement.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.displacement.track.inside[i]<-summary(mean.displacement.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.squareDisplacement.track.inside[i]<-mean(mean.squareDisplacement.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[5]-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.squareDisplacement.track.inside[i]<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.squareDisplacement.track.inside[i]<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  simulated.mean.outreachRatio.track.inside[i]<-mean(mean.outreachRatio.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[5]-summary(mean.outreachRatio.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.outreachRatio.track.inside[i]<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.outreachRatio.track.inside[i]<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.straightness.track.inside[i]<-mean(mean.straightness.track.inside[1:(r-1), i])
  
  
  
  
  
  iqr<-summary(mean.straightness.track.inside[1:(r-1), i])[5]-summary(mean.straightness.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.straightness.track.inside[i]<-summary(mean.straightness.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.straightness.track.inside[i]<-summary(mean.straightness.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  simulated.mean.max.displacement.ratio.track.inside[i]<-mean(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])
  
  
  
  
  iqr<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[5]-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.max.displacement.ratio.track.inside[i]<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.max.displacement.ratio.track.inside[i]<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  
  simulated.mean.meanTurningAngle.track.inside[i]<-mean(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[5]-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.meanTurningAngle.track.inside[i]<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.meanTurningAngle.track.inside[i]<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  simulated.mean.asphericity.track.inside[i]<-mean(mean.simulated.asphericity.inside[1:(r-1), i])
  
  
  
  
  iqr<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[5]-summary(mean.simulated.asphericity.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.asphericity.track.inside[i]<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.asphericity.track.inside[i]<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  ##############################################
  ##############################################    
  
  
  simulated.mean.hurst.track.observation[i]<-mean(mean.hurst.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.hurst.track.observation[1:(r-1), i])[5]-summary(mean.hurst.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.hurst.track.observation[i]<-summary(mean.hurst.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.hurst.track.observation[i]<-summary(mean.hurst.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  simulated.mean.displacement.track.observation[i]<-mean(mean.displacement.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.displacement.track.observation[1:(r-1), i])[5]-summary(mean.displacement.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.displacement.track.observation[i]<-summary(mean.displacement.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.displacement.track.observation[i]<-summary(mean.displacement.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.squareDisplacement.track.observation[i]<-mean(mean.squareDisplacement.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[5]-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.squareDisplacement.track.observation[i]<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.squareDisplacement.track.observation[i]<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  simulated.mean.outreachRatio.track.observation[i]<-mean(mean.outreachRatio.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[5]-summary(mean.outreachRatio.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.outreachRatio.track.observation[i]<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.outreachRatio.track.observation[i]<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.straightness.track.observation[i]<-mean(mean.straightness.track.observation[1:(r-1), i])
  
  
  
  
  
  iqr<-summary(mean.straightness.track.observation[1:(r-1), i])[5]-summary(mean.straightness.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.straightness.track.observation[i]<-summary(mean.straightness.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.straightness.track.observation[i]<-summary(mean.straightness.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  
  
  
  
  
}	   

###################################################################################################################################################################################################################################################




time.msd<-21


counter.time.simulated<-matrix(0, nrow=time.msd)


measure.greg.hurst.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.hurst<-matrix(0, nrow= time.msd, ncol=number.of.cells)



measure.greg.displacement.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.displacement<-matrix(0, nrow= time.msd, ncol=number.of.cells)


measure.greg.max.displacement.ratio.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.max.displacement.ratio<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.squared.displacement.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.squared.displacement<-matrix(0, nrow= time.msd, ncol=number.of.cells)



measure.greg.asphericity.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.asphericity<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.straightness.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.straightness<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.outreachratio.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.outreachratio<-matrix(0, nrow= time.msd, ncol=number.of.cells)



true.radius<-matrix(0, nrow= time.msd, ncol=number.of.cells)






simulated.lower.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.asphericity.track<-matrix(0, nrow= time.msd)




simulated.upper.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.asphericity.track<-matrix(0, nrow= time.msd)



simulated.sd.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.asphericity.track<-matrix(0, nrow= time.msd)





for(i in 1:number.of.cells){
  
  
  
  index.i<-which(data_cell_motility[, 8]==id[i])
  
  
  
  # # # # for(i in 1:length(cluster.1[cluster.1>249])){
  
  # # t<-cluster.1[cluster.1>249]-249; 
  
  # # index.i<-which(data_cell_motility[, 8]==id[t[i]])
  
  
  
  
  l<-length(index.i)
  
  for(k in 2: time.msd){
    
    
    
    measure.greg.asphericity.simulated[k]<-measure.greg.asphericity.simulated[k]+ simulated.mean.asphericity.track.inside[index.i[k]]
    
    separated.sim.based.on.time.asphericity[k, i]<-simulated.mean.asphericity.track.inside[index.i[k]]
    
    
    
    
    # t<-0
    
    
    measure.greg.displacement.simulated[k]<-measure.greg.displacement.simulated[k]+ simulated.mean.displacement.track.inside[index.i[k]]
    
    separated.sim.based.on.time.displacement[k, i]<-simulated.mean.displacement.track.inside[index.i[k]]
    
    
    
    measure.greg.hurst.simulated[k]<-measure.greg.hurst.simulated[k]+ simulated.mean.hurst.track.inside[index.i[k]]
    
    separated.sim.based.on.time.hurst[k, i]<-simulated.mean.hurst.track.inside[index.i[k]]
    
    
    
    
    measure.greg.max.displacement.ratio.simulated[k]<-measure.greg.max.displacement.ratio.simulated[k]+ simulated.mean.max.displacement.ratio.track.inside[index.i[k]]
    
    separated.sim.based.on.time.max.displacement.ratio[k, i]<-simulated.mean.max.displacement.ratio.track.inside[index.i[k]]
    
    
    
    
    
    measure.greg.squared.displacement.simulated[k] <-measure.greg.squared.displacement.simulated[k]+ simulated.mean.squareDisplacement.track.inside[index.i[k]]
    
    separated.sim.based.on.time.squared.displacement[k, i]<-simulated.mean.squareDisplacement.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    measure.greg.straightness.simulated[k] <-measure.greg.straightness.simulated[k] +simulated.mean.straightness.track.inside[index.i[k]]
    
    separated.sim.based.on.time.straightness[k, i]<-simulated.mean.straightness.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    
    measure.greg.outreachratio.simulated[k] <-measure.greg.outreachratio.simulated[k]+simulated.mean.outreachRatio.track.inside[index.i[k]]
    
    separated.sim.based.on.time.outreachratio[k, i]<-simulated.mean.outreachRatio.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    counter.time.simulated[k]<-counter.time.simulated[k]+1
    
    
  }	
  
  
  
  # # # # # # 		            for(k in 30: time.msd){
  
  # # # if(l>=k){
  
  # # # # t<-0
  
  
  # # # measure.greg.displacement.simulated[k]<-measure.greg.displacement.simulated[k]+ simulated.mean.displacement.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.displacement[k, i]<-simulated.mean.displacement.track.inside[index.i[k]]
  
  
  
  
  # # # measure.greg.hurst.simulated[k]<-measure.greg.hurst.simulated[k]+ simulated.mean.hurst.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.hurst[k, i]<-simulated.mean.hurst.track.inside[index.i[k]]
  
  
  
  
  
  # # # measure.greg.squared.displacement.simulated[k] <-measure.greg.squared.displacement.simulated[k]+ simulated.mean.squareDisplacement.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.squared.displacement[k, i]<-simulated.mean.squareDisplacement.track.inside[index.i[k]]
  
  
  
  
  
  
  # # # measure.greg.straightness.simulated[k] <-measure.greg.straightness.simulated[k] +simulated.mean.straightness.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.straightness[k, i]<-simulated.mean.straightness.track.inside[index.i[k]]
  
  
  
  
  
  
  
  # # # measure.greg.outreachratio.simulated[k] <-measure.greg.outreachratio.simulated[k]+simulated.mean.outreachRatio.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.outreachratio[k, i]<-simulated.mean.outreachRatio.track.inside[index.i[k]]
  
  
  # # # counter.time.simulated[k]<-counter.time.simulated[k]+1
  
  
  # # # }             
  
  # # # }		
  
  
}  




for(k in 2:time.msd){
  
  
  t<-which(separated.sim.based.on.time.asphericity[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.asphericity[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.asphericity[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.asphericity.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.asphericity.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.asphericity.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  t<-which(separated.sim.based.on.time.max.displacement.ratio[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.max.displacement.ratio[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.max.displacement.ratio[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.max.displacement.ratio.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.max.displacement.ratio.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.max.displacement.ratio.track[k]<-sd(values)
  
  
  
  ##########################################################################################################################
  
  
  
  
  
  
  t<-which(separated.sim.based.on.time.hurst[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.hurst[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.hurst[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.hurst.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.hurst.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.hurst.track[k]<-sd(values)
  
  
  
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  
  
  t<-which(separated.sim.based.on.time.displacement[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.displacement[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.displacement[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.displacement.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.displacement.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.displacement.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  t<-which(separated.sim.based.on.time.squared.displacement[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.squared.displacement[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.squared.displacement[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.squareDisplacement.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.squareDisplacement.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.squareDisplacement.track[k]<-sd(values)
  
  
  
  
  ##########################################################################################################################
  
  
  
  t<-which(separated.sim.based.on.time.straightness[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.straightness[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.straightness[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.straightness.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.straightness.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.straightness.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  t<-which(separated.sim.based.on.time.outreachratio[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.outreachratio[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.outreachratio[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.outreachRatio.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.outreachRatio.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.outreachRatio.track[k]<-sd(values)
  
  
  
  
  
  
  
}	  	   







time.msd<-21



M.G.S.hurst<-measure.greg.hurst.simulated[3: time.msd]/counter.time.simulated[3: time.msd]




M.G.S.diplacement<-measure.greg.displacement.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.max.diplacement.ratio<-measure.greg.max.displacement.ratio.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.asphericity<-measure.greg.asphericity.simulated[2: time.msd]/counter.time.simulated[2: time.msd]




M.G.S.squared.diplacement<-measure.greg.squared.displacement.simulated[2: time.msd]/counter.time.simulated[2: time.msd]




M.G.S.outreachratio<-measure.greg.outreachratio.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.straightness<-measure.greg.straightness.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



##############################################


setwd("/sbb/scratch/torkashe/Spherical-Coordinates-Motility-Metrics/Asphericity/")



ylab<-c(10, 100, 1000)

ylab.values<-log(ylab)/log(10)

x.lab<-seq(2, (time.msd-1), 1)

xlab.values<-log(x.lab)/log(10)




points(x.lab, M.G.S.asphericity, col="cyan", pch=20, type="l", cex=1.5, lwd=1.5)


# arrows(seq(1, (time.msd-1), 1), M.G.S.diplacement-simulated.sd.fence.displacement.track[2: time.msd], seq(1, (time.msd-1), 1), M.G.S.diplacement+simulated.sd.fence.displacement.track[2: time.msd], length=0.05, angle=90, code=3, col="darkmagenta")



# # # # u<-ecdf(M.G.S.diplacement)

# # plot(u, add=T, col="darkred") 



##################################################################################################################################################################################################################################################

remove(list=ls())


##########################################################################################################################    

setwd("/sbb/scratch/torkashe/Clustering-Results/measures-of-motility-2-clusters-spectral-2-states-NHMM/")

load(".RMData")




simulated.mean.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)



simulated.lower.fence.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)




simulated.upper.fence.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)





simulated.mean.hurst.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.asphericity.track.observation<-matrix(0, nrow= number.of.positions.total)




simulated.lower.fence.hurst.track.observation<-matrix(0, nrow= number.of.positions.total)  	   

simulated.lower.fence.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.asphericity.observation<-matrix(0, nrow= number.of.positions.total)




simulated.upper.fence.hurst.track.observation<-matrix(0, nrow= number.of.positions.total) 	   

simulated.upper.fence.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.asphericity.observation<-matrix(0, nrow= number.of.positions.total)









for(i in 1: number.of.positions.total){
  
  
  
  
  simulated.mean.hurst.track.inside[i]<-mean(mean.hurst.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.hurst.track.inside[1:(r-1), i])[5]-summary(mean.hurst.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.hurst.track.inside[i]<-summary(mean.hurst.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.hurst.track.inside[i]<-summary(mean.hurst.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.displacement.track.inside[i]<-mean(mean.displacement.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.displacement.track.inside[1:(r-1), i])[5]-summary(mean.displacement.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.displacement.track.inside[i]<-summary(mean.displacement.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.displacement.track.inside[i]<-summary(mean.displacement.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.squareDisplacement.track.inside[i]<-mean(mean.squareDisplacement.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[5]-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.squareDisplacement.track.inside[i]<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.squareDisplacement.track.inside[i]<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  simulated.mean.outreachRatio.track.inside[i]<-mean(mean.outreachRatio.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[5]-summary(mean.outreachRatio.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.outreachRatio.track.inside[i]<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.outreachRatio.track.inside[i]<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.straightness.track.inside[i]<-mean(mean.straightness.track.inside[1:(r-1), i])
  
  
  
  
  
  iqr<-summary(mean.straightness.track.inside[1:(r-1), i])[5]-summary(mean.straightness.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.straightness.track.inside[i]<-summary(mean.straightness.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.straightness.track.inside[i]<-summary(mean.straightness.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  simulated.mean.max.displacement.ratio.track.inside[i]<-mean(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])
  
  
  
  
  iqr<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[5]-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.max.displacement.ratio.track.inside[i]<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.max.displacement.ratio.track.inside[i]<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  
  simulated.mean.meanTurningAngle.track.inside[i]<-mean(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[5]-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.meanTurningAngle.track.inside[i]<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.meanTurningAngle.track.inside[i]<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  simulated.mean.asphericity.track.inside[i]<-mean(mean.simulated.asphericity.inside[1:(r-1), i])
  
  
  
  
  iqr<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[5]-summary(mean.simulated.asphericity.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.asphericity.track.inside[i]<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.asphericity.track.inside[i]<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  ##############################################
  ##############################################    
  
  
  simulated.mean.hurst.track.observation[i]<-mean(mean.hurst.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.hurst.track.observation[1:(r-1), i])[5]-summary(mean.hurst.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.hurst.track.observation[i]<-summary(mean.hurst.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.hurst.track.observation[i]<-summary(mean.hurst.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  simulated.mean.displacement.track.observation[i]<-mean(mean.displacement.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.displacement.track.observation[1:(r-1), i])[5]-summary(mean.displacement.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.displacement.track.observation[i]<-summary(mean.displacement.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.displacement.track.observation[i]<-summary(mean.displacement.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.squareDisplacement.track.observation[i]<-mean(mean.squareDisplacement.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[5]-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.squareDisplacement.track.observation[i]<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.squareDisplacement.track.observation[i]<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  simulated.mean.outreachRatio.track.observation[i]<-mean(mean.outreachRatio.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[5]-summary(mean.outreachRatio.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.outreachRatio.track.observation[i]<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.outreachRatio.track.observation[i]<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.straightness.track.observation[i]<-mean(mean.straightness.track.observation[1:(r-1), i])
  
  
  
  
  
  iqr<-summary(mean.straightness.track.observation[1:(r-1), i])[5]-summary(mean.straightness.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.straightness.track.observation[i]<-summary(mean.straightness.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.straightness.track.observation[i]<-summary(mean.straightness.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  
  
  
  
  
}	   

###################################################################################################################################################################################################################################################




time.msd<-21


counter.time.simulated<-matrix(0, nrow=time.msd)


measure.greg.hurst.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.hurst<-matrix(0, nrow= time.msd, ncol=number.of.cells)



measure.greg.displacement.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.displacement<-matrix(0, nrow= time.msd, ncol=number.of.cells)


measure.greg.max.displacement.ratio.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.max.displacement.ratio<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.squared.displacement.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.squared.displacement<-matrix(0, nrow= time.msd, ncol=number.of.cells)



measure.greg.asphericity.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.asphericity<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.straightness.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.straightness<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.outreachratio.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.outreachratio<-matrix(0, nrow= time.msd, ncol=number.of.cells)



true.radius<-matrix(0, nrow= time.msd, ncol=number.of.cells)






simulated.lower.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.asphericity.track<-matrix(0, nrow= time.msd)




simulated.upper.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.asphericity.track<-matrix(0, nrow= time.msd)



simulated.sd.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.asphericity.track<-matrix(0, nrow= time.msd)





for(i in 1:number.of.cells){
  
  
  
  index.i<-which(data_cell_motility[, 8]==id[i])
  
  
  
  # # # # for(i in 1:length(cluster.1[cluster.1>249])){
  
  # # t<-cluster.1[cluster.1>249]-249; 
  
  # # index.i<-which(data_cell_motility[, 8]==id[t[i]])
  
  
  
  
  l<-length(index.i)
  
  for(k in 2: time.msd){
    
    
    
    measure.greg.asphericity.simulated[k]<-measure.greg.asphericity.simulated[k]+ simulated.mean.asphericity.track.inside[index.i[k]]
    
    separated.sim.based.on.time.asphericity[k, i]<-simulated.mean.asphericity.track.inside[index.i[k]]
    
    
    
    
    # t<-0
    
    
    measure.greg.displacement.simulated[k]<-measure.greg.displacement.simulated[k]+ simulated.mean.displacement.track.inside[index.i[k]]
    
    separated.sim.based.on.time.displacement[k, i]<-simulated.mean.displacement.track.inside[index.i[k]]
    
    
    
    measure.greg.hurst.simulated[k]<-measure.greg.hurst.simulated[k]+ simulated.mean.hurst.track.inside[index.i[k]]
    
    separated.sim.based.on.time.hurst[k, i]<-simulated.mean.hurst.track.inside[index.i[k]]
    
    
    
    
    measure.greg.max.displacement.ratio.simulated[k]<-measure.greg.max.displacement.ratio.simulated[k]+ simulated.mean.max.displacement.ratio.track.inside[index.i[k]]
    
    separated.sim.based.on.time.max.displacement.ratio[k, i]<-simulated.mean.max.displacement.ratio.track.inside[index.i[k]]
    
    
    
    
    
    measure.greg.squared.displacement.simulated[k] <-measure.greg.squared.displacement.simulated[k]+ simulated.mean.squareDisplacement.track.inside[index.i[k]]
    
    separated.sim.based.on.time.squared.displacement[k, i]<-simulated.mean.squareDisplacement.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    measure.greg.straightness.simulated[k] <-measure.greg.straightness.simulated[k] +simulated.mean.straightness.track.inside[index.i[k]]
    
    separated.sim.based.on.time.straightness[k, i]<-simulated.mean.straightness.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    
    measure.greg.outreachratio.simulated[k] <-measure.greg.outreachratio.simulated[k]+simulated.mean.outreachRatio.track.inside[index.i[k]]
    
    separated.sim.based.on.time.outreachratio[k, i]<-simulated.mean.outreachRatio.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    counter.time.simulated[k]<-counter.time.simulated[k]+1
    
    
  }	
  
  
  
  # # # # # # 		            for(k in 30: time.msd){
  
  # # # if(l>=k){
  
  # # # # t<-0
  
  
  # # # measure.greg.displacement.simulated[k]<-measure.greg.displacement.simulated[k]+ simulated.mean.displacement.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.displacement[k, i]<-simulated.mean.displacement.track.inside[index.i[k]]
  
  
  
  
  # # # measure.greg.hurst.simulated[k]<-measure.greg.hurst.simulated[k]+ simulated.mean.hurst.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.hurst[k, i]<-simulated.mean.hurst.track.inside[index.i[k]]
  
  
  
  
  
  # # # measure.greg.squared.displacement.simulated[k] <-measure.greg.squared.displacement.simulated[k]+ simulated.mean.squareDisplacement.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.squared.displacement[k, i]<-simulated.mean.squareDisplacement.track.inside[index.i[k]]
  
  
  
  
  
  
  # # # measure.greg.straightness.simulated[k] <-measure.greg.straightness.simulated[k] +simulated.mean.straightness.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.straightness[k, i]<-simulated.mean.straightness.track.inside[index.i[k]]
  
  
  
  
  
  
  
  # # # measure.greg.outreachratio.simulated[k] <-measure.greg.outreachratio.simulated[k]+simulated.mean.outreachRatio.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.outreachratio[k, i]<-simulated.mean.outreachRatio.track.inside[index.i[k]]
  
  
  # # # counter.time.simulated[k]<-counter.time.simulated[k]+1
  
  
  # # # }             
  
  # # # }		
  
  
}  




for(k in 2:time.msd){
  
  
  t<-which(separated.sim.based.on.time.asphericity[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.asphericity[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.asphericity[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.asphericity.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.asphericity.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.asphericity.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  t<-which(separated.sim.based.on.time.max.displacement.ratio[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.max.displacement.ratio[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.max.displacement.ratio[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.max.displacement.ratio.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.max.displacement.ratio.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.max.displacement.ratio.track[k]<-sd(values)
  
  
  
  ##########################################################################################################################
  
  
  
  
  
  
  t<-which(separated.sim.based.on.time.hurst[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.hurst[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.hurst[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.hurst.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.hurst.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.hurst.track[k]<-sd(values)
  
  
  
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  
  
  t<-which(separated.sim.based.on.time.displacement[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.displacement[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.displacement[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.displacement.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.displacement.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.displacement.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  t<-which(separated.sim.based.on.time.squared.displacement[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.squared.displacement[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.squared.displacement[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.squareDisplacement.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.squareDisplacement.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.squareDisplacement.track[k]<-sd(values)
  
  
  
  
  ##########################################################################################################################
  
  
  
  t<-which(separated.sim.based.on.time.straightness[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.straightness[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.straightness[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.straightness.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.straightness.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.straightness.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  t<-which(separated.sim.based.on.time.outreachratio[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.outreachratio[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.outreachratio[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.outreachRatio.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.outreachRatio.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.outreachRatio.track[k]<-sd(values)
  
  
  
  
  
  
  
}	  	   







time.msd<-21



M.G.S.hurst<-measure.greg.hurst.simulated[3: time.msd]/counter.time.simulated[3: time.msd]




M.G.S.diplacement<-measure.greg.displacement.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.max.diplacement.ratio<-measure.greg.max.displacement.ratio.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.asphericity<-measure.greg.asphericity.simulated[2: time.msd]/counter.time.simulated[2: time.msd]




M.G.S.squared.diplacement<-measure.greg.squared.displacement.simulated[2: time.msd]/counter.time.simulated[2: time.msd]




M.G.S.outreachratio<-measure.greg.outreachratio.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.straightness<-measure.greg.straightness.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



##############################################


setwd("/sbb/scratch/torkashe/Spherical-Coordinates-Motility-Metrics/Asphericity/")



ylab<-c(10, 100, 1000)

ylab.values<-log(ylab)/log(10)

x.lab<-seq(2, (time.msd-1), 1)

xlab.values<-log(x.lab)/log(10)




points(x.lab, M.G.S.asphericity, col="cyan", pch=11, type="l", cex=1.5, lwd=1.5, lty=3)


# arrows(seq(1, (time.msd-1), 1), M.G.S.diplacement-simulated.sd.fence.displacement.track[2: time.msd], seq(1, (time.msd-1), 1), M.G.S.diplacement+simulated.sd.fence.displacement.track[2: time.msd], length=0.05, angle=90, code=3, col="darkmagenta")



# # # # u<-ecdf(M.G.S.diplacement)

# # plot(u, add=T, col="darkred") 



###################################################################################################################################################################################################################################################


########  to work with held-out data: in the first step we call the simulated values from the fitted models. In the next steps, we plot them against the values obtained from held-out data.

remove(list=ls())


setwd("/sbb/scratch/torkashe/Clustering-Results/measures-of-motility-2-clusters-ET-2-states-HHMM/")



library(matlib)

library(sde)

library("rstan")

library("MotilityLab")

library(ggplot2)


####  Functions to generate trajectories.

##########   Generating trajectories


new.inc=function(r,phi,theta,x1){ 
  
  fake<-sum(x1^2)
  
  if(fake ==0){return("Error, first increment has zero length")}
  if( x1[3] != 0 ){
    A =cbind(x1,as.matrix(c(1,0,0)),as.matrix(c(0,1,0)))
    A.qr=qr(A)
    Qm=qr.Q(A.qr)
    u=Qm[,2]
    v=Qm[,3]
    c1=r*cos(phi)
    r2=sqrt(r^2-c1^2)
    c2=r2*cos(theta)
    c3=r2*sin(theta)
    x2=c1*x1/sqrt(fake)+c2*u+c3*v	
  }	
  if( (x1[3]==0)& (x1[1]!=0)){
    A =cbind(x1,as.matrix(c(0,0,1)),as.matrix(c(0,1,0)))
    A.qr=qr(A)
    Qm=qr.Q(A.qr)
    u=Qm[,2]
    v=Qm[,3]
    c1=r*cos(phi)
    r2=sqrt(r^2-c1^2)
    c2=r2*cos(theta)
    c3=r2*sin(theta)
    x2=c1*x1/sqrt(fake)+c2*u+c3*v
  }
  if( (x1[3]==0)& (x1[1]==0)){
    A =cbind(x1,as.matrix(c(0,0,1)),as.matrix(c(1,0,0)))
    A.qr=qr(A)
    Qm=qr.Q(A.qr)
    u=Qm[,2]
    v=Qm[,3]
    c1=r*cos(phi)
    r2=sqrt(r^2-c1^2)
    c2=r2*cos(theta)
    c3=r2*sin(theta)
    x2=c1*x1/sqrt(fake)+c2*u+c3*v
  }
  
  x2
}


s.angle=function(x1,x2){ #this calculates the second angle going from increment x1 to x2
  
  fake1<-sum(x1^2)
  
  fake2<-sum(x2^2)
  
  if( t(x1)%*%x2==sqrt(fake1)*sqrt(fake2) | t(x1)%*%x2==-sqrt(fake1)*sqrt(fake2)){return(0)}
  
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




###################################################################################################################
###############################################################  Spherical representation and angular change 


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







##################################################################################################################################################################################################################################################################################################  training Data  #########################################################





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


######################################################
######################################################
######################################################
########   measures obtained from the real data



mean.displacement.track.real<-matrix(0, nrow= number.of.positions.total)

mean.hurst.track.real<-matrix(0, nrow= number.of.positions.total)

mean.squareDisplacement.track.real<-matrix(0, nrow= number.of.positions.total)

mean.outreachRatio.track.real<-matrix(0, nrow= number.of.positions.total)

mean.straightness.track.real<-matrix(0, nrow= number.of.positions.total)

mean.max.displacement.ratio.track.real<-matrix(0, nrow= number.of.positions.total)

mean.meanTurningAngle.track.real<-matrix(0, nrow= number.of.positions.total)

mean.asphericity.track.real<-matrix(0, nrow= number.of.positions.total)

increment.x.total<-increment.y.total<-increment.z.total<-NULL

index.measures<-index.meanTurningAngle <-NULL

theta.cell.total<-matrix(0, nrow= number.of.positions.total)



for(i in 1: number.of.cells){
  
  
  hurst.track<-NULL
  
  displacement.track<-NULL
  
  squareDisplacement.track<-NULL
  
  outreachRatio.track<-NULL
  
  straightness.track<-NULL
  
  max.displacement.ratio.track<-NULL
  
  meanTurningAngle.track<-NULL
  
  asphericity.track<-NULL
  
  
  
  
  index.i<-which(data_cell_motility[, 8]==id[i])
  
  l<-length(index.i)
  
  index.measures<-c(index.measures, index.i[2:l])
  
  index.meanTurningAngle<-c(index.meanTurningAngle, index.i[2:(l-1)])
  
  
  data_i<-matrix(c(data_cell_motility[index.i, 7], data_cell_motility[index.i, 1], data_cell_motility[index.i, 2], data_cell_motility[index.i, 3]), nrow=length(index.i), ncol=4)
  
  
  
  increment.x<-increment.y<-increment.z<-NULL
  
  increment.x[1]<-data_cell_motility[index.i[1], 1]
  
  increment.y[1]<-data_cell_motility[index.i[1], 2]
  
  increment.z[1]<-data_cell_motility[index.i[1], 3]
  
  for(j in 2: length(index.i)){
    
    increment.x[j]<-data_cell_motility[index.i[j], 1]-data_cell_motility[index.i[j-1], 1]
    
    increment.y[j]<-data_cell_motility[index.i[j], 2]-data_cell_motility[index.i[j-1], 2]
    
    increment.z[j]<-data_cell_motility[index.i[j], 3]-data_cell_motility[index.i[j-1], 3]
    
  }
  
  theta.cell<-NULL
  
  
  for(j in 1:(length(index.i)-1)){
    
    
    t.1<-c(increment.x[j], increment.y[j], increment.z[j])
    
    t.2<-c(increment.x[j+1], increment.y[j+1], increment.z[j+1])
    
    theta<-s.angle(t.1, t.2)
    
    
    theta.cell<-c(theta.cell, theta)
    
    
    
  }
  
  theta.cell.total[index.i[1:(l-1)]]<-theta.cell
  
  increment.x.total<-c(increment.x[2:(l-1)], increment.x.total)
  
  increment.y.total<-c(increment.y[2:(l-1)], increment.y.total)
  
  increment.z.total<-c(increment.z[2:(l-1)], increment.z.total)
  
  
  
  for(k in 2:l){
    
    
    if(k>2) hurst.track[k-1]<-hurstExponent(data_i[1:k, ])
    
    displacement.track[k-1]<-displacement(data_i[1:k, ])
    
    squareDisplacement.track[k-1]<-squareDisplacement(data_i[1:k, ])
    
    outreachRatio.track[k-1]<-outreachRatio(data_i[1:k, ])
    
    straightness.track[k-1]<-straightness(data_i[1:k, ])
    
    max.displacement.ratio.track[k-1]<- displacementRatio(data_i[1:k, ])
    
    if(k>2) meanTurningAngle.track[k-1]<-meanTurningAngle(data_i[1:k, ])
    
    asphericity.track[k-1]<-asphericity(data_i[1:k, ])
    
    
  }
  
  mean.hurst.track.real[index.i[3:l]]<-hurst.track[2:(l-1)]
  
  mean.displacement.track.real[index.i[2:l]]<-displacement.track
  
  mean.squareDisplacement.track.real[index.i[2:l]]<-squareDisplacement.track
  
  mean.outreachRatio.track.real[index.i[2:l]]<-outreachRatio.track
  
  mean.straightness.track.real[index.i[2:l]]<-straightness.track
  
  mean.max.displacement.ratio.track.real[index.i[2:l]]<-max.displacement.ratio.track
  
  mean.meanTurningAngle.track.real[index.i[2:(l-1)]]<-meanTurningAngle.track[2:(l-1)]
  
  mean.asphericity.track.real[index.i[2:l]]<-asphericity.track
  
  
}


time.msd<-21




separated.data.based.on.time.asphericity<-matrix(0, nrow= time.msd, ncol=number.of.cells)

measure.greg.asphericity.data <-matrix(0, nrow=time.msd)



separated.data.based.on.time.max.displacement.ratio<-matrix(0, nrow= time.msd, ncol=number.of.cells)

measure.greg.max.displacement.ratio.data <-matrix(0, nrow=time.msd)




separated.data.based.on.time.hurst<-matrix(0, nrow= time.msd, ncol=number.of.cells)

measure.greg.hurst.data <-matrix(0, nrow=time.msd)




separated.data.based.on.time.displacement<-matrix(0, nrow= time.msd, ncol=number.of.cells)

measure.greg.displacement.data <-matrix(0, nrow=time.msd)





separated.data.based.on.time.squared.displacement<-matrix(0, nrow= time.msd, ncol=number.of.cells)

measure.greg.squared.displacement.data <-matrix(0, nrow=time.msd)




separated.data.based.on.time.straightness<-matrix(0, nrow= time.msd, ncol=number.of.cells)

measure.greg.straightness.data <-matrix(0, nrow=time.msd)



separated.data.based.on.time.outreachratio<-matrix(0, nrow= time.msd, ncol=number.of.cells)

measure.greg.outreachratio.data <-matrix(0, nrow=time.msd)



true.radius<-matrix(0, nrow= time.msd, ncol=number.of.cells)


data.sd.fence.displacement.track<-matrix(0, nrow= time.msd)

data.sd.fence.hurst.track<-matrix(0, nrow= time.msd)

data.sd.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

data.sd.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

data.sd.fence.straightness.track<-matrix(0, nrow= time.msd)

data.sd.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

data.sd.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

data.sd.fence.asphericity.track<-matrix(0, nrow= time.msd)







data.lower.fence.displacement.track<-matrix(0, nrow= time.msd)

data.lower.fence.hurst.track<-matrix(0, nrow= time.msd)

data.lower.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

data.lower.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

data.lower.fence.straightness.track<-matrix(0, nrow= time.msd)

data.lower.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

data.lower.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

data.lower.fence.asphericity.track<-matrix(0, nrow= time.msd)





data.upper.fence.displacement.track<-matrix(0, nrow= time.msd)

data.upper.fence.hurst.track<-matrix(0, nrow= time.msd)

data.upper.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

data.upper.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

data.upper.fence.straightness.track<-matrix(0, nrow= time.msd)

data.upper.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

data.upper.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

data.upper.fence.asphericity.track<-matrix(0, nrow= time.msd)







time.msd<-21


counter.time.data<-matrix(0, nrow=time.msd)




for(i in 1:number.of.cells){
  
  
  
  index.i<-which(data_cell_motility[, 8]==id[i])
  
  
  
  # # # # for(i in 1:length(cluster.1[cluster.1>249])){
  
  # # t<-cluster.1[cluster.1>249]-249; 
  
  # # index.i<-which(data_cell_motility[, 8]==id[t[i]])
  
  
  
  
  l<-length(index.i)
  
  for(k in 2: 21){
    
    # t<-0
    
    
    
    
    separated.data.based.on.time.displacement[k, i]<-mean.displacement.track.real[index.i[k]]
    
    measure.greg.displacement.data[k]<-measure.greg.displacement.data[k]+ mean.displacement.track.real[index.i[k]]
    
    
    
    separated.data.based.on.time.hurst[k, i]<-mean.hurst.track.real[index.i[k]]
    
    measure.greg.hurst.data[k]<-measure.greg.hurst.data[k]+ mean.hurst.track.real[index.i[k]]
    
    
    
    
    separated.data.based.on.time.squared.displacement[k, i]<-mean.squareDisplacement.track.real[index.i[k]]
    
    measure.greg.squared.displacement.data[k] <-measure.greg.squared.displacement.data[k]+ mean.squareDisplacement.track.real[index.i[k]]
    
    
    
    
    
    
    
    
    
    
    separated.data.based.on.time.straightness[k, i]<-mean.straightness.track.real[index.i[k]]
    
    measure.greg.straightness.data[k] <-measure.greg.straightness.data[k]+mean.straightness.track.real[index.i[k]]
    
    
    
    
    
    separated.data.based.on.time.asphericity[k, i]<-mean.asphericity.track.real[index.i[k]]
    
    measure.greg.asphericity.data[k] <-measure.greg.asphericity.data[k]+mean.asphericity.track.real[index.i[k]]
    
    
    
    
    
    
    
    separated.data.based.on.time.max.displacement.ratio[k, i]<-mean.max.displacement.ratio.track.real[index.i[k]]
    
    measure.greg.max.displacement.ratio.data[k] <-measure.greg.max.displacement.ratio.data[k]+ mean.max.displacement.ratio.track.real[index.i[k]]
    
    
    
    
    
    separated.data.based.on.time.outreachratio[k, i]<-mean.outreachRatio.track.real[index.i[k]]
    
    measure.greg.outreachratio.data[k] <-measure.greg.outreachratio.data[k]+ mean.outreachRatio.track.real[index.i[k]]
    
    
    counter.time.data[k]<-counter.time.data[k]+1
    
    
  }	
  
  
  
  # # # # # # 		             for(k in 30: time.msd){
  
  # # # if(l>=k){
  
  # # # # t<-0
  
  
  
  
  # # # separated.data.based.on.time.displacement[k, i]<-mean.displacement.track.real[index.i[k]]
  
  # # # measure.greg.displacement.data[k]<-measure.greg.displacement.data[k]+ mean.displacement.track.real[index.i[k]]
  
  
  
  
  
  
  
  
  # # # separated.data.based.on.time.squared.displacement[k, i]<-mean.squareDisplacement.track.real[index.i[k]]
  
  # # # measure.greg.squared.displacement.data[k] <-measure.greg.squared.displacement.data[k]+ mean.squareDisplacement.track.real[index.i[k]]
  
  
  
  
  
  
  
  # # # separated.data.based.on.time.straightness[k, i]<-mean.straightness.track.real[index.i[k]]
  
  # # # measure.greg.straightness.data[k] <-measure.greg.straightness.data[k]+mean.straightness.track.real[index.i[k]]
  
  
  
  
  
  # # # separated.data.based.on.time.outreachratio[k, i]<-mean.outreachRatio.track.real[index.i[k]]
  
  # # # measure.greg.outreachratio.data[k] <-measure.greg.outreachratio.data[k]+ mean.outreachRatio.track.real[index.i[k]]
  
  
  # # # counter.time.data[k]<-counter.time.data[k]+1
  
  
  # # # }             
  
  # # # }		
  
  
}  



for(k in 2:time.msd){
  
  
  
  t<-which(separated.data.based.on.time.displacement[k, ]==0)
  
  if(length(t)>0) values<-separated.data.based.on.time.displacement[k, -t]
  
  if(length(t)==0) values<-separated.data.based.on.time.displacement[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  data.lower.fence.displacement.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  data.upper.fence.displacement.track[k]<-summary(values)[5]+1.5*iqr
  
  data.sd.fence.displacement.track[k]<-sd(values)
  
  
  ########################################################################################################################## 
  
  
  
  
  t<-which(separated.data.based.on.time.hurst[k, ]==0)
  
  if(length(t)>0) values<-separated.data.based.on.time.hurst[k, -t]
  
  if(length(t)==0) values<-separated.data.based.on.time.hurst[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  data.lower.fence.hurst.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  data.upper.fence.hurst.track[k]<-summary(values)[5]+1.5*iqr
  
  data.sd.fence.hurst.track[k]<-sd(values)
  
  
  ########################################################################################################################## 
  
  
  
  
  
  
  t<-which(separated.data.based.on.time.max.displacement.ratio[k, ]==0)
  
  if(length(t)>0) values<-separated.data.based.on.time.max.displacement.ratio[k, -t]
  
  if(length(t)==0) values<-separated.data.based.on.time.max.displacement.ratio[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  data.lower.fence.max.displacement.ratio.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  data.upper.fence.max.displacement.ratio.track[k]<-summary(values)[5]+1.5*iqr
  
  data.sd.fence.max.displacement.ratio.track[k]<-sd(values)
  
  
  ##########################################################################################################################   
  
  
  
  
  
  t<-which(separated.data.based.on.time.asphericity[k, ]==0)
  
  if(length(t)>0) values<-separated.data.based.on.time.asphericity[k, -t]
  
  if(length(t)==0) values<-separated.data.based.on.time.asphericity[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  data.lower.fence.asphericity.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  data.upper.fence.asphericity.track[k]<-summary(values)[5]+1.5*iqr
  
  data.sd.fence.asphericity.track[k]<-sd(values)
  
  
  ##########################################################################################################################   
  
  
  
  
  
  
  
  
  t<-which(separated.data.based.on.time.squared.displacement[k, ]==0)
  
  if(length(t)>0) values<-separated.data.based.on.time.squared.displacement[k, -t]
  
  if(length(t)==0) values<-separated.data.based.on.time.squared.displacement[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  data.lower.fence.squareDisplacement.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  data.upper.fence.squareDisplacement.track[k]<-summary(values)[5]+1.5*iqr
  
  data.sd.fence.squareDisplacement.track[k]<-sd(values)
  
  
  
  ##########################################################################################################################
  
  
  
  
  
  t<-which(separated.data.based.on.time.straightness[k, ]==0)
  
  if(length(t)>0) values<-separated.data.based.on.time.straightness[k, -t]
  
  if(length(t)==0) values<-separated.data.based.on.time.straightness[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  data.lower.fence.straightness.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  data.upper.fence.straightness.track[k]<-summary(values)[5]+1.5*iqr
  
  data.sd.fence.straightness.track[k]<-sd(values)
  
  
  
  ##########################################################################################################################
  
  
  
  
  
  
  t<-which(separated.data.based.on.time.outreachratio[k, ]==0)
  
  if(length(t)>0) values<-separated.data.based.on.time.outreachratio[k, -t]
  
  if(length(t)==0) values<-separated.data.based.on.time.outreachratio[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  data.lower.fence.outreachRatio.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  data.upper.fence.outreachRatio.track[k]<-summary(values)[5]+1.5*iqr
  
  data.sd.fence.outreachRatio.track[k]<-sd(values)
  
  
}	 



time.msd<-21




M.G.D.asphericity<-measure.greg.asphericity.data[2: time.msd]/counter.time.data[2: time.msd]



M.G.D.max.displacement.ratio<-measure.greg.max.displacement.ratio.data[2: time.msd]/counter.time.data[2: time.msd]


M.G.D.hurst<-measure.greg.hurst.data[3: time.msd]/counter.time.data[3: time.msd]



M.G.D.displacement<-measure.greg.displacement.data[2: time.msd]/counter.time.data[2: time.msd]



M.G.D.squared.displacement<-measure.greg.squared.displacement.data[2: time.msd]/counter.time.data[2: time.msd]



M.G.D.outreachratio<-measure.greg.outreachratio.data[2: time.msd]/counter.time.data[2: time.msd]



M.G.D.straightness<-measure.greg.straightness.data[2: time.msd]/counter.time.data[2: time.msd]



##############################################





xlab.values<-seq(1, (time.msd-1), 1)


# plot(M.G.D.outreachratio, col="black",  xlab="Step", ylab="", main="Displacement Ratio", pch=20, cex.lab=2, cex.main=2.5)


# plot(M.G.D.max.displacement.ratio, col="black",  xlab="Step", ylab="", main="Displacement Ratio", pch=20, cex.lab=2, cex.main=2.5, ylim=c(0.5, 1))


# plot(M.G.D.asphericity, col="black",  xlab="Step", ylab="", main="Asphericity", pch=20, cex.lab=2, cex.main=2.5, ylim=c(0.5, 1))

# # # # axis(1, at=seq(0, 3, 3/19), labels= xlab.values)




# # # # arrows(seq(0, 3, 3/19), log(M.G.D.displacement-data.sd.fence.displacement.track[2: time.msd]/sqrt(counter.time.data[2: time.msd]))/log(10), seq(0, 3, 3/19), log(M.G.D.displacement+data.sd.fence.displacement.track[2: time.msd]/sqrt(counter.time.data[2: time.msd]))/log(10), length=0.05, angle=90, code=3, col="black")






# # # # # # # # # # # # 

# # # # par(mfrow=c(2, 1), lty=2)


# # # # 

setwd("/sbb/scratch/torkashe/Spherical-Coordinates-Motility-Metrics/Asphericity/")


ylab<-c(10, 100, 1000)


ylab.values<-log(ylab)/log(10)


x.lab<-seq(2, (time.msd-1), 1)


xlab.values<-log(x.lab)/log(10)


points(x.lab, M.G.D.asphericity, col="black", lty=6, lwd=1.5, pch=20, type="l", cex=1.5)


arrows(x.lab, M.G.D.asphericity-data.sd.fence.asphericity.track[3: time.msd]/sqrt(counter.time.data[3: time.msd]), x.lab, M.G.D.asphericity+ data.sd.fence.asphericity.track[3: time.msd]/sqrt(counter.time.data[3: time.msd]), length=0.05, angle=90, code=3, col="black")

axis(4, at= ylab.values, labels= ylab.values)



################################################################################################################################################################################################################################################


remove(list=ls())		

##########################################################################################################################    

setwd("/sbb/scratch/torkashe/Clustering-Results/no-clustering-spherical-modeling-NHMM")

load(".RMData")




simulated.mean.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)



simulated.lower.fence.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)




simulated.upper.fence.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)





simulated.mean.hurst.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.asphericity.track.observation<-matrix(0, nrow= number.of.positions.total)




simulated.lower.fence.hurst.track.observation<-matrix(0, nrow= number.of.positions.total)  	   

simulated.lower.fence.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.asphericity.observation<-matrix(0, nrow= number.of.positions.total)




simulated.upper.fence.hurst.track.observation<-matrix(0, nrow= number.of.positions.total) 	   

simulated.upper.fence.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.asphericity.observation<-matrix(0, nrow= number.of.positions.total)









for(i in 1: number.of.positions.total){
  
  
  
  
  simulated.mean.hurst.track.inside[i]<-mean(mean.hurst.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.hurst.track.inside[1:(r-1), i])[5]-summary(mean.hurst.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.hurst.track.inside[i]<-summary(mean.hurst.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.hurst.track.inside[i]<-summary(mean.hurst.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.displacement.track.inside[i]<-mean(mean.displacement.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.displacement.track.inside[1:(r-1), i])[5]-summary(mean.displacement.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.displacement.track.inside[i]<-summary(mean.displacement.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.displacement.track.inside[i]<-summary(mean.displacement.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.squareDisplacement.track.inside[i]<-mean(mean.squareDisplacement.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[5]-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.squareDisplacement.track.inside[i]<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.squareDisplacement.track.inside[i]<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  simulated.mean.outreachRatio.track.inside[i]<-mean(mean.outreachRatio.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[5]-summary(mean.outreachRatio.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.outreachRatio.track.inside[i]<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.outreachRatio.track.inside[i]<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.straightness.track.inside[i]<-mean(mean.straightness.track.inside[1:(r-1), i])
  
  
  
  
  
  iqr<-summary(mean.straightness.track.inside[1:(r-1), i])[5]-summary(mean.straightness.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.straightness.track.inside[i]<-summary(mean.straightness.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.straightness.track.inside[i]<-summary(mean.straightness.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  simulated.mean.max.displacement.ratio.track.inside[i]<-mean(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])
  
  
  
  
  iqr<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[5]-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.max.displacement.ratio.track.inside[i]<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.max.displacement.ratio.track.inside[i]<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  
  simulated.mean.meanTurningAngle.track.inside[i]<-mean(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[5]-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.meanTurningAngle.track.inside[i]<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.meanTurningAngle.track.inside[i]<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  simulated.mean.asphericity.track.inside[i]<-mean(mean.simulated.asphericity.inside[1:(r-1), i])
  
  
  
  
  iqr<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[5]-summary(mean.simulated.asphericity.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.asphericity.track.inside[i]<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.asphericity.track.inside[i]<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  ##############################################
  ##############################################    
  
  
  simulated.mean.hurst.track.observation[i]<-mean(mean.hurst.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.hurst.track.observation[1:(r-1), i])[5]-summary(mean.hurst.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.hurst.track.observation[i]<-summary(mean.hurst.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.hurst.track.observation[i]<-summary(mean.hurst.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  simulated.mean.displacement.track.observation[i]<-mean(mean.displacement.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.displacement.track.observation[1:(r-1), i])[5]-summary(mean.displacement.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.displacement.track.observation[i]<-summary(mean.displacement.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.displacement.track.observation[i]<-summary(mean.displacement.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.squareDisplacement.track.observation[i]<-mean(mean.squareDisplacement.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[5]-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.squareDisplacement.track.observation[i]<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.squareDisplacement.track.observation[i]<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  simulated.mean.outreachRatio.track.observation[i]<-mean(mean.outreachRatio.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[5]-summary(mean.outreachRatio.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.outreachRatio.track.observation[i]<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.outreachRatio.track.observation[i]<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.straightness.track.observation[i]<-mean(mean.straightness.track.observation[1:(r-1), i])
  
  
  
  
  
  iqr<-summary(mean.straightness.track.observation[1:(r-1), i])[5]-summary(mean.straightness.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.straightness.track.observation[i]<-summary(mean.straightness.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.straightness.track.observation[i]<-summary(mean.straightness.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  
  
  
  
  
}	   

###################################################################################################################################################################################################################################################




time.msd<-21


counter.time.simulated<-matrix(0, nrow=time.msd)


measure.greg.hurst.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.hurst<-matrix(0, nrow= time.msd, ncol=number.of.cells)



measure.greg.displacement.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.displacement<-matrix(0, nrow= time.msd, ncol=number.of.cells)


measure.greg.max.displacement.ratio.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.max.displacement.ratio<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.squared.displacement.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.squared.displacement<-matrix(0, nrow= time.msd, ncol=number.of.cells)



measure.greg.asphericity.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.asphericity<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.straightness.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.straightness<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.outreachratio.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.outreachratio<-matrix(0, nrow= time.msd, ncol=number.of.cells)



true.radius<-matrix(0, nrow= time.msd, ncol=number.of.cells)






simulated.lower.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.asphericity.track<-matrix(0, nrow= time.msd)




simulated.upper.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.asphericity.track<-matrix(0, nrow= time.msd)



simulated.sd.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.asphericity.track<-matrix(0, nrow= time.msd)





for(i in 1:number.of.cells){
  
  
  
  index.i<-which(data_cell_motility[, 8]==id[i])
  
  
  
  # # # # for(i in 1:length(cluster.1[cluster.1>249])){
  
  # # t<-cluster.1[cluster.1>249]-249; 
  
  # # index.i<-which(data_cell_motility[, 8]==id[t[i]])
  
  
  
  
  l<-length(index.i)
  
  for(k in 2: time.msd){
    
    
    
    measure.greg.asphericity.simulated[k]<-measure.greg.asphericity.simulated[k]+ simulated.mean.asphericity.track.inside[index.i[k]]
    
    separated.sim.based.on.time.asphericity[k, i]<-simulated.mean.asphericity.track.inside[index.i[k]]
    
    
    
    
    # t<-0
    
    
    measure.greg.displacement.simulated[k]<-measure.greg.displacement.simulated[k]+ simulated.mean.displacement.track.inside[index.i[k]]
    
    separated.sim.based.on.time.displacement[k, i]<-simulated.mean.displacement.track.inside[index.i[k]]
    
    
    
    measure.greg.hurst.simulated[k]<-measure.greg.hurst.simulated[k]+ simulated.mean.hurst.track.inside[index.i[k]]
    
    separated.sim.based.on.time.hurst[k, i]<-simulated.mean.hurst.track.inside[index.i[k]]
    
    
    
    
    measure.greg.max.displacement.ratio.simulated[k]<-measure.greg.max.displacement.ratio.simulated[k]+ simulated.mean.max.displacement.ratio.track.inside[index.i[k]]
    
    separated.sim.based.on.time.max.displacement.ratio[k, i]<-simulated.mean.max.displacement.ratio.track.inside[index.i[k]]
    
    
    
    
    
    measure.greg.squared.displacement.simulated[k] <-measure.greg.squared.displacement.simulated[k]+ simulated.mean.squareDisplacement.track.inside[index.i[k]]
    
    separated.sim.based.on.time.squared.displacement[k, i]<-simulated.mean.squareDisplacement.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    measure.greg.straightness.simulated[k] <-measure.greg.straightness.simulated[k] +simulated.mean.straightness.track.inside[index.i[k]]
    
    separated.sim.based.on.time.straightness[k, i]<-simulated.mean.straightness.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    
    measure.greg.outreachratio.simulated[k] <-measure.greg.outreachratio.simulated[k]+simulated.mean.outreachRatio.track.inside[index.i[k]]
    
    separated.sim.based.on.time.outreachratio[k, i]<-simulated.mean.outreachRatio.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    counter.time.simulated[k]<-counter.time.simulated[k]+1
    
    
  }	
  
  
  
  # # # # # # 		            for(k in 30: time.msd){
  
  # # # if(l>=k){
  
  # # # # t<-0
  
  
  # # # measure.greg.displacement.simulated[k]<-measure.greg.displacement.simulated[k]+ simulated.mean.displacement.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.displacement[k, i]<-simulated.mean.displacement.track.inside[index.i[k]]
  
  
  
  
  # # # measure.greg.hurst.simulated[k]<-measure.greg.hurst.simulated[k]+ simulated.mean.hurst.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.hurst[k, i]<-simulated.mean.hurst.track.inside[index.i[k]]
  
  
  
  
  
  # # # measure.greg.squared.displacement.simulated[k] <-measure.greg.squared.displacement.simulated[k]+ simulated.mean.squareDisplacement.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.squared.displacement[k, i]<-simulated.mean.squareDisplacement.track.inside[index.i[k]]
  
  
  
  
  
  
  # # # measure.greg.straightness.simulated[k] <-measure.greg.straightness.simulated[k] +simulated.mean.straightness.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.straightness[k, i]<-simulated.mean.straightness.track.inside[index.i[k]]
  
  
  
  
  
  
  
  # # # measure.greg.outreachratio.simulated[k] <-measure.greg.outreachratio.simulated[k]+simulated.mean.outreachRatio.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.outreachratio[k, i]<-simulated.mean.outreachRatio.track.inside[index.i[k]]
  
  
  # # # counter.time.simulated[k]<-counter.time.simulated[k]+1
  
  
  # # # }             
  
  # # # }		
  
  
}  




for(k in 2:time.msd){
  
  
  t<-which(separated.sim.based.on.time.asphericity[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.asphericity[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.asphericity[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.asphericity.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.asphericity.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.asphericity.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  t<-which(separated.sim.based.on.time.max.displacement.ratio[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.max.displacement.ratio[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.max.displacement.ratio[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.max.displacement.ratio.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.max.displacement.ratio.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.max.displacement.ratio.track[k]<-sd(values)
  
  
  
  ##########################################################################################################################
  
  
  
  
  
  
  t<-which(separated.sim.based.on.time.hurst[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.hurst[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.hurst[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.hurst.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.hurst.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.hurst.track[k]<-sd(values)
  
  
  
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  
  
  t<-which(separated.sim.based.on.time.displacement[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.displacement[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.displacement[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.displacement.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.displacement.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.displacement.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  t<-which(separated.sim.based.on.time.squared.displacement[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.squared.displacement[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.squared.displacement[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.squareDisplacement.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.squareDisplacement.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.squareDisplacement.track[k]<-sd(values)
  
  
  
  
  ##########################################################################################################################
  
  
  
  t<-which(separated.sim.based.on.time.straightness[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.straightness[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.straightness[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.straightness.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.straightness.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.straightness.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  t<-which(separated.sim.based.on.time.outreachratio[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.outreachratio[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.outreachratio[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.outreachRatio.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.outreachRatio.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.outreachRatio.track[k]<-sd(values)
  
  
  
  
  
  
  
}	  	   








time.msd<-21



M.G.S.hurst<-measure.greg.hurst.simulated[3: time.msd]/counter.time.simulated[2: time.msd]




M.G.S.diplacement<-measure.greg.displacement.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.max.diplacement.ratio<-measure.greg.max.displacement.ratio.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.asphericity<-measure.greg.asphericity.simulated[2: time.msd]/counter.time.simulated[2: time.msd]




M.G.S.squared.diplacement<-measure.greg.squared.displacement.simulated[2: time.msd]/counter.time.simulated[2: time.msd]




M.G.S.outreachratio<-measure.greg.outreachratio.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.straightness<-measure.greg.straightness.simulated[2: time.msd]/counter.time.simulated[2: time.msd]





##############################################


setwd("/sbb/scratch/torkashe/Spherical-Coordinates-Motility-Metrics/Asphericity/")


ylab<-c(10, 100, 1000)

ylab.values<-log(ylab)/log(10)

x.lab <-seq(2, (time.msd-1), 1)

xlab.values<-log(x.lab)/log(10)

# points(M.G.S.max.diplacement.ratio, col="violetred", pch=20, type="p")

points(x.lab, M.G.S.asphericity, col="darkblue", pch=20, type="l", cex=1.5, lwd=1.5, lty=4)

# points(xlab , M.G.S.asphericity, col="violetred", pch=20, type="p")



# arrows(seq(1, (time.msd-1), 1), M.G.S.diplacement-simulated.sd.fence.displacement.track[2: time.msd], seq(1, (time.msd-1), 1), M.G.S.diplacement+simulated.sd.fence.displacement.track[2: time.msd], length=0.05, angle=90, code=3, col="darkgreen")



# # # # u<-ecdf(M.G.S.diplacement)

# # plot(u, add=T, col="darkred") 




###################################################################################################################################################################################################################################################


remove(list=ls())		

##########################################################################################################################    

setwd("/sbb/scratch/torkashe/Clustering-Results/cluster-evaluation-2-states-NHMM-Greg")

load(".RMGammaData")




simulated.mean.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.mean.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)



simulated.lower.fence.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)




simulated.upper.fence.hurst.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.displacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.squareDisplacement.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.outreachRatio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.straightness.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.max.displacement.ratio.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.meanTurningAngle.track.inside<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.asphericity.track.inside<-matrix(0, nrow= number.of.positions.total)





simulated.mean.hurst.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.mean.asphericity.track.observation<-matrix(0, nrow= number.of.positions.total)




simulated.lower.fence.hurst.track.observation<-matrix(0, nrow= number.of.positions.total)  	   

simulated.lower.fence.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.lower.fence.asphericity.observation<-matrix(0, nrow= number.of.positions.total)




simulated.upper.fence.hurst.track.observation<-matrix(0, nrow= number.of.positions.total) 	   

simulated.upper.fence.displacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.squareDisplacement.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.outreachRatio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.straightness.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.max.displacement.ratio.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.meanTurningAngle.track.observation<-matrix(0, nrow= number.of.positions.total)

simulated.upper.fence.asphericity.observation<-matrix(0, nrow= number.of.positions.total)









for(i in 1: number.of.positions.total){
  
  
  
  
  simulated.mean.hurst.track.inside[i]<-mean(mean.hurst.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.hurst.track.inside[1:(r-1), i])[5]-summary(mean.hurst.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.hurst.track.inside[i]<-summary(mean.hurst.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.hurst.track.inside[i]<-summary(mean.hurst.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.displacement.track.inside[i]<-mean(mean.displacement.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.displacement.track.inside[1:(r-1), i])[5]-summary(mean.displacement.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.displacement.track.inside[i]<-summary(mean.displacement.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.displacement.track.inside[i]<-summary(mean.displacement.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.squareDisplacement.track.inside[i]<-mean(mean.squareDisplacement.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[5]-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.squareDisplacement.track.inside[i]<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.squareDisplacement.track.inside[i]<-summary(mean.squareDisplacement.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  simulated.mean.outreachRatio.track.inside[i]<-mean(mean.outreachRatio.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[5]-summary(mean.outreachRatio.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.outreachRatio.track.inside[i]<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.outreachRatio.track.inside[i]<-summary(mean.outreachRatio.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.straightness.track.inside[i]<-mean(mean.straightness.track.inside[1:(r-1), i])
  
  
  
  
  
  iqr<-summary(mean.straightness.track.inside[1:(r-1), i])[5]-summary(mean.straightness.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.straightness.track.inside[i]<-summary(mean.straightness.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.straightness.track.inside[i]<-summary(mean.straightness.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  simulated.mean.max.displacement.ratio.track.inside[i]<-mean(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])
  
  
  
  
  iqr<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[5]-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.max.displacement.ratio.track.inside[i]<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.max.displacement.ratio.track.inside[i]<-summary(mean.simulated.max.displacepent.ratio.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  
  simulated.mean.meanTurningAngle.track.inside[i]<-mean(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])
  
  
  
  iqr<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[5]-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.meanTurningAngle.track.inside[i]<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.meanTurningAngle.track.inside[i]<-summary(mean.simulated.meanTurningAngle.track.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  simulated.mean.asphericity.track.inside[i]<-mean(mean.simulated.asphericity.inside[1:(r-1), i])
  
  
  
  
  iqr<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[5]-summary(mean.simulated.asphericity.inside[1:(r-1), i])[2]
  
  simulated.lower.fence.asphericity.track.inside[i]<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.asphericity.track.inside[i]<-summary(mean.simulated.asphericity.inside[1:(r-1), i])[5]+1.5*iqr
  
  
  
  ##############################################
  ##############################################    
  
  
  simulated.mean.hurst.track.observation[i]<-mean(mean.hurst.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.hurst.track.observation[1:(r-1), i])[5]-summary(mean.hurst.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.hurst.track.observation[i]<-summary(mean.hurst.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.hurst.track.observation[i]<-summary(mean.hurst.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  simulated.mean.displacement.track.observation[i]<-mean(mean.displacement.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.displacement.track.observation[1:(r-1), i])[5]-summary(mean.displacement.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.displacement.track.observation[i]<-summary(mean.displacement.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.displacement.track.observation[i]<-summary(mean.displacement.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.squareDisplacement.track.observation[i]<-mean(mean.squareDisplacement.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[5]-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.squareDisplacement.track.observation[i]<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.squareDisplacement.track.observation[i]<-summary(mean.squareDisplacement.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  simulated.mean.outreachRatio.track.observation[i]<-mean(mean.outreachRatio.track.observation[1:(r-1), i])
  
  
  
  iqr<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[5]-summary(mean.outreachRatio.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.outreachRatio.track.observation[i]<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.outreachRatio.track.observation[i]<-summary(mean.outreachRatio.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  simulated.mean.straightness.track.observation[i]<-mean(mean.straightness.track.observation[1:(r-1), i])
  
  
  
  
  
  iqr<-summary(mean.straightness.track.observation[1:(r-1), i])[5]-summary(mean.straightness.track.observation[1:(r-1), i])[2]
  
  simulated.lower.fence.straightness.track.observation[i]<-summary(mean.straightness.track.observation[1:(r-1), i])[2]-1.5*iqr
  
  simulated.upper.fence.straightness.track.observation[i]<-summary(mean.straightness.track.observation[1:(r-1), i])[5]+1.5*iqr
  
  
  
  
  
  
  
  
  
  
  
}	   

###################################################################################################################################################################################################################################################




time.msd<-21


counter.time.simulated<-matrix(0, nrow=time.msd)


measure.greg.hurst.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.hurst<-matrix(0, nrow= time.msd, ncol=number.of.cells)



measure.greg.displacement.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.displacement<-matrix(0, nrow= time.msd, ncol=number.of.cells)


measure.greg.max.displacement.ratio.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.max.displacement.ratio<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.squared.displacement.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.squared.displacement<-matrix(0, nrow= time.msd, ncol=number.of.cells)



measure.greg.asphericity.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.asphericity<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.straightness.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.straightness<-matrix(0, nrow= time.msd, ncol=number.of.cells)




measure.greg.outreachratio.simulated <-matrix(0, nrow=time.msd)

separated.sim.based.on.time.outreachratio<-matrix(0, nrow= time.msd, ncol=number.of.cells)



true.radius<-matrix(0, nrow= time.msd, ncol=number.of.cells)






simulated.lower.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.lower.fence.asphericity.track<-matrix(0, nrow= time.msd)




simulated.upper.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.upper.fence.asphericity.track<-matrix(0, nrow= time.msd)



simulated.sd.fence.displacement.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.hurst.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.squareDisplacement.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.outreachRatio.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.straightness.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.max.displacement.ratio.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.meanTurningAngle.track<-matrix(0, nrow= time.msd)

simulated.sd.fence.asphericity.track<-matrix(0, nrow= time.msd)





for(i in 1:number.of.cells){
  
  
  
  index.i<-which(data_cell_motility[, 8]==id[i])
  
  
  
  # # # # for(i in 1:length(cluster.1[cluster.1>249])){
  
  # # t<-cluster.1[cluster.1>249]-249; 
  
  # # index.i<-which(data_cell_motility[, 8]==id[t[i]])
  
  
  
  
  l<-length(index.i)
  
  for(k in 2: time.msd){
    
    
    
    measure.greg.asphericity.simulated[k]<-measure.greg.asphericity.simulated[k]+ simulated.mean.asphericity.track.inside[index.i[k]]
    
    separated.sim.based.on.time.asphericity[k, i]<-simulated.mean.asphericity.track.inside[index.i[k]]
    
    
    
    
    # t<-0
    
    
    measure.greg.displacement.simulated[k]<-measure.greg.displacement.simulated[k]+ simulated.mean.displacement.track.inside[index.i[k]]
    
    separated.sim.based.on.time.displacement[k, i]<-simulated.mean.displacement.track.inside[index.i[k]]
    
    
    
    measure.greg.hurst.simulated[k]<-measure.greg.hurst.simulated[k]+ simulated.mean.hurst.track.inside[index.i[k]]
    
    separated.sim.based.on.time.hurst[k, i]<-simulated.mean.hurst.track.inside[index.i[k]]
    
    
    
    
    measure.greg.max.displacement.ratio.simulated[k]<-measure.greg.max.displacement.ratio.simulated[k]+ simulated.mean.max.displacement.ratio.track.inside[index.i[k]]
    
    separated.sim.based.on.time.max.displacement.ratio[k, i]<-simulated.mean.max.displacement.ratio.track.inside[index.i[k]]
    
    
    
    
    
    measure.greg.squared.displacement.simulated[k] <-measure.greg.squared.displacement.simulated[k]+ simulated.mean.squareDisplacement.track.inside[index.i[k]]
    
    separated.sim.based.on.time.squared.displacement[k, i]<-simulated.mean.squareDisplacement.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    measure.greg.straightness.simulated[k] <-measure.greg.straightness.simulated[k] +simulated.mean.straightness.track.inside[index.i[k]]
    
    separated.sim.based.on.time.straightness[k, i]<-simulated.mean.straightness.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    
    measure.greg.outreachratio.simulated[k] <-measure.greg.outreachratio.simulated[k]+simulated.mean.outreachRatio.track.inside[index.i[k]]
    
    separated.sim.based.on.time.outreachratio[k, i]<-simulated.mean.outreachRatio.track.inside[index.i[k]]
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    counter.time.simulated[k]<-counter.time.simulated[k]+1
    
    
  }	
  
  
  
  # # # # # # 		            for(k in 30: time.msd){
  
  # # # if(l>=k){
  
  # # # # t<-0
  
  
  # # # measure.greg.displacement.simulated[k]<-measure.greg.displacement.simulated[k]+ simulated.mean.displacement.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.displacement[k, i]<-simulated.mean.displacement.track.inside[index.i[k]]
  
  
  
  
  # # # measure.greg.hurst.simulated[k]<-measure.greg.hurst.simulated[k]+ simulated.mean.hurst.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.hurst[k, i]<-simulated.mean.hurst.track.inside[index.i[k]]
  
  
  
  
  
  # # # measure.greg.squared.displacement.simulated[k] <-measure.greg.squared.displacement.simulated[k]+ simulated.mean.squareDisplacement.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.squared.displacement[k, i]<-simulated.mean.squareDisplacement.track.inside[index.i[k]]
  
  
  
  
  
  
  # # # measure.greg.straightness.simulated[k] <-measure.greg.straightness.simulated[k] +simulated.mean.straightness.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.straightness[k, i]<-simulated.mean.straightness.track.inside[index.i[k]]
  
  
  
  
  
  
  
  # # # measure.greg.outreachratio.simulated[k] <-measure.greg.outreachratio.simulated[k]+simulated.mean.outreachRatio.track.inside[index.i[k]]
  
  # # # separated.sim.based.on.time.outreachratio[k, i]<-simulated.mean.outreachRatio.track.inside[index.i[k]]
  
  
  # # # counter.time.simulated[k]<-counter.time.simulated[k]+1
  
  
  # # # }             
  
  # # # }		
  
  
}  




for(k in 2:time.msd){
  
  
  t<-which(separated.sim.based.on.time.asphericity[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.asphericity[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.asphericity[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.asphericity.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.asphericity.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.asphericity.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  t<-which(separated.sim.based.on.time.max.displacement.ratio[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.max.displacement.ratio[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.max.displacement.ratio[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.max.displacement.ratio.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.max.displacement.ratio.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.max.displacement.ratio.track[k]<-sd(values)
  
  
  
  ##########################################################################################################################
  
  
  
  
  
  
  t<-which(separated.sim.based.on.time.hurst[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.hurst[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.hurst[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.hurst.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.hurst.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.hurst.track[k]<-sd(values)
  
  
  
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  
  
  t<-which(separated.sim.based.on.time.displacement[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.displacement[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.displacement[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.displacement.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.displacement.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.displacement.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  t<-which(separated.sim.based.on.time.squared.displacement[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.squared.displacement[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.squared.displacement[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.squareDisplacement.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.squareDisplacement.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.squareDisplacement.track[k]<-sd(values)
  
  
  
  
  ##########################################################################################################################
  
  
  
  t<-which(separated.sim.based.on.time.straightness[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.straightness[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.straightness[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.straightness.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.straightness.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.straightness.track[k]<-sd(values)
  
  
  
  
  
  ##########################################################################################################################
  
  
  
  
  t<-which(separated.sim.based.on.time.outreachratio[k, ]==0)
  
  if(length(t)>0) values<-separated.sim.based.on.time.outreachratio[k, -t]
  
  if(length(t)==0) values<-separated.sim.based.on.time.outreachratio[k, ]
  
  
  
  iqr<-summary(values)[5]-summary(values)[2]
  
  simulated.lower.fence.outreachRatio.track[k]<-max(c(0, summary(values)[2]-1.5*iqr))
  
  simulated.upper.fence.outreachRatio.track[k]<-summary(values)[5]+1.5*iqr
  
  simulated.sd.fence.outreachRatio.track[k]<-sd(values)
  
  
  
  
  
  
  
}	  	   







time.msd<-21



M.G.S.hurst<-measure.greg.hurst.simulated[3: time.msd]/counter.time.simulated[3: time.msd]




M.G.S.diplacement<-measure.greg.displacement.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.max.diplacement.ratio<-measure.greg.max.displacement.ratio.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.asphericity<-measure.greg.asphericity.simulated[2: time.msd]/counter.time.simulated[2: time.msd]




M.G.S.squared.diplacement<-measure.greg.squared.displacement.simulated[2: time.msd]/counter.time.simulated[2: time.msd]




M.G.S.outreachratio<-measure.greg.outreachratio.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



M.G.S.straightness<-measure.greg.straightness.simulated[2: time.msd]/counter.time.simulated[2: time.msd]



##############################################

setwd("/sbb/scratch/torkashe/Spherical-Coordinates-Motility-Metrics/Asphericity/")



ylab<-c(10, 100, 1000)

ylab.values<-log(ylab)/log(10)

x.lab<-seq(2, (time.msd-1), 1)

xlab.values<-log(x.lab)/log(10)

# points(M.G.S.max.diplacement.ratio, col="violetred", pch=20, type="p")

points(x.lab, M.G.S.asphericity, col="black", pch=18, type="l", cex=1.5, lwd=1.5)

# points(xlab , M.G.S.asphericity, col="violetred", pch=20, type="p")



# arrows(seq(1, (time.msd-1), 1), M.G.S.diplacement-simulated.sd.fence.displacement.track[2: time.msd], seq(1, (time.msd-1), 1), M.G.S.diplacement+simulated.sd.fence.displacement.track[2: time.msd], length=0.05, angle=90, code=3, col="darkgreen")



# # # # u<-ecdf(M.G.S.diplacement)

# # plot(u, add=T, col="darkred") 






################################################################################################################################################################################################################################################




setwd("/sbb/scratch/torkashe/Spherical-Coordinates-Motility-Metrics/Asphericity")

dev.print(pdf, 'Rplot-Asphericity.pdf', width=10, height=6, paper='special')


dev.copy(jpeg, filename="Rplot-Asphericity.jpeg");

dev.off()







