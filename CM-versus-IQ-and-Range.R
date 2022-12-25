####  Testing Compactness Measure   ####


       remove(list=ls())


       sigma.test<-c(0.1, 0.5, 0.8, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5)
       
       
       df<-c(1, 2, 3, 4, 5, 6)


       B<-B.boot<-10000
       
       n<-c(30, 35, 40, 45, 50, 55, 60)
       
       R.sim<-150
       
       
       range.classifier <- matrix(0, nrow=R.sim, ncol=length(n))
       
       inter.quantile.classifier <- matrix(0, nrow=R.sim, ncol=length(n))

       cm.classifier <- matrix(0, nrow=R.sim, ncol=length(n))
       
       
      
    for(i in 1: length(df)){   	
    
    
    # for(i in 1: length(sigma.test)){  
       
       
       for(k in 1: length(n)){
       	
       	# for(k in 1: 2){
       
			     
			      for(r in 1: R.sim){
			      	
			       	
			       	
			       	   t.sample <- rt(n[k], df[i])
			       	   
			       	   # t.sample <- rnorm(n[k], 0, sigma.test[k])
			       	   
			       	   
			       	   summary.t.sample<-summary(t.sample)
			       	
			       	
			       	   normal.samples = matrix( rnorm(n[k]*B, mean(t.sample), sd(t.sample)), B, n[k])
			       	   
			       	   
			       	   normal.statistics = apply(normal.samples, 1, summary)
			       	   
			       	   
			       	   first.quantile<-normal.statistics[2, ]
                    
                    
	                   third.quantile<-normal.statistics[5, ]
	                    
	                   
	                   min.boot<-normal.statistics[1, ]
	                    
	                   
	                   max.boot<-normal.statistics[6, ]
	                    
	                   
	                   cm_normal<-(third.quantile-first.quantile)/(max.boot-min.boot)
	                   
	                   
	                   range_normal<-max.boot-min.boot
	                   
	                   
	                   iq_normal<-third.quantile-first.quantile
	                   
	                   
	                   interval.iq = quantile(iq_normal, c(0.025, 0.975))
	                   
	                   
	                   iq.t.sample = summary.t.sample[5]-summary.t.sample[2]
	                   
	                   
	                   interval.range = quantile(range_normal, c(0.025, 0.975))	 

	                   
	                   range.t.sample = summary.t.sample[6]-summary.t.sample[1]
	                   
	                   
	                   interval.cm = quantile(cm_normal, c(0.025, 0.975))
	                   
	                   
	                   cm.t.sample = iq.t.sample/range.t.sample
	                   
	                   
	                   
	                   if(iq.t.sample > interval.iq[2]){
	                   	
	                   		                   	     
	                   	     inter.quantile.classifier[r, k] <- 1
	                   	
	                   	
	                   }
	                   
	                   
	                   
	                   
	                   if(range.t.sample > interval.range[2]){
	                   	
	                   	
	                   	     range.classifier[r, k] <- 1
	                   	
	                   	
	                   }                  
	                   
	                   
	                   
	                   
	                   if(cm.t.sample < interval.cm[1]){
	                   	
	                   	
	                   	     cm.classifier[r, k]<-1
	                   	
	                   }
	                   
			       	   
			       	
			       	
			       	
			       	
			       }


       
       }
       



       average.range.classifier <- matrix(0, nrow=length(n))
       
       average.inter.quantile.classifier <- matrix(0, nrow=length(n))

       average.cm.classifier <- matrix(0, nrow=length(n))
       
       
       
       for(k in 1:length(n)){
       	
       	
       	average.range.classifier[k] <- mean(range.classifier[, k])
       	
       	
       	average.inter.quantile.classifier[k] <- mean(inter.quantile.classifier[, k])
       	
       	
       	average.cm.classifier[k] <- mean(cm.classifier[, k])
       	
       	
       }
       
       
       print("############################  df  ###############################")
       
       
       print(df)
       
       
       print("###################################################################")

       
       
       print("average.range.classifier")
       
       
       print(average.range.classifier)
       
       
       
       print("average.inter.quantile.classifier")
       
       
       print(average.inter.quantile.classifier)
       
       
       
       print("average.cm.classifier")
       
       print(average.cm.classifier)
       
       
       print("###################################################################")
       
       
       
       
       
       }
       
       
       
       

       