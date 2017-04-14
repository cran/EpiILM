epiBR0<-  function(x=NULL,y=NULL,contact=NULL,alpha,beta,spark=NULL,infperiod,Sformula=NULL,tmax,niter){
          
####### checks  #######
   lambda <- infperiod

   ns <- length(alpha)
   ni <- length(beta)
   

   if (is.null(spark)){
       spark=0
   }
   
   if (!is.null(x)){
       n <- length(x)
       if ((length(y) !=n) || (length(x) !=n))
       stop('epiBR0: Length of x or y is not compatible ')
   }
   
   if(!is.null(contact)){
       n <- dim(contact)[1]
   }
   
   if (length(lambda) != n){
        stop('epiBR0: Length of lambda is not compatible.')
   }
     
#### formula for susceptibility function

    if(!is.null(Sformula)){
        covmat <- model.matrix(Sformula)
        
        if ((ncol(covmat) == length(all.vars(Sformula))) & (ns != length(all.vars(Sformula))))
        stop('epiBR0: Check Sformula (no intercept term) and the dimension of alpha')
        
        if ((ncol(covmat) > length(all.vars(Sformula))) & (ns != ncol(covmat)))
        stop('epiBR0: Check Sformula (intercept term) and the dimension of alpha')
        
    }
    else{
        if (ns==1)
        covmat=matrix(1.0,nrow=n,ncol=ns)
        if (ns >1)
        stop('epiBR0: Please specify covariate')
    }

    if (!is.null(contact)){
        if (length(contact)/(n*n) !=ni)
        stop('epiBR0:  Dimension of beta  and the number of contact networks are not matching')
        network <- array(contact,c(n,n,ni))
    }


    val=0.0

    if (!is.null(contact)) {
        tmp <- .Fortran("rconsir",
        n = as.integer(n),tmax = as.integer(tmax),
        ns= as.integer(ns),ni=as.integer(ni),
        lambda=as.integer(lambda),
        alpha=as.numeric(alpha),beta=as.numeric(beta),spark=as.numeric(spark),
        covmat=as.vector(covmat),network=as.vector(network),
        sim=as.integer(niter),val= as.double(val),countinf=as.integer(rep(0,niter)))
       
	result1 <- list(BasicR0 =tmp$val, simulated_BR0 = tmp$countinf)
    }
    else {
        tmp <- .Fortran("rxysir",
        n = as.integer(n),tmax = as.integer(tmax),
        ns= as.integer(ns),ni=as.integer(ni),
        alpha=as.numeric(alpha), beta=as.numeric(beta),spark=as.numeric(spark),
        covmat=as.vector(covmat),
        lambda=as.integer(lambda),x=as.double(x),y=as.double(y),
        sim=as.integer(niter),val=as.double(val),countinf=as.integer(rep(0,niter)))
        
	result1 <- list(BasicR0 =tmp$val, simulated_BR0 = tmp$countinf)

       }
    return(result1)
}



