

epidata<-  function(type,n,tmin=NULL,tmax,alpha,beta,spark=NULL,
Sformula=NULL,x=NULL,y=NULL,inftime=NULL,infperiod=NULL,contact=NULL){
          
####### checks  #######

   tau <- inftime
   lambda <- infperiod

   if (is.null(type) || !(type %in% c("SI", "SIR"))) {
           stop("Specify type as \"SI\" or \"SIR\" ", call. = FALSE)
       }

   ns <- length(alpha)
   ni <- length(beta)
   
   if (is.null(tmin)){
       tmin=1
   }
   if (is.null(spark)){
       spark=0
   }
   
   
   if (!is.null(x)){
       if ((length(y) !=n) || (length(x) !=n))
       stop('epidata: Length of x or y is not compatible ')
   }
   
   A <- as.integer((n-1)*runif(1)+1)
          
   if (!is.null(tau)){
       if ((length(tau) !=n))
       stop('epidata: Length of tau is not compatible ')
   }
   else{
##### Initialize infectious state
       tau <- rep(0,n)
       tau[A] <- tmin
   }

    if (is.null(lambda) && type=="SIR") {
        stop(' epidata: Specify removal distance,lambda ')
    }
    if (!is.null(lambda)){
        if (length(lambda) != n)
        stop('epidata: Length of lambda is not compatible.')
        if(type=="SI")
        stop('epidata: Type must be "SIR".')
        
        remt <- rep(0,n)
        remt[A] <- tau[A]+lambda[A]
        
    }
    
#### formula for susceptibility function

    if(!is.null(Sformula)){
        covmat <- model.matrix(Sformula)
        
        if ((ncol(covmat) == length(all.vars(Sformula))) & (ns != length(all.vars(Sformula))))
        stop('epidata: Check Sformula (no intercept term) and the dimension of alpha')
        
        if ((ncol(covmat) > length(all.vars(Sformula))) & (ns != ncol(covmat)))
        stop('epidata: Check Sformula (intercept term) and the dimension of alpha')
        
    }
    else{
        if (ns==1)
        covmat=matrix(1.0,nrow=n,ncol=ns)
        if (ns >1)
        stop('epidata: Please specify covariate')
    }

    
#### Purely Spatial #####
            
    if (type=="SI"){
        
        tmp <- .Fortran("dataxy",
        x=as.double(x),y=as.double(y),
        n = as.integer(n),tmin = as.integer(tmin),tmax = as.integer(tmax),
        ns= as.integer(ns),ni=as.integer(ni),
        alpha=as.numeric(alpha),beta=as.numeric(beta),spark=as.numeric(spark),
        covmat=as.vector(covmat),
        tau=as.integer(tau) )
        
        result1 <- list(inftime=tmp$tau)
    }

    if (type=="SIR"){
        
        tmp <- .Fortran("dataxysir",
        n = as.integer(n),tmin = as.integer(tmin),tmax = as.integer(tmax),
        ns= as.integer(ns),ni=as.integer(ni),
        alpha=as.numeric(alpha), beta=as.numeric(beta),spark=as.numeric(spark),
        covmat=as.vector(covmat),
        lambda=as.integer(lambda),x=as.double(x),y=as.double(y),
        tau=as.integer(tau), remt=as.integer(remt))
        
        result1 <- list(inftime=tmp$tau, removaltime=tmp$remt )
        
    }

#### Contact networks #####

    if (!is.null(contact)){
        if (length(contact)/(n*n) !=ni)
        stop('epidata:  Dimension of beta  and the number of contact networks are not matching')
        network <- array(contact,c(n,n,ni))
    }

    if ((type=="SI") && !is.null(contact)) {

        tmp <- .Fortran("datacon",
        n = as.integer(n),tmin = as.integer(tmin),tmax = as.integer(tmax),
        ns= as.integer(ns),ni=as.integer(ni),
        alpha=as.numeric(alpha),beta=as.numeric(beta),spark=as.numeric(spark),
        covmat=as.vector(covmat),network=as.vector(network),
        tau=as.integer(tau) )

        result1 <- list(inftime=tmp$tau)
    }


    if ((type=="SIR") && !is.null(contact)) {
        
        tmp <- .Fortran("dataconsir",
        n = as.integer(n),tmin = as.integer(tmin),tmax = as.integer(tmax),
        ns= as.integer(ns),ni=as.integer(ni),
        lambda=as.integer(lambda),
        alpha=as.numeric(alpha),beta=as.numeric(beta),spark=as.numeric(spark),
        covmat=as.vector(covmat),network=as.vector(network),
        tau=as.integer(tau),remt=as.integer(remt) )
        
        result1 <- list(inftime=tmp$tau, removaltime=tmp$remt)
    }

    return(result1)
    }



