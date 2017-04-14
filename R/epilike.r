

epilike<- function(type,x=NULL,y=NULL,inftime,infperiod=NULL,tmin=NULL,tmax,alpha,beta,
           spark=NULL,Sformula=NULL,contact=NULL){
    
####### checks  #######
 
 tau <- inftime
 lambda <- infperiod
 
    if (is.null(type) || !(type %in% c("SI", "SIR"))) {
        stop("elike: Specify type as \"SI\" or \"SIR\" ", call. = FALSE)
    }
    
    ns <- length(alpha)
    ni <- length(beta)
    
    if (!is.vector(tau))
    stop('epilike: tau is not a vector')
    n <- length(tau)
    
    if (!is.null(x)){
        if ((length(y) !=n) || (length(x) !=n))
        stop('epilike: Length of x/y is not compatible ')
    }
    
    if (is.null(lambda) && type=="SIR") {
        stop(' epilike: Specify removal distance,lambda ')
    }
    if (!is.null(lambda)){
        if (length(lambda) != n)
        stop('epilike: Length of lambda is not compatible.')
        if(type=="SI")
        stop('epilike: Type must be "SIR".')
    }
    
    if (is.null(tmin)){
        tmin=1
    }
    if (is.null(spark)){
        spark=0
    }
    
#### formula for susceptibility function

    if(!is.null(Sformula)){
        covmat <- model.matrix(Sformula)
        
        if ((ncol(covmat) == length(all.vars(Sformula))) & (ns != length(all.vars(Sformula))))
        stop('epilike: Check Sformula (no intercept term) and the dimension of alpha')
        
        if ((ncol(covmat) > length(all.vars(Sformula))) & (ns != ncol(covmat)))
        stop('epilike: Check Sformula (intercept term) and the dimension of alpha')
    }
    else{
        if (ns==1)
        covmat=matrix(1.0,nrow=n,ncol=ns)
        if (ns >1)
        stop('epilike: Please specify covariate')
    }
    

    val=0.0
    
#### Purely Spatial #####

    if(type=="SI"){
        tmp1 <- .Fortran("like",
        x=as.numeric(x),y=as.numeric(y),tau=as.integer(tau),
        n=as.integer(n),tmin=as.integer(tmin),tmax=as.integer(tmax),
        ns=as.integer(ns),ni=as.integer(ni),
        alpha=as.double(alpha), beta=as.double(beta),spark=as.double(spark),
        covmat= as.vector(covmat),val=as.double(val))
    }
    if(type=="SIR"){
        tmp1 <- .Fortran("likesir",
        x=as.numeric(x),y=as.numeric(y),tau=as.integer(tau),
        lambda=as.integer(lambda),n=as.integer(n),tmin=as.integer(tmin),
        tmax=as.integer(tmax),
        ns=as.integer(ns),ni=as.integer(ni),
        alpha=as.double(alpha), beta=as.double(beta),spark=as.double(spark),
        covmat= as.vector(covmat),val=as.double(val))
    }
    
#### Contact networks #####
    
    
    if (!is.null(contact)){
        if (length(contact)/(n*n) !=ni)
        stop('epilike:  Dimension of beta  and the number of contact networks are not matching')
        network <- array(contact,c(n,n,ni))
    }
    
    
    if ((type=="SI") && !is.null(contact)) {
        tmp1 <- .Fortran("likecon",
        tau=as.integer(tau),n = as.integer(n),
        ns= as.integer(ns),ni=as.integer(ni),
        tmin=as.integer(tmin),tmax = as.integer(tmax),
        alpha=as.numeric(alpha),beta=as.numeric(beta),spark=as.double(spark),
        covmat=as.vector(covmat),network=as.vector(network),
        val=as.double(val) )
    }
    
    if ((type=="SIR") && !is.null(contact)) {
        tmp1 <- .Fortran("likeconsir",
        tau=as.integer(tau),lambda=as.integer(lambda),
        n = as.integer(n),ns= as.integer(ns),
        ni=as.integer(ni),tmin=as.integer(tmin),tmax = as.integer(tmax),
        alpha=as.numeric(alpha),beta=as.numeric(beta),spark=as.double(spark),
        covmat=as.vector(covmat),network=as.vector(network),
        val=as.double(val) )
    }

    
    return(tmp1$val)
}

