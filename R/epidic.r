

epidic<- function(burnin,niter,LLchain,LLpostmean){
    
    loglike <- LLchain[burnin:niter]
    D.bar <- -2*mean(loglike)

    D.hat <- -2* LLpostmean

    pD <- D.bar - D.hat
    dic <- pD+D.bar

return(dic=dic)
}

