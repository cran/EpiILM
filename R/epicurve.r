

epicurve<- function(type,plottype,inftime,removaltime=NULL,tmin=NULL,timepoints=NULL){
    
tau <- inftime
if (is.null(type) || !(type %in% c("SI", "SIR"))) {
        stop("epicurve: Specify type as \"SI\" or \"SIR\" ", call. = FALSE)
    }

if (is.null(plottype) || !(plottype %in% c("complete","susceptible","totalinfect","newinfect"))) {
    stop("epicurve: Specify plottype as \"complete\" , \"susceptible\",\"totalinfect\" or \"newinfect\"", call. = FALSE)
}

n <- length(tau)
 
if (is.null(removaltime) && type=="SIR") {
     stop(' epicurve: Specify removal time ')
 }
if (!is.null(removaltime)){
     if (length(removaltime) != n)
     stop('epicurve: Length of removaltime is not compatible.')
     if(type=="SI")
     stop('epicurve: Type must be "SIR".')
 }

if (is.null(tmin)){
    tmin=1
}


totalinf <- rep(0)
sus <- rep(0)
newinf <- rep(0)
removed <- rep(0)

tmax <- max(tau)
time <- rep(tmin:tmax)

if (type=="SI"){

    for (i in tmin:tmax){
        newinf[i] <- length(tau[tau==i])
        xc <- subset(tau,tau<=i & tau !=0)
        totalinf[i] <- length(xc)
         sus[i] <- n-totalinf[i]
    }
    
    if (tmin>1){
        newinf <- newinf[tmin:tmax]
        totalinf <- totalinf[tmin:tmax]
        sus <- sus[tmin:tmax]
    }
    
    if (plottype=="complete"){
        plot (time,sus,xlim=c(tmin, tmax),ylim=c(1,n),main="Epidemic Curves", ylab="Number of individuals ",xlab="time",type="l",lwd=2,cex=0.5,xaxt="n")
        lines(time, totalinf, col="red",lwd=2)
        axis(1, at=1:tmax)
        legend("topright",inset=.001,cex=0.8,bty = "n",
        c("Infected","Susceptible"),
        horiz=TRUE,lty=c(1,1),lwd=c(2,2),col=c("red","black"))
        
    }
    if ((plottype=="complete") & (!is.null(timepoints))){
        plot (time,sus,xlim=c(timepoints[1], timepoints[2]),ylim=c(1,n),main="Epidemic Curves", ylab="Number of individuals",xlab="time",type="l",lwd=2,cex=0.5,xaxt="n")
        lines(time, totalinf, col="red",lwd=2)
        axis(1, at=1:timepoints[2])
        legend("topright",inset=.001,cex=0.8,bty = "n",
        c("Infected","Susceptible"),
        horiz=TRUE,lty=c(1,1),lwd=c(2,2),col=c("red","black"))
        
    }
}

if (type=="SIR"){
    
    dat <- data.frame(tau,removaltime)
    
    for (i in tmin:tmax){
        xcc <- subset(dat,tau <=i & tau !=0 & i <removaltime)
        totalinf[i] <- length(xcc$tau)
    }
    for (i in tmin:tmax){
        newinf[i] <- length(tau[tau==i])
    }
    for (i in tmin:tmax){
        xcc <- subset(dat, tau<=i & tau !=0)
        xc <- subset(xcc, i>= removaltime)
        removed[i] <- length(xc$tau)
        sus[i] <- n- length(xcc$tau)
    }
    
    if (tmin>1){
        newinf <- newinf[tmin:tmax]
        totalinf <- totalinf[tmin:tmax]
        removed <- removed[tmin:tmax]
        sus <- sus[tmin:tmax]
    }
    
    
    if (plottype=="complete"){
        plot (time,sus,xlim=c(tmin, tmax),ylim=c(1,n),main="Epidemic Curve",ylab="Number of individuals" ,xlab="time",type="l",lwd=2,cex=0.5,xaxt="n")
        lines(time, totalinf, col="red",lwd=2)
        lines(time, removed, col="blue",lwd=2)
        axis(1, at=1:tmax)
        legend("topright",inset=.001,cex=0.8,bty = "n",
        c("Infected","Susceptible","Removed"),
        horiz=TRUE,lty=c(1,1,1),lwd=c(2,2,2),col=c("red","black","blue"))
        
    }
    if ((plottype=="complete") & (!is.null(timepoints))){
        plot (time,sus,xlim=c(timepoints[1], timepoints[2]),ylim=c(1,n),main="Epidemic Curve", ylab="Number of individuals",xlab="time",type="l",lwd=2,cex=0.5,xaxt="n")
        lines(time, totalinf, col="red",lwd=2)
        lines(time, removed, col="blue",lwd=2)
        axis(1, at=1:timepoints[2])
        legend("topright",inset=.001,bty = "n",
        c("Infected","Susceptible","Removed"),
        horiz=TRUE,lty=c(1,1,1),lwd=c(2,2,2),col=c("red","black","blue"))
        
    }
}


if (plottype=="totalinfect"){
    plot (time,totalinf,xlim=c(min(time), max(time)),ylim=c(0,(max(totalinf)+1)),main="Epidemic Curve", ylab="Number of infected individuals",xlab="time",type="b",pch=20,lwd=2,xaxt="n")
    axis(1, at=1:max(time))
}

if (plottype=="newinfect"){
    plot (time,newinf,xlim=c(min(time), max(time)),ylim=c(0,(max(newinf)+1)),main="Epidemic Curve",
    ylab="Number of new infections ",xlab="time",type="b",pch=20,lwd=2,xaxt="n")
    axis(1, at=1:max(time))
}

if (plottype=="susceptible"){
    plot (time,sus,xlim=c(min(time),max(time)),ylim=c(0,(max(sus)+1)),main="Epidemic Curve",
    ylab="Number of susceptibles  ",xlab="time",type="b",pch=20,lwd=2,xaxt="n")
    axis(1, at=1:max(time))
}


if (!is.null(timepoints)){
    
    if (plottype=="totalinfect"){
        plot (time,totalinf,xlim=c(timepoints[1], timepoints[2]),ylim=c(1,n),main="Epidemic Curve",
        ylab="Number of infected individuals ",xlab="time",type="b",pch=20,lwd=2,xaxt="n")
        axis(1, at=1:timepoints[2])
    }
    
    if (plottype=="newinfect"){
        plot (time,newinf,xlim=c(timepoints[1], timepoints[2]),ylim=c(0,(max(newinf)+1)),main="Epidemic Curve",
        ylab="Number of new infections ",xlab="time",type="b",pch=20,lwd=2,xaxt="n")
        axis(1, at=1:timepoints[2])
    }

if (plottype=="susceptible"){
    plot (time,sus,xlim=c(timepoints[1], timepoints[2]),ylim=c(1,n),main="Epidemic curve",
    ylab="Number of susceptibles  ",xlab="time",type="b",pch=20,lwd=2,xaxt="n")
    axis(1, at=1:timepoints[2])
    
}
}

}




