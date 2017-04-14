
epispatial<- function(type,x,y,inftime,removaltime=NULL,time=NULL,tmin=NULL){
    
tau <- inftime
if (is.null(type) || !(type %in% c("SI", "SIR"))) {
        stop("epiplot: Specify type as \"SI\" or \"SIR\" ", call. = FALSE)
    }
n <- length(tau)
 
if ((length(y) !=n) || (length(x) !=n))
     stop('epispatial: Length of x or y is not compatible ')
   
 
if (is.null(removaltime) && type=="SIR") {
     stop(' epispatial: Specify removal time ')
 }
if (!is.null(removaltime)){
     if (length(removaltime) != n)
     stop('epispatial: Length of removaltime is not compatible.')
     if(type=="SI")
     stop('epispatial: Type must be "SIR".')
 }

if (is.null(tmin)){
    tmin=1
}

old.par<-par(mfrow=c(3,3))
par(cex = 0.5)
par(xpd=T)

if (type=="SI"){
    
    dat <- data.frame(x,y,tau)
    
    if (is.null(time)){
        for(i in tmin:max(tau)){
            xcc<-subset(dat, tau<=i & tau!=0)
            plot(x,y,xlim=c(min(x),max(x)), ylim=c(min(y),max(y)),panel.first=grid(),sub= paste("time ",i))
            points(xcc$x,xcc$y,pch=16,col="red")
            legend(x=min(x),y=max(y)+2,c("Infected","Susceptible"),col=c("red","black"),pch=c(16,21),bty = "n",horiz=TRUE)
        }
    }
    if (!is.null(time)){
            for(i in 1:length(time)){
                xcc<-subset(dat, tau<=time[i] & tau!=0)
                plot(x,y,xlim=c(min(x),max(x)), ylim=c(min(y),max(y)),panel.first=grid(),sub= paste("time ",time[i]))
                points(xcc$x,xcc$y,pch=16,col="red")
                legend(x=min(x),y=max(y)+2,c("Infected","Susceptible"),col=c("red","black"),pch=c(16,21),bty = "n",horiz=TRUE)
            }
        
    }

    
}


if (type=="SIR"){
    dat <- data.frame(x,y,tau,removaltime)
    if (is.null(time)){
        for(i in tmin:max(tau)){
            xcc<-subset(dat, tau<=i & tau!=0)
            plot(x,y,xlim=c(min(x),max(x)), ylim=c(min(y),max(y)),pch=21,panel.first=grid(),sub= paste("time ",i))
            
            xred <- subset(xcc, i < xcc$removaltime)
            points(xred$x,xred$y,pch=16,col="red")
            xblue <- subset(xcc, i >= xcc$removaltime)
            points(xblue$x,xblue$y,pch=4,col="blue")
            
            legend(x=min(x),y=max(y)+2,c("Infected","Susceptible","Removed"),col=c("red","black","blue"),pch=c(16,21,4),bty = "n",horiz=TRUE)
        }
    }
    if (!is.null(time)){
            for(i in 1:length(time)){
                xcc<-subset(dat, tau<=time[i] & tau!=0)
                plot(x,y,xlim=c(min(x),max(x)), ylim=c(min(y),max(y)),pch=21,panel.first=grid(),sub= paste("time ",time[i]))
                
                xred <- subset(xcc, time[i] < xcc$removaltime)
                points(xred$x,xred$y,pch=16,col="red")
                xblue <- subset(xcc, time[i] >= xcc$removaltime)
                points(xblue$x,xblue$y,pch=4,col="blue")
                
                legend(x=min(x),y=max(y)+2,c("Infected","Susceptible","Removed"),col=c("red","black","blue"),pch=c(16,21,4),bty = "n",horiz=TRUE)
            }
    }
        
}

par(old.par)
}

