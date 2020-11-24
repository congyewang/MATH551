#hess evaluates the variance-covariance matrix for the mle, based
#on observed info. y is the mle, x is the data, and f is minus the
# log likelihood

hess<-function(f,y,...){

	ep<-0.0001
	eps<- ep*y
	n <- length(y)
	m <- matrix(0,ncol=n, nrow=n)
	for(i in 1:n){
		for( j in 1:n){
		y1<-y
		y1[i]<-y1[i]+eps[i]
		y1[j]<-y1[j]+eps[j]
		y2<-y
		y2[i]<-y2[i]+eps[i]
		y2[j]<-y2[j]-eps[j]
		y3<-y
		y3[i]<-y3[i]-eps[i]
		y3[j]<-y3[j]+eps[j]
		y4<-y
		y4[i]<-y4[i]-eps[i]
		y4[j]<-y4[j]-eps[j]
	m[i,j]<-(f(y1,...)-f(y2,...)-f(y3,...)+f(y4,...))/(4*eps[i]*eps[j])
	}
	}

	return(m)
      }

##function for calculating the derivative at mle - g is function, y is mle
##... is any other arguments needed for g
deriv <- function(g,y,...)
  {
    n<-length(y)
    ep<-0.0001
    eps<- ep*y

    d <- rep(0,n)
    
    for(i in 1:n){
      yt<-y
      yt[i]<-yt[i]+eps[i]
      d[i] <- (g(yt,...)-g(y,...))/eps[i]
    }

    return(d)
    
  }
