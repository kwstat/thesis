Q_function(nu, mu, w, muk)
{
	temp <- 0
	for(i in 1:length(w)) {
		temp <- temp + (w[i] - mu)^2/(nu + (w[i] - muk)^2)
	}
	temp <- temp * (nu + 1)
	return(temp)
}


#-------------------------------------------------------------------------

nu=.05
n=4, observed data = -20,1,2,3

Here is the whole loglikelihood:
  
LogLike_function(mu,w){
  L _ 4*log(gamma(.5025)) -2*log(pi/20) -4*log(gamma(.025))
  for(i in 1:length(w)){
    L _ L - .5025 * log(1 +20*(w[i]-mu)^2)
  }
  return(L)
}

The loglikelihood is proportional to:

LogLike_function(mu,w){
  L _ 0
  for(i in 1:length(w)){
    L _ L - log(1 +20*(w[i]-mu)^2)
  }
  return(L)
}


postscript(file="uni_t.ps")
ps.options(horizontal=F)
mu_-250:150/10
Lmu_LogLike(mu,w)
plot(mu,Lmu,type='l',xlab="",ylab="")
dev.off()
