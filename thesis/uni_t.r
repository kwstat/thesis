# uni_t.r

# nu=.05
# n=4, observed data = -20,1,2,3

# Here is the whole loglikelihood:
  
# LogLike <- function(mu,w){
#   L <- 4*log(gamma(.5025)) -2*log(pi/20) -4*log(gamma(.025))
#   for(i in 1:length(w)){
#     L <- L - .5025 * log(1 +20*(w[i]-mu)^2)
#   }
#   return(L)
# }

# The loglikelihood is proportional to:

LogLike <- function(mu,w){
  L <- 0
  for(i in 1:length(w)){
    L <- L - log(1 +20*(w[i]-mu)^2)
  }
  return(L)
}



setwd("x:/kw/Research/Thesis/")
pdf("uni_t.pdf")
#postscript("uni_t.ps",width=5,height=5,horizontal=FALSE)
w <- c(-20,1,2,3)  # Observed data
mu <- -250:150/10
Lmu <- LogLike(mu,w)
plot(mu,Lmu,type='l',xlab="",ylab="")
dev.off()
