
# mult_qfun.r

# Create a graph with different q functions for different values of
# p_k

setwd("x:/kw/Research/Thesis/")
pdf("mult_qfun.pdf")
#postscript("mult_qfun.ps",width=5,height=5,horizontal=FALSE)

p <- 1:99/100

pk <- .1
qp <- (125*pk/(2+pk) +34)*log(p) + 38*log(1-p)
plot(p,qp,type='l',xlab="",ylab="")
pk <- .2
qp <- (125*pk/(2+pk) +34)*log(p) + 38*log(1-p)
lines(p,qp,lty=2)
pk <- .3
qp <- (125*pk/(2+pk) +34)*log(p) + 38*log(1-p)
lines(p,qp,lty=2)
pk <- .4
qp <- (125*pk/(2+pk) +34)*log(p) + 38*log(1-p)
lines(p,qp,lty=2)
pk <- .5
qp <- (125*pk/(2+pk) +34)*log(p) + 38*log(1-p)
lines(p,qp,lty=2)
pk <- .6
qp <- (125*pk/(2+pk) +34)*log(p) + 38*log(1-p)
lines(p,qp,lty=2)
pk <- .7
qp <- (125*pk/(2+pk) +34)*log(p) + 38*log(1-p)
lines(p,qp,lty=2)
pk <- .8
qp <- (125*pk/(2+pk) +34)*log(p) + 38*log(1-p)
lines(p,qp,lty=2)
pk <- .9
qp <- (125*pk/(2+pk) +34)*log(p) + 38*log(1-p)
lines(p,qp,lty=1)

title("",xlab="p",ylab="q(p|.)")

dev.off()




