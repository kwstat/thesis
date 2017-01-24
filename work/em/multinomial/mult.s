
#-------------------------------------------------------------------------
# Create a graph with different q functions for different values of
# p_k

postscript(file="/home/kwright/I/EM/Multinomial/mult.qfun.ps")
ps.options(horizontal=F)
p_1:99/100

pk_.1
qp_(125*pk/(2+pk) +34)*log(p) + 38*log(1-p)
plot(p,qp,type='l',xlab="",ylab="")
pk_.2
qp_(125*pk/(2+pk) +34)*log(p) + 38*log(1-p)
lines(p,qp,lty=2)
pk_.3
qp_(125*pk/(2+pk) +34)*log(p) + 38*log(1-p)
lines(p,qp,lty=2)
pk_.4
qp_(125*pk/(2+pk) +34)*log(p) + 38*log(1-p)
lines(p,qp,lty=2)
pk_.5
qp_(125*pk/(2+pk) +34)*log(p) + 38*log(1-p)
lines(p,qp,lty=2)
pk_.6
qp_(125*pk/(2+pk) +34)*log(p) + 38*log(1-p)
lines(p,qp,lty=2)
pk_.7
qp_(125*pk/(2+pk) +34)*log(p) + 38*log(1-p)
lines(p,qp,lty=2)
pk_.8
qp_(125*pk/(2+pk) +34)*log(p) + 38*log(1-p)
lines(p,qp,lty=2)
pk_.9
qp_(125*pk/(2+pk) +34)*log(p) + 38*log(1-p)
lines(p,qp,lty=1)

#legend(.4,-160,legend=c("Lower"),lty=1)
title("",xlab="p",ylab="q(p|.)")

dev.off()




