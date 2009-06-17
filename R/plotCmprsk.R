plotCmprsk <- function(x, plottype=1, dc=TRUE, p1, p2, filename, main=" ", 
ylim=c(-1,1.5), xlim, legloc=1, xlab="Mark", ylab="Estimated Risk", ...){


V1 <- x[[1]]$V1
V2 <- x[[1]]$V2
VEoverallPH <- x[[1]]$VEhatPH
VEoverallCI <- x[[1]]$VEhatCI
matrix2 <- x[[2]] 
matrix3 <- x[[3]]
if(missing(p1)) {pval1 <- x[[1]]$pval13} else {pval1 <- p1}
if(missing(p2)) {pval2 <- x[[1]]$pvalnp2} else {pval2 <- p2}
if(missing(xlim)) {xlim <- c(V1,V2)} else {}
if(missing(filename)) {stop("Please enter a value for filename.")}
if(!exists("type")){type <- "n"}
if(length(legloc)>1) {legloc <- c(xlim[1], ylim[2])} else {legloc <- "topright"}


if(dc==FALSE){
mark <- matrix2[,3]
VE <- matrix2[,8]
varVE <- matrix2[,9]
lowlim <- matrix2[,10]
uplim <- matrix2[,11]
} else {
mark <- matrix3[,3]
VE <- matrix3[,8]
varVE <- matrix3[,9]
lowlim <- matrix3[,10]
uplim <- matrix3[,11]
}

inds <- VE > -10 & varVE < 10
mark <- mark[inds]
VE <- VE[inds]
lowlim <- lowlim[inds]
uplim  <- uplim[inds]
VE <- ifelse(VE<1,VE,1)
uplim <- ifelse(uplim<1,uplim,1)
keep <- mark >= V1 & mark <= V2

## Plots of VE
if(plottype==1) {
postscript(paste(filename,".ps",sep=""),horizontal=T)
par(las=1)
plot(mark[keep],VE[keep],main=main,ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab,type,...)
lines(mark[keep],VE[keep],lwd=3,col="red")
lines(mark[keep],lowlim[keep],lwd=3,lty=2,col="red")
lines(mark[keep],uplim[keep],lwd=3,lty=2,col="red")
abline(h=0,lty=2,lwd=2,col="green")
abline(h=VEoverallPH,lty=1,lwd=2,col="blue")
if(dc==FALSE){
legend(legloc,legend=c("Estimated VE^c(t,v)","95% CIs","Overall VE^c(t)",paste("P-value for H0^0: VE(t,v) = 0: ",round(pval1,4)),paste("P-value for H0: VE(t,v) = VE(t): ",round(pval2,4))),lty=c(1,3,1,0,0), lwd=c(2,3,2,2,20),
col=c("red","red","blue","green","green"))
}else{
legend(legloc,legend=c("Estimated VE^dc(t,v)","95% CIs","Overall VE^dc(t)"),lty=c(1,3,1),lwd=c(2,3,2),col=c("red","red","blue"))}
dev.off()

} else {

## Plots of RR
postscript(paste(filename,".ps",sep=""),horizontal=T)
par(las=1)
plot(mark[keep],1-VE[keep],main=main,ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab,type,...)
lines(mark[keep],1-VE[keep],lwd=3,col="red")
lines(mark[keep],1-lowlim[keep],lwd=3,lty=2,col="red")
lines(mark[keep],1-uplim[keep],lwd=3,lty=2,col="red")
abline(h=1,lty=2,lwd=2,col="green")
abline(h=1-VEoverallPH,lty=1,lwd=2,col="blue")
if(dc==FALSE){
legend(legloc,legend=c("Estimated RR^c(t,v)","95% CIs","Overall RR^c(t)", paste("P-value for H0^0: RR(t,v) = 1: ",round(pval1,4)),paste("P-value for H0: RR(t,v) = RR(t): ", round(pval2,4))), lty=c(1,3,1,0,0),lwd=c(2,3,2), col=c("red","red","blue"))} else { 
legend(legloc,legend=c("Estimated RR^dc(t,v)","95% CIs","Overall RR^dc(t)"),lty=c(1,3,1),col=c("red","red","blue"),lwd=c(2,3,2))}
dev.off()
}
}

