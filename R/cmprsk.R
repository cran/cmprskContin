cmprsk <- function(gp, ftime, ftype, mark, nboot=5000, ngrid=25, ngridv=25, 
T1=0, T2=0, UT1=0, UT2=0, ttanal=0, BAND1=0, BAND2=0, TAILSL=1, TAILSU=1, 
V1=0, V2=0, UV1=0, UV2=0, BANDV1=0, BANDV2=0, BANDVLOW=0, BANDVUP=0, TAILSV=1) {

## Checks
if(!exists("gp") | !exists("ftime") | !exists("ftype") | !exists("mark")) {stop("All of gp, ftime, ftype and mark are required inputs.")}
if(any(is.na(c(gp,ftime,ftype,mark)))==TRUE) {stop("Missing mark values should be coded as 99. Missing data for gp, ftime and ftype are not allowed. Please see the documention for details.")}
if(T1>0 & UT1>T1) {stop("The value of UT1 must be less than or equal to T1.")}
if(UT2>0 & UT2<T2) {stop("The value of UT2 must be greater than or equal to T2.")} 
if(V1>0 &  UV1>V1) {stop("The value of UV1 must be less than or equal to V1.")}
if(UV2>0 &  UV2<V2) {stop("The value of UV2 must be greater than or equal to V2.")} 

## Inputs for Fortran:
maxn <- length(gp)
maxm <- length(which(ftype==1)) 
keep <- ftype==0 | (ftype==1 & mark >=0 & mark <= 1)
gp <- gp[keep]
ftime <- ftime[keep]
ftype <- ftype[keep]
mark <- mark[keep]
df <- data.frame(gp, ftime, ftype, mark) 
names(df) <- c("Vx","time","delta","mark")
vac <- subset(df, df$Vx==1) 
placebo <- subset(df, df$Vx!=1) 
time <- c(vac$time, placebo$time)
censor <- c(vac$delta, placebo$delta)
cause <- c(vac$mark, placebo$mark)
gpsub <- c(length(vac$Vx), length(placebo$Vx))
tsub <- length(vac$Vx) + length(placebo$Vx)
dim1 <- ngridv*tsub
dim2 <- ngridv
dim3 <- tsub
doubIN <- c(T1,T2,UT1,UT2,ttanal,BAND1,BAND2,BANDV1,BANDV2,BANDVLOW,BANDVUP,V1,V2,UV1,UV2) 
intIN <- c(nboot,ngrid,ngridv,TAILSL,TAILSU,TAILSV)
mark <- round(mark, digits=4)
time <- round(time, digits=4)

## Fortran call
ans = .Fortran('cmprskContin', 
	as.integer(dim1), as.integer(dim2), as.integer(dim3), as.integer(maxn), as.integer(maxm),
	as.integer(gpsub), as.double(time), as.integer(censor), as.double(cause), as.integer(intIN), as.double(doubIN), 
	OUTV3=integer(3), OUTV30=double(30), OUTM1=integer(ngridv), OUTM2=double(ngridv), 
	OUTM_VEc1=double(ngridv), OUTM_VEc2=double(ngridv), OUTM_VEc3=double(ngridv), OUTM_VEc4=double(ngridv), 
	OUTM_VEc5=double(ngridv), OUTM_VEc6=double(ngridv), OUTM_VEc7=double(ngridv), OUTM_VEc8=double(ngridv), 
	OUTM_VEc9=double(ngridv), OUTM_VEdc1=double(ngridv), OUTM_VEdc2=double(ngridv), OUTM_VEdc3=double(ngridv), 
	OUTM_VEdc4=double(ngridv), OUTM_VEdc5=double(ngridv), OUTM_VEdc6=double(ngridv), OUTM_VEdc7=double(ngridv), 
	OUTM_VEdc8=double(ngridv), OUTM_VEdc9=double(ngridv), OUTtm1=double(dim1), OUTtm2=double(dim1), 
	OUTtm3=double(dim1), OUTtm4=double(dim1), OUTtm5=double(dim1), OUTtm6=double(dim1), OUTtm7=double(dim1), 
	OUTtm8=double(dim1), OUTtm9=double(dim1), PACKAGE="cmprskContin")

## Test statistics and p-values
STATvec <- as.list(c(ans$OUTV3, ans$OUTV30))
names(STATvec) <- c("nsamp1","nsamp2","nboot","T1","T2","ttanal","AvgV1","AvgV2","BAND1","BAND2","BANDV1","BANDV2",
"V1","V2","VEhatCI","VEhatPH","LogRankZ","U11","U12","U13","U14","pval11","pval12","pval13","pval14","Unp1","Unp2",
"pvalnp1","pvalnp2","Usp1","Usp2","pvalsp1","pvalsp2")

## Data for estimating VE^c(ttanal,v) vs v:
VECmat <- cbind(ans$OUTM1,ans$OUTM2,ans$OUTM_VEc1,ans$OUTM_VEc2,ans$OUTM_VEc3,ans$OUTM_VEc4,ans$OUTM_VEc5,ans$OUTM_VEc6,ans$OUTM_VEc7,ans$OUTM_VEc8,ans$OUTM_VEc9)
colnames(VECmat) <- c("index","ttanal","mark","F1","F2","varF1","varF2","VEC","varVEC","CIlow","CIhigh")

## Data for estimating VE^dc(ttanal,v) vs v:
VEDCmat <- cbind(ans$OUTM1,ans$OUTM2,ans$OUTM_VEdc1,ans$OUTM_VEdc2,ans$OUTM_VEdc3,ans$OUTM_VEdc4,ans$OUTM_VEdc5,ans$OUTM_VEdc6,ans$OUTM_VEdc7,ans$OUTM_VEdc8,ans$OUTM_VEdc9)
colnames(VEDCmat) <- c("index","ttanal","mark","F1dc","F2dc","varF1dc","varF2dc","VEDC","varVEDC","CIlow","CIhigh")

## Data for VE^c at all timepoints between T1 and T2, not just at ttanal:  
timeMark <- cbind(ans$OUTtm1,ans$OUTtm2,ans$OUTtm3,ans$OUTtm4,ans$OUTtm5,ans$OUTtm6,ans$OUTtm7,ans$OUTtm8,ans$OUTtm9)
# remove unused portion of matrix
ind <- which(timeMark[,1]==0)
timeMark <- timeMark[-ind,]
colnames(timeMark) <- c("nsamp","eventtime","mark","F1","F2","varF1","varF2","VEC","varVEC")

outlist <- list(STATvec,VECmat,VEDCmat,timeMark)

outlist
}
