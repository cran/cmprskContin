cc	cmprskContin (6/1/09)
cc
cc   Peter Gilbert, Ian McKeague, Yanqing Sun 2007
cc
cc   Computes the test statistics and p-values,
cc   and point and confidence interval estimates of VE^c(t,v)
cc   and VE^dc(t,v) for a fixed v, as described in Gilbert, McKeague,
cc   and Sun (2008, Biostatistics)
cc
cc
cc------------------
cc	INPUTS:
cc------------------
cc	dim1IN = mark gridpoints multiplied by total subjects
cc	dim2IN = number of mark gridpoints 
cc	dim3IN = total subjects
cc	maxn = total subjects
cc	maxm = number of failures
cc	sampleIN = vector of length 2, # subs in group 1 and # subs in group 2
cc	timeIN = vector of right-censored failure times
cc	censorIN = vector indicating infected (1) or not infected(0)
cc	causeIN = vector of continuous mark value
cc	intIN = vector of length 6 containing the following:
cc		nboot = number of bootstrap replicates used for computing p-values
cc		ngrid = number of failure time gridpoints 
cc		ngridv = number of mark gridpoints
cc		TAILSL = left tail correction for kernel estimation of hazard function (1=correction,2=no correction) 
cc		TAILSU = right tail correction for kernel estimation of hazard function (1=correction,2=no correction) 
cc		TAILSV = upper and lower tail correction for smoothing in the mark (1=smoothing, 2=no smoothing)
cc	doubIN = vector of length 15 containing the following:
cc		T1 = smallest failure time used to compute test statistics
cc		T2 = largest failure time used to compute test statistics
cc		UT1 = lower envelope for smoothing the failure time (UT1<T1)
cc		UT2 = upper envelope for smoothing the failure time (T2<UT2) 
cc		ttanal = failure time at which VE^c and VE^dc are evaluated
cc		BAND1 = bandwidth for kernel estimation of group 1 hazard 
cc		BAND2 = bandwidth for kernel estimation of group 2 hazard 
cc		BANDV1 = bandwidth for smoothing the group 1 mark
cc		BANDV2 = bandwidth for smoothing the group 2 mark
cc		BANDVLOW = lower range of the mark values to search over when optimizing BANDV1 and BANDV2
cc		BANDVUP = upper range of the mark values to search over when optimizing BANDV1 and BANDV2
cc		V1 = left limit for smoothing in the mark
cc		V2 = right limit for smoothing in the mark
cc		UV1 = left envelope for smoothing in the mark
cc		UV2 = right envelope for smoothing in the mark
cc	
cc	
cc------------------
cc	OUTPUTS:
cc------------------
cc	OUTV3 = vector of length 3: 
cc		nsamp1 = # sub in group1 
cc		nsamp2 = # sub in group 2 
cc		nboot = number of bootstrap replicates used for computing p-values
cc
cc	OUTV30 = vector of length 30:
cc		ttanal = failure time at which VE^c and VE^dc were evaluated 
cc		AvgV1 = average mark among infected group 1 subjects
cc		AvgV2 = average mark among infected group 2 subjects
cc		BAND1 = bandwidth for kernel estimation of group 1 hazard
cc		BAND2 = bandwidth for kernel estimation of group 2 hazard
cc		BANDV1 = bandwidth for smoothing the group 1 mark
cc		BANDV2 = bandwidth for smoothing the group 2 mark
cc		V1 = smallest mark used for estimation 
cc		V2 = largest mark used for estimation 
cc		VEhatCI = cumulative incidence estimate of overall VE 
cc		VEhatPH = proportional hazard estimate of overall VE 
cc		LogRankZ = log-rank statistic for comparing overall hazard rates 
cc		U11 = nonparametric test statistic for 1-sided testing of H_0^0: VE(t,v)=0 
cc		U12 = nonparametric test statistic for 1-sided testing of H_0^0: VE(t,v)=0 
cc		U12 = nonparametric test statistic for 2-sided testing of H_0^0: VE(t,v)=0 
cc		U12 = nonparametric test statistic for 2-sided testing of H_0^0: VE(t,v)=0 
cc		pval11 = p-value corresponding to U11 
cc		pval12 = p-value corresponding to U12 
cc		pval13 = p-value corresponding to U13 
cc		pval14 = p-value corresponding to U14 
cc		Unp1 = nonparametric test statistic for 1-sided testing of H_0: VE(t,v)=VE(t) 
cc		Unp2 = nonparametric test statistic for 2-sided testing of H_0: VE(t,v)=VE(t) 
cc		pvalnp1 = p-value corresponding to Unp1 
cc		pvalnp2 = p-value corresponding to Unp2 
cc		Usp1 = semiparametric test statistic for 1-sided testing of H_0: VE(t,v)=VE(t)  
cc		Usp2 = semiparametric test statistic for 2-sided testing of H_0: VE(t,v)=VE(t)  
cc		pvalsp1 = p-value corresponding to Usp1 
cc		pvalsp2 = p-value corresponding to Usp2 
cc
cc	All vectors below are of length ngridv
cc	OUTM1 = vector of the index of the analysis time (ttanal) 
cc	OUTM2 = vector of ttanal 
cc	OUTM_VEc1 = vector of mark values 
cc	OUTM_VEc2 = vector of estimated cumulative incidence function [F1(ttanal,v)] for group 1 
cc	OUTM_VEc3 = vector of estimated cumulative incidence function for [(F2(ttanal,v)] group 2   
cc	OUTM_VEc4 = vector of estimated variance of F1(ttanal,v) 
cc	OUTM_VEc5 = vector of estimated variance of F2(ttanal,v) 
cc	OUTM_VEc6 = vector of estimated cumulative vaccine efficacy (VEC) 
cc	OUTM_VEc7 = vector of estimated variance of VEC 
cc	OUTM_VEc8 = vector of lower 95% CI for VEC 
cc	OUTM_VEc9 = vector of upper 95% CI for VEC 
cc
cc	OUTM_VEdc1 = vector of mark values 
cc	OUTM_VEdc2 = vector of estimated doubly cumulative incidence function (F1dc) for group 1 
cc	OUTM_VEdc3 = vector of estimated doubly cumulative incidence function for (F2dc) group 2
cc	OUTM_VEdc4 =  vector of estimated variance of F1dc 
cc	OUTM_VEdc5 =  vector of estimated variance of F2dc 
cc	OUTM_VEdc6 = vector of estimated doubly cumulative vaccine efficacy (VEDC) 
cc	OUTM_VEdc7 = vector of estimated variance of VEDC 
cc	OUTM_VEdc8 = vector of lower 95% CI for VEDC 
cc	OUTM_VEdc9 = vector of upper 95% CI for VEDC 
cc
cc	OUTtm1 = vector of sample index 
cc	OUTtm2 = vector of observed event times 
cc	OUTtm3 = vector of mark values 
cc	OUTtm4 = vector of estimated cumulative incidence function [(F1(t,v) for group 1 
cc	OUTtm5 = vector of estimated cumulative incidence function [(F2(t,v) for group 2 
cc	OUTtm6 = vector of estimated variance of F1(t,v) 
cc	OUTtm7 = vector of estimated variance of F2(t,v)
cc	OUTtm8 = vector of estimated cumulative vaccine efficacy 
cc	OUTtm9 = vector of estimated variance of cumulative vaccine efficacy 
cc	
cc	
cc-----------------------------------------
cc	FUNCTIONS AND SUBROUTINES:
cc-----------------------------------------
cc	ran1 : generates a single random number
cc	rstart : initializes the table for the f(97,33,-mod 1.)
cc generator and the values for the arithmetic sequence.
cc Note that a BLOCK DATA initializes the table, so it is not
cc imperative to call this subroutine. The seeds used for
cc generating the default table is 12,34,56,78.
cc	rnor : returns a standard normal variate using the Ziggurat method
cc	uni : combines, combines, with subtraction mod 1,
cc an f(97,33,-mod 1) generator with the element c in the arithmetic
cc sequence generated by c=c-cd mod(16777213./16777216.), period 2**24-3.
cc period of combined generator is (2**97-1)(2**24-3)2**23, about 2**144.
cc	iuni : combines, with subtraction mod 2**24,
cc an f(97,33,-mod 2**24) generator with the element c in the arithmetic
cc sequence generated by c=c-cd mod(16777213), period 2**24-3.
cc period of combined generator is (2**97-1)(2**24-3)2**23, about 2**144.
cc	ivni : combines, with subtraction mod 2**24,
cc an f(97,33,-mod 2**24) generator with the element c in the arithmetic
cc sequence generated by c=c-cd mod(16777213), period 2**24-3.
cc period of combined generator is (2**97-1)(2**24-3)2**23, about 2**144.
cc	indexx : indexes the array arrin, outputs the array indx 
cc	EPAN : Epanechnikov kernel function
cc	EPANL : Epanechnikov lower tail kernel function
cc	EPANU : Epanechnikov upper tail kernel function
cc	alphdblp : Computes second derivative bandwidth for kernel smoothing. Details in Andersen, P.K., Borgan, O., Gill, R.D., and Keiding (1993)
cc	estp : Computes maximum likelihood estimates for Cox model
cc	GAUSSJ : Gaussian elimination with pivoting
c
cc----------------------------------------------------------------------

	SUBROUTINE cmprskContin(
     +dim1IN,dim2IN,dim3IN,maxn,maxm,
     +sampleIN,timeIN,censorIN,causeIN,intIN,doubIN,
     +OUTV3,OUTV30,OUTM1,OUTM2,
     +OUTM_VEc1,OUTM_VEc2,OUTM_VEc3,OUTM_VEc4,OUTM_VEc5,
     +OUTM_VEc6,OUTM_VEc7,OUTM_VEc8,OUTM_VEc9,OUTM_VEdc1,
     +OUTM_VEdc2,OUTM_VEdc3,OUTM_VEdc4,OUTM_VEdc5,OUTM_VEdc6,
     +OUTM_VEdc7,OUTM_VEdc8,OUTM_VEdc9,OUTtm1,OUTtm2,
     +OUTtm3,OUTtm4,OUTtm5,OUTtm6,OUTtm7,OUTtm8,OUTtm9)

	integer mycount,mystart,ran1flag 
	integer dim1IN,dim2IN,dim3IN,intIN(6) 
	integer maxn,maxm 
	integer sampleIN(2),censorIN(dim3IN)
	double precision causeIN(dim3IN),timeIN(dim3IN)
	double precision doubIN(15)
	integer OUTV3(3)
	double precision OUTV30(30)
	integer OUTM1(dim2IN)
	double precision OUTM2(dim2IN)
	double precision OUTM_VEc1(dim2IN),OUTM_VEc2(dim2IN)
	double precision OUTM_VEc3(dim2IN),OUTM_VEc4(dim2IN)
	double precision OUTM_VEc5(dim2IN),OUTM_VEc6(dim2IN)
	double precision OUTM_VEc7(dim2IN),OUTM_VEc8(dim2IN)
	double precision OUTM_VEc9(dim2IN)
	double precision OUTM_VEdc1(dim2IN),OUTM_VEdc2(dim2IN)
	double precision OUTM_VEdc3(dim2IN),OUTM_VEdc4(dim2IN)
	double precision OUTM_VEdc5(dim2IN),OUTM_VEdc6(dim2IN)
	double precision OUTM_VEdc7(dim2IN),OUTM_VEdc8(dim2IN)
	double precision OUTM_VEdc9(dim2IN)
	double precision OUTtm1(dim1IN),OUTtm2(dim1IN),OUTtm3(dim1IN)
	double precision OUTtm4(dim1IN),OUTtm5(dim1IN),OUTtm6(dim1IN)
	double precision OUTtm7(dim1IN),OUTtm8(dim1IN),OUTtm9(dim1IN) 

      integer seeda,seedb,seedc,seedd,nsamp1,nsamp2,nsamp,indx(maxn)
      integer  nv1,nv2,nv,ntau,nboot,nrun,javg,oldind,irun,ngrid
      integer  indstm1(maxn),indstm2(maxn),indsv1(maxn),indsv2(maxn)
      integer  ind,den,T1ind1,T1ind2,nsamp1tr,nsamp2tr
      integer indtt1,indtt2,indtt,ngridv,indsCV1(maxn),ihalf
      real stmppr1(0:maxn,maxn),stmppr2(0:maxn,maxn)
      real KM01,KMC1,KMT1(0:maxn),KM02,KMC2,KMT2(0:maxn)
      real T1,T2,UT1,UT2,BAND1,BAND2,TAILSL,TAILSU,theta1,theta2,beta1
      real beta2,censor(maxn),cause(maxn),ocause(maxn)
	real time1(maxn),time2(maxn),time(maxn)
	real censor1(maxn),cause1(maxn),ocause1(maxn)
      real censor2(maxn),cause2(maxn),ocause2(maxn)
      real copy(maxn),hatl1(maxn,maxm),hatl2(maxn,maxm)
      real centim,wt,wt1,wt2,temp1,temp2,ttu,sprocess(0:maxn,maxm)
      real process(0:maxn,maxm),tmpproc1(0:maxn,maxm),num
      real tmproc1m(0:maxn,maxm),tmproc2m(0:maxn,maxm)
      real stmpproc1(0:maxn,maxm),stmpproc2(0:maxn,maxm)
      real processm(0:maxn,maxm),Hnm1(maxn),Hnm2(maxn)
      real tmpproc2(0:maxn,maxm),ee1(maxn),ee2(maxn),group(maxn)
      real NA1(0:maxn),NA2(0:maxn),JNA,vavg1,vavg2,delta11,delta12
      real gtime1(maxn),gcensor1(maxn),gcause1(maxn)
      real gtime2(maxn),gcensor2(maxn),gcause2(maxn)
      real rannums(maxn),cumW1(maxn),cumW2(maxn),dat(5)
      real h1hat1(0:maxn,maxm),h1hat2(0:maxn,maxm),h1hat3(0:maxn,maxm)
      real h1hat4(0:maxn,maxm),h2hat1(0:maxn,maxm),h2hat2(0:maxn,maxm)
      real h2hat20(0:maxn,maxn),h2hat10(0:maxn,maxn)
      real h2hat3(0:maxn,maxm),h2hat4(0:maxn,maxm)
      real h1hat5(0:maxn,maxm),h2hat5(0:maxn,maxm)
      real nrej1,nrej2,u1,u2,u3,u1obs,u2obs,u3obs,mucen,Y12(maxn)
      real u4,u4obs,pvalue1,pvalue2,pvalue3,pvalue4,u5,u5obs
      real Y21(maxn),nrej3,nrej4,nrej5,vv,vvnext,pvalue5
      real nrej52,nrej53,nrej54,nrej52ph,nrej53ph,nrej54ph
      real W1(maxn),W2(maxn),Lambda2pr2(0:maxn,maxm),obscause(maxn)
      real ahat(maxn),bhat(maxn),EE21(maxn),Lambda2pr1(0:maxn,maxm)
      real ttemp1,ttemp2,ttemp3,ttemp4,temp3,temp4,temp5,ttemp5
      real F1T2,F2T2,VE,VEoverall,u1mean,u2mean,u3mean,ttanal
      real u4mean,weightv(maxm)
      real F1vectTT(maxm),F2vectTT(maxm),VEvectTT(maxm),vVEvcTT(maxm)
      real F1vectTTcum(maxm),F2vectTTcum(maxm),VEvectTTcum(maxm)
      real vVEvcTTcum(maxm)
      real F1vect1(0:maxn,maxm),F2vect2(0:maxn,maxm)
      real F1vectCV1(0:maxn,maxm),F2vectCV1(0:maxn,maxm)
      real F1vectCV2(0:maxn,maxm),F2vectCV2(0:maxn,maxm)
      real F1vect1cum(0:maxn,maxm),F2vect2cum(0:maxn,maxm)
      real F1vect(0:maxn,maxm),F2vect(0:maxn,maxm),vVEvect(0:maxn,maxm)
      real F1vectcum(0:maxn,maxm),F2vectcum(0:maxn,maxm)
      real vVEvectcum(0:maxn,maxm)
      real vF1vect1(0:maxn,maxm),vF2vect2(0:maxn,maxm)
      real vF1vect1cum(0:maxn,maxm),vF2vect2cum(0:maxn,maxm)
      real vF1vect(0:maxn,maxm),vF2vect(0:maxn,maxm)
      real vF1vectcum(0:maxn,maxm),vF2vectcum(0:maxn,maxm)
      real vF1vcTT(maxm),vF2vcTT(maxm)
      real vF1vcTTcum(maxm),vF2vcTTcum(maxm)
      real VEvect(0:maxn,maxm),VEvectcum(0:maxn,maxm)
      real gprocess(0:maxn,maxm)
      real proc1(0:maxn,maxm),proc2(0:maxn,maxm)
      real gproc1(0:maxn,maxm),gproc2(0:maxn,maxm)
      real gtmpproc1(0:maxn,maxm),gtmpproc2(0:maxn,maxm)
      real gtmppr1(0:maxn,maxn),gtmppr2(0:maxn,maxn)
      real obsproc(0:maxn,maxm),btproc1(0:maxn,maxm)
      real btproc2(0:maxn,maxm),btproc3(0:maxn,maxm)
      real btproc4(0:maxn,maxm),btproc5(0:maxn,maxm)
      real btproc6(0:maxn,maxm),btproc7(0:maxn,maxm)
      real btproc8(0:maxn,maxm),BANDV1,BANDV2
      real BANDVLOW,BANDVUP,MISEOLD
      real TBANDV1,TBANDV2
      real band111,band112,k2K,intK2,delt11,delt12
      real delt21,delt22,alpy1,alpy2,intalpdp1,intalpdp2
      real V1,V2,lowlim(maxn),uplim(maxn)
      real lowlimcum(maxn),uplimcum(maxn)
      real UV1,UV2,TAILSV,MISE
      real h2hat1ph(0:maxn,maxn),h2hat2ph(0:maxn,maxn)

      real nrej6ph,nrej7ph
      integer NDEAD,NMIS,NP
      real Jnbeta,Unbeta,CL,GS,BAND22
      real s0(maxn),s1(1,maxn),s2(1,1),var
      real BZ(maxn,maxn),U12,U11,beta,covar(1,maxn,maxn)
      real U1beta(maxn),U2beta(maxn),U1betapc2(maxn)
      real U2betapc2(maxn)
      real cumW1U1b(maxn),cumW2U2b(maxn)
      real stproc1ph(0:maxn,maxn)
      real stpr1ph(0:maxn,maxn)
      real stproc2ph(0:maxn,maxn),delta11ph,delta12ph
      real stpr2ph(0:maxn,maxn),u1ph,u4ph,u5ph
      real u6,u7,u6ph,u7ph
      real nrej6,nrej7,pvalue6,pvalue7
      real pvalue52,pvalue53,pvalue54
      real pvalue52ph,pvalue53ph,pvalue54ph
      real sprocph(0:maxn,maxn),u2ph,u3ph
      real pvalue3ph,nrej2ph,nrej3ph  
      real pvalue1ph,pvalue4ph,pvalue5ph
      real nrej1ph,nrej4ph,nrej5ph
      real nrej20,nrej30u20,u30 
      real nrej10,nrej60,nrej70,nrej80
      real u2meanph,u3meanph,u2obsph,u3obsph
      real u1meanph,u4meanph,u5meanph
      real u6meanph,u7meanph,u6mean,u7mean
      real u1obsph,u4obsph,u5obsph
      real u6obsph,u7obsph,u6obs,u7obs
      real u2obs0,u3obs0,u2mean0,u3mean0
      real gprocph(0:maxn,maxn),gproc0(0:maxn,maxn)
      real gtmppr10(0:maxn,maxn)
      real gtmppr2ph(0:maxn,maxn),gtmppr20(0:maxn,maxn)
      real gtproc10(0:maxn,maxn)
      real gtproc2ph(0:maxn,maxn),gtproc20(0:maxn,maxn)
      real delta110,delta120,temp10
      real stproc10(0:maxn,maxn)
      real stpr10(0:maxn,maxn),stpr20(0:maxn,maxn)
      real stproc20(0:maxn,maxn)
      real sproc0(0:maxn,maxn),VEPH
      real u40,u50,u4obs0,u5obs0,nrej40,nrej50
      real u1obs0,u6obs0,u7obs0,u8obs0
      real pvalue6ph,pvalue7ph
      real u8,u9,u8ph,u9ph,u8obs,u9obs,u8obsph,u9obsph
      real u52obs,u53obs,u54obs,u52obsph,u53obsph,u54obsph
      real nrej8,nrej9,nrej8ph,nrej9ph
      real pvalue8,pvalue8ph,pvalue9,pvalue9ph
      real u8mean,u9mean,u8meanph,u9meanph,ZVE
      real pvalue20,pvalue30,pvalue40,pvalue50
      real pvalue10,pvalue60,pvalue70,pvalue80


cc Flag for ran1 function:
cc Set flag to 1 once at start of program
cc When value is 1, iy and iv in ran1 will be set to zero
	ran1flag=1

cc----------------------------------------------------------------------
cc Relationship between variables in the code and in the manuscript: 
cc----------------------------------------------------------------------
cc   In the paper, U^r_1 corresponds to u6 in the program (1-sided alt)
cc                 U^r_2 corresponds to u53 in the program (1-sided alt)
cc                 U^r_3 corresponds to u54 in the program (2-sided alt)
cc                 U^0_1 corresponds to u10 in the program (1-sided alt)
cc                 U^0_2 corresponds to u80 in the program (2-sided alt)


cc u8, u8ph are the 1-sided integrated versions of u5, u5ph
cc u9, u9ph are the 2-sided integrated versions of u5, u5ph

cc "0" at the end of the variable means it is for the
cc test of VE(t,v)=0

cc F1vect1 is F_1(t,v) defined for the group 1 failure times and marks
cc F2vect2 is F_2(t,v) defined for the group 2 failure times and marks
cc F1vect is F_1(t,v) defined for the pooled failure times and marks
cc F2vect is F_2(t,v) defined for the pooled failure times and marks

cc F1vect1cum, F2vect2cum, F1vectcum, F2vectcum are defined similarly
cc for cumulative VE based on Pr(T<=t,V<=v|k)

cc indtt1 is the largest i such that time1(i) .le. T2 and is an event
cc indtt2 is the largest i such that time2(i) .le. T2 and is an event
cc indtt is the largest i such that time(i) .le. T2 and is an event

cc Use ahat(s) = 1/EE2(iicount) and
cc bhat(s) = EE1(iicount)/(EE2(iicount)**2)

cc process is for the 2-sided alternative

cc Hnmk is the weight process for the test statistic
cc U1 for detecting the monotone alternative, calculated
cc at the group k failure times, k=1,2 
cc ttanal is the time at which VE(ttanal,v)=1-F_1(t,v)/F_2(t,v)

cc--------------------------------------
cc Read in data 
cc--------------------------------------
cc Vaccine group was put in vectors first followed by placebo
	 item1=sampleIN(1)
	 item2=sampleIN(2)
	 nsamp1=sampleIN(1)
	 nsamp2=sampleIN(2)
	 nsamp1tr=nsamp1
	 nsamp2tr=nsamp2
	 nsamp=nsamp1+nsamp2

	 time1(1:nsamp1)=timeIN(1:nsamp1)
	 time2(1:nsamp2)=timeIN(nsamp1+1:nsamp)
	 censor1(1:nsamp1)=censorIN(1:nsamp1)
	 censor2(1:nsamp2)=censorIN(nsamp1+1:nsamp)
	 cause1(1:nsamp1)=causeIN(1:nsamp1)
	 cause2(1:nsamp2)=causeIN(nsamp1+1:nsamp)

	 nboot=intIN(1)
	 ngrid=intIN(2)
	 ngridv=intIN(3)
	 TAILSL=intIN(4)
	 TAILSU=intIN(5)
	 TAILSV=intIN(6)

	 T1=doubIN(1)
	 T2=doubIN(2)
	 UT1=doubIN(3)
	 UT2=doubIN(4)
	 ttanal=doubIN(5)
	 BAND1=doubIN(6)
	 BAND2=doubIN(7)
	 BANDV1=doubIN(8)
	 BANDV2=doubIN(9)
	 BANDVLOW=doubIN(10)
	 BANDVUP=doubIN(11)
	 V1=doubIN(12)
	 V2=doubIN(13)
	 UV1=doubIN(14)
	 UV2=doubIN(15)


	
	 OUTV3(1:3)=(/nsamp1,nsamp2,nboot/)
	 
cc ngrid is the number of time points used in the grid for t
cc ngridv is the number of mark points used in the grid for v
cc      nboot:  number of simulated test processes used 
cc              to find P-value (nboot=1000 gives good results).
cc      ntau:  number of observations < tau 
cc                      (tau is not needed as input)
cc              (ntau is used in test statistics u1, u2, u3, u4)
cc      mucen: mean of censoring     

	 ntau=nsamp
 
      seeda=24
      seedb=54
      seedc=23
      seedd=56
      call rstart(seeda,seedb,seedc,seedd)

      ttu=0.0

      num=0.
      den=0
      do i=1,nsamp1
      time(i)=time1(i)
      censor(i)=censor1(i)
      cause(i)=cause1(i)
      group(i)=1
      if (cause1(i) .lt. 2.) then
      num=num+cause1(i)
      den=den+1
      end if
      end do

      num=0.
      den=0
      do i=1,nsamp2    
      time(i+nsamp1)=time2(i)
      censor(i+nsamp1)=censor2(i)
      cause(i+nsamp1)=cause2(i)  
      if (cause2(i) .lt. 2.) then
      num=num+cause2(i)
      den=den+1
      end if
      group(i+nsamp1)=2
      end do


cc
cc     Order the data by time 
cc      
      call indexx(nsamp1,time1,indx)
      do 20 i=1,nsamp1
         copy(i)=time1(i)
 20   continue
      do 25 i=1,nsamp1
         time1(i)=copy(indx(i))
 25   continue
      do 30 i=1,nsamp1
         copy(i)=censor1(i)
 30   continue
      do 35 i=1,nsamp1
         censor1(i)=copy(indx(i))
 35   continue
      do 40 i=1,nsamp1
         copy(i)=cause1(i)
 40   continue
      do 45 i=1,nsamp1
         cause1(i)=copy(indx(i))
 45   continue

            
      call indexx(nsamp2,time2,indx)
      do 50 i=1,nsamp2
         copy(i)=time2(i)
 50   continue
      do 55 i=1,nsamp2
         time2(i)=copy(indx(i))
 55   continue
      do 60 i=1,nsamp2
         copy(i)=censor2(i)
 60   continue
      do 65 i=1,nsamp2
         censor2(i)=copy(indx(i))
 65   continue
      do 70 i=1,nsamp2
         copy(i)=cause2(i)
 70   continue
      do 75 i=1,nsamp2
         cause2(i)=copy(indx(i))
 75   continue


      call indexx(nsamp,time,indx)
      do 77 i=1,nsamp
         copy(i)=time(i)
 77   continue

      icount1=1
      icount2=1
      do 79 i=1,nsamp
         time(i)=copy(indx(i))
         if (indx(i) .le. nsamp1) then
         indstm1(icount1)=i
         icount1=icount1+1
	 do l=1,nsamp
	    covar(1,i,l)=1.
	 end do
         else
         indstm2(icount2)=i
         icount2=icount2+1
	 do l=1,nsamp
	    covar(1,i,l)=0.
         end do
         endif
 79   continue
      do 81 i=1,nsamp
         copy(i)=censor(i)
 81   continue
      do 82 i=1,nsamp
         censor(i)=copy(indx(i))
 82   continue
      do 83 i=1,nsamp
         copy(i)=cause(i)
 83   continue
      do 84 i=1,nsamp
         cause(i)=copy(indx(i))
 84   continue
      do 85 i=1,nsamp
         copy(i)=group(i)
 85   continue
      do 86 i=1,nsamp
         group(i)=copy(indx(i))
 86   continue


cc--------------------------------------------------------
cc   Define the time intervals for kernel estimation 
cc--------------------------------------------------------
cc Define T1 as the smallest failure time such that there are 2 failure times
cc smaller than T1 in each group
cc
cc Define T2 as the largest failure time such that there are 2 failure times
cc larger than T2 in each group
cc
cc Let T1ind1 be the smallest index for time1 that
cc is >= T1:
      
cc    Compute the default T1 and T2 values if they equal 0
	if (T1.eq.0) then
      ind=1
      do 90 i=1,nsamp1
      if (censor1(i) .gt. 0.1) then
      ind=ind+1
      endif
      if ((ind .eq. 2) .and. (censor1(i) .gt. 0.1)) then
      T1=time1(i)
      endif
 90   continue

      ind=1
      do 91 i=1,nsamp2
      if (censor2(i) .gt. 0.1) then
      ind=ind+1
      endif
      if ((ind .eq. 2) .and. (T1 .lt. time2(i)) 
     #    .and. (censor2(i) .gt. 0.1)) then
      T1=time2(i)
      endif
 91   continue
      endif

      if (T2.eq.0.) then
      ind=1
      do 92 i=1,nsamp1
      ii=nsamp1+1-i
      if (censor1(ii) .gt. 0.1) then
      ind=ind+1
      endif
      if ((ind .eq. 2) .and. (censor1(ii) .gt. 0.1)) then
      T2=time1(ii)
      endif
 92   continue

      ind=1
      do 93 i=1,nsamp2
      ii=nsamp2+1-i
      if (censor2(ii) .gt. 0.1) then
      ind=ind+1
      end if
      if ((ind .eq. 2) .and. (censor2(ii) .gt. 0.1)) then
      if (T2 .gt. time2(ii)) then
      T2=time2(ii)
      endif
      endif
 93   continue
      endif

      ind=1
      do 94 i=1,nsamp1
      if ((ind .eq. 1) .and. (time1(i) .ge. T1)) then
      T1ind1=i
      ind=ind+1
      endif
 94   continue

      ind=1
      do 95 i=1,nsamp2
      if ((ind .eq. 1) .and. (time2(i) .ge. T1)) then
      T1ind2=i
      ind=ind+1
      endif
 95   continue

cc--------------------------------------
cc Defaults for UT1, UT2,ttanal
cc--------------------------------------
 
cc Set default UT1 and ensure UT1<T1
      if (UT1 .eq. 0.) then
      UT1=T1/2.
      endif
      if (UT1 .gt. T1) then
      UT1=T1
      endif
      
cc Set default UT2 and ensure UT2<T2
	if (UT2 .eq. 0.) then
      UT2=min(time1(nsamp1),time2(nsamp2))
      endif
      if (UT2 .lt. T2) then
      UT2=T2
      endif

cc Set default for ttanal
      if (ttanal .eq. 0.) then
      ttanal=T2
      endif

cc Count percentage of uncensored before analysis time-point
	do 96 i=1,ntau
      if (time1(i) .le. ttanal) then
      ttu=ttu+censor1(i)
      end if
 96   continue

      do 97 i=1,ntau
      if (time2(i) .le. ttanal) then
      ttu=ttu+censor2(i)
      end if
 97   continue

cc    Set the bandwidth to 1/4th of the follow-up period 
      BAND=(T2-T1)/float(4)


      do i=1,nsamp1
        if (time1(i) .le. T2 .and. censor1(i) .gt. 0.1) then
        indtt1=i
        end if
      end do

      do i=1,nsamp2
        if (time2(i) .le. T2 .and. censor2(i) .gt. 0.1) then
        indtt2=i
        end if
      end do

      do i=1,nsamp
        if (time(i) .le. T2 .and. censor(i) .gt. 0.1) then
        indtt=i
        end if
      end do


cc-------------------------------------------
cc    Number of failures, average mark
cc-------------------------------------------

cc    Order the causes of failure 
      call indexx(nsamp1,cause1,indx)
      do i=1,nsamp1
        ocause1(i)=cause1(indx(i))
      end do

cc   nv1 = number of causes of failure < or = 1
cc         for Group 1.
cc       = number of uncensored observations
cc         for Group 1.
      nv1=0
      vavg1=0.
      do 100 i=1,nsamp1
         if (ocause1(i).le.1.) then
            nv1=nv1+1
            vavg1=vavg1+ocause1(i)
         endif
 100  continue
      vavg1=vavg1/float(nv1)


      call indexx(nsamp2,cause2,indx)
      do 101 i=1,nsamp2
         ocause2(i)=cause2(indx(i))
 101   continue

cc   nv2 = number of causes of failure < or = 1
cc         for Group 2.
cc       = number of uncensored observations
cc         for Group 2.
      nv2=0
      vavg2=0.
      do 102 i=1,nsamp2
         if (ocause2(i).le.1.) then
            nv2=nv2+1
            vavg2=vavg2+ocause2(i)
         endif
 102  continue
      vavg2=vavg2/float(nv2)

	 
	 OUTV30(1:5)=(/T1,T2,ttanal,vavg1,vavg2/)

      call indexx(nsamp,cause,indx)
      icount1=1
      icount2=1
      do 103 i=1,nsamp
         ocause(i)=cause(indx(i))
         obscause(i)=ocause(i)
         if (ocause(i) .le. 1.) then
         if (group(indx(i)) .eq. 1) then
         indsv1(icount1)=i
         icount1=icount1+1
         else
         indsv2(icount2)=i
         icount2=icount2+1
         endif
         endif
 103  continue


cc   nv = number of causes of failure < or = 1
cc       in both groups
cc      = number of uncensored observations in both groups
      nv=0
      do 104 i=1,nsamp
         if (ocause(i).le.1.) then
            nv=nv+1
         endif
 104  continue


cc---------------------------------------------------------------------
cc Setting defaults of V1,V2,UV1,UV2,BANDVUP,BANDVLOW if not specified
cc---------------------------------------------------------------------
cc Set V1 to be the larger of the smallest
cc observed marks for groups 1 and 2:
       if (V1 .eq. 0.) then 
	 V1=ocause1(1)
       temp=ocause2(1)
       if (temp .gt. V1) then
	 V1=temp
       endif
       endif

cc Set V2 to be the smaller of the two largest
cc observed marks for groups 1 and 2:

       if (V2 .eq. 0.) then 
	 V2=ocause1(nv1)
       temp=ocause2(nv2)
       if (temp .lt. V2) then
       V2=temp
       endif
       endif

cc Set default for UV1 and ensure UV1>V1
	 if (UV1 .eq. 0.) then
       UV1=V1
	 endif
	 if (UV1 .gt. V1) then
       UV1=V1
	 endif


cc Set default for UV2 and ensure UV2>V2
	 if (UV2 .eq. 0.) then
       UV2=V2
	 endif
	 if (UV2 .lt. V2) then
       UV2=V2
	 endif

cc      Range of bandwidths to optimize over with 2-fold CV
      if (BANDVUP .eq. 0.) then
      BANDVUP=V2/2.
      endif
      if (BANDVLOW .eq. 0.) then
      BANDVLOW=V1+(2.*BANDVUP-V1)*.05
      endif

cc-------------------------------------
cc Nonparametric estimation
cc-------------------------------------

cc  compute: 1. increments in the
cc              nonparametric estimators 
cc              hatl1(i,j) and hatl2(i,j)
cc           2. the overall nonparametric
cc              kernel hazard estimates 
cc              ee1 and ee2 (\hat \lambda_1
cc              and \hat \lambda_2)
cc           3. the test process process(i,j)
cc

         NP=1
         do kl=1,10

         call estp(nsamp,NP,time,covar,censor,s0,s1,s2,beta,
     $    BZ,CL,GS,U12,U11,var,NDEAD,NMIS,maxn)
         end do

cc See if reject H0: VE = 0 using standard Cox model
         VEPH=1.-exp(beta)
	   ZVE=beta/sqrt(var)


cc Compute Jnbeta, U1betai, and U2betai in the paper

cc NOTE: If problems, remove the restriction that time is in [T1,T2]
      temp1=0.
      do 108 i=1,nsamp1
      U1beta(i)=0.
      U1betapc2(i)=0.


       if (censor1(i) .gt. 0.1) then 
        y1j=0.
        y2j=0.
        do 114 k=1,nsamp1
          if (time1(k).ge.time1(i)) then
            y1j=y1j+1.
          end if
 114    continue

        do 115 k=1,nsamp2
          if (time2(k).ge.time1(i)) then
            y2j=y2j+1.
          end if
 115    continue
	temp1=temp1+(y1j*y2j*exp(beta))/
     $        ((y1j*exp(beta) + y2j)**2.) 
        U1beta(i)=y2j/(y1j*exp(beta) + y2j)
      end if


 108  continue

      temp2=0.
      do 118 i=1,nsamp2
        U2beta(i)=0.
        U2betapc2(i)=0.

cc Piece 2
       if (censor2(i) .gt. 0.1) then
        y1j=0.
        y2j=0.
        do 124 k=1,nsamp1
          if (time1(k).ge.time2(i)) then
            y1j=y1j+1.
          end if
 124    continue
        do 125 k=1,nsamp2
          if (time2(k).ge.time2(i)) then
            y2j=y2j+1.
          end if
 125    continue
	temp2=temp2+(y1j*y2j*exp(beta))/
     $          ((y1j*exp(beta) + y2j)**2.) 
        U2beta(i)=(y1j*exp(beta))/(y1j*exp(beta) + y2j)
      end if
 118  continue

      Jnbeta=temp1+temp2


      do 126 j=1,nv
         do 126 i=1,nsamp
         process(i,j)=0.
 126  continue


      do 127 j=1,nv+1
            tmpproc1(0,j)=0.
	    tmpproc2(0,j)=0.
 127  continue

       do 200 j=1,nv1
         do 200 i=1,nsamp1
           if (cause1(i).le.ocause1(j)) then
              hatl1(i,j)=censor1(i)/float(nsamp1-i+1)
           else
              hatl1(i,j)=0.
           endif
 200  continue


       ind=0
       do 220 j=1,nv2
         do 220 i=1,nsamp2
           if (cause2(i).le.ocause2(j)) then
              hatl2(i,j)=censor2(i)/float(nsamp2-i+1)
              if (ocause2(j) .gt. vavg2 .and. ind .eq. 0) then
                javg=j
                ind=ind+1
              end if
           else
              hatl2(i,j)=0.
           endif

 220  continue


        if (hatl2(1,javg) .gt. 0.) then
        Hnm2(1)=float(1)/hatl2(1,javg)
        else 
        Hnm2(1)=0.
        end if
      do 221 i=2,nsamp2
        if (hatl2(i,javg) .gt. 0.) then
        Hnm2(i)=float(1)/hatl2(i,javg)
        else
        Hnm2(i)=Hnm2(i-1)
        end if
 221  continue


cc For now, put Hnm2 = 1 always
      do 222 i=1,nsamp2
        Hnm2(i)=1.
 222  continue


cc Now compute Hnm1 (at the failure times for group 1)
cc Difficulty: Hnm2 is only defined at the group 2 failure times

cc first for i=1:
           ind=0
           do 225 j2=1,nv2
             do 225 i2=1,nsamp2
             if (cause2(i2).le.ocause2(j2)) then
             if (ind .eq. 0) then
             ind=ind+1
             if (time1(1) .ge. time2(i2)
     #       .and. hatl2(i2,javg) .gt. 0.) then
             Hnm1(1)=float(1)/hatl2(i2,javg)
             else
             Hnm1(1)=0.
             end if
             end if
             end if
 225  continue

         do 226 i=2,nsamp1
           oldind=1
           ind=0
           do 227 i2=1,nsamp2
           do 227 j2=1,nv2
             if (cause2(i2).le.ocause2(j2)) then
             if (time2(i2) .gt. time1(i)) then
             if (ind .eq. 0) then
             ind=ind+1
             if (hatl2(oldind,javg) .gt. 0.) then
             Hnm1(i)=float(1)/hatl2(oldind,javg)
             else
             Hnm1(i)=Hnm1(i-1)
             end if
             end if
             end if
             oldind=i2
             end if
 227  continue  
 226  continue

cc For now just set Hnm1 = 1 always

      do 228 i=1,nsamp1
        Hnm1(i)=1.
 228  continue

cc--------------------------------------------------
cc  Nonparametric kernel estimates 
cc--------------------------------------------------

cc First compute the Nelson-Aalen estimates

      KMC1=1.0
      KMC2=1.0
      KM01=1.0
      KM02=1.0
      KMT1(0)=1.
      KMT2(0)=1.
      NA1(0)=0.0
      NA2(0)=0.0

      DO 232 i=1,nsamp1  

cc Compute group-specific Kaplan-Meier curves:
        if (censor1(i) .lt. 0.1) then
          KMC1=KMC1*(1.0-(1.0/float(nsamp1-i+1)))
        else
          KM01=KM01*(1.0-(1.0/float(nsamp1-i+1)))
        end if
          KMT1(i)=KM01 
          
        if (censor1(i) .gt. 0.1) then
           JNA=1.0/float(nsamp1-i+1)
        else
           JNA=0.0
        endif
         NA1(i)=NA1(i-1)+JNA
        if (time1(i) .le. T2) then
c        F1T2=1.0-exp(-NA1(i))
        end if
        if (censor1(i) .gt. 0.1) then
        endif
 232  CONTINUE 

      DO 237 I=1,nsamp2  
        if (censor2(i) .lt. 0.1) then
          KMC2=KMC2*(1.0-(1.0/float(nsamp2-i+1)))
        else
          KM02=KM02*(1.0-(1.0/float(nsamp2-i+1)))
        end if
          KMT2(i)=KM02
        
        if (censor2(i) .gt. 0.1) then
           JNA=1.0/float(nsamp2-i+1)
        else
           JNA=0.0
        endif
         NA2(i)=NA2(i-1)+JNA
        if (time2(i) .le. T2) then
c        F2T2=1.0-exp(-NA2(i))
        end if
        if (censor2(i) .gt. 0.1) then
        endif
 237  CONTINUE 

cc----------------------------------
cc Compute the optimal bandwidths
cc----------------------------------
      BAND111=(T2-T1)/4.
      BAND112=(T2-T1)/4.

      k2K=0.20
      intK2=0.6

      delt11=0.
      delt12=0.
      delt21=0.
      delt22=0.


      do 238 i=1,nsamp1
      if (censor1(i).gt. 0.01 .and. time1(i).le.T1) then
      delt11=delt11+(1./float(nsamp1-i+1))**2.
      end if
      if (censor1(i).gt. 0.01 .and. time1(i).le.T2) then
      delt12=delt12+(1./float(nsamp1-i+1))**2.
      end if
 238  continue

      do 239 i=1,nsamp2
      if (censor2(i).gt. 0.01 .and. time2(i).le.T1) then
      delt21=delt21+(1./float(nsamp2-i+1))**2.
      end if
      if (censor2(i).gt. 0.01 .and. time2(i).le.T2) then
      delt22=delt22+(1./float(nsamp2-i+1))**2.
      end if
 239  continue

      alpy1=nsamp1*(delt12-delt11)
      alpy2=nsamp2*(delt22-delt21)

      intalpdp1=alphdblp(time1(T1ind1),band111,time1,censor1,na1,
     #          nsamp1tr,T1ind1,maxn)**2.

      ind=T1ind1
      do 240 i=1,nsamp1
      if ((time1(i).ge.T1) .and. (time1(i).le.T2)) then
      if (censor1(i).gt. 0.01) then
      intalpdp1=intalpdp1+(alphdblp(time1(i),band111,time1,censor1,
     #        na1,nsamp1tr,T1ind1,maxn)**2.)*(time1(i)-time1(ind))
      ind=i
      end if
      end if
 240  continue


      intalpdp2=alphdblp(time2(T1ind2),band112,time2,censor2,na2,
     #          nsamp2tr,T1ind2,maxn)**2.
	ind=T1ind2
      do 241 i=(T1ind2+1),nsamp2
      if ((time2(i).ge.T1) .and. (time2(i).le.T2)) then
      if (censor2(i) .gt. 0.01) then
      intalpdp2=intalpdp2+(alphdblp(time2(i),band112,time2,censor2,
     #       na2,nsamp2tr,T1ind2,maxn)**2.)*(time2(i)-time2(ind))

      ind=i
      end if
      end if
 241  continue

	if (BAND1 .eq. 0.) then
      BAND1=(1./(k2K**0.4))*((intK2*alpy1)**0.2)*
     #      (1./(intalpdp1**0.2))*(1./(nsamp1**0.2))
      endif
      if (BAND2 .eq. 0.) then
      BAND2=(1./(k2K**0.4))*((intK2*alpy2)**0.2)*
     #      (1./(intalpdp2**0.2))*(1./(nsamp2**0.2))
      endif
	 
	OUTV30(6:7)=(/BAND1,BAND2/) 

      VEoverall=1.0-(1.0-KMT1(indtt1))/(1.0-KMT2(indtt2))

cc Compute EE1, accommodate dealing with the tail problem:
      DO 252 I=1,nsamp2
       EE1(I)=0.
       IF (time2(I) .GE. UT1 .and. time2(I) .LE. UT2) THEN
       IF ((time2(I) .GE. (UT1+BAND1)) .AND.
     #     (time2(I) .LE. (UT2-BAND1))) THEN

       DO 253 L=1,nsamp1tr
       IF (ABS(censor1(L)).LT.0.00000001) GOTO 253

         EE1(I)=EE1(I)+EPAN((time2(I)-time1(L))/BAND1)*
     #                 ((NA1(L)-NA1(L-1))/BAND1)
 253   continue

      ELSE IF (time2(I) .LT. (UT1+BAND1)) THEN

      IF (TAILSL .EQ. 2) THEN

       DO 254 L=1,nsamp1tr                                                                                       
       IF (ABS(censor1(L)).LT.0.00000001) GOTO 254
         EE1(I)=EE1(I)+EPAN((time2(I)-time1(L))/BAND1)*
     #                 ((NA1(L)-NA1(L-1))/BAND1)
 254   continue
      ELSE
       DO 255 L=1,nsamp1tr
       IF (ABS(censor1(L)).LT.0.00000001) GOTO 255

         EE1(I)=EE1(I)+EPANL((time2(I)-time1(L))/BAND1,
     #   (time2(I)-UT1)/BAND1)*((NA1(L)-NA1(L-1))/BAND1)
 255  continue
      ENDIF

      ELSE IF (time2(I) .gt. (UT2-BAND1)) THEN

      IF (TAILSU .EQ. 2) THEN
       DO 256 L=1,nsamp1tr
       IF (ABS(censor1(L)).LT.0.00000001) GOTO 256

         EE1(I)=EE1(I)+EPAN((time2(I)-time1(L))/BAND1)*
     #                 ((NA1(L)-NA1(L-1))/BAND1)
 256   continue
      ELSE
       DO 257 L=1,nsamp1tr
       IF (ABS(censor1(L)).LT.0.00000001) GOTO 257

         EE1(I)=EE1(I)+EPANU((time2(I)-time1(L))/BAND1,
     #   -(UT2-time2(I))/BAND1)*((NA1(L)-NA1(L-1))/BAND1)
 257   continue
      ENDIF
      ENDIF
      ENDIF
 
 252  CONTINUE

cc Compute EE2, accommodate dealing with the tail problem:

      DO 262 I=1,nsamp2
       EE2(I)=0.

       IF (time2(I) .GE. UT1 .AND. time2(I) .LE. UT2) THEN
       IF ((time2(I) .GE. (UT1+BAND2)) .AND.
     #     (time2(I) .LE. (UT2-BAND2))) THEN

       DO 263 L=1,nsamp2tr
       IF (ABS(censor2(L)).LT.0.00000001) GOTO 263

         EE2(I)=EE2(I)+EPAN((time2(I)-time2(L))/BAND2)*
     #                 ((NA2(L)-NA2(L-1))/BAND2)
 263   continue

      ELSE IF (time2(I) .LT. (UT1+BAND2)) THEN

      IF (TAILSL .EQ. 2) THEN

       DO 264 L=1,nsamp2tr                                                                                      
       IF (ABS(censor2(L)).LT.0.00000001) GOTO 264
         EE2(I)=EE2(I)+EPAN((time2(I)-time2(L))/BAND2)*
     #                 ((NA2(L)-NA2(L-1))/BAND2)
 264   continue

      ELSE

       DO 265 L=1,nsamp2tr                                                                                      
       IF (ABS(censor2(L)).LT.0.00000001) GOTO 265
         EE2(I)=EE2(I)+EPANL((time2(I)-time2(L))/BAND2,
     #           (time2(I)-UT1)/BAND2)*((NA2(L)-NA2(L-1))/BAND2)
 265   continue
      ENDIF

      ELSE IF (time2(I) .gt. (UT2-BAND2)) THEN

      IF (TAILSU .EQ. 2) THEN

       DO 266 L=1,nsamp2tr                                                                                      
       IF (ABS(censor2(L)).LT.0.00000001) GOTO 266

         EE2(I)=EE2(I)+EPAN((time2(I)-time2(L))/BAND2)*
     #                 ((NA2(L)-NA2(L-1))/BAND2)
 266   continue
      ELSE

       DO 267 L=1,nsamp2tr                                                                                      
       IF (ABS(censor2(L)).LT.0.00000001) GOTO 267

         EE2(I)=EE2(I)+EPANU((time2(I)-time2(L))/BAND2,
     #   -(UT2-time2(I))/BAND2)*((NA2(L)-NA2(L-1))/BAND2)
 267   continue
      ENDIF
      ENDIF
      ENDIF
 262  CONTINUE


cc Also calculate EE2 at the failure times of group 1:

cc Deal with the tail issue:
      DO 272 I=1,nsamp1
       EE21(I)=0.

       IF (time1(I) .GE. UT1 .AND. time1(I) .LE. UT2) THEN
       IF ((time1(I) .GE. (UT1+BAND2)) .AND.
     #     (time1(I) .LE. (UT2-BAND2))) THEN

       DO 273 L=1,nsamp2tr
       IF (ABS(censor2(L)).LT.0.00000001) GOTO 273

         EE21(I)=EE21(I)+EPAN((time1(I)-time2(L))/BAND2)*
     #                 ((NA2(L)-NA2(L-1))/BAND2)
 273   continue

      ELSE IF (time1(I) .LT. (UT1+BAND2)) THEN

      IF (TAILSL .EQ. 2) THEN

       DO 274 L=1,nsamp2tr
       IF (ABS(censor2(L)).LT.0.00000001) GOTO 274
         EE21(I)=EE21(I)+EPAN((time1(I)-time2(L))/BAND2)*
     #                 ((NA2(L)-NA2(L-1))/BAND2)
 274   continue

      ELSE

       DO 275 L=1,nsamp2tr
       IF (ABS(censor2(L)).LT.0.00000001) GOTO 275
         EE21(I)=EE21(I)+EPANL((time1(I)-time2(L))/BAND2,
     #   (time1(I)-UT1)/BAND2)*((NA2(L)-NA2(L-1))/BAND2)
 275   continue

      ENDIF

      ELSE IF (time1(I) .GT. (UT2-BAND2)) THEN

      IF (TAILSU .EQ. 2) THEN
       DO 276 L=1,nsamp2tr
       IF (ABS(censor2(L)).LT.0.00000001) GOTO 276
         EE21(I)=EE21(I)+EPAN((time1(I)-time2(L))/BAND2)*
     #                 ((NA2(L)-NA2(L-1))/BAND2)
 276   continue
      ELSE


       DO 277 L=1,nsamp2tr
       IF (ABS(censor2(L)).LT.0.00000001) GOTO 277
         EE21(I)=EE21(I)+EPANU((time1(I)-time2(L))/BAND2,
     #   -(UT2-time1(I))/BAND2)*((NA2(L)-NA2(L-1))/BAND2)
 277  continue
      ENDIF
      ENDIF
      ENDIF
 272  CONTINUE

cc----------------------------------------------------------------
cc Compute 2-Fold Cross-validation MISE for the bandwidth, group 1
cc----------------------------------------------------------------
cc BANDV1
      if (BANDV1 .eq. 0.) then
      MISEOLD=1.
      do 278 il=1,20
      BANDV1=BANDVLOW+float(il-1)*(BANDVUP-BANDVLOW)/
     #                float(20-1)
cc Randomly split the data in half:
      do i=1,nsamp1
      rannums(i)=ran1(seeda,ran1flag)
      ran1flag=0 
	indsCV1(i)=0
      end do

cc      Order the indices 1, ..., n1 by the random #s
      call indexx(nsamp1,rannums,indx)

      ihalf=nsamp1/2
      do i=1,ihalf
      indsCV1(indx(i))=1
      end do

cc indsCV1 is of length nsamp1, with 1s indicating the first half of
cc the group 1 sample and 0s indicating the second half

cc Compute estimates on first half:
       do 279 j=1,ngridv
       vv=V1+float(j-1)*(V2-V1)/float(ngridv-1)
         do 280 I=1,nsamp1
         F1vectCV1(I,j)=0.
         do 281 L=1,nsamp1tr
       IF ((ABS(censor1(L)).LT.0.00000001) .or.
     #     (time1(L) .gt. time1(I)) .or.
     #     (indsCV1(L) .eq. 0)) GOTO 281
cc indsCV1(L) = 0 means that the lth point is not in the first sample

       IF ((cause1(L) .ge. UV1) .and.
     #     (cause1(L) .le. UV2)) THEN
       IF ((vv .GE. (UV1+BANDV1)) .AND.
     #     (vv .LE. (UV2-BANDV1))) THEN
        F1vectCV1(I,j)=F1vectCV1(I,j)+EPAN((vv-cause1(L))/BANDV1)*
     #               (KMT1(L)/float(nsamp1-L+1))/BANDV1
        ELSE IF (vv .LT. (UV1+BANDV1)) THEN

        IF (TAILSV .eq. 1) THEN
cc Deal with the tail problem
        F1vectCV1(I,j)=F1vectCV1(I,j)+
     #  EPANL((vv-cause1(L))/BANDV1,(vv-UV1)/BANDV1)*
     #               (KMT1(L)/float(nsamp1-L+1))/BANDV1
        ELSE
        F1vectCV1(I,j)=F1vectCV1(I,j)+EPAN((vv-cause1(L))/BANDV1)*
     #               (KMT1(L)/float(nsamp1-L+1))/BANDV1
        ENDIF

        ELSE IF (vv .GT. (UV2-BANDV1)) THEN

        IF (TAILSV .eq. 1) THEN
        F1vectCV1(I,j)=F1vectCV1(I,j)+
     #  EPANU((vv-cause1(L))/BANDV1,-(UV2-vv)/BANDV1)*
     #               (KMT1(L)/float(nsamp1-L+1))/BANDV1
        ELSE
        F1vectCV1(I,j)=F1vectCV1(I,j)+EPAN((vv-cause1(L))/BANDV1)*
     #               (KMT1(L)/float(nsamp1-L+1))/BANDV1
        ENDIF

        ENDIF 
        ENDIF
 281  CONTINUE
 280  CONTINUE
 279  CONTINUE

cc Have F1vectCV1 for the first half

cc Now compute F1vectCV2 for the second half

       do 282 j=1,ngridv
       vv=V1+float(j-1)*(V2-V1)/float(ngridv-1)
         do 283 I=1,nsamp1
         F1vectCV2(I,j)=0.
         do 284 L=1,nsamp1tr
       IF ((ABS(censor1(L)).LT.0.00000001) .or.
     #     (time1(L) .gt. time1(I)) .or.
     #     (indsCV1(L) .eq. 1)) GOTO 284
cc indsCV1(L) = 0 means that the lth point is not in the first sample

       IF ((cause1(L) .ge. UV1) .and.
     #     (cause1(L) .le. UV2)) THEN
       IF ((vv .GE. (UV1+BANDV1)) .AND.
     #     (vv .LE. (UV2-BANDV1))) THEN
        F1vectCV2(I,j)=F1vectCV2(I,j)+EPAN((vv-cause1(L))/BANDV1)*
     #               (KMT1(L)/float(nsamp1-L+1))/BANDV1
        ELSE IF (vv .LT. (UV1+BANDV1)) THEN

        IF (TAILSV .eq. 1) THEN
cc Deal with the tail problem
        F1vectCV2(I,j)=F1vectCV2(I,j)+
     #  EPANL((vv-cause1(L))/BANDV1,(vv-UV1)/BANDV1)*
     #               (KMT1(L)/float(nsamp1-L+1))/BANDV1
        ELSE
        F1vectCV2(I,j)=F1vectCV2(I,j)+EPAN((vv-cause1(L))/BANDV1)*
     #               (KMT1(L)/float(nsamp1-L+1))/BANDV1
        ENDIF

        ELSE IF (vv .GT. (UV2-BANDV1)) THEN

        IF (TAILSV .eq. 1) THEN
        F1vectCV2(I,j)=F1vectCV2(I,j)+
     #  EPANU((vv-cause1(L))/BANDV1,-(UV2-vv)/BANDV1)*
     #               (KMT1(L)/float(nsamp1-L+1))/BANDV1
        ELSE
        F1vectCV2(I,j)=F1vectCV2(I,j)+EPAN((vv-cause1(L))/BANDV1)*
     #               (KMT1(L)/float(nsamp1-L+1))/BANDV1
        ENDIF

        ENDIF 
        ENDIF
 284  CONTINUE
 283  CONTINUE
 282  CONTINUE

cc Compute MISE for CV sample 1 predicting sample 2

      MISE=0.
      do 285 j=1,ngridv
       do 285 I=1,nsamp1tr
        IF (indsCV1(I) .eq. 0) THEN
        MISE=MISE + (F1vectCV1(I,j) - F1vectCV2(I,j))**2.
        ENDIF
 285  CONTINUE

cc Vice versa, Compute MISE for CV sample 2 predicting sample 1

      do 286 j=1,ngridv
        do 286  I=1,nsamp1tr
        IF (indsCV1(I) .eq. 1) THEN
        MISE=MISE + (F1vectCV2(I,j) - F1vectCV1(I,j))**2.
        ENDIF
 286  CONTINUE

      MISE=MISE/(ngridv*nsamp1tr)
      MISE=MISE/2.
      IF (MISE .lt. MISEOLD) THEN
	TBANDV1=BANDV1
      END IF 
      MISEOLD=MISE
278   continue

      BANDV1=TBANDV1
      endif

cc Now compute the estimates

cc adjust for the tail problem

       do 287 j=1,ngridv
       vv=V1+float(j-1)*(V2-V1)/float(ngridv-1)
         do 288 I=1,nsamp1
         F1vect1(I,j)=0.
         vF1vect1(I,j)=0.
	 F1vect1cum(I,j)=0.
	 vF1vect1cum(I,j)=0.
         do 289 L=1,nsamp1tr
       IF ((ABS(censor1(L)).LT.0.00000001) .or.
     #     (time1(L) .gt. time1(I))) GOTO 289
       IF (cause1(L) .le. vv) THEN
        F1vect1cum(I,j)=F1vect1cum(I,j)+
     #                  (KMT1(L)/float(nsamp1-L+1))
        vF1vect1cum(I,j)=vF1vect1cum(I,j)+
     #   ((KMT1(L)/float(nsamp1-L+1)))**2.
       END IF

       IF ((cause1(L) .ge. UV1) .and.
     #     (cause1(L) .le. UV2)) THEN
       IF ((vv .GE. (UV1+BANDV1)) .AND.
     #     (vv .LE. (UV2-BANDV1))) THEN
        F1vect1(I,j)=F1vect1(I,j)+EPAN((vv-cause1(L))/BANDV1)*
     #               (KMT1(L)/float(nsamp1-L+1))/BANDV1
        vF1vect1(I,j)=vF1vect1(I,j)+(EPAN((vv-cause1(L))/
     #  BANDV1)*((KMT1(L)/float(nsamp1-L+1))/BANDV1))**2.

        ELSE IF (vv .LT. (UV1+BANDV1)) THEN

        IF (TAILSV .eq. 1) THEN
cc Deal with the tail problem
        F1vect1(I,j)=F1vect1(I,j)+
     #  EPANL((vv-cause1(L))/BANDV1,(vv-UV1)/BANDV1)*
     #               (KMT1(L)/float(nsamp1-L+1))/BANDV1
        vF1vect1(I,j)=vF1vect1(I,j)+
     #  (EPANL((vv-cause1(L))/BANDV1,(vv-UV1)/BANDV1)*
     #                ((KMT1(L)/float(nsamp1-L+1))/BANDV1))**2.
        ELSE
        F1vect1(I,j)=F1vect1(I,j)+EPAN((vv-cause1(L))/BANDV1)*
     #               (KMT1(L)/float(nsamp1-L+1))/BANDV1
        vF1vect1(I,j)=vF1vect1(I,j)+(EPAN((vv-cause1(L))/
     #  BANDV1)*((KMT1(L)/float(nsamp1-L+1))/BANDV1))**2.
        ENDIF

        ELSE IF (vv .GT. (UV2-BANDV1)) THEN

        IF (TAILSV .eq. 1) THEN
        F1vect1(I,j)=F1vect1(I,j)+
     #  EPANU((vv-cause1(L))/BANDV1,-(UV2-vv)/BANDV1)*
     #               (KMT1(L)/float(nsamp1-L+1))/BANDV1
        vF1vect1(I,j)=vF1vect1(I,j)+
     #  (EPANU((vv-cause1(L))/BANDV1,-(UV2-vv)/
     #  BANDV1)*((KMT1(L)/float(nsamp1-L+1))/BANDV1))**2.
        ELSE
        F1vect1(I,j)=F1vect1(I,j)+EPAN((vv-cause1(L))/BANDV1)*
     #               (KMT1(L)/float(nsamp1-L+1))/BANDV1
        vF1vect1(I,j)=vF1vect1(I,j)+(EPAN((vv-cause1(L))/
     #  BANDV1)*((KMT1(L)/float(nsamp1-L+1))/BANDV1))**2.
        ENDIF

        ENDIF 
        ENDIF
 289  CONTINUE
 288  CONTINUE
 287  CONTINUE

cc------------------------------------------------------------------
cc Compute 2-Fold Cross-validation MISE for the bandwidth, group 2
cc------------------------------------------------------------------
cc BANDV2

      if (BANDV2 .eq. 0.) then
      MISEOLD=1.
      do 293 il=1,20
      BANDV2=BANDVLOW+float(il-1)*(BANDVUP-BANDVLOW)/
     #                float(20-1)

cc Randomly split the data in half:
      do i=1,nsamp2
      rannums(i)=ran1(seeda,ran1flag)
      indsCV1(i)=0
      end do

cc      Order the indices 1, ..., n1 by the random #s
      call indexx(nsamp2,rannums,indx)

      ihalf=nsamp2/2.
      do i=1,ihalf
      indsCV1(i)=indx(i)
      end do

cc indsCV1 is of length nsamp2, with 1s indicating the first half of
cc the group 1 sample and 0s indicating the second half

cc Compute estimates on first half:
       do 294 j=1,ngridv
       vv=V1+float(j-1)*(V2-V1)/float(ngridv-1)
         do 295 I=1,nsamp2
         F2vectCV1(I,j)=0.
         do 296 L=1,nsamp2tr
       IF ((ABS(censor2(L)).LT.0.00000001) .or.
     #     (time2(L) .gt. time2(I)) .or.
     #     (indsCV1(L) .eq. 0)) GOTO 296
cc indsCV1(L) = 0 means that the lth point is not in the first sample

       IF ((cause2(L) .ge. UV1) .and.
     #     (cause2(L) .le. UV2)) THEN
       IF ((vv .GE. (UV1+BANDV2)) .AND.
     #     (vv .LE. (UV2-BANDV2))) THEN
        F2vectCV1(I,j)=F2vectCV1(I,j)+EPAN((vv-cause2(L))/BANDV2)*
     #               (KMT2(L)/float(nsamp2-L+1))/BANDV2
        ELSE IF (vv .LT. (UV1+BANDV2)) THEN

        IF (TAILSV .eq. 1) THEN
cc Deal with the tail problem
        F2vectCV1(I,j)=F2vectCV1(I,j)+
     #  EPANL((vv-cause2(L))/BANDV2,(vv-UV1)/BANDV2)*
     #               (KMT2(L)/float(nsamp2-L+1))/BANDV2
        ELSE
        F2vectCV1(I,j)=F2vectCV1(I,j)+EPAN((vv-cause2(L))/BANDV2)*
     #               (KMT2(L)/float(nsamp2-L+1))/BANDV2
        ENDIF

        ELSE IF (vv .GT. (UV2-BANDV2)) THEN

        IF (TAILSV .eq. 1) THEN
        F2vectCV1(I,j)=F2vectCV1(I,j)+
     #  EPANU((vv-cause2(L))/BANDV2,-(UV2-vv)/BANDV2)*
     #               (KMT2(L)/float(nsamp2-L+1))/BANDV2
        ELSE
        F2vectCV1(I,j)=F2vectCV1(I,j)+EPAN((vv-cause2(L))/BANDV2)*
     #               (KMT2(L)/float(nsamp2-L+1))/BANDV2
        ENDIF

        ENDIF 
        ENDIF
 296  CONTINUE
 295  CONTINUE
 294  CONTINUE

cc Have F2vectCV1 for the first half
cc Now compute F2vectCV2 for the second half

       do 297 j=1,ngridv
       vv=V1+float(j-1)*(V2-V1)/float(ngridv-1)
         do 298 I=1,nsamp2
         F2vectCV2(I,j)=0.
         do 299 L=1,nsamp2tr
       IF ((ABS(censor2(L)).LT.0.00000001) .or.
     #     (time2(L) .gt. time2(I)) .or.
     #     (indsCV1(L) .eq. 1)) GOTO 299
cc indsCV1(L) = 0 means that the lth point is not in the first sample

       IF ((cause2(L) .ge. UV1) .and.
     #     (cause2(L) .le. UV2)) THEN
       IF ((vv .GE. (UV1+BANDV2)) .AND.
     #     (vv .LE. (UV2-BANDV2))) THEN
        F2vectCV2(I,j)=F2vectCV2(I,j)+EPAN((vv-cause2(L))/BANDV2)*
     #               (KMT2(L)/float(nsamp2-L+1))/BANDV2
        
        ELSE IF (vv .LT. (UV1+BANDV2)) THEN

        IF (TAILSV .eq. 1) THEN
cc Deal with the tail problem
        F2vectCV2(I,j)=F2vectCV2(I,j)+
     #  EPANL((vv-cause2(L))/BANDV2,(vv-UV1)/BANDV2)*
     #               (KMT2(L)/float(nsamp2-L+1))/BANDV2
        ELSE
        F2vectCV2(I,j)=F2vectCV2(I,j)+EPAN((vv-cause2(L))/BANDV2)*
     #               (KMT2(L)/float(nsamp2-L+1))/BANDV2
        ENDIF

        ELSE IF (vv .GT. (UV2-BANDV2)) THEN

        IF (TAILSV .eq. 1) THEN
        F2vectCV2(I,j)=F2vectCV2(I,j)+
     #  EPANU((vv-cause2(L))/BANDV2,-(UV2-vv)/BANDV2)*
     #               (KMT2(L)/float(nsamp2-L+1))/BANDV2
        ELSE
        F2vectCV2(I,j)=F2vectCV2(I,j)+EPAN((vv-cause2(L))/BANDV2)*
     #               (KMT2(L)/float(nsamp2-L+1))/BANDV2
        ENDIF

        ENDIF 
        ENDIF
 299  CONTINUE
 298  CONTINUE
 297  CONTINUE

cc Compute MISE for CV sample 1 predicting sample 2

      MISE=0.
      do 300 j=1,ngridv
       vv=V1+float(j-1)*(V2-V1)/float(ngridv-1)
       do 300 I=1,nsamp2tr
        IF (indsCV1(I) .eq. 0) THEN
        MISE=MISE + (F2vectCV1(I,j) - F2vectCV2(I,j))**2.
        ENDIF
 300  CONTINUE

cc Vice versa, Compute MISE for CV sample 2 predicting sample 1

      do 301 j=1,ngridv
       vv=V1+float(j-1)*(V2-V1)/float(ngridv-1)
         do 301 I=1,nsamp2tr
        IF (indsCV1(I) .eq. 1) THEN
        MISE=MISE + (F2vectCV2(I,j) - F2vectCV1(I,j))**2.
        ENDIF
 301  CONTINUE

      MISE=MISE/(ngridv*nsamp2tr)
      MISE=MISE/2.

      IF (MISE .lt. MISEOLD) THEN
	TBANDV2=BANDV2
      END IF 
      MISEOLD=MISE
 293  CONTINUE

      BANDV2=TBANDV2
      endif
	 
	 OUTV30(8:11)=(/BANDV1,BANDV2,V1,V2/)  

cc Now compute the estimates

cc Calculate F2vect (F_2(t,v)), at the group 2 failure times
       do 302 j=1,ngridv
       vv=V1+float(j-1)*(V2-V1)/float(ngridv-1)
         do 303 I=1,nsamp2
         F2vect2(I,j)=0.
         vF2vect2(I,j)=0.
	 F2vect2cum(I,j)=0.
	 vF2vect2cum(I,j)=0.
         do 304 L=1,nsamp2tr
       IF ((ABS(censor2(L)).LT.0.00000001) .or.
     #     (time2(L) .gt. time2(I))) GOTO 304

       IF (cause2(L) .le. vv) THEN
        F2vect2cum(I,j)=F2vect2cum(I,j)+
     #                  (KMT2(L)/float(nsamp2-L+1))
        vF2vect2cum(I,j)=vF2vect2cum(I,j)+
     #   ((KMT2(L)/float(nsamp2-L+1)))**2.
       END IF

       IF ((cause2(L) .ge. UV1) .and.
     #     (cause2(L) .le. UV2)) THEN
       IF ((vv .GE. (UV1+BANDV2)) .AND.
     #     (vv .LE. (UV2-BANDV2))) THEN
        F2vect2(I,j)=F2vect2(I,j)+EPAN((vv-cause2(L))/BANDV2)*
     #               (KMT2(L)/float(nsamp2-L+1))/BANDV2
        
        vF2vect2(I,j)=vF2vect2(I,j)+(EPAN((vv-cause2(L))/
     #  BANDV2)*((KMT2(L)/float(nsamp2-L+1))/BANDV2))**2.

        ELSE IF (vv .LT. (UV1+BANDV2)) THEN

        IF (TAILSV .eq. 1) THEN
cc Deal with the tail problem
        F2vect2(I,j)=F2vect2(I,j)+
     #  EPANL((vv-cause2(L))/BANDV2,(vv-UV1)/BANDV2)*
     #               (KMT2(L)/float(nsamp2-L+1))/BANDV2
        vF2vect2(I,j)=vF2vect2(I,j)+
     #  (EPANL((vv-cause2(L))/BANDV2,(vv-UV1)/BANDV2)*
     #                ((KMT2(L)/float(nsamp2-L+1))/BANDV2))**2.
        ELSE
        F2vect2(I,j)=F2vect2(I,j)+EPAN((vv-cause2(L))/BANDV2)*
     #               (KMT2(L)/float(nsamp2-L+1))/BANDV2
        
        vF2vect2(I,j)=vF2vect2(I,j)+(EPAN((vv-cause2(L))/
     #  BANDV2)*((KMT2(L)/float(nsamp2-L+1))/BANDV2))**2.
        ENDIF

        ELSE IF (vv .GT. (UV2-BANDV2)) THEN

        IF (TAILSV .eq. 1) THEN
        F2vect2(I,j)=F2vect2(I,j)+
     #  EPANU((vv-cause2(L))/BANDV2,-(UV2-vv)/BANDV2)*
     #               (KMT2(L)/float(nsamp2-L+1))/BANDV2
        vF2vect2(I,j)=vF2vect2(I,j)+
     #  (EPANU((vv-cause2(L))/BANDV2,-(UV2-vv)/
     #  BANDV2)*((KMT2(L)/float(nsamp2-L+1))/BANDV2))**2.
    
        ELSE
        F2vect2(I,j)=F2vect2(I,j)+EPAN((vv-cause2(L))/BANDV2)*
     #               (KMT2(L)/float(nsamp2-L+1))/BANDV2

        vF2vect2(I,j)=vF2vect2(I,j)+(EPAN((vv-cause2(L))/
     #  BANDV2)*((KMT2(L)/float(nsamp2-L+1))/BANDV2))**2.
        ENDIF

        ENDIF 
        ENDIF
 304  CONTINUE
 303  CONTINUE
 302  CONTINUE

cc----------------------------------------------
cc Compute the weights for the test process
cc----------------------------------------------

cc Down-weight the right-tail;
cc Set it to lambda2hat(t)[(Y_1(t)/n_1)*(Y_2(t)/n_2)]^{1/2}
cc as suggested by Yanqing:

cc Compute Y1 at the group 2 failure times:
      do 305 i=1,nsamp2
        Y12(i)=0.
      do 305 ii=1,nsamp1
        if (time1(ii) .ge. time2(i)) then
          Y12(i)=Y12(i)+1.
        end if
 305  continue

cc Suggested by Yanqing:      
      do i=1,nsamp2
        Hnm2(i)=(((Y12(i)/nsamp1)*(float(nsamp2-i+1)/nsamp2))**0.5)
cc     #          EE2(i)
      end do

cc Down-weight the right-tail:
cc Compute Y2 at the group 1 failure times

      do 306 i=1,nsamp1
        Y21(i)=0.
      do 306 ii=1,nsamp2
        if (time2(ii) .ge. time1(i)) then
          Y21(i)=Y21(i)+1.
        end if
 306  continue
                                                                                                                                    
cc Set it to lambda2hat(t)*((Y_1(t)/n_1)*(Y_2(t)/n_2))^{1/2}
cc as suggested by Yanqing:
      do i=1,nsamp1
        Hnm1(i)=(((float(nsamp1-i+1)/nsamp1)*(Y21(i)/nsamp2))**0.5)
cc     #          EE21(i)
      end do

      wt=sqrt((float(nsamp1)*float(nsamp2))/
     #        float(nsamp))
      wt1=sqrt(float(nsamp2)/(float(nsamp1)*float(nsamp)))
      wt2=sqrt(float(nsamp1)/(float(nsamp2)*float(nsamp)))

cc-----------------------------------------------------
cc Weights for the test statistics
cc-----------------------------------------------------
cc weightv is a vector of the number of failures with
cc marks between the grid values

c	Change ngridv-1 to ngridv: 

c      do 307 j=1,(ngridv-1)
       do 307 j=1,(ngridv)
        weightv(j)=0.
        vv=V1+float(j-1)*(V2-V1)/float(ngridv-1)
        vvnext=V1+float(j)*(V2-V1)/float(ngridv-1)
         do 308 I=1,nsamp1
       IF ((ABS(censor1(I)).gt.0.5) .and.
     #     (cause1(I) .ge. vv) .and.
     #     (cause1(I) .lt. vvnext)) THEN
           weightv(j)=weightv(j)+1.

       END IF
 308     continue
         do 309 I=1,nsamp2
       IF ((ABS(censor2(I)).gt.0.5) .and.
     #     (cause2(I) .ge. vv) .and.
     #     (cause2(I) .lt. vvnext)) THEN
           weightv(j)=weightv(j)+1.
       END IF
 309    continue
        weightv(j)=weightv(j)/nv
 307  continue        
      
cc-----------------------------------------------------
cc Cumulative incidence functions 
cc-----------------------------------------------------
cc For sample 1: 
cc indstm1 is such that time[indstm1] = time1
cc indsv1 is such that ocause[indsv1] = ocause1

      do 314 j=1,ngridv
         icount=1
         iicount=1
         F1vect(0,j)=0.
         F1vectcum(0,j)=0.
         vF1vect(0,j)=0.
         vF1vectcum(0,j)=0.
         do 316 i=1,nsamp
         
	 if (i .eq. indstm1(iicount)) then 
         F1vect(i,j)=F1vect1(iicount,j)
	 vF1vect(i,j)=vF1vect1(iicount,j)
	 F1vectcum(i,j)=F1vect1cum(iicount,j)
	 vF1vectcum(i,j)=vF1vect1cum(iicount,j)
         iicount=iicount+1
                  
         else
         F1vect(i,j)=F1vect(i-1,j)
         vF1vect(i,j)=vF1vect(i-1,j)
	 F1vectcum(i,j)=F1vectcum(i-1,j)
	 vF1vectcum(i,j)=vF1vectcum(i-1,j)
         end if         

         icount=icount+1         
 316  continue
 314  continue

cc For sample 2:
cc indstm2 is such that time[indstm2] = time2
cc indsv1 is such that ocause[indsv2] = ocause2

cc Calculate Lambda2pr1, at the group 1 failure times
cc simplify, ignore the tail problem:
      do 320 j=1,nv2
         do 322 I=1,nsamp1
         Lambda2pr1(I,j)=0.
         do 324 L=1,nsamp2tr
       IF (ABS(censor2(L)).LT.0.00000001) GOTO 324

       IF (cause2(L) .le. ocause2(j)) THEN

       IF (time2(L) .GE. T1 .AND. time2(L) .LE. T2) THEN
         Lambda2pr1(I,j)=Lambda2pr1(I,j)
     #   +EPAN((time1(I)-time2(L))/BAND2)*((NA2(L)-NA2(L-1))/BAND2)
c     #   +EPAN((time1(I)-time2(L))/BAND2)/(BAND2*(float(nsamp2)-L+1))
 
       ENDIF
       ENDIF
 324   CONTINUE
 322   CONTINUE
 320   CONTINUE

cc Calculate Lambda2pr2, at the group 2 failure times
cc simplify, ignore the tail problem:
      do 325 j=1,nv2
         do 330 I=1,nsamp2
         Lambda2pr2(I,j)=0.
         do 347 L=1,nsamp2tr
       IF (ABS(censor2(L)).LT.0.00000001) GOTO 347

       IF (cause2(L) .le. ocause2(j)) THEN

       IF (time2(L) .GE. T1 .AND. time2(L) .LE. T2) THEN
         Lambda2pr2(I,j)=Lambda2pr2(I,j)
     #   +EPAN((time2(I)-time2(L))/BAND2)*((NA2(L)-NA2(L-1))/BAND2)
c     #   +EPAN((time2(I)-time2(L))/BAND2)/(BAND2*(float(nsamp2)-L+1))
       ENDIF
       ENDIF
 347   CONTINUE
 330   CONTINUE
 325   CONTINUE

      do 381 j=1,ngridv
         icount=1
         iicount=1
         F2vect(0,j)=0.
         F2vectcum(0,j)=0.
         vF2vect(0,j)=0.
         vF2vectcum(0,j)=0.
         do 383 i=1,nsamp
        
	 if (i .eq. indstm2(iicount)) then 
         F2vect(i,j)=F2vect2(iicount,j)
	 
       vF2vect(i,j)=vF2vect2(iicount,j)
	 F2vectcum(i,j)=F2vect2cum(iicount,j)
	 vF2vectcum(i,j)=vF2vect2cum(iicount,j)
         iicount=iicount+1
                  
         else
         F2vect(i,j)=F2vect(i-1,j)
         vF2vect(i,j)=vF2vect(i-1,j)
	 F2vectcum(i,j)=F2vectcum(i-1,j)
	 vF2vectcum(i,j)=vF2vectcum(i-1,j)
         end if         

         icount=icount+1         
 383  continue
 381  continue


      do 392 j=1,nv1
        gtmpproc1(0,j)=0.
	gtproc10(0,j)=0.
 392  continue

      do 396 j=1,nv1
      ii=T1ind1-1
      do 396 l=1,ngrid
      tt=T1+float(l-1)*(T2-T1)/float(ngrid-1)
      temp1=0.
      i=ii+1
 397  if ((time1(i) .ge. time1(ii+1)) .and. (time1(i) .le. tt)) then     
      temp1=temp1+wt*Hnm1(i)*hatl1(i,j)
      i=i+1
      goto 397
      else
      ii=i-1
      end if
      gtmpproc1(l,j)=gtmpproc1(l-1,j)+temp1                
      gtproc10(l,j)=gtproc10(l-1,j)+temp1

 396  continue

      do 399 jj=1,ngridv
      vv=V1+float(jj-1)*(V2-V1)/float(ngridv-1)
      do 400 j=1,nv1
      if (ocause1(j) .le. vv) then
      do 401 l=1,ngrid
      gtmppr1(l,jj)=gtmpproc1(l,j)
      gtmppr10(l,jj)=gtproc10(l,j)
 401  continue
      end if
 400  continue
 399  continue


      do 404 j=1,nv2
        gtmpproc2(0,j)=0.
	gtproc2ph(0,j)=0.
	gtproc20(0,j)=0.
 404  continue

      do 405 j=1,nv2
      ii=T1ind2-1
      do 405 l=1,ngrid
      tt=T1+float(l-1)*(T2-T1)/float(ngrid-1)
      temp1=0.
      temp10=0.
      temp1ph=0.
      i=ii+1
 406  if ((time2(i) .ge. time2(ii+1)) .and. (time2(i) .le. tt)) then     
      temp1=temp1+wt*Hnm2(i)*(EE1(i)/EE2(i))*hatl2(i,j)
      temp1ph=temp1ph+wt*Hnm2(i)*exp(beta)*hatl2(i,j)
      temp10=temp10+wt*Hnm2(i)*hatl2(i,j)

      i=i+1
      goto 406
      else
      ii=i-1
      end if
      gtmpproc2(l,j)=gtmpproc2(l-1,j)+temp1                
      gtproc2ph(l,j)=gtproc2ph(l-1,j)+temp1ph
      gtproc20(l,j)=gtproc20(l-1,j)+temp10
 405  continue

      do 407 jj=1,ngridv
      vv=V1+float(jj-1)*(V2-V1)/float(ngridv-1)
      do 408 j=1,nv2
      if (ocause2(j) .le. vv) then
      do 409 l=1,ngrid
        gtmppr2(l,jj)=gtmpproc2(l,j)
	gtmppr2ph(l,jj)=gtproc2ph(l,j)
	gtmppr20(l,jj)=gtproc20(l,j)
 409  continue
      end if
 408  continue
 407  continue
       do 550 j=1,ngridv
       do 550 i=1,ngrid
         gprocess(i,j)=gtmppr1(i,j)-gtmppr2(i,j)
         obsproc(i,j)=gprocess(i,j)
	 gprocph(i,j)=gtmppr1(i,j)-gtmppr2ph(i,j)
	 gproc0(i,j)=gtmppr20(i,j)-gtmppr10(i,j)
 550   continue

       u1=0.
       u1ph=0.
       do 600 ii=1,(ngrid-1)
         do 603 ij=(ii+1),ngrid
         do 605 i=1,(ngridv-1)
             do 610 j=(i+1),ngridv
               delta11=abs(gprocess(ij,j)-gprocess(ij,i)
     #                - gprocess(ii,j)+gprocess(ii,i))
               if (delta11.ge.u1) then
                   u1=delta11
               endif          
             delta11ph=abs(gprocph(ij,j)-gprocph(ij,i)
     #                - gprocph(ii,j)+gprocph(ii,i))
               if (delta11ph.ge.u1ph) then
                   u1ph=delta11ph
               endif          
 610  continue
 605  continue
 603  continue
 600  continue

      u2=0.
      u2ph=0.
      u6=0.
      u6ph=0.
      do 620 ii=1,(ngrid-1)
      do 623 ij=(ii+1),ngrid
      do 625 i=1,ngridv
               delta12=gprocess(ii,i)-gprocess(ij,i)
               if(delta12.ge.u2) then
                   u2=delta12
               endif
               delta12ph=gprocph(ii,i)-gprocph(ij,i)
               if(delta12ph.ge.u2ph) then
                   u2ph=delta12ph
               endif
         u6=u6+weightv(i)*(gprocess(ii,i)-gprocess(ij,i))
         u6ph=u6ph+weightv(i)*(gprocph(ii,i)-gprocph(ij,i))
 625  continue
 623  continue
 620  continue


cc-----------------------------------------------------
cc Test statistics U11, U12, U13, U14
cc-----------------------------------------------------


      u3=0.
      u3ph=0.
      u7=0.
      u7ph=0.
         do 627 ij=1,ngrid
         do 629 i=1,(ngridv-1)
             do 631 j=(i+1),ngridv
               delta11=abs(gprocess(ij,j)-gprocess(ij,i))
               if (delta11.ge.u3) then
                   u3=delta11
               endif
               delta11ph=abs(gprocph(ij,j)-gprocph(ij,i))
               if (delta11ph.ge.u3ph) then
                   u3ph=delta11ph
               endif
         u7=u7+(((float(i)-float(j))/float(ngridv))**2.)*
     #         (gprocess(ij,j)-gprocess(ij,i))**2.
         u7ph=u7ph+(((float(i)-float(j))/float(ngridv))**2.)*
     #         (gprocph(ij,j)-gprocph(ij,i))**2.
 631  continue
 629  continue
 627  continue


      u10=0.
      u20=0.
      u30=0.
      u40=0.
      u50=0.
      u60=0.
      u70=0.
      u80=0.
         do 633 ij=1,ngrid
         do 635 i=1,ngridv
               delta110=gproc0(ij,i)
               if (delta110.ge.u20) then
                   u20=delta110
               endif
               delta120=abs(gproc0(ij,i))
               if (delta120.ge.u30) then
                   u30=delta120
               endif
               if (ij .eq. ngrid) then
	       u40=u40+gproc0(ij,i)*weightv(i)
	       u50=u50+(gproc0(ij,i)**2.)*weightv(i)

               end if
	       if (ij.eq.ngrid) then
	       if (delta110.ge.u70) then
		  u70=delta110
               end if
	       if (delta120.ge.u80) then
		  u80=delta120
               end if
	       end if
 635  continue
 633  continue
      u10=gproc0(ngrid,ngridv)
      u60=abs(gproc0(ngrid,ngridv))

      u4=0.
      u4ph=0.
         do 637 ij=1,ngrid
         do 639 i=1,(ngridv-1)
               u4=u4+weightv(i)*(gprocess(ij,i+1)-gprocess(ij,i))
         u4ph=u4ph+weightv(i)*(gprocph(ij,i+1)-gprocph(ij,i))
 639  continue
 637  continue
      u4=u4/float(ngrid)
      u4ph=u4ph/float(ngrid)

cc-----------------------------------------------------
cc Test statistics Unp1, Unp2, Usp1, Usp2 
cc-----------------------------------------------------

      u5=0.
      u52=0.
      u53=0.
      u54=0.
      u5ph=0.
      u52ph=0.
      u53ph=0.
      u54ph=0.
      u8=0.
      u8ph=0.
      u9=0.
      u9ph=0.
         do 641 ii=1,(ngrid-1)
         do 643 ij=(ii+1),ngrid
         do 645 i=1,(ngridv-1)
             do 647 j=(i+1),ngridv
         delta11=gprocess(ij,i)+gprocess(ij,j)-2.*gprocess(ij,(i+j)/2)
     #        -(gprocess(ii,i)+gprocess(ii,j)-2.*gprocess(ii,(i+j)/2))
               if (delta11.ge.u5) then
                   u5=delta11
               endif
               if (abs(delta11).ge.u52) then
                   u52=abs(delta11)
               endif
         u8=u8+weightv(i)*delta11
         u9=u9+(((float(i)-float(j))/float(ngridv))**2.)*
     #         (delta11**2.)
               if (ii.eq.1 .and. ij.eq.ngrid) then
         delta11=gprocess(ij,i)+gprocess(ij,j)-2.*gprocess(ij,(i+j)/2)
               if (delta11.ge.u53) then
                   u53=delta11
               end if
               if (abs(delta11).ge.u54) then
                   u54=abs(delta11)
               end if
               delta13=abs(gprocess(ij,i))
               if (delta13.ge.u11) then
                  u11=delta13
               end if
               end if
         delta12=gprocph(ij,i)+gprocph(ij,j)-2.*gprocph(ij,(i+j)/2)
     #        -(gprocph(ii,i)+gprocph(ii,j)-2.*gprocph(ii,(i+j)/2))
               if (delta12.ge.u5ph) then
                   u5ph=delta12
               endif
               if (abs(delta12).ge.u52ph) then
                   u52ph=abs(delta12)
               endif
         u8ph=u8ph+weightv(i)*delta12
         u9ph=u9ph+(((float(i)-float(j))/float(ngridv))**2.)*
     #         (delta12**2.)
               if (ii.eq.1 .and. ij.eq.ngrid) then
         delta12=gprocph(ij,i)+gprocph(ij,j)-2.*gprocph(ij,(i+j)/2)
               if (delta12.ge.u53ph) then
                   u53ph=delta12
               end if
               if (abs(delta12).ge.u54ph) then
                   u54ph=abs(delta12)
               end if
               delta13=abs(gprocph(ij,i))
               if (delta13.ge.u11ph) then
                  u11ph=delta13
               end if
               end if
 647  continue
 645  continue
 643  continue
 641  continue


      u1obs=u1
      u2obs=u2
      u3obs=u3
      u4obs=u4
      u5obs=u5
      u52obs=u52
      u53obs=u53
      u54obs=u54
      u6obs=u6
      u7obs=u7
      u8obs=u8
      u9obs=u9
      u1obs0=u10
      u2obs0=u20
      u3obs0=u30
      u4obs0=u40
      u5obs0=u50
      u6obs0=u60
      u7obs0=u70
      u8obs0=u80
      u1obsph=u1ph
      u2obsph=u2ph
      u3obsph=u3ph
      u4obsph=u4ph
      u5obsph=u5ph
      u52obsph=u52ph
      u53obsph=u53ph
      u54obsph=u54ph
      u6obsph=u6ph
      u7obsph=u7ph
      u8obsph=u8ph
      u9obsph=u9ph

cc----------------------------------------------------------------------
cc Mark-specific vaccine efficacy on a grid for the mark V:
cc----------------------------------------------------------------------
      do 730 j=1,ngridv
      do 735 i=1,nsamp
cc A convention is used to prevent dividing by 0
	if (F2vect(i,j).eq.0.) then
	VEvect(i,j)=0.
	vVEvect(i,j)=1000.
	else
        VEvect(i,j)=1.-(F1vect(i,j)/F2vect(i,j))
cc Use the delta method to compute the variance estimate:
        vVEvect(i,j)=(vF1vect(i,j)/(F2vect(i,j)**2.)) +
     #  (vF2vect(i,j)*(F1vect(i,j)**2.))/(F2vect(i,j)**4.)
	endif

	if (F2vectcum(i,j).eq.0.) then
	VEvectcum(i,j)=0.
	vVEvectcum(i,j)=1000.
	else
	VEvectcum(i,j)=1.-(F1vectcum(i,j)/F2vectcum(i,j))
        vVEvectcum(i,j)=(vF1vectcum(i,j)/(F2vectcum(i,j)**2.)) +
     #  (vF2vectcum(i,j)*(F1vectcum(i,j)**2.))/
     #  (F2vectcum(i,j)**4.)
	endif
 735  continue
	  F1vectTT(j)=F1vect(indtt,j)
        F2vectTT(j)=F2vect(indtt,j)
        vF1vcTT(j)=vF1vect(indtt,j)
        vF2vcTT(j)=vF2vect(indtt,j)

        F1vectTTcum(j)=F1vectcum(indtt,j)
        F2vectTTcum(j)=F2vectcum(indtt,j)
        vF1vcTTcum(j)=vF1vectcum(indtt,j)
        vF2vcTTcum(j)=vF2vectcum(indtt,j)

        VEvectTT(j)=VEvect(indtt,j)
        vVEvcTT(j)=vVEvect(indtt,j)

        VEvectTTcum(j)=VEvectcum(indtt,j)
        vVEvcTTcum(j)=vVEvectcum(indtt,j)

cc The confidence intervals are computed by transforming symmetric
cc limits about log(F1vectTT(j)/F2vectTT(j))
        lowlim(j)=1.0 - (1.0-VEvectTT(j))*exp(1.96*sqrt(
     #  (vF1vcTT(j)/(F1vectTT(j)**2.)) +
     #      (vF2vcTT(j)/(F2vectTT(j)**2.))))
        uplim(j)=1.0 - (1.0-VEvectTT(j))*exp(-1.96*sqrt(
     #  (vF1vcTT(j)/(F1vectTT(j)**2.)) +
     #      (vF2vcTT(j)/(F2vectTT(j)**2.))))

        lowlimcum(j)=1.0 - (1.0-VEvectTTcum(j))*exp(1.96*sqrt(
     #  (vF1vcTTcum(j)/(F1vectTTcum(j)**2.)) +
     #      (vF2vcTTcum(j)/(F2vectTTcum(j)**2.))))
        uplimcum(j)=1.0 - (1.0-VEvectTTcum(j))*exp(-1.96*sqrt(
     #  (vF1vcTTcum(j)/(F1vectTTcum(j)**2.)) +
     #      (vF2vcTTcum(j)/(F2vectTTcum(j)**2.))))

 730  continue

cc   Loop to find P-values from simulated cause of failure data
cc   Note: nv1, nv2, Hnm, EE1, EE2 do not change here.
      
      nrej1=0.
      nrej2=0.
      nrej3=0.
      nrej4=0.
      nrej5=0.
      nrej52=0.
      nrej53=0.
      nrej54=0.
      nrej6=0.
      nrej7=0.
      nrej8=0.
      nrej9=0.
      nrej1ph=0.
      nrej2ph=0.
      nrej3ph=0.
      nrej4ph=0.
      nrej5ph=0.
      nrej52ph=0.
      nrej53ph=0.
      nrej54ph=0.
      nrej6ph=0.
      nrej7ph=0. 
      nrej8ph=0.
      nrej9ph=0.
      nrej10=0.
      nrej20=0.
      nrej30=0.
      nrej40=0.
      nrej50=0.
      nrej60=0.
      nrej70=0.
      nrej80=0.

      do 1000 iboot=1,nboot

cc   For group 1:

      do 2112 j=1,nv
         do 2112 i=1,ngrid
         sprocess(i,j)=0.
 2112 continue

      do 2120 j=1,nv+1
            stmpproc1(0,j)=0.
            stmpproc2(0,j)=0.
	    stproc1ph(0,j)=0.
	    stproc2ph(0,j)=0.
	    stproc10(0,j)=0.
	    stproc20(0,j)=0.
 2120 continue

      do 2121 i=1,ngrid
      do 2123 j=1,nv+1
            stmpproc1(i,j)=0.
            stmpproc2(i,j)=0.
	    stproc1ph(i,j)=0.
	    stproc2ph(i,j)=0.
	    stproc10(i,j)=0.
	    stproc20(i,j)=0.
 2123 continue
 2121 continue

      cumW1(nsamp1)=0.
      cumW1U1b(nsamp1)=0.
      cumW1(nsamp1+1)=0.
      cumW1U1b(nsamp1+1)=0.
      do 2140 i=1,nsamp1
        ii=nsamp1-i+1
        W1(ii)=rnor()
        cumW1(ii)=cumW1(ii+1)+W1(ii)
	cumW1U1b(ii)=cumW1U1b(ii+1)+W1(ii)*U1beta(ii)
 2140  continue


      cumW2(nsamp2)=0.
      cumW2U2b(nsamp2)=0.
      cumW2(nsamp2+1)=0.
      cumW2U2b(nsamp2+1)=0.
      do 2142 i=1,nsamp2
        ii=nsamp2-i+1
        W2(ii)=rnor()
        cumW2(ii)=cumW2(ii+1)+W2(ii) 
	cumW2U2b(ii)=cumW2U2b(ii+1)+W2(ii)*U2beta(ii)
 2142  continue

      do 2143 j=1,nv1
        h1hat1(0,j)=0.
        h1hat2(0,j)=0.
        h1hat3(0,j)=0.
        h1hat4(0,j)=0.
 2143  continue

cc Start i as the index of the smallest time1 >= T1
      do 2232 j=1,nv1
      ii=T1ind1-1
      do 2233 l=1,ngrid
      tt=T1+float(l-1)*(T2-T1)/float(ngrid-1)
      temp1=0.
      temp3=0.
      i=ii+1
 2234 if ((time1(i) .ge. time1(ii+1)) .and. (time1(i) .le. tt)) then
      ttemp1=0.
      ttemp3=0.
      if (censor1(i) .gt. 0.1) then
      if (cause1(i) .le. ocause1(j)) then
        ttemp1=wt1*Hnm1(i)*float(nsamp1)/float(nsamp1-i+1)
      end if
         ahat(i)=1/EE21(i)
         ttemp3=wt1*(Hnm1(i)*(float(nsamp1)/float(nsamp1-i+1)))*
     #   ahat(i)*Lambda2pr1(i,j)
      end if
      temp1=temp1 + ttemp1*W1(i)
      temp3=temp3 + ttemp3*W1(i)
      i=i+1
      goto 2234
      else
      ii=i-1
      end if
         h1hat1(l,j)=h1hat1(l-1,j)+temp1
         h1hat3(l,j)=h1hat3(l-1,j)+temp3
         if (j .eq. 10) then
         end if
 2233 continue
 2232 continue

cc Start k as the index of the smallest time1 >= T1
      do 2252 j=1,nv1
      ii=T1ind1-1
      do 2252 l=1,ngrid
      tt=T1+float(l-1)*(T2-T1)/float(ngrid-1)

      temp2=0.
      temp4=0.
      k=ii+1
 2255 if ((time1(k) .ge. time1(ii+1)) .and. (time1(k) .le. tt)) then
       ttemp2=0.
       ttemp4=0.
         if ((censor1(k) .gt. 0.1)) then
	 ahat(k)=1/EE21(k)
         if (cause1(k) .le. ocause1(j)) then 
           ttemp2=wt1*(Hnm1(k)*(float(nsamp1))/ 
     #     (float(nsamp1-k+1)**2.))*cumW1(k)

         end if
         ttemp4=wt1*(Hnm1(k)*(float(nsamp1))/
     #              (float(nsamp1-k+1)**2.))*
     #      ahat(k)*Lambda2pr1(k,j)*cumW1(k)
	 end if
      temp2=temp2 + ttemp2
      temp4=temp4 + ttemp4
      k=k+1
      goto 2255
      else
      ii=k-1
      end if
         h1hat2(l,j)=h1hat2(l-1,j)+temp2
         h1hat4(l,j)=h1hat4(l-1,j)+temp4
         if (j .eq. 10) then
         end if
 2252 continue

      do 2300 j=1,nv2
        h2hat1(0,j)=0.
	  h2hat10(0,j)=0.
        h2hat2(0,j)=0.
	  h2hat20(0,j)=0.
        h2hat3(0,j)=0.
        h2hat4(0,j)=0.
        h2hat1ph(0,j)=0.
        h2hat2ph(0,j)=0.
 2300  continue




cc For Group 2

      do 2335 j=1,nv2
      ii=T1ind2-1
      do 2335 l=1,ngrid
      tt=T1+float(l-1)*(T2-T1)/float(ngrid-1)
      temp1=0.
      temp10=0.
      temp1ph=0.
      temp3=0.
      i=ii+1
 2336 if ((time2(i) .ge. time2(ii+1)) .and. (time2(i) .le. tt)) then
      ttemp1=0.
      ttemp10=0.
      ttemp1ph=0.
      ttemp3=0.
      if (censor2(i) .gt. 0.1) then
      if (cause2(i) .le. ocause2(j)) then
	ttemp1=wt2*Hnm2(i)*(float(nsamp2)/float(nsamp2-i+1)) 
     #       *(EE1(i)/EE2(i))
	ttemp10=wt2*Hnm2(i)*(float(nsamp2)/float(nsamp2-i+1)) 
      ttemp1ph=wt2*Hnm2(i)*(float(nsamp2)/float(nsamp2-i+1))
     #       *(exp(beta))
      end if
         bhat(i)=EE1(i)/(EE2(i)**2.)      
         ttemp3=wt2*(Hnm2(i)*(float(nsamp2)/
     #   float(nsamp2-i+1)))*bhat(i)*Lambda2pr2(i,j)
      end if
      temp1=temp1 + ttemp1*W2(i)
      temp10=temp10 + ttemp10*W2(i)
      temp1ph=temp1ph + ttemp1ph*W2(i)
      temp3=temp3 + ttemp3*W2(i)
      i=i+1
      goto 2336
      else
      ii=i-1
      end if
         h2hat1(l,j)=h2hat1(l-1,j)+temp1
	   h2hat10(l,j)=h2hat10(l-1,j)+temp10
         h2hat1ph(l,j)=h2hat1ph(l-1,j)+temp1ph
         h2hat3(l,j)=h2hat3(l-1,j)+temp3
 2335 continue

cc Start k as the index of the smallest time2 >= T1

      do 2352 j=1,nv2
      ii=T1ind2-1
      do 2352 l=1,ngrid
      tt=T1+float(l-1)*(T2-T1)/float(ngrid-1)

      temp2=0.
      temp20=0.
      temp2ph=0.
      temp4=0.
      k=ii+1
 2355 if ((time2(k) .ge. time2(ii+1)) .and. (time2(k) .le. tt)) then
       ttemp2=0.
       ttemp20=0.
       ttemp2ph=0.
       ttemp4=0.
         if ((censor2(k) .gt. 0.1)) then
	 bhat(k)=EE1(k)/(EE2(k)**2.)
         if (cause2(k) .le. ocause2(j)) then
           ttemp2=wt2*((EE1(k)/EE2(k))*(Hnm2(k)*
     #     (float(nsamp2))/(float(nsamp2-k+1)**2.)))*
     #     cumW2(k)
           ttemp20=wt2*Hnm2(k)*
     #     (float(nsamp2)/(float(nsamp2-k+1)**2.))*
     #     cumW2(k)
           ttemp2ph=wt2*((exp(beta))*(Hnm2(k)*
     #     (float(nsamp2))/(float(nsamp2-k+1)**2.)))*
     #     cumW2(k)
         end if
           ttemp4=wt2*(Hnm2(k)*((float(nsamp2))/
     # (float(nsamp2-k+1)**2.)))*bhat(k)*Lambda2pr2(k,j)*
     #     cumW2(k)
	 end if
      temp2=temp2 + ttemp2
      temp20=temp20 + ttemp20
      temp2ph=temp2ph + ttemp2ph
      temp4=temp4 + ttemp4
      k=k+1
      goto 2355
      else
      ii=k-1
      end if
         h2hat2(l,j)=h2hat2(l-1,j)+temp2
	 h2hat20(l,j)=h2hat20(l-1,j)+temp20
         h2hat2ph(l,j)=h2hat2ph(l-1,j)+temp2ph
         h2hat4(l,j)=h2hat4(l-1,j)+temp4
 2352 continue


cc New, simpler way, on the grid of mark values V:
      
      do 2400 i=1,ngrid
      do 2402 j=1,nv1
         stmpproc1(i,j)=h1hat1(i,j)
     #   -h1hat2(i,j)-h1hat3(i,j)+h1hat4(i,j)
	 stproc1ph(i,j)=h1hat1(i,j)-h1hat2(i,j)
	 stproc10(i,j)=h1hat2(i,j)-h1hat1(i,j)
 2402 continue
 2400 continue

      do 2409 jj=1,ngridv
      vv=V1+float(jj-1)*(V2-V1)/float(ngridv-1)
      do 2410 j=1,nv1
      if (ocause1(j) .le. vv) then
      do 2411 l=1,ngrid
      stmppr1(l,jj)=stmpproc1(l,j)
      stpr1ph(l,jj)=stproc1ph(l,j)
      stpr10(l,jj)=stproc10(l,j)
 2411 continue
      end if
 2410 continue
 2409 continue

cc Add the 3rd piece of the PH test process

      do 2413 l=1,ngrid
      do 2415 j=1,ngridv
	 stpr1ph(l,j)=stpr1ph(l,j)-(cumW1U1b(1)*
     #   (exp(beta)/(Jnbeta))*gtmppr20(l,j))
 2415 continue
 2413 continue

      do 2420 i=1,ngrid
      do 2422 j=1,nv2
         stmpproc2(i,j)=h2hat1(i,j)
     #   -h2hat2(i,j)-h2hat3(i,j)+h2hat4(i,j)
	 stproc2ph(i,j)=h2hat1ph(i,j)-h2hat2ph(i,j)
	 stproc20(i,j)=h2hat20(i,j)-h2hat10(i,j)
 2422 continue
 2420 continue
	
      do 2429 jj=1,ngridv
      vv=V1+float(jj-1)*(V2-V1)/float(ngridv-1)
      do 2430 j=1,nv2
      if (ocause2(j) .le. vv) then
      do 2431 l=1,ngrid
      stmppr2(l,jj)=stmpproc2(l,j)
      stpr2ph(l,jj)=stproc2ph(l,j)
      stpr20(l,jj)=stproc20(l,j)
 2431 continue
      end if
 2430 continue
 2429 continue

cc Add the 3rd piece of the PH test process

      do 2433 l=1,ngrid
      do 2435 j=1,ngridv
	 stpr2ph(l,j)=stpr2ph(l,j)-(cumW2U2b(1)*
     #   (exp(beta)/(Jnbeta))*gtmppr20(l,j))
 2435 continue
 2433 continue

       do 2550 j=1,ngridv
       do 2550 i=1,ngrid
         sprocess(i,j)=stmppr1(i,j)-stmppr2(i,j)
	 sprocph(i,j)=stpr1ph(i,j)-stpr2ph(i,j)
	 sproc0(i,j)=stpr20(i,j)-stpr10(i,j)
 2550  continue

cc
cc   compute u1-u7, u1ph-u7ph, u20-u50
cc
       u1=0.
       u1ph=0.
       do 2600 ii=1,(ngrid-1)
         do 2603 ij=(ii+1),ngrid
         do 2605 i=1,(ngridv-1)
             do 2610 j=(i+1),ngridv
               delta11=abs(sprocess(ij,j)-sprocess(ij,i)
     #                - sprocess(ii,j)+sprocess(ii,i))
               if (delta11.ge.u1) then
                   u1=delta11
               endif          
               delta11=abs(sprocph(ij,j)-sprocph(ij,i)
     #                - sprocph(ii,j)+sprocph(ii,i))
               if (delta11.ge.u1ph) then
                   u1ph=delta11
               endif          
 2610 continue
 2605 continue
 2603 continue
 2600 continue

      u2=0.
      u2ph=0.
      u6=0.
      u6ph=0.
      do 2620 ii=1,(ngrid-1)
      do 2625 ij=(ii+1),ngrid
      do 2626 i=1,ngridv
               delta12=sprocess(ii,i)-sprocess(ij,i)
               if(delta12.ge.u2) then
                   u2=delta12
               endif
               delta12ph=sprocph(ii,i)-sprocph(ij,i)
               if(delta12ph.ge.u2ph) then
                   u2ph=delta12ph
               endif
         u6=u6+weightv(i)*(sprocess(ii,i)-sprocess(ij,i))
         u6ph=u6ph+weightv(i)*(sprocph(ii,i)-sprocph(ij,i))
 2626 continue
 2625 continue
 2620 continue

       u3=0.
       u3ph=0.
       u7=0.
       u7ph=0.
         do 2627 ij=1,ngrid
         do 2629 i=1,(ngridv-1)
             do 2631 j=(i+1),ngridv
               delta11=abs(sprocess(ij,j)-sprocess(ij,i))
               if (delta11.ge.u3) then
                   u3=delta11
               endif
	       delta11ph=abs(sprocph(ij,j)-sprocph(ij,i))
               if (delta11ph.ge.u3ph) then
                   u3ph=delta11ph
               endif
         u7=u7+(((float(i)-float(j))/float(ngridv))**2.)*
     #         (sprocess(ij,j)-sprocess(ij,i))**2.
         u7ph=u7ph+(((float(i)-float(j))/float(ngridv))**2.)*
     #         (sprocph(ij,j)-sprocph(ij,i))**2.
 2631 continue
 2629 continue
 2627 continue 

      u10=0.
      u20=0.
      u30=0.
      u40=0.
      u50=0.
      u60=0.
      u70=0.
      u80=0.

         do 2633 ij=1,ngrid
         do 2635 i=1,ngridv
               delta110=sproc0(ij,i)
               if (delta110.ge.u20) then
                   u20=delta110
               endif
               delta120=abs(sproc0(ij,i))
               if (delta120.ge.u30) then
                   u30=delta120
               endif
               if (ij .eq. ngrid) then
	       u40=u40+sproc0(ij,i)*weightv(i)
	       u50=u50+(sproc0(ij,i)**2.)*weightv(i)

               end if
	       if (ij.eq.ngrid) then
	       if (delta110.ge.u70) then
		  u70=delta110
               end if
	       if (delta120.ge.u80) then
		  u80=delta120
               end if
	       end if
 2635   continue
 2633   continue
	u10=sproc0(ngrid,ngridv)
	u60=abs(sproc0(ngrid,ngridv))
       

       u4=0.
       u4ph=0.
         do 2637 ij=1,ngrid
         do 2639 i=1,(ngridv-1)
               u4=u4+weightv(i)*(sprocess(ij,i+1)-sprocess(ij,i))
         u4ph=u4ph+weightv(i)*(sprocph(ij,i+1)-sprocph(ij,i))
 2639 continue
 2637 continue
      u4=u4/float(ngrid)
      u4ph=u4ph/float(ngrid)

      u5=0.
      u52=0.
      u53=0.
      u54=0.
      u5ph=0.
      u52ph=0.
      u53ph=0.
      u54ph=0.
      u8=0.
      u8ph=0.
      u9=0.
      u9ph=0.
         do 2641 ii=1,(ngrid-1)
         do 2643 ij=(ii+1),ngrid
         do 2645 i=1,(ngridv-1)
             do 2647 j=(i+1),ngridv
         delta11=sprocess(ij,i)+sprocess(ij,j)-2.*sprocess(ij,(i+j)/2)
     #        -(sprocess(ii,i)+sprocess(ii,j)-2.*sprocess(ii,(i+j)/2))
               if (delta11.ge.u5) then
                   u5=delta11
               endif
               if (abs(delta11).ge.u52) then
                   u52=abs(delta11)
               endif
         u8=u8+weightv(i)*delta11
         u9=u9+(((float(i)-float(j))/float(ngridv))**2.)*
     #         (delta11**2.)
               if (ii.eq.1 .and. ij.eq.ngrid) then
         delta11=sprocess(ij,i)+sprocess(ij,j)-2.*sprocess(ij,(i+j)/2)
               if (delta11.ge.u53) then
                   u53=delta11
               end if
               if (abs(delta11).ge.u54) then
                   u54=abs(delta11)
               end if
               delta13=abs(sprocess(ij,i))
               if (delta13.ge.u11) then
                  u11=delta13
               end if
               end if
         delta12=sprocph(ij,i)+sprocph(ij,j)-2.*sprocph(ij,(i+j)/2)
     #        -(sprocph(ii,i)+sprocph(ii,j)-2.*sprocph(ii,(i+j)/2))
               if (delta12.ge.u5ph) then
                   u5ph=delta12
               endif
               if (abs(delta12).ge.u52ph) then
                   u52ph=abs(delta12)
               endif
         u8ph=u8ph+weightv(i)*delta12
         u9ph=u9ph+(((float(i)-float(j))/float(ngridv))**2.)*
     #         (delta12**2.)
               if (ii.eq.1 .and. ij.eq.ngrid) then
         delta12=sprocph(ij,i)+sprocph(ij,j)-2.*sprocph(ij,(i+j)/2)
               if (delta12.ge.u53ph) then
                   u53ph=delta12
               end if
               if (abs(delta12).ge.u54ph) then
                   u54ph=abs(delta12)
               end if
               delta13=abs(sprocph(ij,i))
               if (delta13.ge.u11ph) then
                  u11ph=delta13
               end if
               end if
 2647  continue
 2645  continue
 2643  continue
 2641  continue

      if (u1.gt.u1obs) then
         nrej1=nrej1+1.
      endif    
      if (u2.gt.u2obs) then
         nrej2=nrej2+1.
      endif  
      if (u3.gt.u3obs) then
         nrej3=nrej3+1.
      endif
      if (u4.gt.u4obs) then
         nrej4=nrej4+1.
      endif
      if (u5.gt.u5obs) then
         nrej5=nrej5+1.
      endif
      if (u52.gt.u52obs) then
         nrej52=nrej52+1.
      endif
      if (u53.gt.u53obs) then
         nrej53=nrej53+1.
      endif
      if (u54.gt.u54obs) then
         nrej54=nrej54+1.
      endif
      if (u6.gt.u6obs) then
         nrej6=nrej6+1.
      endif
      if (u7.gt.u7obs) then
         nrej7=nrej7+1.
      endif
      if (u8.gt.u8obs) then
         nrej8=nrej8+1.
      endif
      if (u9.gt.u9obs) then
         nrej9=nrej9+1.
      endif

      if (u1ph.gt.u1obsph) then
         nrej1ph=nrej1ph+1.
      endif  
      if (u2ph.gt.u2obsph) then
         nrej2ph=nrej2ph+1.
      endif  
      if (u3ph.gt.u3obsph) then
         nrej3ph=nrej3ph+1.
      endif
      if (u4ph.gt.u4obsph) then
         nrej4ph=nrej4ph+1.
      endif  
      if (u5ph.gt.u5obsph) then
         nrej5ph=nrej5ph+1.
      endif
      if (u52ph.gt.u52obsph) then
         nrej52ph=nrej52ph+1.
      endif
      if (u53ph.gt.u53obsph) then
         nrej53ph=nrej53ph+1.
      endif
      if (u54ph.gt.u54obsph) then
         nrej54ph=nrej54ph+1.
      endif
      if (u6ph.gt.u6obsph) then
         nrej6ph=nrej6ph+1.
      endif
      if (u7ph.gt.u7obsph) then
         nrej7ph=nrej7ph+1.
      endif
      if (u8ph.gt.u8obsph) then
         nrej8ph=nrej8ph+1.
      endif
      if (u9ph.gt.u9obsph) then
         nrej9ph=nrej9ph+1.
      endif


      if (u10.gt.u1obs0) then
         nrej10=nrej10+1.
      endif  
      if (u20.gt.u2obs0) then
         nrej20=nrej20+1.
      endif  
      if (u30.gt.u3obs0) then
         nrej30=nrej30+1.
      endif
      if (u40.gt.u4obs0) then
         nrej40=nrej40+1.
      endif  
      if (u50.gt.u5obs0) then
         nrej50=nrej50+1.
      endif
      if (u60.gt.u6obs0) then
         nrej60=nrej60+1.
      endif
      if (u70.gt.u7obs0) then
         nrej70=nrej70+1.
      endif
      if (u80.gt.u8obs0) then
         nrej80=nrej80+1.
      endif
 1000  continue

      pvalue1=nrej1/float(nboot)
      pvalue2=nrej2/float(nboot)
      pvalue3=nrej3/float(nboot)
      pvalue4=nrej4/float(nboot)
      pvalue5=nrej5/float(nboot)
      pvalue52=nrej52/float(nboot)
      pvalue53=nrej53/float(nboot)
      pvalue54=nrej54/float(nboot)
      pvalue6=nrej6/float(nboot)
      pvalue7=nrej7/float(nboot) 
      pvalue8=nrej8/float(nboot)
      pvalue9=nrej9/float(nboot)
      pvalue1ph=nrej1ph/float(nboot)
      pvalue2ph=nrej2ph/float(nboot)
      pvalue3ph=nrej3ph/float(nboot)
      pvalue4ph=nrej4ph/float(nboot)
      pvalue5ph=nrej5ph/float(nboot)
      pvalue52ph=nrej52ph/float(nboot)
      pvalue53ph=nrej53ph/float(nboot)
      pvalue54ph=nrej54ph/float(nboot)
      pvalue6ph=nrej6ph/float(nboot)
      pvalue7ph=nrej7ph/float(nboot)
      pvalue8ph=nrej8ph/float(nboot)
      pvalue9ph=nrej9ph/float(nboot)
      pvalue10=nrej10/float(nboot)
      pvalue20=nrej20/float(nboot)
      pvalue30=nrej30/float(nboot)
      pvalue40=nrej40/float(nboot)
      pvalue50=nrej50/float(nboot)
      pvalue60=nrej60/float(nboot)
      pvalue70=nrej70/float(nboot)
      pvalue80=nrej80/float(nboot)

	 OUTV30(12:14)=(/VEoverall,VEPH,ZVE/)   

cc---------------------------------------------------
cc Write the test statistics described in the paper:
cc---------------------------------------------------

	 OUTV30(15:30)=(/u1obs0,u4obs0,u6obs0,u5obs0, 
     + pvalue10,pvalue40,pvalue60,pvalue50,u5obs,u52obs, 
     + pvalue5,pvalue52,u5obsph,u52obsph,pvalue5ph,pvalue52ph/)

cc------------------------------------------------------
cc Write out data at all event times between T1 and T2 
cc------------------------------------------------------
      mycount=0 
	do 8600 i=1,nsamp
      if (time(i) .ge. T1 .and. time(i) .le. T2 .and. 
     #    censor(i) .gt. 0.1) then
	mycount=mycount+1
        if (mycount .eq. 1.) then
	  mystart=1
	  else
	  mystart= (mycount-1)*ngridv + 1
	  end if
	do 8610 j=0,(ngridv-1)
      vv=V1+float(j)*(V2-V1)/float(ngridv-1)
	OUTtm1(mystart+j)=i
	OUTtm2(mystart+j)=time(i)
	OUTtm3(mystart+j)=vv
	OUTtm4(mystart+j)=F1vect(i,(j+1))
	OUTtm5(mystart+j)=F2vect(i,(j+1))
	OUTtm6(mystart+j)=vF1vect(i,(j+1))
	OUTtm7(mystart+j)=vF2vect(i,(j+1))
	OUTtm8(mystart+j)=VEvect(i,(j+1))
	OUTtm9(mystart+j)=vVEvect(i,(j+1))

 8610 continue    
      endif
 8600 continue

cc------------------------------------------------------
cc Write out data for estimating VE^c(ttanal,v) vs v: 
cc------------------------------------------------------

	do 8700 j=1,ngridv
      vv=V1+float(j-1)*(V2-V1)/float(ngridv-1)
	OUTM1(j)=indtt
	OUTM2(j)=time(indtt)
	OUTM_VEc1(j)=vv
	OUTM_VEc2(j)=F1vectTT(j)
	OUTM_VEc3(j)=F2vectTT(j)
	OUTM_VEc4(j)=vF1vcTT(j)
	OUTM_VEc5(j)=vF2vcTT(j)
	OUTM_VEc6(j)=VEvectTT(j)
	OUTM_VEc7(j)=vVEvcTT(j)
	OUTM_VEc8(j)=lowlim(j)
	OUTM_VEc9(j)=uplim(j)

 8700 continue

cc------------------------------------------------------
cc Write out data for estimating VE^dc(ttanal,v) vs v: 
cc------------------------------------------------------
      do 8720 j=1,ngridv
      vv=V1+float(j-1)*(V2-V1)/float(ngridv-1)
	OUTM_VEdc1(j)=vv
	OUTM_VEdc2(j)=F1vectTTcum(j)
	OUTM_VEdc3(j)=F2vectTTcum(j)
	OUTM_VEdc4(j)=vF1vcTTcum(j)
	OUTM_VEdc5(j)=vF2vcTTcum(j)
	OUTM_VEdc6(j)=VEvectTTcum(j)
	OUTM_VEdc7(j)=vVEvcTTcum(j)
	OUTM_VEdc8(j)=lowlimcum(j)
	OUTM_VEdc9(j)=uplimcum(j)

 8720 continue

       end
     

cc-----------------------------
cc Functions and subroutines
cc-----------------------------
      FUNCTION ran1(idum,flagin)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
	INTEGER flagin
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
c      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/

      if (flagin.eq.1) then 
	iv = 0  
      iy = 0
      endif

      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END


      SUBROUTINE RSTART(i,j,k,l)
cc This subroutine initializes the table for the f(97,33,-mod 1.)
cc generator and the values for the arithmetic sequence.
cc
cc For our convenience on the Sun, this routine seeds the IVNI generator 
cc also by setting a sign bit.
cc
cc For speed, the tables are maintained as integers and the result is 
cc floated.
cc
cc Note that a BLOCK DATA initializes the table, so it is not
cc imperative to call this subroutine. The seeds used for
cc generating the default table is 12,34,56,78.
      external unidata
      integer u(97),c,vnisign,ip,jp
      common /unidat/ u,c,vnisign,ip,jp
      c=362436
      do 2 ii=1,97
         is=0
         do 3 jj=1,24
           m=mod(i*mod(j*k,179),179)
           i=j
           j=k
           k=m
           l=mod(53*l+1,169)
           is = is+is
    3      if(mod(l*m,64) .ge. 32) is=is+1
    2    u(ii)=is
       vnisign = and(l*m,1)
       ip = 97
       jp = 33
       return
       end


      function rnor()
cc Returns a standard Normal variate using the Ziggurat method.
cc G. Marsaglia & Wai Tsang, Siam J. Sci. Stat. Comput. Vol 5
cc June 1984, pp 349-359.
cc
cc April 5, 1990. Made a correction to the tail part so that 
cc Log of a zero UNI is avoided. 
      real v(0:256)
      DATA AA,B,C/25.74023263217, .3194187376868, 25.96701302767/
      DATA RMAX/5.960464478e-8/
      DATA C1,C2,PC,XN/.9674937559244, 1.035680693510,
     +                  .00489575835,      3.289869847629/
C
      data (v(j), j=0, 29)/
     + .24698083002457, .30686964662795, .35153978938895,
     + .38807372507154, .41944111113146, .44719524398379,
     + .47226002099856, .49523381630958, .51652849544015,
     + .53644082893859, .55519243431656, .57295362330211,
     + .58985839761905, .60601428557041, .62150902494419,
     + .63641523999656, .65079379879381, .66469627686500,
     + .67816680018702, .69124344747922, .70395933341549,
     + .71634345674087, .72842137244342, .74021573037902,
     + .75174671122491, .76303238257455, .77408899225176,
     + .78493121178444, .79557233995228, .80602447408300/
C
      data (v(j), j=30, 59)/
     + .81629865509377, .82640499100558, .83635276268797,
     + .84615051484361, .85580613466031, .86532692010181,
     + .87471963944859, .88399058341281, .89314561092262,
     + .90219018948527, .91112943088938, .91996812288381,
     + .92871075737076, .93736155556822, .94592449052850,
     + .95440330734290, .96280154131524, .97112253434701,
     + .97936944974390, .98754528562459, .99565288708942,
     +1.00369495728548,1.01167406748786,1.01959266630145,
     +1.02745308807506,1.03525756060863,1.04300821222465,
     +1.05070707826663,1.05835610708045,1.06595716552808/
C
      data (v(j), j=60, 89)/
     +1.07351204407750,1.08102246150823,1.08849006926751,
     +1.09591645550830,1.10330314883748,1.11065162179932,
     +1.11796329411693,1.12523953571236,1.13248166952351,
     +1.13969097413490,1.14686868623717,1.15401600292913,
     +1.16113408387481,1.16822405332690,1.17528700202672,
     +1.18232398899035,1.18933604318933,1.19632416513391,
     +1.20328932836597,1.21023248086826,1.21715454639590,
     +1.22405642573584,1.23093899789935,1.23780312125216,
     +1.24464963458674,1.25147935814061,1.25829309456444,
     +1.26509162984343,1.27187573417500,1.27864616280591/
C
      data (v(j), j=90, 119)/
     +1.28540365683140,1.29214894395899,1.29888273923931,
     +1.30560574576607,1.31231865534737,1.31902214915020,
     +1.32571689831988,1.33240356457626,1.33908280078811,
     +1.34575525152740,1.35242155360456,1.35908233658634,
     +1.36573822329735,1.37238983030648,1.37903776839932,
     +1.38568264303767,1.39232505480705,1.39896559985334,
     +1.40560487030916,1.41224345471127,1.41888193840939,
     +1.42552090396760,1.43216093155888,1.43880259935355,
     +1.44544648390252,1.45209316051569,1.45874320363660,
     +1.46539718721367,1.47205568506879,1.47871927126394/
C
      data (v(j), j=120, 149)/
     +1.48538852046644,1.49206400831339,1.49874631177596,
     +1.50543600952421,1.51213368229290,1.51883991324899,
     +1.52555528836149,1.53228039677414,1.53901583118162,
     +1.54576218820999,1.55252006880182,1.55929007860681,
     +1.56607282837854,1.57286893437799,1.57967901878452,
     +1.58650371011513,1.59334364365258,1.60019946188339,
     +1.60707181494618,1.61396136109153,1.62086876715397,
     +1.62779470903712,1.63473987221291,1.64170495223585,
     +1.64869065527338,1.65569769865342,1.66272681143024,
     +1.66977873496983,1.67685422355612,1.68395404501935/
C
      data (v(j), j=150, 179)/
     +1.69107898138789,1.69822982956533,1.70540740203402,
     +1.71261252758712,1.71984605209074,1.72710883927814,
     +1.73440177157803,1.74172575097914,1.74908169993340,
     +1.75647056230006,1.76389330433365,1.77135091571826,
     +1.77884441065150,1.78637482898114,1.79394323739797,
     +1.80155073068867,1.80919843305257,1.81688749948677,
     +1.82461911724403,1.83239450736871,1.84021492631583,
     +1.84808166765936,1.85599606389581,1.86395948835000,
     +1.87197335719041,1.88003913156195,1.88815831984492,
     +1.89633248004955,1.90456322235616,1.91285221181234/
C
      data (v(j), j=180, 209)/
     +1.92120117119893,1.92961188407826,1.93808619803888,
     +1.94662602815264,1.95523336066136,1.96391025691183,
     +1.97265885756012,1.98148138706769,1.99038015851462,
     +1.99935757875730,2.00841615396118,2.01755849554212,
     +2.02678732655353,2.03610548856055,2.04551594904689,
     +2.05502180940521,2.06462631356762,2.07433285733931,
     +2.08414499850583,2.09406646779286,2.10410118076694,
     +2.11425325077645,2.12452700304473,2.13492699004145,
     +2.14545800827495,2.15612511666698,2.16693365669349,
     +2.17788927450032,2.18899794523217,2.20026599984763/
C
      data (v(j), j=210, 239)/
     +2.21170015473319,2.22330754447619,2.23509575821221,
     +2.24707288002746,2.25924753397443,2.27162893435095,
     +2.28422694200254,2.29705212753971,2.31011584252007,
     +2.32343029983670,2.33700866478640,2.35086515857565,
     +2.36501517636968,2.37947542241971,2.39426406533604,
     +2.40940091723926,2.42490764135740,2.44080799369296,
     +2.45712810572956,2.47389681687590,2.49114606758128,
     +2.50891136697751,2.52723235275195,2.54615346608274,
     +2.56572477136742,2.58600295987340,2.60705258940068,
     +2.62894763017268,2.65177341290207,2.67562911211061/
C
      data (v(j), j=240, 256)/
     +2.70063095234582,2.72691640673304,2.75464978267448,
     +2.78402978650249,2.81529997720240,2.84876355022267,
     +2.88480481063120,2.92392135126341,2.96677409038864,
     +3.01426863112874,3.06769501549787,3.12898502985483,
     +3.20123093035120,3.28986984762935,3.28986984762935,
     +3.28986984762935,3.28986984762935/
C     
      i=ivni()
      j=and(iabs(i),255)
      rnor=i*rmax*v(j+1)
      if (abs(rnor).le.v(j)) return
      X = (ABS(rnor)-V(J))/(V(J+1)-V(J))
      Y=UNI()
      S=X+Y
      IF (S .GT. C2) GO TO 11
      IF (S .LE. C1) RETURN
      IF (Y .GT. C-AA*EXP(-.5*(B-B*X)**2)) GO TO 11
      IF (EXP(-.5*V(J+1)**2)+Y*PC/V(J+1) .LE. EXP(-.5*rnor**2))
     +RETURN
C      ----------------TAIL PART----------------------

22      cconst = 0.0
2       x = uni()
        if (x .ge. 0.5) then
          x = .3039633925703*(-ALOG(x)+cconst)
        else
          cconst = cconst + 0.693147180559945
          goto 2
        endif

        cconst = 0.0
3       y = uni()
        if (y .ge. 0.5) then
          y = 2.*(-ALOG(y)+cconst)
        else
          cconst = cconst + 0.693147180559945
          goto 3
        endif
        IF (y .LE. x**2) GO TO 22

33      rnor = SIGN(XN-X,rnor)
        RETURN
11    rnor = SIGN(B-B*X,rnor)
      RETURN
      END


       function uni()
cc The uni function sub-program combines, with subtraction mod 1,
cc an f(97,33,-mod 1) generator with the element c in the arithmetic
cc sequence generated by c=c-cd mod(16777213./16777216.), period 2**24-3.
cc period of combined generator is (2**97-1)(2**24-3)2**23, about 2**144.
       external unidata
       integer u(97),c,vnisign,i,j
       common /unidat/ u,c,vnisign,i,j
       iu=and(u(i)-u(j),16777215)
       u(i)=iu
       i=i-1
       if(i.eq.0) i=97
       j=j-1
       if(j.eq.0) j=97
       c=c-7654321
       if(c.lt.0) c=and(c+16777213,16777215)
       uni=and(iu-c,16777215)/16777216.0
       vnisign = and(iu,32)
       return
       end


       function iuni()
cc The iuni function sub-program combines, with subtraction mod 2**24,
cc an f(97,33,-mod 2**24) generator with the element c in the arithmetic
cc sequence generated by c=c-cd mod(16777213), period 2**24-3.
cc period of combined generator is (2**97-1)(2**24-3)2**23, about 2**144.
       external unidata	
       integer u(97),c,vnisign,i,j
       common /unidat/ u,c,vnisign,i,j
       iuni=and(u(i)-u(j),16777215)
       u(i)=iuni
       vnisign = and(iuni,32)
       i=i-1
       if(i.eq.0) i=97
       j=j-1
       if(j.eq.0) j=97
       c=c-7654321
       if(c.lt.0) c=and(c+16777213,16777215)
       iuni=and(iuni-c,16777215)
       return
       end


       function ivni()
cc The ivni function sub-program combines, with subtraction mod 2**24,
cc an f(97,33,-mod 2**24) generator with the element c in the arithmetic
cc sequence generated by c=c-cd mod(16777213), period 2**24-3.
cc period of combined generator is (2**97-1)(2**24-3)2**23, about 2**144.
       external unidata
       integer u(97),c,vnisign,i,j
       common /unidat/ u,c,vnisign,i,j
       ivni=and(u(i)-u(j),16777215)
       u(i)=ivni
       vnisign = and(ivni,32)
       i=i-1
       if(i.eq.0) i=97
       j=j-1
       if(j.eq.0) j=97
       c=c-7654321
       if(c.lt.0) c=and(c+16777213,16777215)
       ivni=and(ivni-c,16777215)
       if(vnisign .ne. 0) ivni = -ivni
       return
       end


      block data unidata
cc Initialized values in COMMON block for UNI, VNI, IUNI, IVNI
cc and RSTART.
      integer u(97),c,vnisign,ip,jp
      common /unidat/ u,c,vnisign,ip,jp
      data c,vnisign,ip,jp /362436,0,97,33/
      data (u(j), j=1,97)/
     +13697435, 3833429,12353926, 2287754, 3468638, 1232959, 8059805,
     +10745739, 4236676, 2095136, 1349346, 3672867,14563641,15473517,
     + 9897259, 2207061,  929657, 8109095, 5246947, 1066111, 8460236,
     +13162386,  501474,10402355,  352505, 2104170,12045925, 4350943,
     +13996856, 9897761, 6626452,15057436, 3168599,14038489, 8550848,
     + 5242835,13296102,11969002,   95246, 5917978, 8555838,13557738,
     + 1526088,11197237,15721125,14247931,  897046,15537441,16645456,
     +16279884, 1289925,14032128,10641039, 9961793, 2737638, 5073398,
     + 5231619, 2007688,15753584,12368695,12926325,10522018, 8692194,
     + 8531802,14755384,  276334, 9157821,  989353, 6093627,15866666,
     + 9532882, 3434034,  710155,  672726,12734991,13809842, 4832132,
     + 9753458,11325486,12137466, 3617374, 4913050, 9978642,12740205,
     +15754026, 4928136, 8545553,12893795, 8164497,12420478, 8192378,
     + 2028808, 1183983, 3474722,15616920,16298670,14606645/
      end


      subroutine indexx(n,arrin,indx)
cc  indexes the array arrin, outputs the array indx
      dimension arrin(n),indx(n)
      do 11 j=1,n
         indx(j)=j
 11   continue
      l=n/2+1
      ir=n
 10   continue
         if(l.gt.1)then
            l=l-1
            indxt=indx(l)
            q=arrin(indxt)
         else
            indxt=indx(ir)
            q=arrin(indxt)
            indx(ir)=indx(1)
            ir=ir-1 
            if(ir.eq.1)then
               indx(1)=indxt
               return
            endif
         endif
         i=l
         j=l+l
 20      if(j.le.ir)then
            if(j.lt.ir)then
               if(arrin(indx(j)).lt.arrin(indx(j+1)))j=j+1
            endif
            if(q.lt.arrin(indx(j)))then 
               indx(i)=indx(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
         go to 20
         endif
         indx(i)=indxt
      go to 10
      end 


      real FUNCTION EPAN(X)
cc  Epanechnikov kernel function
      real X
      IF (ABS(X) .LE. 1.) THEN
      EPAN=(3./4.)*(1. - X**2.)
      ELSE
      EPAN=0.
      ENDIF
      RETURN
      END


      real FUNCTION EPANL(X,Q)
cc  Epanechnikov lower tail kernel function
      real X,Q,GAMMAQ,ALPHAQ,BETAQ

      GAMMAQ=(0.75)*(((0.13333)+(Q**3./3.)-(Q**5./5.))*
     &         ((0.666667)+Q-(Q**3./3.)) - (0.0625)*
     &         ((1. - Q**2.)**4.))

      GAMMAQ=1./GAMMAQ

      ALPHAQ=((0.13333)+(Q**3./3.)-(Q**5./5.))*GAMMAQ

      BETAQ=((1. - Q**2.)**2.)*(GAMMAQ/4.)

      IF (X .LE. Q .AND. X .GE. -1.) THEN
      EPANL=(3./4.)*(1. - X**2.)*(ALPHAQ+BETAQ*X)
      ELSE
      EPANL=0.
      ENDIF
      RETURN
      END


      real FUNCTION EPANU(X,Q)
cc  Epanechnikov upper tail kernel function
      real X,Q,GAMMAQ,ALPHAQ,BETAQ,TEMPQ
      TEMPQ=0.13333-(Q**3./3.)+(Q**5./5.)

      GAMMAQ=1./((0.75)*(TEMPQ*((0.66667)-Q+(Q**3./3.)) -
     &       ((1. - Q**2.)**4.)/16.))

      ALPHAQ=TEMPQ*GAMMAQ
      BETAQ=-((1. - Q**2.)**2.)*(GAMMAQ/4.)
      IF (X .GE. Q .AND. X .LE. 1.) THEN
      EPANU=(3./4.)*(1. - X**2.)*(ALPHAQ+BETAQ*X)
      ELSE
      EPANU=0.
      ENDIF
      RETURN
      END


      real function alphdblp(t,band,time,censor,na,
     #              nsamp,ind,maxn)
cc Computes second derivative bandwidth for kernel estimation
      integer nsamp,ind,maxn
      real t,band,time(maxn),censor(maxn),na(maxn),ans
      real x,y

      ans=0.
      x=(t-time(ind))/band
      y=(-105./16.)*(1.-6.*(x**2.)+5.*(x**4.))
      if ((x.le.1.) .and. (x.ge.-1.)) then
      ans=y*(na(ind)-0)
      end if
c      ans=K1dlbpr((t-time(T1ind1))/band)*(na(ind)-0)
      do 10 j=(ind+1),nsamp
      if (abs(censor(j)).lt.0.01) goto 10
      x=(t-time(j))/band
      y=(-105./16.)*(1.-6.*(x**2.)+5.*(x**4.))
      if ((x.le.1.) .and. (x.ge.-1.)) then
c      ans=ans+K1dblpr((t-time(j))/band)*(na(j)-na(j-1))
      ans=ans+y*(na(j)-na(j-1))
      end if
  10  continue
      alphdblp=(ans/(band**3.))
      return
      end


      SUBROUTINE estp(N,NP,X,Z,DELTA,S0,S1,S2,BETA,BZ,C,GS,
     &       U,U2,F,NDEAD,NMIS,mxnsmp)
cc Maximum likelihood estimate for Cox model
      INTEGER mxnsmp
      parameter(mxncov=1)
      REAL X(mxnsmp),DELTA(mxnsmp),Z(mxncov,mxnsmp,mxnsmp)
      REAL S0(mxnsmp),S1(mxncov,mxnsmp),S2(mxncov,mxncov)
      REAL U(mxncov),U2(mxncov),F(mxncov,mxncov)
      REAL BZ(mxnsmp,mxnsmp),BETA(mxncov)

      C=0.0
      NDEAD=0
      NMIS=0
      GS=0.0
      DO 80 J=1,NP
         U(J)=0.0
            DO 60 K=1,NP
               F(J,K)=0.0
   60       CONTINUE
   80 CONTINUE
      DO 200 I=1,N
         IF (ABS(DELTA(I)).GT.1.E-10) NDEAD=NDEAD+1
         IF ((ABS(X(I)).LT.1.E-10).AND.(ABS(DELTA(I)).LT.1.E-10))
     $      NMIS=NMIS+1
            DO 180 K=1,N
              BZ(I,K)=0.0
              DO 170 J=1,NP
                 BZ(I,K)=BZ(I,K)+BETA(J)*Z(J,I,K)
  170         CONTINUE
  180    CONTINUE
  200 CONTINUE
      DO 400 I=1,N
         S0(I)=0.0
         DO 220 J=1,NP
            S1(J,I)=0.0
            DO 210 K=1,NP
               S2(J,K)=0.0
  210       CONTINUE
  220    CONTINUE
         DO 300 L=1,N
            IF (X(L).GE.X(I)) THEN
               S0(I)=S0(I)+EXP(BZ(L,I))
               DO 280 J=1,NP
                  S1(J,I)=S1(J,I)+EXP(BZ(L,I))*Z(J,L,I)
                  DO 260 K=1,NP
                     S2(J,K)=S2(J,K)+EXP(BZ(L,I))*Z(J,L,I)*Z(K,L,I)
  260             CONTINUE
  280          CONTINUE
            ENDIF
  300    CONTINUE
  350    C=C+DELTA(I)*(BZ(I,I)-LOG(S0(I)))
         DO 380 J=1,NP
            U(J)=U(J)+DELTA(I)*(Z(J,I,I)-S1(J,I)/S0(I))
            U2(J)=U(J)
            DO 360 K=1,NP
               F(J,K)=F(J,K)+DELTA(I)*(S2(J,K)/S0(I)
     $                -S1(J,I)*S1(K,I)/S0(I)**2.)
  360       CONTINUE
  380    CONTINUE
  400 CONTINUE
      CALL GAUSSJ(F,NP,NP,U,1,1)
      DO 800 J=1,NP
         BETA(J)=BETA(J)+U(J)
         DO 700 K=1,NP
            GS=GS+U2(J)*U2(K)*F(J,K)
  700    CONTINUE
  800 CONTINUE
      RETURN
      END
      


	SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
      PARAMETER (NMAX=200)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
	call rexit("Singular matrix in Cox")
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.) call rexit("Singular matrix in Cox") 
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END
	
	END 
