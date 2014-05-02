library(marked)
###################################################################################
# Simulation set 1 - set run_sims=TRUE to run the simulations
# variation in se(Phi) as a function of proportion permanently marked 
###################################################################################
# set of proportions of the population that is permanently marked
sf_set=c(0.02,0.05,.1,.2,.33,.5,.75,.99)
# set parameters used to simulate data
surv=0.9
pcap=0.4
beta1=-2
beta3=1
# define function for simulation set 1
sim.tl=function(size,nrep=1)
{
# save std error of survival and the beta for survival estimate
se=matrix(NA,ncol=length(sf_set),nrow=nrep)
beta.Surv=matrix(NA,ncol=length(sf_set),nrow=nrep)
# loop over each value of permanently marked
k=0
for(sf in sf_set)
{
  k=k+1
  cat("\n sf = ",sf)
  # loop over each replicate
  for(j in 1:nrep)
  {
  cat("\n j = ",j)
  # simulate data with dependent double tag loss using specified parameters with 5 cohorts
  # and 6 capture occasions (first release and 5 recapture occasions)  
  simdata=simHMM(data=data.frame(ch=c("++,0,0,0,0,0","0,++,0,0,0,0","0,0,++,0,0,0","0,0,0,++,0,0","0,0,0,0,++,+-"),
            freq=rep(size/5,5)),model="hmmcjs2tl",model.parameters=list(tau=list(formula=~I(tag1+tag2)+tag1:tag2)),
            initial=list(Phi=log(surv/(1-surv)),p=log(pcap/(1-pcap)),tau=c(beta1,beta3)))
  # select random set to be permanently marked and assign perm=yes for that set
  which.perm=sample(1:nrow(simdata),floor(sf*size),replace=F)
  simdata$perm="no"
  simdata$perm[which.perm]="yes"
  simdata$perm=factor(simdata$perm)
  # for those not permanently marked set all -- observations (both tags lost) to 0
  simdata$ch[simdata$perm=="no"]=gsub("--","0",simdata$ch[simdata$perm=="no"])
  # process data using perm for groups
  dp=process.data(simdata,model="hmmcjs2tl",groups="perm")
  # create design data
  ddl=make.design.data(dp)
  # set p=0 for state when both tags are lost and not permanently marked
  ddl$p$fix[ddl$p$perm=="no"&ddl$p$tag1==1&ddl$p$tag2==1]=0
  # fit model
  mod=crm(dp,ddl,model.parameters=list(tau=list(formula=~I(tag1+tag2)+tag1:tag2)),
            initial=list(Phi=2,p=-.3,tau=c(-2,1)),hessian=TRUE)
  # save std error and estimate
  se[j,k]=mod$results$beta.vcv[1,1]
  beta.Surv[j,k]=mod$results$beta$Phi
  }
}
return(list(se=se,beta=beta.Surv))
}
# run 100 simulations using n=500,1000,2500 and 5000
if(run_sims)
{
nrep=100
se.500=sim.tl(500,nrep=nrep)
se.1000=sim.tl(1000,nrep=nrep)
se.2500=sim.tl(2500,nrep=nrep)
se.5000=sim.tl(5000,nrep=nrep)
}
# define some functions to produce the figures used here and other simulation sets
plotsim.est=function(x,xnames,xlab="Fraction of population with permanent marks",ylab="Estimated survival",main="",link="logit")
{
   means=colMeans(x)
   if(link=="logit")
   {
      boxplot(plogis(x),names=xnames,xlab=xlab,ylab=ylab,main=main,outline=FALSE)
      points(1:length(means)+.1,plogis(means))
   }else {
      boxplot(x,names=xnames,xlab=xlab,ylab=ylab,main=main,outline=FALSE)
      points(1:length(means)+.1,means)
	}
   abline(h=surv)
   stderr=sqrt(apply(x,2,var)/nrow(x))
   if(link=="logit")
   {
      lcl=plogis(means-1.96*stderr)
      ucl=plogis(means+1.96*stderr)
	} else {
      lcl=means-1.96*stderr
      ucl=means+1.96*stderr
	}
	
   for(i in 1:length(lcl))
     lines(x=c(i+.1,i+.1),y=c(lcl[i],ucl[i]),lty=3)
}
plotsim.se=function(se,beta,xnames,main="",ylim=NULL,xlab="Fraction of population with permanent marks",ylab="Estimated std error of survival")
{
   boxplot(sqrt(se),names=xnames,xlab=xlab,ylab=ylab,main=main,ylim=ylim,outline=FALSE)
   stdev=sqrt(apply(beta,2,var))
   points(1:ncol(beta),stdev,pch=16)  
   invisible()
 }
library(Hmisc)
plotsim.cic=function(se,beta,xnames,xlab="Fraction of population with permanent marks",ylab="Confidence interval coverage",main="")
{
  y=colMeans(plogis(beta-1.96*sqrt(se))>surv| plogis(beta+1.96*sqrt(se))>surv )
  se.y=sqrt(y*(1-y)/nrow(beta))
  errbar(x=xnames,y=y,y-1.96*se.y,y+1.96*se.y,ylab="Confidence interval coverage",xlab="Fraction of population with permanent marks")
  title(main=main,outer=FALSE)
  abline(h=0.95)
  invisible()
 }

# produce figures 1-3 
pdf("Fig1estimate.pdf")
par(mfrow=c(2,2))
plotsim.est(se.500$beta,sf_set,main="n=500")
plotsim.est(se.1000$beta,sf_set,main="n=1000")
plotsim.est(se.2500$beta,sf_set,main="n=2500")
plotsim.est(se.5000$beta,sf_set,main="n=5000")
dev.off()

pdf("Fig2se.pdf")
par(mfrow=c(2,2))
plotsim.se(se.500$se,se.500$beta,sf_set,main="n=500",ylim=c(0,1))
plotsim.se(se.1000$se,se.1000$beta,sf_set,main="n=1000")
plotsim.se(se.2500$se,se.2500$beta,sf_set,main="n=2500")
plotsim.se(se.5000$se,se.5000$beta,sf_set,main="n=5000")
dev.off()

pdf("Fig3cic.pdf")
par(mfrow=c(2,2))
plotsim.cic(se.500$se,se.500$beta,sf_set,main="n=500")
plotsim.cic(se.1000$se,se.1000$beta,sf_set,main="n=1000")
plotsim.cic(se.2500$se,se.2500$beta,sf_set,main="n=2500")
plotsim.cic(se.5000$se,se.5000$beta,sf_set,main="n=5000")
dev.off()
