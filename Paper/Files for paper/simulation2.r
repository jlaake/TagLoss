library(marked)
###################################################################################
# Simulation set 2 - set run_sims=TRUE to run the simulations
# variation in se(Phi) as a function of proportion permanently marked and p
###################################################################################
# set of capture probabilities
p_set=c(0.05,.1,.2,.4,.5)
# parameters other than p to simulate the data
surv=0.9
beta1=-2
beta3=1
# define function for simulation set 2
sim.tlp=function(size,sf,nrep=1)
{
# save std error of survival and the beta for survival estimate
se=matrix(NA,ncol=length(p_set),nrow=nrep)
beta.Surv=matrix(NA,ncol=length(p_set),nrow=nrep)
# loop over each value of capture probability
k=0
for(p in p_set)
{
  cat("\n p = ",p)
  k=k+1
  # loop over each of the replicates
  for(j in 1:nrep)
  {
  cat("\n j = ",j)
  # simulate data with dependent double tag loss using specified parameters with 5 cohorts
  # and 6 capture occasions (first release and 5 recapture occasions)  
  simdata=simHMM(data=data.frame(ch=c("++,0,0,0,0,0","0,++,0,0,0,0","0,0,++,0,0,0","0,0,0,++,0,0","0,0,0,0,++,+-"),freq=rep(size/5,5)),
                 model="hmmcjs2tl",model.parameters=list(tau=list(formula=~I(tag1+tag2)+tag1:tag2)),
                 initial=list(Phi=log(surv/(1-surv)),p=log(p/(1-p)),tau=c(beta1,beta3)))
  # select random set to be permanently marked using argument sf and assign perm=yes for that set
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
  mod=crm(dp,ddl,model.parameters=list(tau=list(formula=~I(tag1+tag2)+tag1:tag2)),initial=list(Phi=2,p=.3,tau=c(-2,1)),hessian=TRUE)
  # save std error and estimate for survival
  se[j,k]=mod$results$beta.vcv[1,1]
  beta.Surv[j,k]=mod$results$beta$Phi
  }
}
return(list(se=se,beta=beta.Surv))
}
# run 100 simulations using sf=0.1,0.2,0.5,0.9
if(run_sims)
{
nrep=100
se.1=sim.tlp(1000,sf=0.1,nrep=nrep)
se.2=sim.tlp(1000,sf=0.2,nrep=nrep)
se.5=sim.tlp(1000,sf=0.5,nrep=nrep)
se.9=sim.tlp(1000,sf=0.9,nrep=nrep)
}
# produce figures 4-6
pdf("Fig4estimate_byp.pdf")
par(mfrow=c(2,2))
plotsim.est(se.1$beta,p_set,main="10% permanent marks",xlab="Capture probability")
plotsim.est(se.2$beta,p_set,main="20% permanent marks",xlab="Capture probability")
plotsim.est(se.5$beta,p_set,main="50% permanent marks",xlab="Capture probability")
plotsim.est(se.9$beta,p_set,main="90% permanent marks",xlab="Capture probability")
dev.off()
pdf("Fig5se_byp.pdf")
par(mfrow=c(2,2))
plotsim.se(se.1$se,se.1$beta,p_set,main="10% permanent marks",xlab="Capture probability")
plotsim.se(se.2$se,se.2$beta,p_set,main="20% permanent marks",xlab="Capture probability")
plotsim.se(se.5$se,se.5$beta,p_set,main="50% permanent marks",xlab="Capture probability")
plotsim.se(se.9$se,se.9$beta,p_set,main="90% permanent marks",xlab="Capture probability")
dev.off()
pdf("Fig6cic_byp.pdf")
par(mfrow=c(2,2))
plotsim.cic(se.1$se,se.1$beta,p_set,main="10% permanent marks",xlab="Capture probability")
plotsim.cic(se.2$se,se.2$beta,p_set,main="20% permanent marks",xlab="Capture probability")
plotsim.cic(se.5$se,se.5$beta,p_set,main="50% permanent marks",xlab="Capture probability")
plotsim.cic(se.9$se,se.9$beta,p_set,main="90% permanent marks",xlab="Capture probability")
dev.off()
