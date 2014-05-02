library(marked)
###################################################################################
# Simulation set 3 - set run_sims=TRUE to run the simulations
# Bias due to assumed independence with different levels of underlying dependence
################################################################################### 
# set of dependence values for beta3
dep_set=c(0,log(2),log(4),log(10))
# parameters other than dependence to simulate the data
surv=0.9
pcap=0.4
beta1=-2
# define function for simulation set 3
sim.dep=function(size,nrep=1,sf)
{
# save betas for tag loss estimate for fitted dependence model
beta.tau1=matrix(NA,ncol=length(dep_set),nrow=nrep)
beta.tau2=matrix(NA,ncol=length(dep_set),nrow=nrep)
# save beta for tag loss estimate for fitted independence model
beta.tau3=matrix(NA,ncol=length(dep_set),nrow=nrep)
# save beta for survival estimate for fitted dependence model
beta.Surv=matrix(NA,ncol=length(dep_set),nrow=nrep)
# save beta for survival estimate for fitted independence model
beta.Surv1=matrix(NA,ncol=length(dep_set),nrow=nrep)
# loop over each value of dependence (dep=0 is independence model)
k=0
for(dep in dep_set)
{
  k=k+1
  cat("\n dep = ",dep)
  # loop over each of the simulation replicates
  for(j in 1:nrep)
  {
  cat("\n j = ",j)
  # simulate data with dependent double tag loss using specified parameters with 5 cohorts
  # and 6 capture occasions (first release and 5 recapture occasions)  
  simdata=simHMM(data=data.frame(ch=c("++,0,0,0,0,0","0,++,0,0,0,0","0,0,++,0,0,0","0,0,0,++,0,0","0,0,0,0,++,+-"),freq=rep(size/5,5)),
                 model="hmmcjs2tl",model.parameters=list(tau=list(formula=~I(tag1+tag2)+tag1:tag2)),
                 initial=list(Phi=log(surv/(1-surv)),p=log(pcap/(1-pcap)),tau=c(beta1,dep)))
  # select random set to be permanently marked using argument sf and assign perm=yes for that set
  simdata$perm="no"
  which.perm=sample(1:nrow(simdata),floor(sf*size),replace=F)
  simdata$perm[which.perm]="yes"
  simdata$perm=factor(simdata$perm,levels=c("yes","no"))
  # for those not permanently marked set all -- observations (both tags lost) to 0
  simdata$ch[simdata$perm=="no"]=gsub("--","0",simdata$ch[simdata$perm=="no"])
  # process data using perm for groups
  dp=process.data(simdata,model="hmmcjs2tl",groups="perm")
  # create design data
  ddl=make.design.data(dp)
  # set p=0 for state when both tags are lost and not permanently marked
  ddl$p$fix[ddl$p$perm=="no"&ddl$p$tag1==1&ddl$p$tag2==1]=0
  # fit model assuming dependence and a model assuming independence of tag loss events
  mod=crm(dp,ddl,model.parameters=list(tau=list(formula=~I(tag1+tag2)+tag1:tag2)),initial=list(Phi=2,p=-.4,tau=c(-2,dep)),hessian=TRUE)
  mod1=crm(dp,ddl,model.parameters=list(tau=list(formula=~I(tag1+tag2))),initial=list(Phi=2,p=-.4,tau=c(-2)),hessian=TRUE)
  # save betas
  beta.Surv[j,k]=mod$results$beta$Phi 
  beta.tau1[j,k]=mod$results$beta$tau[1]
  beta.tau2[j,k]=mod$results$beta$tau[2]
  beta.tau3[j,k]=mod1$results$beta$tau[1]
  beta.Surv1[j,k]=mod1$results$beta$Phi 
  }
}
return(list(beta.tau1=beta.tau1,beta.tau2=beta.tau2,beta.tau3=beta.tau3,beta.S=beta.Surv,beta.S1=beta.Surv1))
}
# run 100 replicate simulations using sf=0.01, 0.5, 0.99
if(run_sims)
{
nrep=100
se.dep1=sim.dep(5000,nrep=nrep,sf=0.01)
se.dep5=sim.dep(5000,nrep=nrep,sf=0.5)
se.dep99=sim.dep(5000,nrep=nrep,sf=0.99)
}
# produce figures 7-9
pdf("Fig7estimate_dep1.pdf")
surv=0.9
par(mfrow=c(2,3))
plotsim.est(se.dep1$beta.S,formatC(dep_set,digits=3),xlab="Log dependence measure",main="1% permanent marks")
plotsim.est(se.dep5$beta.S,formatC(dep_set,digits=3),xlab="Log dependence measure",main="50% permanent marks")
plotsim.est(se.dep99$beta.S,formatC(dep_set,digits=3),xlab="Log dependence measure",main="99% permanent marks")
plotsim.est(se.dep1$beta.S1,formatC(dep_set,digits=3),xlab="Log dependence measure")
plotsim.est(se.dep5$beta.S1,formatC(dep_set,digits=3),xlab="Log dependence measure")
plotsim.est(se.dep99$beta.S1,formatC(dep_set,digits=3),xlab="Log dependence measure")
dev.off()
pdf("Fig8estimate_dep2.pdf")
surv=plogis(-2)
par(mfrow=c(2,3))
plotsim.est(se.dep1$beta.tau1,formatC(dep_set,digits=3),xlab="Log dependence measure",main="1% permanent marks",ylab="Tag loss estimate")
plotsim.est(se.dep5$beta.tau1,formatC(dep_set,digits=3),xlab="Log dependence measure",main="50% permanent marks",ylab="Tag loss estimate")
plotsim.est(se.dep99$beta.tau1,formatC(dep_set,digits=3),xlab="Log dependence measure",main="99% permanent marks",ylab="Tag loss estimate")
surv=dep_set
plotsim.est(se.dep1$beta.tau2,formatC(dep_set,digits=3),xlab="Log dependence measure",link="identity",ylab="Tag loss dependence estimate")
plotsim.est(se.dep5$beta.tau2,formatC(dep_set,digits=3),xlab="Log dependence measure",link="identity",ylab="Tag loss dependence estimate")
plotsim.est(se.dep99$beta.tau2,formatC(dep_set,digits=3),xlab="Log dependence measure",link="identity",ylab="Tag loss dependence estimate")
dev.off()
pdf("Fig9estimate_dep3.pdf")
surv=plogis(-2)
par(mfrow=c(2,3))
plotsim.est(se.dep1$beta.tau1,formatC(dep_set,digits=3),xlab="Log dependence measure",main="1% permanent marks",ylab="Tag loss estimate")
plotsim.est(se.dep5$beta.tau1,formatC(dep_set,digits=3),xlab="Log dependence measure",main="50% permanent marks",ylab="Tag loss estimate")
plotsim.est(se.dep99$beta.tau1,formatC(dep_set,digits=3),xlab="Log dependence measure",main="99% permanent marks",ylab="Tag loss estimate")
plotsim.est(se.dep1$beta.tau3,formatC(dep_set,digits=3),xlab="Log dependence measure",ylab="Tag loss estimate")
plotsim.est(se.dep5$beta.tau3,formatC(dep_set,digits=3),xlab="Log dependence measure",ylab="Tag loss estimate")
plotsim.est(se.dep99$beta.tau3,formatC(dep_set,digits=3),xlab="Log dependence measure",ylab="Tag loss estimate")
dev.off()
