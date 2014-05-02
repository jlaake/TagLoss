library(marked)
bear=read.table("bear.txt",header=TRUE)
bear$ch=as.character(bear$ch)
run_examples=TRUE

##########################################
# function to fit double tag loss models
##########################################
do_bear=function()
{
#Phi models
Phi.1=list(formula=~sex)
Phi.2=list(formula=~ageclass)
Phi.3=list(formula=~ageclass1)
Phi.4=list(formula=~ageclass+sex)
Phi.5=list(formula=~ageclass1+sex)
Phi.6=list(formula=~ageclass1*sex)
#p models
p.1=list(formula=~time)
p.2=list(formula=~time+sex)
p.3=list(formula=~time+sex+ageclass)
#tau models
# independent tag loss - constant rate across tags and age
tau.1=list(formula=~I(tag1+tag2))
# independent tag loss - constant rate across tags but varying with tag age
tau.2=list(formula=~I(tag1+tag2)+tagAge:I(tag1+tag2))
# dependent tag loss but constant for tags
tau.3=list(formula=~I(tag2+tag2)+tag1:tag2)
# dependent tag loss - constant rate across tags but varying with tag age
tau.4=list(formula=~I(tag1+tag2)+tagAge:I(tag1+tag2)+tag1:tag2)
# dependent tag loss - constant rate across tags but varying with tag age and dependency varies with tagAge but initially independent
tau.5=list(formula=~I(tag1+tag2)+tagAge:I(tag1+tag2)+tag1:tag2:tagAge)
# dependent tag loss - constant rate across tags but varying with tag age and dependency varies with tagAge
tau.6=list(formula=~I(tag1+tag2)+tagAge:I(tag1+tag2)+ tag1:tag2 + tag1:tag2:tagAge)
# + sex
tau.7=list(formula=~I(tag1+tag2)+I(tag1+tag2):male)
tau.8=list(formula=~I(tag1+tag2)+tagAge:I(tag1+tag2)+I(tag1+tag2):male + male:tagAge:I(tag1+tag2))
tau.9=list(formula=~I(tag2+tag2)+I(tag1+tag2):male+tag1:tag2+male:tag1:tag2)
tau.10=list(formula=~I(tag1+tag2)+male:I(tag1+tag2)+tagAge:I(tag1+tag2)+male:tagAge:I(tag1+tag2)+tag1:tag2+male:tag1:tag2)
tau.11=list(formula=~I(tag1+tag2)+male:I(tag1+tag2)+tagAge:I(tag1+tag2)+male:tagAge:I(tag1+tag2)+tag1:tag2:tagAge+male:tag1:tag2:tagAge)
# create list of models which will be 120 (6*2*10)
cml=create.model.list(c("Phi","tau","p"))
# run each model
return(crm.wrapper(cml,data=dp,ddl=ddl,initial=mod0,external=FALSE,hessian=TRUE))
}

do_bear0=function()
{
#Phi models
Phi.01=list(formula=~sex)
Phi.02=list(formula=~ageclass)
Phi.03=list(formula=~ageclass1)
Phi.04=list(formula=~ageclass+sex)
Phi.05=list(formula=~ageclass1+sex)
Phi.06=list(formula=~ageclass1*sex)
#p models
p.01=list(formula=~time)
p.02=list(formula=~time+sex)
p.03=list(formula=~time+sex+ageclass)
#tau models
# independent tag loss - constant rate across tags and age
tau.01=list(formula=~I(tag1+tag2))
# independent tag loss - sex-specific rate across tags and age
tau.02=list(formula=~I(tag1+tag2)+male:I(tag1+tag2))
# independent tag loss - constant rate across tags but varying with tag age
tau.03=list(formula=~I(tag1+tag2)+tagAge:I(tag1+tag2))
# independent tag loss - sex-specific rate across tags and varying with tag age*sex
tau.04=list(formula=~I(tag1+tag2)+tagAge:I(tag1+tag2)+male:I(tag1+tag2)+male:tagAge:I(tag1+tag2))
# create list of models which will be 60 (5*2*6)
cml=create.model.list(c("Phi","tau","p"))
# run each model
return(crm.wrapper(cml,data=dp,ddl=ddl,initial=mod0,external=FALSE,hessian=TRUE))
}

###########################################
# function to fit single tag loss models
###########################################
do_bear_single=function()
{
#Phi models
Phi.s1=list(formula=~sex)
Phi.s2=list(formula=~ageclass)
Phi.s3=list(formula=~ageclass1)
Phi.s4=list(formula=~ageclass+sex)
Phi.s5=list(formula=~ageclass1+sex)
Phi.s6=list(formula=~ageclass1*sex)
#p models
p.s1=list(formula=~time)
p.s2=list(formula=~time+sex)
p.s3=list(formula=~time+sex+ageclass)
#tau models
# constant rate 
tau.s1=list(formula=~1)
# varying with tag age
tau.s2=list(formula=~tagAge)
# varying with sex
tau.s3=list(formula=~sex)
# create list of models which will be 60 (5*2*6)
cml=create.model.list(c("Phi","tau","p"))
# run each model
return(crm.wrapper(cml,data=dp,ddl=ddl,external=FALSE,hessian=TRUE))
}
# function to fit single tag loss models

do_bear_single0=function()
{
#Phi models
Phi.s01=list(formula=~sex)
Phi.s02=list(formula=~ageclass)
Phi.s03=list(formula=~ageclass1)
Phi.s04=list(formula=~ageclass+sex)
Phi.s05=list(formula=~ageclass1+sex)
Phi.s06=list(formula=~ageclass1*sex)
#p models
p.s01=list(formula=~time)
p.s02=list(formula=~time+sex)
p.s03=list(formula=~time+sex+ageclass)
# create list of models which will be 60 (5*2*6)
cml=create.model.list(c("Phi","p"))
# run each model
return(crm.wrapper(cml,data=dp,ddl=ddl,external=FALSE,hessian=TRUE))
}

if(run_examples)
{
################### sampled double tag example code ###########################
# process data frame for model using HMMcjs and 2tl (double tag loss); sex are groups and time begins at 2002
sample.bear=bear
sample.bear$ch[sample.bear$perm=="no"]=gsub("--","0",sample.bear$ch[sample.bear$perm=="no"])
dp=process.data(sample.bear,model="hmmcjs2tl",groups=c("sex","perm"),begin.time=2002)
# create default design data
ddl=make.design.data(dp)
# fix p=0 when both tags missing and not permanently marked
ddl$p$fix[ddl$p$perm=="no"&ddl$p$tag1==1&ddl$p$tag2==1]=0
# create ageclass variables for survival
ddl$Phi$ageclass=cut(ddl$Phi$Age,c(-10,1,2,30),labels=c("Yearling","2yr","3+yrs"))
ddl$Phi$ageclass1=cut(ddl$Phi$Age,c(-10,1,30),labels=c("Yearling","2+yrs"))
# create tagAge variable for tau; tag1 and tag2 variables are created automatically
ddl$tau$tagAge=ddl$tau$Time-ddl$tau$Cohort
ddl$tau$male=ifelse(ddl$tau$sex=="M",1,0)
ddl$p$ageclass=cut(ddl$p$Age,c(-10,2,30),labels=c("2yr","3+yrs"))
mod0=crm(dp,ddl)
doubletag_results=do_bear()

################### 0 perm marked double tag example code ###########################
# process data frame for model using HMMcjs and 2tl (double tag loss); sex are groups and time begins at 2002
sample.bear=bear
sample.bear$perm="no"
sample.bear$perm=factor(rep("no",nrow(sample.bear)),levels=c("no","yes"))
sample.bear$ch[sample.bear$perm=="no"]=gsub("--","0",sample.bear$ch[sample.bear$perm=="no"])
dp=process.data(sample.bear,model="hmmcjs2tl",groups=c("sex","perm"),begin.time=2002)
# create default design data
ddl=make.design.data(dp)
# fix p=0 when both tags missing and not permanently marked
ddl$p$fix[ddl$p$perm=="no"&ddl$p$tag1==1&ddl$p$tag2==1]=0
# create ageclass variables for survival
ddl$Phi$ageclass=cut(ddl$Phi$Age,c(-10,1,2,30),labels=c("Yearling","2yr","3+yrs"))
ddl$Phi$ageclass1=cut(ddl$Phi$Age,c(-10,1,30),labels=c("Yearling","2+yrs"))
# create tagAge variable for tau; tag1 and tag2 variables are created automatically
ddl$tau$tagAge=ddl$tau$Time-ddl$tau$Cohort
ddl$tau$male=ifelse(ddl$tau$sex=="M",1,0)
ddl$p$ageclass=cut(ddl$p$Age,c(-10,2,30),labels=c("2yr","3+yrs"))
mod0=crm(dp,ddl)
doubletag_results0=do_bear0()

################## single tag example code ############################
# single tag data
single.bear=bear
single.bear$ch=sapply(lapply(strsplit(single.bear$ch,","), function (x) ifelse(nchar(x)==1,x,substr(x,1,1))),paste,collapse=",")
single.bear$ch[single.bear$perm=="no"]=gsub("-","0",single.bear$ch[single.bear$perm=="no"])
# process data frame for model using HMMcjs and 1tl (double tag loss); sex are groups and time begins at 2002
dp=process.data(single.bear,model="hmmcjs1tl",groups=c("sex","perm"),begin.time=2002)
# create default design data
ddl=make.design.data(dp)
ddl$p$fix[ddl$p$stratum=="0"&ddl$p$perm=="no"]=0
# create ageclass variables for survival
ddl$Phi$ageclass=cut(ddl$Phi$Age,c(-10,1,2,30),labels=c("Yearling","2yr","3+yrs"))
ddl$Phi$ageclass1=cut(ddl$Phi$Age,c(-10,1,30),labels=c("Yearling","2+yrs"))
# create tagAge variable for tau
ddl$tau$tagAge=ddl$tau$Time-ddl$tau$Cohort
ddl$p$ageclass=cut(ddl$p$Age,c(-10,2,30),labels=c("2yr","3+yrs"))
mod0=crm(dp,ddl)
singletag_results=do_bear_single()

################## single tag example code - no permanent marks ############################
# single tag data
single.bear=bear
single.bear$perm="no"
single.bear$perm=factor(rep("no",nrow(single.bear)),levels=c("no","yes"))
single.bear$ch=sapply(lapply(strsplit(single.bear$ch,","), function (x) ifelse(nchar(x)==1,x,substr(x,1,1))),paste,collapse=",")
single.bear$ch[single.bear$perm=="no"]=gsub("-","0",single.bear$ch[single.bear$perm=="no"])
single.bear$ch[single.bear$perm=="no"]=gsub("\\+","1",single.bear$ch[single.bear$perm=="no"])
# process data frame for model using HMMcjs and 1tl (double tag loss); sex are groups and time begins at 2002
dp=process.data(single.bear,model="hmmcjs",groups=c("sex","perm"),begin.time=2002)
# create default design data
ddl=make.design.data(dp)
# create ageclass variables for survival
ddl$Phi$ageclass=cut(ddl$Phi$Age,c(-10,1,2,30),labels=c("Yearling","2yr","3+yrs"))
ddl$Phi$ageclass1=cut(ddl$Phi$Age,c(-10,1,30),labels=c("Yearling","2+yrs"))
# create tagAge variable for tau
ddl$tau$tagAge=ddl$tau$Time-ddl$tau$Cohort
ddl$p$ageclass=cut(ddl$p$Age,c(-10,2,30),labels=c("2yr","3+yrs"))
mod0=crm(dp,ddl)
singletag_results0=do_bear_single0()
}


