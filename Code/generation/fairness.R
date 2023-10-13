##**********************************************#
## Analyse SPARRA fairness measures          ####
##**********************************************#
## James Liley, 2022-23
##


##**********************************************#
## Directories, packages and scripts         ####
##**********************************************#

# Working directory and location (windows or linux machine)

# Scripts
source("fairness_functions.R")
source("auxiliary.R")

# Packages
if (location=="Windows") {
  library(pROC)
  library(PRROC)
  library(xgboost)
  library(glmnet)
  library(Metrics)
  library(data.table)
  library(cvAUC)
  library(latex2exp)
  library(dplyr)
} else {
  library(h2o)
  library(tidyverse)
  library(fst)
}

# Directories


# In text outputs, everything after this marker should be readable by R.
sink_marker="\n\n\n***************\n\n"
options(width=1e4)


##**********************************************#
## Generate all_pred_eth_geog if nonexistent ####
##**********************************************#

all_pred_file=paste0(out_dir,"all_pred_ethn.RDS")

# Add ethnicity and geography information, if not already done
if (!file.exists(all_pred_file)) {
  all_pred0=readRDS(paste0("all_pred.RDS")) # SPARRAv4 output
  
  # Merge ethnicity data
  eth=read.csv("Ethnicity file.csv.gz")
  all_pred0=merge(all_pred0,eth,by.x="id",by.y="UNIQUE_STUDY_ID",all.x=TRUE,all.y=FALSE)
  
  # Merge geographic data
  post1=read.csv("postcodes_20170401_part1.csv.gz")
  post2=read.csv("postcodes_20170401_part2.csv.gz")
  post3=read.csv("postcodes_20170401_part3.csv.gz")
  post=rbind(post1,post2,post3)[,c("UNIQUE_STUDY_ID","islands_2022","UR2_2016")]
  all_pred0=merge(all_pred0,post,by.x="id",by.y="UNIQUE_STUDY_ID",all.x=TRUE,all.y=FALSE)
  
  # Save
  saveRDS(all_pred0,file=all_pred_file)
}


##**********************************************#
## Read in data matrix                       ####
##**********************************************#


all_pred=readRDS(all_pred_file)


## Compress to only usable rows
dat0=all_pred[which((all_pred$time==min(all_pred$time)) &
                    (is.finite(all_pred$v3 + all_pred$super))),]


## Rename 'super' to 'v4'
colnames(dat0)[which(colnames(dat0)=="super")]="v4"



##**********************************************#
## Settings                                  ####
##**********************************************#

save_plots=TRUE # Save plots to PDF 

alpha=0.05 # Type 1 error for confidence intervals

ncut=20 # N cutoffs
cutoffs=seq(0.01,0.9,length=ncut) # Cutoffs

# Specification of fairness type. The ith row name means we compare
# P(Yhat=prob_yhat,Y=prob_Y,sex=M*prob_g | Yhat=cond_yhat,Y=cond_Y,sex=M*cond_g)
# P(Yhat=prob_yhat,Y=prob_Y,sex=F*prob_g | Yhat=cond_yhat,Y=cond_Y,sex=F*cond_g)
# with the corresponding term admitted if the column value is NA.
specs=rbind( # Yhat  Y     Group  |  Yhat  Y     Group 
  fpp  =  c(   1,    0,    NA,       NA,   NA,   1),    # False positive parity: compare P(Yhat=1,Y=0|group)
  fdrp =  c(   NA,   0,    NA,       1,    NA,   1),    # False discovery rate parity: compare P(Y=0|Yhat=1,group)
  fprp =  c(   1,    NA,   NA,       NA,   0,    1),    # False positive rate parity: compare P(Yhat=1,Y=0|group)
  rp   =  c(   1,    NA,   NA,       NA,   1,    1),    # Recall parity: compare P(Yhat=1|Y=1,group)
  irp  =  c(   0,    NA,   NA,       NA,   0,    1),    # Inverse recall parity: compare P(Yhat=0|Y=0,group)
  fnp  =  c(   0,    1,    NA,       NA,   NA,   1),    # False negative parity: compare P(Yhat=0,Y=1|group)
  forp =  c(   NA,   1,    NA,       0,    NA,   1),    # False positive parity: compare P(Yhat=1,Y=0|group)
  fnrp =  c(   0,    NA,   NA,       NA,   1,    1)     # False negative rate parity: compare P(Yhat=0|Y=1,group)
)
colnames(specs)=c("prob_yhat","prob_y","prob_g","cond_yhat","cond_y","cond_g")




##**********************************************#
## Define groupings                          ####
##**********************************************#


# List of groupings to consider. We compare the two groups in each list element: for instance, male and female
# Potentially we could extend this code to cover more than two groups.
groupings=list(
  Sex=list(Male=dat0$id[which(dat0$sexM==TRUE)],Female=dat0$id[which(dat0$sexM==FALSE)]), # Compare between male/female
  SIMD=list(Most_deprived=dat0$id[which(dat0$simd %in% c(1,2))],Least_deprived=dat0$id[which(dat0$simd %in% c(9,10))]), # Compare between more or less deprived
  Age=list(O65=dat0$id[which(dat0$age>65)], U25=dat0$id[which(dat0$age<25)]), # Compare older and younger people
  Ethnicity=list(White=dat0$id[which(dat0$ethnic_group_desc_short=="White")],Nonwhite=dat0$id[which(dat0$ethnic_group_desc_short!="White")]),
  Urban_rural=list(Urban=dat0$id[which(dat0$UR2_2016==1)],Rural=dat0$id[which(dat0$UR2_2016==2)]),
  Mainland_island=list(Mainland=dat0$id[which(dat0$islands_2022==0)],Island=dat0$id[which(dat0$islands_2022>0)])
)

categories=list(all=dat0$id)
cnames=c("all")
ind=2
for (i in 1:length(groupings)) {
  for (ig in 1:length(groupings[[i]])) {
    categories[[ind]]=groupings[[i]][[ig]]
    cnames[ind]=names(groupings[[i]])[ig]
    ind=ind+1
  }
}
names(categories)=cnames





##**********************************************#
## Dataset details                           ####
##**********************************************#


# Categories including those above, but also 'neither' and 'NA' for each groupiing
allcats=list(all=dat0$id)
cnames=c("all")
ind=2
for (i in 1:length(groupings)) {
  for (ig in 1:length(groupings[[i]])) {
    allcats[[ind]]=groupings[[i]][[ig]]
    cnames[ind]=paste0(names(groupings)[i],"_",names(groupings[[i]])[ig])
    ind=ind+1
  }
}
names(allcats)=cnames

# IDs of individuals for which this grouping category is missing
allcats[['Sex_missing']]=dat0$id[which(!is.finite(dat0$sexM))]
allcats[['Age_missing']]=dat0$id[which(!is.finite(dat0$age))]
allcats[['SIMD_missing']]=dat0$id[which(!is.finite(dat0$simd))]
allcats[['Ethnicity_missing']]=dat0$id[which(is.na(dat0$ethnic_group_desc_short))]
allcats[['Urban_rural_missing']]=dat0$id[which(!is.finite(dat0$UR2_2016))]
allcats[['Mainland_island_missing']]=dat0$id[which(!is.finite(dat0$islands_2022))]

# IDs of individuals for which this grouping category takes a value not in categories
allcats[['Sex_other']]=setdiff(dat0$id[which(is.finite(dat0$sexM))], unlist(groupings['Sex']))
allcats[['Age_other']]=setdiff(dat0$id[which(is.finite(dat0$age))], unlist(groupings['Age']))
allcats[['SIMD_other']]=setdiff(dat0$id[which(is.finite(dat0$simd))], unlist(groupings['SIMD']))
allcats[['Ethnicity_other']]=setdiff(dat0$id[which(!is.na(dat0$ethnic_group_desc_short))], unlist(groupings['Ethnicity']))
allcats[['Urban_rural_other']]=setdiff(dat0$id[which(is.finite(dat0$UR2_2016))], unlist(groupings['Urban_rural']))
allcats[['Mainland_island_other']]=setdiff(dat0$id[which(is.finite(dat0$islands_2022))], unlist(groupings['Mainland_island']))


bigtab=c()
xmean=function(x) sum(x,na.rm=T)/length(x)
for (i in 1:length(allcats)) {
  dat1=dat0[match(allcats[[i]],dat0$id),]
  if (dim(dat1)[1]>0) {
    atab=c(
      dim(dat1)[1], # N
      xmean(dat1$sexM), # proportion male
      xmean(dat1$age), # mean age
      sd(dat1$age), # SD age
      xmean(dat1$simd<=2), # proportion SIMD 1,2
      xmean(dat1$simd>2 & dat1$simd<=4), # proportion SIMD 3.4
      xmean(dat1$simd>4 & dat1$simd<=6), # proportion SIMD 5,6
      xmean(dat1$simd>6 & dat1$simd<=8), # proportion SIMD 7,8
      xmean(dat1$simd>8), # proportion SIMD 9,10
      xmean(dat1$ethnic_group_desc_short=="White"), # proportion white
      xmean(dat1$ethnic_group_desc_short!="White"), # proportion non-white
      xmean(dat1$UR2_2016==1), # proportion urban
      xmean(dat1$UR2_2016==2), # proportion rural
      xmean(dat1$islands_2022==0), # proportion mainland
      xmean(dat1$islands_2022>0) # island
    )
    bigtab=rbind(bigtab,atab)
  } 
}
rownames(bigtab)=names(allcats)[which(unlist(lapply(allcats,length))>0)]
colnames(bigtab)=c("N",
                   "sex_prop_male",
                   "age_mean","age_sd",
                   paste0("simd_quintile",1:5),
                   "ethnicity_white",
                   "ethnicity_nonwhite",
                   "Urbanrural_urban",
                   "Urbanrural_rural",
                   "Mainlandisland_mainland",
                   "Mainlandisland_island")

# Privacy check
nbigtab=bigtab # Will become numerators
for (i in 2:dim(nbigtab)[2]) nbigtab[,i]=nbigtab[,i]*nbigtab[,1]
bigtab[which(nbigtab<5)]=0 # Whenever a fraction was calculated from fewer than 5 cases, set it to 0.


# Adjoin mean and median scores
scoretab=c()
for (i in 1:length(allcats)) {
  dat1=dat0[match(allcats[[i]],dat0$id),]
  if (dim(dat1)[1]>0) {
    atab=c(
      mean(dat1$v3), # Mean v3
      median(dat1$v3), # Median v3
      sd(dat1$v3), # Standard deviation v3
      mean(dat1$v4), # Mean v3
      median(dat1$v4), # Median v3
      sd(dat1$v4) # Standard deviation v3
    )
    scoretab=rbind(scoretab,atab)
  } 
}
colnames(scoretab)=c("mean_v3",
                   "median_v3",
                   "sd_v3",
                   "mean_v4",
                   "median_v4",
                   "sd_v4")
bigtab=cbind(bigtab,scoretab)
     
write.table(bigtab,file=paste0(plot_dir,"overview_table.csv"),sep=",")
write.table(nbigtab,file=paste0(plot_dir,"overview_table_numerators.csv"),sep=",")






##**********************************************#
## Initialise text output files              ####
##**********************************************#


for (g in 1:length(groupings)) {
  out_text_file=paste0(plot_dir,names(groupings)[g],"/summary_",names(groupings)[g],".txt")
  file.create(out_text_file)
}



##**********************************************#
## Performance (ROC, PRC, CAL) by group      ####
##**********************************************#

for (version in c(3,4)) {
  
  # If version==3, exclude individuals who died during prediction period
  if (version==3){
    wsub=which((dat0$target==0)|dat0$reason %in% c("B","E"))
    dat=dat0[wsub,]
    score=dat$v3 
  } else {
    dat=dat0
    score=dat$v4
  }
  target=dat$target
  
  # Overall calibration and ROC
  roc_overall=getroc(target,score)
  prc_overall=getprc(target,score)
  cal_overall=plotcal(target,score,kernel=FALSE,n=10,plot=FALSE)  
  
  for (g in 1:length(groupings)) {
   
    # Where to save plots and write to file
    file_prefix=paste0(names(groupings)[g],"/")
    out_text_file=paste0(plot_dir,names(groupings)[g],"/summary_",names(groupings)[g],".txt")

        
    # Define groups (male/female, old/young etc)
    groupA=which(dat$id %in% groupings[[g]][[1]])
    groupB=which(dat$id %in% groupings[[g]][[2]])
    
    
    ### ROC curves
    rocA=getroc(target[groupA],score[groupA])
    rocB=getroc(target[groupB],score[groupB])

    xname=paste0(plot_dir,file_prefix,"overall_roc_v",version,"_",names(groupings)[g])
    if (save_plots) pdf(paste0(xname,".pdf"),width=3,height=3.5)
    
    aucT=signif(roc_overall$auc,digits=2); seT=signif(roc_overall$se,digits=2)
    aucA=signif(rocA$auc,digits=2); seA=signif(rocA$se,digits=2);
    aucB=signif(rocB$auc,digits=2); seB=signif(rocA$se,digits=2);
    labs=c(paste0("Overall: ",aucT," (",seT,")"),
           paste0(names(groupings[[g]])[1],": ",aucA," (",seA,")"),
           paste0(names(groupings[[g]])[2],": ",aucB," (",seB,")"))
    xcol=c("black","red","blue")
    roc_2panel(list(roc_overall,rocA,rocB),labels=labs,col=xcol,title="AUROC (SE)",text.col=xcol,title.col="black")
    if (save_plots) dev.off()
    
    # Object to save in file to reproduce plot
    obj=list(roc_overall=roc_overall,rocA=rocA,rocB=rocB,labs=labs,xcol=xcol)
    obj_name=paste0("roc_v",version,"_",names(groupings)[g])
    
    # Save
    sink(out_text_file,append=TRUE)
    cat("\n\n")
    cat(obj_name,"="); dput(obj); 
    cat("\n\n")
    sink()

    
    ### PR curves
    prcA=getprc(target[groupA],score[groupA])
    prcB=getprc(target[groupB],score[groupB])
    
    xname=paste0(plot_dir,file_prefix,"overall_prc_v",version,"_",names(groupings)[g])
    if (save_plots) pdf(paste0(xname,".pdf"),width=3,height=3.5)
    
    aucT=signif(prc_overall$auc,digits=2); seT=signif(prc_overall$se,digits=2)
    aucA=signif(prcA$auc,digits=2); seA=signif(prcA$se,digits=2);
    aucB=signif(prcB$auc,digits=2); seB=signif(prcA$se,digits=2);
    labs=c(paste0("Overall: ",aucT," (",seT,")"),
           paste0(names(groupings[[g]])[1],": ",aucA," (",seA,")"),
           paste0(names(groupings[[g]])[2],": ",aucB," (",seB,")"))
    xcol=c("black","red","blue")
    prc_2panel(list(prc_overall,prcA,prcB),labels=labs,col=xcol,title="AUPRC (SE)",text.col=xcol,title.col="black")
    if (save_plots) dev.off()
    
    # Object to save in file to reproduce plot
    obj=list(prc_overall=prc_overall,prcA=prcA,prcB=prcB,labs=labs,xcol=xcol)
    obj_name=paste0("prc_v",version,"_",names(groupings)[g])
    
    # Save
    sink(out_text_file,append=TRUE)
    cat("\n\n")
    cat(obj_name,"="); dput(obj); 
    cat("\n\n")
    sink()
    
    
    
    ### Calibration curves
    calA=plotcal(target[groupA],score[groupA],kernel=FALSE,n=10,plot=FALSE) 
    calB=plotcal(target[groupB],score[groupB],kernel=FALSE,n=10,plot=FALSE) 
    
    xname=paste0(plot_dir,file_prefix,"overall_cal_v",version,"_",names(groupings)[g])
    if (save_plots) pdf(paste0(xname,".pdf"),width=3,height=3.5)
    labs=c("Overall",names(groupings[[g]]))
    xcol=c("black","red","blue")
    cal_2panel(list(cal_overall,calA,calB),labels=labs,col=xcol,text.col=xcol)
    if (save_plots)    dev.off()
    
    # Object to save in file to reproduce plot
    obj=list(cal_overall=cal_overall,calA=calA,calB=calB,labs=labs,xcol=xcol)
    obj_name=paste0("cal_v",version,"_",names(groupings)[g])
    
    # Save
    sink(out_text_file,append=TRUE)
    cat("\n\n")
    cat(obj_name,"="); dput(obj); 
    cat("\n\n")
    sink()
    
    
    
  }
}



##**********************************************#
## Main fairness metrics: DP, FPP, FORP etc  ####
##**********************************************#


for (g in 1:length(groupings)) {
  
  file_prefix=paste0(names(groupings)[g],"/")
  out_text_file=paste0(plot_dir,names(groupings)[g],"/summary_",names(groupings)[g],".txt")
  

  # Loop over v3 and v4
  for (version in c(3,4)) {
    
    # Loop over categories
    subcats=which(!(names(categories) %in% names(groupings[[g]])))
    for (cc in subcats) {
      
      # If version==3, exclude individuals who died during prediction period
      if (version==3){
        wsub=which((dat0$target==0)|dat0$reason %in% c("B","E"))
        dat=dat0[intersect(which(dat0$id %in% categories[[cc]]),wsub),]
        groupA=which(dat$id %in% groupings[[g]][[1]])
        groupB=which(dat$id %in% groupings[[g]][[2]])
        score=dat$v3 
      } else {
        dat=dat0[(which(dat0$id %in% categories[[cc]])),]
        groupA=which(dat$id %in% groupings[[g]][[1]])
        groupB=which(dat$id %in% groupings[[g]][[2]])
        score=dat$v4
      }
      
      target=dat$target
      

      ## Demographic parity                        ####

      ## Demographic parity: are the distributions of scores the same for men and women?
      #ks.test(score[male],score[female]) # no
      
      pname="dp"
      xname=paste0(plot_dir,file_prefix,"v",version,"_",pname,"_",names(categories)[cc])
      if (save_plots) pdf(paste0(xname,".pdf"),width=3,height=3.5)
      
      xx=seq(log(0.001),log(0.99),length=200)

      yA=0*xx; yB=0*xx; cA=0*xx; cB=0*xx
      nA=length(groupA); nB=length(groupB); xci=-qnorm(alpha/2)
      for (xi in 1:length(xx)) {
        xA=sum(score[groupA] <= exp(xx[xi]))
        xB=sum(score[groupB] <= exp(xx[xi]))
        # Privacy / disclosure control
        xA[which(xA < 5)]=0
        xB[which(xB < 5)]=0
        
        yA[xi]=xA/nA; yB[xi]=xB/nB
        cA[xi]=xci*sqrt(yA[xi]*(1-yA[xi])/nA); 
        cB[xi]=xci*sqrt(yB[xi]*(1-yB[xi])/nB)
      }
      cdif=sqrt(cA^2 + cB^2)

      ## Plot
      par(mar=c(1,4,0.1,0.1))
      layout(mat=rbind(matrix(1,4,4),matrix(2,2,4)))
      xrange=range(log(score)); if (!is.finite(sum(xrange))) xrange=c(0,1)
      
      plot(0,xlim=xrange,ylim=c(0,1),xaxt="n",ylab="Cumul. prob",xlab="")
      polygon(c(xx,rev(xx)),c(yA + cA,rev(yA-cA)),col=rgb(0,0,0,alpha=0.5),border=NA)
      polygon(c(xx,rev(xx)),c(yB + cB,rev(yB-cB)),col=rgb(1,0,0,alpha=0.5),border=NA)
      lines(xx,yA,col="black")
      lines(xx,yB,col="red")
      legend("bottomright",names(groupings[[g]]),col=c("black","red"),lty=1)

      par(mar=c(4,4,0.1,0.1))
      yrr=range(yA-yB,na.rm=T); if (!is.finite(sum(yrr))) yrr=c(0,1)
      plot(0,xlim=xrange,ylim=yrr,type="n",
           xlab="Score",ylab=expression(paste(Delta,"prob (black-red)")),xaxt="n")
      polygon(c(xx,rev(xx)),c(yA-yB+cdif,rev(yA-yB-cdif)),col=rgb(0.5,0.5,0.5,alpha=0.5),border=NA)
      lines(xx,yA-yB)
      abline(h=0,lty=2)
      pz=pretty(xrange,n=10)
      axis(1,at=pz,labels=round(100*exp(pz)))
      
      if (save_plots) dev.off()

      # Object to save in file to reproduce plot
      obj=list(xx=xx,yA=yA,yB=yB,cA=cA,cB=cB,xrange=xrange)
      obj_name=paste0(pname,"_v",version,"_",names(groupings)[g],"_",names(categories)[cc])

      # Save
      sink(out_text_file,append=TRUE)
      cat("\n\n")
      cat(obj_name,"="); dput(obj); 
      cat("\n\n")
      sink()
      
      
            
      
      
      
      
      
      ## Other main metrics                        ####
      
      for (s in 1:dim(specs)[1]) {
        
        p_AA=rep(NA,ncut); p_BB=rep(NA,ncut)
        n_AA=rep(NA,ncut); n_BB=rep(NA,ncut)
        for (i in 1:ncut) {
          cond_AA=1:dim(dat)[1]
          cond_BB=1:dim(dat)[1]
          if (is.finite(specs[s,4])) {
            if (specs[s,4]==1) yhatc=which(score>=cutoffs[i]) else yhatc=which(score < cutoffs[i]) 
            cond_AA=intersect(cond_AA,yhatc)
            cond_BB=intersect(cond_BB,yhatc)
          }
          if (is.finite(specs[s,5])) {
            yc=which(target==specs[s,5])
            cond_AA=intersect(cond_AA,yc)
            cond_BB=intersect(cond_BB,yc)
          }
          if (is.finite(specs[s,6])) {
            cond_AA=intersect(cond_AA,groupA)
            cond_BB=intersect(cond_BB,groupB)
          }
          
          prob_AA=cond_AA
          prob_BB=cond_BB
          if (is.finite(specs[s,1])) {
            if (specs[s,1]==1) yhatc=which(score>=cutoffs[i]) else yhatc=which(score < cutoffs[i]) 
            prob_AA=intersect(prob_AA,yhatc)
            prob_BB=intersect(prob_BB,yhatc)
          }
          if (is.finite(specs[s,2])) {
            yc=which(target==specs[s,2])
            prob_AA=intersect(prob_AA,yc)
            prob_BB=intersect(prob_BB,yc)
          }
          if (is.finite(specs[s,3])) {
            prob_AA=intersect(prob_AA,groupA)
            prob_BB=intersect(prob_BB,groupB)
          }
          
          ## Privacy/DISCLOSURE CONTROL: important
          if (length(prob_AA)<5) prob_AA=c()
          if (length(prob_BB)<5) prob_BB=c()
          
          p_AA[i]=length(prob_AA)/length(cond_AA)
          p_BB[i]=length(prob_BB)/length(cond_BB)
          n_AA[i]=length(cond_AA)
          n_BB[i]=length(cond_BB)
        }
        
        ci_AA=-qnorm(alpha/2)*sqrt(p_AA*(1-p_AA)/n_AA)  
        ci_BB=-qnorm(alpha/2)*sqrt(p_BB*(1-p_BB)/n_BB)  
        
        
        pname=rownames(specs)[s]
        xname=paste0(plot_dir,file_prefix,"v",version,"_",pname,"_",names(categories)[cc])
        
        if (save_plots) pdf(paste0(xname,".pdf"),width=3,height=3.5)
        
        ## Plot
        par(mar=c(1,4,0.1,0.1))
        layout(mat=rbind(matrix(1,4,4),matrix(2,2,4)))
        
        yrr=range(c(p_AA,p_BB),na.rm=T); if (!is.finite(sum(yrr))) yrr=c(-1,1)
        plot(0,xlim=c(0,1),ylim=yrr,type="n",xlab="",xaxt="n",ylab="Prob.")
        lines(cutoffs,p_AA,col="black")
        lines(cutoffs,p_BB,col="red")
        polygon(c(cutoffs,rev(cutoffs)),c(p_AA + ci_AA,rev(p_AA- ci_AA)),col=rgb(0.5,0.5,0.5,alpha=0.5),border=NA)
        polygon(c(cutoffs,rev(cutoffs)),c(p_BB + ci_BB,rev(p_BB- ci_BB)),col=rgb(1,0,0,alpha=0.5),border=NA)
        legend("topright",names(groupings[[g]]),col=c("black","red"),lty=1)
        polygon(c(0.4,0.6,0.6,0.4),c(-10,-10,10,10),border=NA,col=rgb(1,0,0,alpha=0.2))
        
        ci_dif=sqrt(ci_AA^2 + ci_BB^2)
        par(mar=c(4,4,0.1,0.1))
        yrr=range(p_AA-p_BB,na.rm=T); if (!is.finite(sum(yrr))) yrr=c(-1,1)
        plot(0,xlim=c(0,1),ylim=yrr,type="n",
             xlab="Cutoffs",ylab=expression(paste(Delta,"prob (black-red)")),xaxt="n")
        polygon(c(cutoffs,rev(cutoffs)),c(p_AA-p_BB + ci_dif,rev(p_AA-p_BB- ci_dif)),col=rgb(0.5,0.5,0.5,alpha=0.5),border=NA)
        lines(cutoffs,p_AA-p_BB)
        abline(h=0,lty=2)
        pz=pretty(range(cutoffs),n=10)
        axis(1,at=pz,labels=pz)
        polygon(c(0.4,0.6,0.6,0.4),c(-10,-10,10,10),border=NA,col=rgb(1,0,0,alpha=0.2))
        
        if (save_plots) dev.off()

        # Object to save in file to reproduce plot
        obj=list(cutoffs=cutoffs,p_AA=p_AA,p_BB=p_BB,ci_AA=ci_AA,ci_BB=ci_BB)
        obj_name=paste0(pname,"_v",version,"_",names(groupings)[g],"_",names(categories)[cc])
        
        # Save
        sink(out_text_file,append=TRUE)
        cat("\n\n")
        cat(obj_name,"="); dput(obj); 
        cat("\n\n")
        sink()
        
      }
      
    }
    
  }
}



  
##**********************************************#
## Adjusted FORP/FDRP                        ####
##**********************************************#


# Number of bootstrap resamples
nboot=100

for (version in c(3,4)) {
  
  # If version==3, exclude individuals who died during prediction period
  if (version==3){
    wsub=which((dat0$target==0)|dat0$reason %in% c("B","E"))
    dat=dat0[wsub,]
    score=dat$v3 
  } else {
    dat=dat0
    score=dat$v4
  }
  target=dat$target
  
  for (g in 1:length(groupings)) {
    
    # Where to save, or write context to single file 
    file_prefix=paste0(names(groupings)[g],"/")
    out_text_file=paste0(plot_dir,names(groupings)[g],"/summary_",names(groupings)[g],".txt")
    
    
    # Define groups (male/female, old/young etc)
    groupA=dat$id[which(dat$id %in% groupings[[g]][[1]])]
    groupB=dat$id[which(dat$id %in% groupings[[g]][[2]])]
    
    # Adjust for these categories...
    adj_cats=c("sexM","age","simd")
    
    #... unless they are part of the group definition
    if (names(groupings)[g]=="Sex") adj_cats=setdiff(adj_cats,"sexM")
    if (names(groupings)[g]=="SIMD") adj_cats=setdiff(adj_cats,"simd")
    if (names(groupings)[g]=="Age") adj_cats=setdiff(adj_cats,"age")
    
    if (version==3) {
      dsub=dat[,c(adj_cats,"id","v3","target")]
      colnames(dsub)=c(adj_cats,"id","score","target")
    }
    if (version==4) {
      dsub=dat[,c(adj_cats,"id","v4","target")]
      colnames(dsub)=c(adj_cats,"id","score","target")
    }
    comm=paste0("dsub=dsub[with(dsub,order(",paste(adj_cats,collapse=","),",score)),]")
    eval(parse(text=comm))
    xcat=rep("",dim(dsub)[1]); for (i in 1:length(adj_cats)) xcat=paste0(xcat,dsub[,adj_cats[i]],"_")
    dsub$cat=xcat
    
    dt=table(dsub$cat)[unique(dsub$cat)]
    cdt=c(0,cumsum(dt))
    dsub$group=rep(NA,dim(dsub)[1])
    dsub$group[which(dsub$id %in% groupA)]="A"
    dsub$group[which(dsub$id %in% groupB)]="B"
    
    
    
    # Matrix of counts of A=a|Yhat<c and A=a|Yhat>= c for cutoffs c
    ucat=unique(dsub$cat)
    dt_forp=matrix(0,length(ucat),length(cutoffs))
    dt_fdrp=matrix(0,length(ucat),length(cutoffs))
    for (i in 1:length(cutoffs)) {
      dt_forp[,i]=table(dsub$cat[which(dsub$score<cutoffs[i])])[ucat]
      dt_fdrp[,i]=table(dsub$cat[which(dsub$score>=cutoffs[i])])[ucat]
    }
    dt_forp[which(is.na(dt_forp))]=0
    dt_fdrp[which(is.na(dt_fdrp))]=0
    
    
    # Arrays of counts of A=a|Yhat<c and A=a|Yhat>= c for cutoffs c for bootstrap samples
    b_dt_forp=array(0,dim=c(length(ucat),length(cutoffs),nboot))
    b_dt_fdrp=array(0,dim=c(length(ucat),length(cutoffs),nboot))
    b_n=array(0,dim=c(length(ucat),nboot))
    b_A=array(0,dim=c(length(ucat),nboot))
    b_B=array(0,dim=c(length(ucat),nboot))
    for (j in 1:nboot) {
      dsub_b=dsub[sample(dim(dsub)[1],rep=T),]
      b_n[,j]=table(dsub_b$cat)[ucat]
      b_A[,j]=table(dsub_b$cat[which(dsub_b$group=="A")])[ucat]
      b_B[,j]=table(dsub_b$cat[which(dsub_b$group=="B")])[ucat]
      for (i in 1:length(cutoffs)) {
        b_dt_forp[,i,j]=table(dsub_b$cat[which(dsub_b$score<cutoffs[i])])[ucat]
        b_dt_fdrp[,i,j]=b_n[,j]-b_dt_forp[,i,j]
      }
      print(j)
    }
    b_n[which(is.na(b_n))]=0
    b_A[which(is.na(b_A))]=0
    b_B[which(is.na(b_B))]=0
    b_dt_forp[which(is.na(b_dt_forp))]=0
    b_dt_fdrp[which(is.na(b_dt_fdrp))]=0
    
    
    # FORP/FDRP for each category
    forp_A=matrix(NA,length(dt),length(cutoffs))
    forp_B=matrix(NA,length(dt),length(cutoffs))
    fdrp_A=matrix(NA,length(dt),length(cutoffs))
    fdrp_B=matrix(NA,length(dt),length(cutoffs))
    
    # Denominators for each category
    nforp_A=matrix(0,length(dt),length(cutoffs))
    nforp_B=matrix(0,length(dt),length(cutoffs))
    nfdrp_A=matrix(0,length(dt),length(cutoffs))
    nfdrp_B=matrix(0,length(dt),length(cutoffs))
    
    b_forp_A=array(NA,dim=c(length(dt),length(cutoffs),nboot))
    b_fdrp_A=array(NA,dim=c(length(dt),length(cutoffs),nboot))
    b_forp_B=array(NA,dim=c(length(dt),length(cutoffs),nboot))
    b_fdrp_B=array(NA,dim=c(length(dt),length(cutoffs),nboot))
    
    
    for (i in 1:length(dt)) {
      if ((i %% 100)==0) print(i)
      sub=dsub[(1+cdt[i]):(cdt[i+1]),]
      
      subA=sub[which(sub$group=="A"),]
      if (length(unique(subA$score))>2) {
        nA1=round(ecdf(subA$score)(cutoffs)*length(subA$score))
        nA2=round((1-ecdf(subA$score)(cutoffs))*length(subA$score))
        w1=which(nA1>0); w2=which(nA2>0)
        forp_A[i,w1]=cumsum(subA$target)[nA1[w1]]/nA1[w1]
        fdrp_A[i,w2]=cumsum(rev(1-subA$target))[nA2[w2]]/nA2[w2]
      }
      for (j in 1:nboot) {
        subAc=subA[sort(sample(dim(subA)[1],b_A[i,j],rep=T)),]
        if (length(unique(subAc$score))>2) {
          nA1=round(ecdf(subAc$score)(cutoffs)*length(subAc$score))
          nA2=round((1-ecdf(subAc$score)(cutoffs))*length(subAc$score))
          w1=which(nA1>0); w2=which(nA2>0)
          b_forp_A[i,w1,j]=cumsum(subAc$target)[nA1[w1]]/nA1[w1]
          b_fdrp_A[i,w2,j]=cumsum(rev(1-subAc$target))[nA2[w2]]/nA2[w2]
        }
      }
      
      subB=sub[which(sub$group=="B"),]
      if (length(unique(subB$score))>2) {
        nB1=round(ecdf(subB$score)(cutoffs)*length(subB$score))
        nB2=round((1-ecdf(subB$score)(cutoffs))*length(subB$score))
        w1=which(nB1>0); w2=which(nB2>0)
        forp_B[i,w1]=cumsum(subB$target)[nB1[w1]]/nB1[w1]
        fdrp_B[i,w2]=cumsum(rev(1-subB$target))[nB2[w2]]/nB2[w2]
      }
      for (j in 1:nboot) {
        subBc=subB[sort(sample(dim(subB)[1],b_B[i,j],rep=T)),]
        if (length(unique(subBc$score))>2) {
          nB1=round(ecdf(subBc$score)(cutoffs)*length(subBc$score))
          nB2=round((1-ecdf(subBc$score)(cutoffs))*length(subBc$score))
          w1=which(nB1>0); w2=which(nB2>0)
          b_forp_B[i,w1,j]=cumsum(subBc$target)[nB1[w1]]/nB1[w1]
          b_fdrp_B[i,w2,j]=cumsum(rev(1-subBc$target))[nB2[w2]]/nB2[w2]
        }
      }
    }
    
    FORP_A1=rep(0,length(cutoffs))
    FORP_B1=rep(0,length(cutoffs))
    FDRP_A1=rep(0,length(cutoffs))
    FDRP_B1=rep(0,length(cutoffs))
    b_FORP_A1=matrix(0,length(cutoffs),nboot)
    b_FDRP_A1=matrix(0,length(cutoffs),nboot)
    b_FORP_B1=matrix(0,length(cutoffs),nboot)
    b_FDRP_B1=matrix(0,length(cutoffs),nboot)
    
    for (i in 1:length(cutoffs)) {
      FORP_A1[i]=sum(forp_A[,i]*dt_forp[,i],na.rm=T) 
      FORP_B1[i]=sum(forp_B[,i]*dt_forp[,i],na.rm=T) 
      FDRP_A1[i]=sum(fdrp_A[,i]*dt_fdrp[,i],na.rm=T) 
      FDRP_B1[i]=sum(fdrp_B[,i]*dt_fdrp[,i],na.rm=T) 
      for (j in 1:nboot) {
        b_FORP_A1[i,j]=sum(b_forp_A[,i,j]*b_dt_forp[,i,j],na.rm=T) 
        b_FORP_B1[i,j]=sum(b_forp_B[,i,j]*b_dt_forp[,i,j],na.rm=T) 
        b_FDRP_A1[i,j]=sum(b_fdrp_A[,i,j]*b_dt_fdrp[,i,j],na.rm=T) 
        b_FDRP_B1[i,j]=sum(b_fdrp_B[,i,j]*b_dt_fdrp[,i,j],na.rm=T) 
      }
    }
    
    # Privacy/disclosure control
    FORP_A1[which(FORP_A1<5)]=0
    FORP_B1[which(FORP_B1<5)]=0
    FDRP_A1[which(FDRP_A1<5)]=0
    FDRP_B1[which(FDRP_B1<5)]=0
    
    # Normalise
    FORP_A=rep(0,length(cutoffs))
    FORP_B=rep(0,length(cutoffs))
    FDRP_A=rep(0,length(cutoffs))
    FDRP_B=rep(0,length(cutoffs))
    b_FORP_A=matrix(0,length(cutoffs),nboot)
    b_FORP_B=matrix(0,length(cutoffs),nboot)
    b_FDRP_A=matrix(0,length(cutoffs),nboot)
    b_FDRP_B=matrix(0,length(cutoffs),nboot)
    for (i in 1:length(cutoffs)) {
      FORP_A[i]=FORP_A1[i]/sum(dt_forp[which(is.finite(forp_A[,i])),i])
      FORP_B[i]=FORP_B1[i]/sum(dt_forp[which(is.finite(forp_B[,i])),i])
      FDRP_A[i]=FDRP_A1[i]/sum(dt_fdrp[which(is.finite(fdrp_A[,i])),i])
      FDRP_B[i]=FDRP_B1[i]/sum(dt_fdrp[which(is.finite(fdrp_B[,i])),i])
      for (j in 1:nboot) {
        b_FORP_A[i,j]=b_FORP_A1[i,j]/sum(b_dt_forp[which(is.finite(b_forp_A[,i,j])),i,j])
        b_FORP_B[i,j]=b_FORP_B1[i,j]/sum(b_dt_forp[which(is.finite(b_forp_B[,i,j])),i,j])
        b_FDRP_A[i,j]=b_FDRP_A1[i,j]/sum(b_dt_fdrp[which(is.finite(b_fdrp_A[,i,j])),i,j])
        b_FDRP_B[i,j]=b_FDRP_B1[i,j]/sum(b_dt_fdrp[which(is.finite(b_fdrp_B[,i,j])),i,j])
      }
    }
    
    
    # Standard errors and confidence intervals are approximate, and take weights as fixed
    qx=function(x,alpha) if (any(is.finite(x))) quantile(x,alpha,na.rm=T) else NA
    ci_FORP_A=apply(b_FORP_A,1,function(x) qx(x,1-alpha)-qx(x,alpha))/2
    ci_FORP_B=apply(b_FORP_B,1,function(x) qx(x,1-alpha)-qx(x,alpha))/2
    ci_FDRP_A=apply(b_FDRP_A,1,function(x) qx(x,1-alpha)-qx(x,alpha))/2
    ci_FDRP_B=apply(b_FDRP_B,1,function(x) qx(x,1-alpha)-qx(x,alpha))/2
    
    for (ptype in 1:2) {
      
      if (ptype==1) {  
        xname=paste0(plot_dir,names(groupings)[g],"/v",version,"_FORP_adjusted_",names(groupings)[g])
        p_AA=FORP_A; p_BB=FORP_B
        ci_AA=ci_FORP_A; ci_BB=ci_FORP_B
        obj_name=paste0("forp_adjusted_v",version,"_",names(groupings)[g])
      } else {
        xname=paste0(plot_dir,names(groupings)[g],"/v",version,"_FDRP_adjusted_",names(groupings)[g])
        p_AA=FDRP_A; p_BB=FDRP_B
        ci_AA=ci_FDRP_A; ci_BB=ci_FDRP_B
        obj_name=paste0("fdrp_adjusted_v",version,"_",names(groupings)[g])
      }
      
      
      
      
      if (save_plots) pdf(paste0(xname,".pdf"),width=3,height=3.5)
      
      
      ## Plot
      par(mar=c(1,4,0.1,0.1))
      layout(mat=rbind(matrix(1,4,4),matrix(2,2,4)))
      
      yrr=range(c(p_AA,p_BB),na.rm=T); if (!is.finite(sum(yrr))) yrr=c(-1,1)
      plot(0,xlim=c(0,1),ylim=yrr,type="n",xlab="",xaxt="n",ylab="Prob.")
      lines(cutoffs,p_AA,col="black")
      lines(cutoffs,p_BB,col="red")
      polygon(c(cutoffs,rev(cutoffs)),c(p_AA + ci_AA,rev(p_AA- ci_AA)),col=rgb(0.5,0.5,0.5,alpha=0.5),border=NA)
      polygon(c(cutoffs,rev(cutoffs)),c(p_BB + ci_BB,rev(p_BB- ci_BB)),col=rgb(1,0,0,alpha=0.5),border=NA)
      legend("topright",names(groupings[[g]]),col=c("black","red"),lty=1)
      polygon(c(0.4,0.6,0.6,0.4),c(-10,-10,10,10),border=NA,col=rgb(1,0,0,alpha=0.2))
      
      ci_dif=sqrt(ci_AA^2 + ci_BB^2)
      par(mar=c(4,4,0.1,0.1))
      yrr=range(p_AA-p_BB,na.rm=T); if (!is.finite(sum(yrr))) yrr=c(-1,1)
      plot(0,xlim=c(0,1),ylim=yrr,type="n",
           xlab="Cutoffs",ylab=expression(paste(Delta,"prob (black-red)")),xaxt="n")
      polygon(c(cutoffs,rev(cutoffs)),c(p_AA-p_BB + ci_dif,rev(p_AA-p_BB- ci_dif)),col=rgb(0.5,0.5,0.5,alpha=0.5),border=NA)
      lines(cutoffs,p_AA-p_BB)
      abline(h=0,lty=2)
      pz=pretty(range(cutoffs),n=10)
      axis(1,at=pz,labels=pz)
      polygon(c(0.4,0.6,0.6,0.4),c(-10,-10,10,10),border=NA,col=rgb(1,0,0,alpha=0.2))
      
      if (save_plots) dev.off()
      
      # Object to save in file to reproduce plot
      obj=list(cutoffs=cutoffs,p_AA=p_AA,p_BB=p_BB,ci_AA=ci_AA,ci_BB=ci_BB)
      
      # Save
      sink(out_text_file,append=TRUE)
      cat("\n\n")
      cat(obj_name,"="); dput(obj); 
      cat("\n\n")
      sink()
      
      
    }
    
  }
  
}









##**********************************************#
## Counterfactual fairness                   ####
##**********************************************#



##' counterfactual_yhat
##' 
##' Estimation of counterfactual quantities by resampling.
##' 
##' Counterfactual fairness is with respect to the causal graph:
##'
##'
##'   G <----------U
##'   |         /   |
##'   |     /       |
##'   | /           |
##'   VV            V
##'   X -------->  Yhat
##'   
##' where
##'
##' G=group (usually sensitive attribute);  Yhat=outcome; X=set of variables 
##'  through which G can act on Yhat, U=set of background variables;
##'
##' We want the counterfactual Yhat_{g' <- G} |X=x,G=g (or alternatively 
##'  Yhat_{g' <- G} |G=g).
##' 
##' This can be interpreted as:
##'  `the distribution of values of Yhat amongst indivdiduals whose values of U 
##'   are distributed as though they were in group G=g (and, optionally, had 
##'   values X=x, but whose value of G is g'`
##' 
##' Essentially, comparison of the counterfactual quantity above to the 
##'  conditional Yhat|G=g isolates the difference in Yhat due to the effect of 
##'  G on Yhat through X, removing any effect due to different distributions of 
##'  U due to different values of G.
##' 
##' To estimate Y'=Yhat_{g' <- G} | G=g, we need to
##'  1. Compute U'~(U|G=g)
##'  2. Compute the distribution X' as X'~(X|U~U',G=g')
##'  3. Sample Y'~(Yhat|X~X',U~U')
##' 
##' To estimate Y'=Yhat_{g' <- G} |X=x, G=g, we need to
##'  1. Compute U'~(U|G=g,X=x)
##'  2. Compute the distribution X' as X'~(X|U~U',G=g')
##'  3. Sample Y'~(Yhat|X~X',U~U')
##'  
##' This function approximates this samplying procedure as follows
##'  1. Look at individuals with G=g (and optionally X=x)
##'  2. Find the values of U for these individuals
##'  3. Find a second set of individuals with the same values of U but for whom G=g'
##'  4. Return the indices of these individuals
##'  
##' The values of Yhat for these individuals constitute a sample from the 
##'  desired counterfactual. 
##'  
##' @param dat data frame containing variables in X, Yhat, G, excl. Variables in U are assumed to be colnames(dat)\(X U Yhat U G U excl)
##' @param X set of variables and values on which to 'condition'; essentially we assume that any causal pathway from G to Yhat is through X
##' @param x values of variables X on which to condition; can be null in which case we use marginal distribution of X
##' @param G grouping variable, usually a sensitive attribute
##' @param g conditioned value of g
##' @param gdash counterfactual value of g
##' @param excl variable names to exclude from U
##' @param n number of samples; if NULL return all
##' @export 
##' @return indices representing sample(s) from counterfactual  Yhat_{g' <- G} |X=x,G=g
counterfactual_yhat=function(dat,X,x=NULL,G,g,gdash,excl=NULL,n=NULL) {
  
  U=setdiff(colnames(dat),c(X,G,excl))
  dat$id=1:dim(dat)[1]
  odat=dat[do.call(order,dat[,U]),]
  
  datXG=odat[,c("id",X,G)]
  w=1:dim(dat)[1]
  if (!is.null(x)) for (i in 1:length(x)) w=intersect(w,which(datXG[,X]==x[i]))
  w=intersect(w,which(datXG[,dim(datXG)[2]]==g))
  
  w2=which(datXG[[G]]==gdash)
  
  ax=approx(w2,1:length(w2),w,rule=2)$y
  wz1=w2[floor(ax)]
  wz2=w2[ceiling(ax)]
  
  w12=which(abs(wz2-w)<abs(wz1-w))
  wz=wz1; wz[w12]=wz2[w12]
  
  if (!is.null(n)) out=sample(datXG$id[wz],n) else out=datXG$id[wz]
  
  return(out)
}


# Read full dataset
fullmatrix=paste0("James/Analysis/Data/full_model/full_data_matrix.RDS")
nhs=readRDS(fullmatrix)
dat=nhs[which((nhs$time==dat0$time[1]) & (nhs$id %in% dat0$id)),]
rm(list=c("nhs"))

# Basic processing
dat=dat[match(dat0$id,dat$id),]
dat$v3=dat0$v3
dat$v4=dat0$v4
names(dat)[which(names(dat)=="SIMD_DECILE_2016_SCT")]="simd"


# General variables not to go in U
notU=c("id","v3"," v4","time","target","reason","cv","age","sexM","simd","v3score",grep("topic",colnames(dat),val=T))




for (g in 1:length(groupings)) {
  print(names(groupings)[g])
  
  # Define groups (male/female, old/young etc)
  W=rep(NA,dim(dat)[1])
  W[which(dat$id %in% groupings[[g]][[1]])]=1
  W[which(dat$id %in% groupings[[g]][[2]])]=2
  dat$W=W
  
  X0=c("age","sexM","simd")
  if (names(groupings)[g]=="Sex") X0=setdiff(X0,"sexM")
  if (names(groupings)[g]=="Age") X0=setdiff(X0,"age")
  if (names(groupings)[g]=="SIMD") X0=setdiff(X0,"simd")
  
  Yc2=dat$id[counterfactual_yhat(dat,
                          X=X0,
                          x=NULL,
                          G="W",
                          g=1,
                          gdash=2,
                          excl=notU,
                          n=NULL)]
  Yc1=dat$id[which(dat$W==1)]
  
  
  
  d1=match(Yc1,dat$id)
  d2=match(Yc2,dat$id)
  
  for (version in c(3,4)) {
    
    # Where to save, or write context to single file 
    file_prefix=paste0(names(groupings)[g],"/")
    out_text_file=paste0(plot_dir,names(groupings)[g],"/summary_",names(groupings)[g],".txt")
    

    xname=paste0(plot_dir,file_prefix,"v",version,"_counterfactual_",names(groupings)[g])
    if (save_plots) pdf(paste0(xname,".pdf"),width=3,height=3.5)
    
    if (version==3) {
      yhat1=dat$v3[d1]
      yhat2=dat$v3[d2]
    }
    if (version==4) {
      yhat1=dat$v4[d1]
      yhat2=dat$v4[d2]
    }
    
    xx=seq(log(0.01),log(0.99),length=200)
    
    yA=0*xx; yB=0*xx; cA=0*xx; cB=0*xx
    nA=length(yhat1); nB=length(yhat2); xci=-qnorm(alpha/2)
    for (xi in 1:length(xx)) {
      xA=sum(yhat1 <= exp(xx[xi]))
      xB=sum(yhat2 <= exp(xx[xi]))
      
      # Privacy/disclosure control
      xA[which(xA<5)]=0
      xB[which(xB<5)]=0
      
      yA[xi]=xA/nA; yB[xi]=xB/nB
      cA[xi]=xci*sqrt(yA[xi]*(1-yA[xi])/nA); 
      cB[xi]=xci*sqrt(yB[xi]*(1-yB[xi])/nB)
    }
    xrange=range(xx)
    if (!is.finite(sum(xrange))) xrange=c(-1,1)
    
  
    ## Plot
    par(mar=c(1,4,0.1,0.1))
    layout(mat=rbind(matrix(1,4,4),matrix(2,2,4)))
    
    plot(0,xlim=xrange,ylim=c(0,1),xaxt="n",ylab="Cumul. prob",xlab="")
    polygon(c(xx,rev(xx)),c(yA + cA,rev(yA-cA)),col=rgb(0,0,0,alpha=0.5),border=NA)
    polygon(c(xx,rev(xx)),c(yB + cB,rev(yB-cB)),col=rgb(1,0,0,alpha=0.5),border=NA)
    lines(xx,yA,col="black")
    lines(xx,yB,col="red")
    legend("bottomright",names(groupings[[g]]),col=c("black","red"),lty=1)
    
    par(mar=c(4,4,0.1,0.1))
    cdif=sqrt(cA^2 + cB^2)
    yrr=range(yA-yB,na.rm=T); if (!is.finite(sum(yrr))) yrr=c(-1,1)
    plot(0,xlim=xrange,ylim=yrr,type="n",
         xlab="Score",ylab=expression(paste(Delta,"prob (black-red)")),xaxt="n")
    polygon(c(xx,rev(xx)),c(yA-yB+cdif,rev(yA-yB-cdif)),col=rgb(0.5,0.5,0.5,alpha=0.5),border=NA)
    lines(xx,yA-yB)
    abline(h=0,lty=2)
    pz=pretty(xrange,n=10)
    axis(1,at=pz,labels=round(100*exp(pz)))
    
    if (save_plots) dev.off()
    
    # Object to save in file to reproduce plot
    obj=list(xx=xx,yA=yA,yB=yB,cA=cA,cB=cB,xrange=xrange)
    obj_name=paste0("counterfactual_dp_v",version,"_",names(groupings)[g])
    
    # Save
    sink(out_text_file,append=TRUE)
    cat("\n\n")
    cat(obj_name,"="); dput(obj); 
    cat("\n\n")
    sink()
    
    
    
  }
}




##**********************************************#
## FOR decomposition                         ####
##**********************************************#

# Read deaths table
death_tab=haven:::read_spss("../../../1718-0370/Linked Data/Final_deaths_extract_incl_UniqueStudyID_v2.zsav")

# Get icd10 codes, including those for death
icd10=dat0$main_condition
w=which(dat0$reason=="D")
icd10[w]=death_tab$PRIMARY_CAUSE_OF_DEATH[match(dat0$id[w],death_tab$UNIQUE_STUDY_ID)]

# Convert ICD10 codes to causes of admission
icd2=as.factor(substring(icd10,1,3))
causes=icd2sick(levels(icd2))
levels(icd2)=causes
icd2=as.character(icd2)

# Differentiate deaths from other admissions
wd = which(dat0$reason=="D")
icd2[wd]=paste0("Died of ",icd2[wd])

# Include NAs
icd2[which(dat0$target==1 & is.na(icd2))]="Not recorded"
icd2[which(icd2=="Died of NA")]="Died of unrecorded"

# Convert back to factor for easy tabulation
icd2=as.factor(icd2)

# Construct a giant table of the distribution of admission causes in various cohorts
dtype=data.frame()
dcolnames=levels(icd2)
dcolnames=c(sort(grep("Died",dcolnames,inv=T,val=T)),
            sort(grep("Died",dcolnames,val=T)))
dnames=c()

# Score cutoffs
score_cutoffs=seq(0,1,length=21)

for (i in 1:length(categories)) {
  for (scut in 1:(length(score_cutoffs)-1)) {
    w3=which(dat0$id %in% categories[[i]] & 
               dat0$v3>= score_cutoffs[scut] & 
               dat0$v3 < score_cutoffs[scut+1])
    w4=which(dat0$id %in% categories[[i]] & 
               dat0$v4>= score_cutoffs[scut] & 
               dat0$v4 < score_cutoffs[scut+1])
    dtype=rbind(dtype,table(icd2[w3])[dcolnames])
    dtype=rbind(dtype,table(icd2[w4])[dcolnames])
    dnames=c(dnames,
             paste0("v3_",names(categories)[i],"_q",scut),
             paste0("v4_",names(categories)[i],"_q",scut))
    # print(c(i,scut))
  }
}
rownames(dtype)=dnames
colnames(dtype)=dcolnames


# Disclosure control
for (i in 1:dim(dtype)[2]) {
  sub=dtype[,i]
  sub[which(sub<5)]=0
  dtype[,i]=sub
}

# Write to file
write.table(dtype,file=paste0(plot_dir,"decomp_table.csv"),row.names=TRUE,col.names=TRUE,sep=",")









