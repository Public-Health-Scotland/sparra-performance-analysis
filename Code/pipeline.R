##**********************************************#
## Analyse SPARRA fairness measures          ####
##**********************************************#
## James Liley, 2022-23
##

# Run this after setting working directory to the home directory of the repository.


##**********************************************#
## Directories, packages and scripts         ####
##**********************************************#


# Scripts
# Data - includes plotting data and explanations
sfiles=list.files("./Data/raw_data/",pattern="summary_*",full.names=TRUE)
for (s in sfiles) source(s) # Main analysis data
decomp=read.csv("./Data/raw_data/decomp_table.csv") # Table of decomposition data
otab0=read.csv("./Data/raw_data/overview_table.csv") # Table of demographic data


# Directories
out_dir="./Figures/All/" # Main directory to which to write figure data
table_dir="./Tables/"

# Packages
library(SPARRAfairness)

# Options
save_plots=TRUE # Save plots to PDF
highlight_value=0.1 # Highlight this cutoff (none)




##**********************************************#
## Specifications                            ####
##**********************************************#


groupings=list(
  Sex=list(Male=NULL,Female=NULL),
  SIMD=list(Most_deprived=NULL,Least_deprived = NULL),
  Age=list(O65=NULL, U25 = NULL),
  Ethnicity = list(White=NULL,Nonwhite=NULL),
  Mainland_island = list(Mainland=NULL,Island=NULL),
  Urban_rural = list(Urban=NULL,Rural=NULL)
)

plot_details=list(
  "dp"                = list(lpos=c(1,0),     yrange_upper=c(0,1),    yrange_lower=c(-0.3,0.3)),
  "counterfactual_dp" = list(lpos=c(1,0),     yrange_upper=c(0,1),    yrange_lower=c(-0.3,0.3)),
  "fdrp"              = list(lpos=c(1,0.7),   yrange_upper=c(0,1),    yrange_lower=c(-0.3,0.3)),
  "fdrp_adjusted"     = list(lpos=c(1,0.7),   yrange_upper=c(0,1),    yrange_lower=c(-0.3,0.3)),
  "fnp"               = list(lpos=c(1,0),     yrange_upper=c(0,0.2),  yrange_lower=c(-0.1,0.1)),
  "fnrp"              = list(lpos=c(1,0),     yrange_upper=c(0,1),    yrange_lower=c(-0.2,0.2)),
  "forp"              = list(lpos=c(0.3,0.7), yrange_upper=c(0,0.15), yrange_lower=c(-0.15,0.15)),
  "forp_adjusted"     = list(lpos=c(0.3,0.7), yrange_upper=c(0,0.15), yrange_lower=c(-0.15,0.15)),
  "fpp"               = list(lpos=c(1,0.7),   yrange_upper=c(0,1),    yrange_lower=c(-0.3,0.3)),
  "fprp"              = list(lpos=c(1,0.7),   yrange_upper=c(0,1),    yrange_lower=c(-0.5,0.5)),
  "irp"               = list(lpos=c(1,0),     yrange_upper=c(0,1),    yrange_lower=c(-0.3,0.3)),
  "rp"                = list(lpos=c(1,0.7),   yrange_upper=c(0,1),    yrange_lower=c(-0.5,0.5)),
  "roc"               = list(lpos=c(1,0),     yrange_upper=c(0,1),    yrange_lower=c(-0.2,0.2)),
  "prc"               = list(lpos=c(1,0.7),   yrange_upper=c(0,1),    yrange_lower=c(-0.3,0.3)),
  "cal"               = list(lpos=c(1,0),     yrange_upper=c(0,1),    yrange_lower=c(-0.1,0.1))
)


# Lower Yrange should be dependent on data?
for (i in 1:length(plot_details)) plot_details[[i]]$yrange_lower=NULL

# Group fairness metrics
metrics=c("fpp", "fdrp", "fprp", "rp", "irp", "fnp", "forp", "fnrp")

# List of categories
categories=list(all=c(0))
cnames=c("all")
ind=2
for (i in 1:length(groupings)) {
  for (ig in 1:length(groupings[[i]])) {
    categories[[ind]]=c(0)
    cnames[ind]=names(groupings[[i]])[ig]
    ind=ind+1
  }
}
names(categories)=cnames

print("Completed setup")

##**********************************************#
## Performance (ROC, PRC, CAL) by group      ####
##**********************************************#

for (version in c(3,4)) {

  for (g in 1:length(groupings)) {

    # ROC ####
    obj_name=paste0("roc_v",version,"_",names(groupings)[g])
    xname=paste0(out_dir,obj_name,".pdf")

    pdf(xname,width=3,height=3.5)

    obj=get(obj_name)
    with(obj,{
      pname="roc"
      aucT=signif(roc_overall$auc,digits=2); seT=signif(roc_overall$se,digits=2)
      aucA=signif(rocA$auc,digits=2); seA=signif(rocA$se,digits=2);
      aucB=signif(rocB$auc,digits=2); seB=signif(rocA$se,digits=2);
      labs=c(paste0("Overall: ",aucT," (",seT,")"),
             paste0(names(groupings[[g]])[1],": ",aucA," (",seA,")"),
             paste0(names(groupings[[g]])[2],": ",aucB," (",seB,")"))
      xcol=phs_colours(c("phs-blue","phs-magenta","phs-purple"))
      roc_2panel(list(roc_overall,rocA,rocB),
                 labels = labs,
                 col=xcol,
                 highlight=highlight_value,
                 #yrange=plot_details[[input$pname]]$yrange,
                 #lpos=plot_details[[input$pname]]$lpos,
                 yrange_lower=plot_details[[pname]]$yrange_lower,
                 legend_title="AUROC (SE)")
    })
    dev.off()



    # PR curves ####
    obj_name=paste0("prc_v",version,"_",names(groupings)[g])
    xname=paste0(out_dir,obj_name,".pdf")

    pdf(xname,width=3,height=3.5)

    obj=get(obj_name)
    with(obj,{
      pname="prc"
      aucT=signif(prc_overall$auc,digits=2); seT=signif(prc_overall$se,digits=2)
      aucA=signif(prcA$auc,digits=2); seA=signif(prcA$se,digits=2);
      aucB=signif(prcB$auc,digits=2); seB=signif(prcA$se,digits=2);
      labs=c(paste0("Overall: ",aucT," (",seT,")"),
             paste0(names(groupings[[g]])[1],": ",aucA," (",seA,")"),
             paste0(names(groupings[[g]])[2],": ",aucB," (",seB,")"))
      xcol=phs_colours(c("phs-blue","phs-magenta","phs-purple"))
      prc_2panel(list(prc_overall,prcA,prcB),
                 labels = labs,
                 col=xcol,
                 highlight=highlight_value,
                 #yrange=plot_details[[input$pname]]$yrange,
                 #lpos=plot_details[[input$pname]]$lpos,
                 yrange_lower=plot_details[[pname]]$yrange_lower,
                 legend_title="AUPRC (SE)")
    })
    dev.off()


    # Calibration curves ####
    obj_name=paste0("cal_v",version,"_",names(groupings)[g])
    xname=paste0(out_dir,obj_name,".pdf")

    pdf(xname,width=3,height=3.5)

    obj=get(obj_name)
    with(obj,{
      pname="cal"
      labs=c("Overall",names(groupings[[g]]))
      xcol=phs_colours(c("phs-blue","phs-magenta","phs-purple"))
      cal_2panel(list(cal_overall,calA,calB),
                 labels = labs,
                 col=xcol,
                 highlight=highlight_value,
                 #yrange=plot_details[[input$pname]]$yrange,
                 #lpos=plot_details[[input$pname]]$lpos,
                 yrange_lower=plot_details[[pname]]$yrange_lower,
                 legend_title="Group")

    })
    dev.off()
  }
}

print("Completed ROC/PRC/CAL")


##**********************************************#
## Main fairness metrics: DP, FPP, FORP etc  ####
##**********************************************#


for (g in 1:length(groupings)) {

  # Loop over v3 and v4
  for (version in c(3,4)) {

    # Loop over categories
    subcats=which(!(names(categories) %in% names(groupings[[g]])))
    for (cc in subcats) {


      ## Demographic parity                        ####

      ## Demographic parity: are the distributions of scores the same for men and women?
      #ks.test(score[male],score[female]) # no

      pname="dp"

      obj_name=paste0(pname,"_v",version,"_",names(groupings)[g],"_",names(categories)[cc])
      obj=get(obj_name)
      obj_list=list(
        list(x=exp(obj$xx),y=obj$yA,ci=obj$cA),
        list(x=exp(obj$xx),y=obj$yB,ci=obj$cB)
      )

      xname=paste0(out_dir,obj_name,".pdf")
      pdf(xname,width=3,height=3.5)

      groupmetric_2panel(obj_list,
                         labels=names(groupings[[g]]),
                         col=phs_colours(c("phs-blue","phs-magenta")),
                         ci_col=phs_colours(c("phs-blue","phs-magenta")),
                         yrange = plot_details[[pname]]$yrange,
                         lpos = plot_details[[pname]]$lpos,
                         yrange_lower = plot_details[[pname]]$yrange_lower,
                         highlight=highlight_value,
                         logscale=TRUE)

      dev.off()


      ## Other main metrics                        ####

      for (s in 1:length(metrics)) {

        pname=metrics[s]
        obj_name=paste0(pname,"_v",version,"_",names(groupings)[g],"_",names(categories)[cc])
        obj=get(obj_name)
        xname=paste0(out_dir,obj_name,".pdf")
        obj_list=list(
          list(x=obj$cutoffs,y=obj$p_AA,ci=obj$ci_AA),
          list(x=obj$cutoffs,y=obj$p_BB,ci=obj$ci_BB)
        )

        pdf(xname,width=3,height=3.5)

        groupmetric_2panel(obj_list,
                           labels= names(groupings[[g]]),
                           col=phs_colours(c("phs-blue","phs-magenta")),
                           ci_col=phs_colours(c("phs-blue","phs-magenta")),
                           yrange = plot_details[[pname]]$yrange,
                           lpos = plot_details[[pname]]$lpos,
                           yrange_lower = plot_details[[pname]]$yrange_lower,
                           highlight=highlight_value)

         dev.off()

      }

    }

  }
}

print("Completed main fairness plots")



##**********************************************#
## Adjusted FORP/FDRP                        ####
##**********************************************#


for (version in c(3,4)) {


  for (g in 1:length(groupings)) {


    for (ptype in 1:2) {

      if (ptype==1) {
        obj_name=paste0("forp_adjusted_v",version,"_",names(groupings)[g])
        pname="forp_adjusted"
      } else {
        obj_name=paste0("fdrp_adjusted_v",version,"_",names(groupings)[g])
        pname="fdrp_adjusted"
      }

      xname=paste0(out_dir,obj_name,".pdf")
      obj=get(obj_name)
      obj_list=list(
        list(x=obj$cutoffs,y=obj$p_AA,ci=obj$ci_AA),
        list(x=obj$cutoffs,y=obj$p_BB,ci=obj$ci_BB)
      )

      pdf(xname,width=3,height=3.5)
      groupmetric_2panel(obj_list,
                         labels=names(groupings[[g]]),
                         col=phs_colours(c("phs-blue","phs-magenta")),
                         ci_col=phs_colours(c("phs-blue","phs-magenta")),
                         yrange = plot_details[[pname]]$yrange,
                         lpos = plot_details[[pname]]$lpos,
                         yrange_lower = plot_details[[pname]]$yrange_lower,
                         highlight=highlight_value)
      dev.off()

    }

  }

}

print("Completed adjusted plots")


##**********************************************#
## Counterfactual fairness                   ####
##**********************************************#

for (g in 1:length(groupings)) {
  for (version in c(3,4)) {

    pname="counterfactual_dp"
    obj_name=paste0("counterfactual_dp_v",version,"_",names(groupings)[[g]])
    obj=get(obj_name)
    xname=paste0(out_dir,obj_name,".pdf")
    obj_list=list(
      list(x=exp(obj$xx),y=obj$yA,ci=obj$cA),
      list(x=exp(obj$xx),y=obj$yB,ci=obj$cB)
    )

    pdf(xname,width=3,height=3.5)

    groupmetric_2panel(obj_list,
                       labels=names(groupings[[g]]),
                       col=phs_colours(c("phs-blue","phs-magenta")),
                       ci_col=phs_colours(c("phs-blue","phs-magenta")),
                       yrange = plot_details[[pname]]$yrange,
                       lpos = plot_details[[pname]]$lpos,
                       yrange_lower = plot_details[[pname]]$yrange_lower,
                       highlight=highlight_value,
                       logscale=TRUE)

    dev.off()


  }
}

print("Completed counterfactual plots")


##**********************************************#
## FORP decomposition                        ####
##**********************************************#

cutoff=10/100 # 10% risk score threshold

for (version in c(3,4)) {
  
  for (g in 1:length(groupings)) {
    
    obj_name=paste0("forp_decomposition_v",version,"_",names(groupings)[[g]])
    xname=paste0(out_dir,obj_name,".pdf")
    
    
    names_group1=paste0("v",version,"_",names(groupings[[g]])[1],"_q",1:20)
    names_group2=paste0("v",version,"_",names(groupings[[g]])[2],"_q",1:20)
    decomp1=decomposition_matrix[names_group1,]
    decomp2=decomposition_matrix[names_group2,]
    
    
    pdf(xname,width=6,height=6)
    plot_decomp(decomp1,
                decomp2,
                threshold=cutoff,
                labels=names(groupings[[g]]))
    dev.off()
    
    
    for (subgroup in names(groupings[[g]])) {
      
      obj_name=paste0("for_breakdown_v",version,"_",subgroup)
      xname=paste0(out_dir,obj_name,".pdf")
      pdf(xname,width=3,height=3.5)
      names_group=paste0("v",version,"_",subgroup,"_q",1:20)
      decomp=decomposition_matrix[names_group,]
      for_breakdown(decomp,group = subgroup,threshold = cutoff,ylimit=c(-0.065,0.065),ldiff=0.003)
      dev.off()
      
    }
    
  }
  
  obj_name=paste0("forp_decomposition_v",version,"_all")
  xname=paste0(out_dir,obj_name,".pdf")
  
  
  names_group1=paste0("v",version,"_all_q",1:20)
  names_group2=paste0("v",version,"_all_q",1:20)
  decomp1=decomposition_matrix[names_group1,]
  decomp2=decomposition_matrix[names_group2,]
  
  
  pdf(xname,width=6,height=6)
  plot_decomp(decomp1,
              decomp2,
              threshold=cutoff,
              labels=names(groupings[[g]]))
  dev.off()
  
  
  obj_name=paste0("for_breakdown_v",version,"_all")
  xname=paste0(out_dir,obj_name,".pdf")
  pdf(xname,width=3,height=3.5)
  names_group=paste0("v",version,"_all_q",1:20)
  decomp=decomposition_matrix[names_group,]
  for_breakdown(decomp,group = subgroup,threshold = cutoff,ylimit=c(-0.065,0.065),ldiff=0.003)
  dev.off()
  
}


print("Completed decomposition")


##**********************************************#
## Organise figures for manuscript           ####
##**********************************************#

gg=c("all",names(groupings))
cc=c("all",names(categories))

# Calibration
calnames=paste0("cal_v3_",gg,".pdf")
for (i in 1:length(calnames)) {
  file.copy(paste0(out_dir,calnames[i]),paste0("Figures/paper_figures/Calibration/",calnames[i]),overwrite=TRUE)
}

# Counterfactuals
cfnames=paste0("counterfactual_dp_v3_",gg,".pdf")
for (i in 1:length(cfnames)) {
  file.copy(paste0(out_dir,cfnames[i]),paste0("Figures/paper_figures/Counterfactual/",cfnames[i]),overwrite=TRUE)
}

# Decomposition
dnames=c(
  paste0("for_breakdown_v3_",cc,".pdf"),
  paste0("forp_decomposition_v3_",gg,".pdf"))
for (i in 1:length(dnames)) {
  file.copy(paste0(out_dir,dnames[i]),paste0("Figures/paper_figures/Decomposition/",dnames[i]),overwrite=TRUE)
}

# DP
dpnames=c(
  paste0("dp_v3_",gg,"_all.pdf"),
  "dp_v3_Age_Female.pdf",
  "dp_v3_SIMD_Rural.pdf")
for (i in 1:length(dpnames)) {
  file.copy(paste0(out_dir,dpnames[i]),paste0("Figures/paper_figures/Demographic_parity/",dpnames[i]),overwrite=TRUE)
}

# FDRP
fdnames=c(
  paste0("fdrp_v3_",gg,"_all.pdf"),
  paste0("fdrp_adjusted_v3_",gg,".pdf"))
for (i in 1:length(fdnames)) {
  file.copy(paste0(out_dir,fdnames[i]),paste0("Figures/paper_figures/FDRP/",fdnames[i]),overwrite=TRUE)
}

# FORP
fonames=c(
  paste0("forp_v3_",gg,"_all.pdf"),
  paste0("forp_adjusted_v3_",gg,".pdf"))
for (i in 1:length(fonames)) {
  file.copy(paste0(out_dir,fonames[i]),paste0("Figures/paper_figures/FORP/",fonames[i]),overwrite=TRUE)
}

# ROC
rocnames=paste0("roc_v3_",gg,".pdf")
for (i in 1:length(rocnames)) {
  file.copy(paste0(out_dir,rocnames[i]),paste0("Figures/paper_figures/ROC/",rocnames[i]),overwrite=TRUE)
}




##**********************************************#
## Generate main table                       ####
##**********************************************#

# Writes demographic table as latex table.
otab=otab0

# Fix problem with integer overflow
for (i in 1:dim(otab)[2]) otab[which(otab[,i]<1e-5),i]=0

# Decimal places
for (i in 1:dim(otab)[2]) {
  if ("numeric" %in% class(otab[,i])) otab[,i]=signif(otab[,i],digits=3)
}

# Means and SDs for age
age=otab[,3]; age_sd=otab[,4]
otab[,3]=paste0(otab[,3]," (",otab[,4],")")
otab=otab[,-c(4,6:8)]

# Rename rows
rownames(otab)=c("All", "Male", "Female", "Most dep.", "Least dep.", 
                 "Over 65", "Under 25", "White eth.", "Nonwhite eth.", 
                 "Urban", "Rural", "Mainland", 
                 "Island", "Missing eth.", "Post. missing", 
                 "Post2. missing", "Other age", "Mid dep.")
ovec=c(1:5,18,6:7,17,8:9,14,10:11,15,12:13)
otab=otab[ovec,]
age=age[ovec]; age_sd=age_sd[ovec]

# Highlight differences
rdif=0.05 # Highlight this absolute percentage difference
sedif=5 # at this many SEs

low_text="{\\color{blue} " # do this to low values
high_text="{\\color{red} " # do this to high values

N=otab$N

# Age
age_se=age_sd/sqrt(N)
w=which(age>age[1]*(1+rdif) & abs(age-age[1])/age_se>5)
z=which(age<age[1]*(1-rdif) & abs(age-age[1])/age_se>5)
otab[w,3]=paste0(high_text,otab[w,3],"}")
otab[z,3]=paste0(low_text,otab[z,3],"}")

# Everything else
for (i in c(2,4:dim(otab)[2])) {
  xse=sqrt(otab[,i]*(100-otab[,i])/(100^2 * N))
  w=which(otab[,i]>otab[1,i]*(1+rdif) & abs(otab[,i]-otab[1,i])/xse>5)
  z=which(otab[,i]<otab[1,i]*(1-rdif) & abs(otab[,i]-otab[1,i])/xse>5)
  otab[w,i]=paste0(high_text,otab[w,i],"}")
  otab[z,i]=paste0(low_text,otab[z,i],"}")
}

otab=cbind(rownames(otab),otab)
otabx=apply(otab,1,function(x) paste(paste(x,collapse=" & ")," \\\\"))
names(otabx)=NULL
x="\\hline"
otab_lines=c(otabx[1],x,otabx[2:3],x,otabx[4:6],x,
             otabx[7:9],x,otabx[10:12],x,otabx[13:15],x,otabx[16:17])

con=file("Tables/demtab.txt")
writeLines(otab_lines, con)
close(con)


print("Generated demographic table")