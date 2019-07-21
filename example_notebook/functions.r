
library(ggplot2)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(broom)
library(gplots)
library(reshape2)
library(openxlsx)
library(heatmap.plus)
library(ggrepel)
library(xtable)
library(gee)
library(geepack)
library(igraph)

library(effsize)
library(ggpubr)

parseSummary.GEEpack <- function(mod,odsig=3,logit=FALSE){
    #sample size
    observation_count = length(mod$id)
    group_count = length(unique(mod$id))
    #test statistics
    summary=coef(summary(mod))
    # residuals
    res = resid(mod)
    #confidence intervals
    if(logit){
        CI=paste(   '(',   signif( exp(summary$Estimate)  - 1.96 * exp(summary$Estimate) * summary[,2] ,odsig),' - ',
                 signif( exp(summary$Estimate)  + 1.96 * exp(summary$Estimate) * summary[,2] ,odsig), ')',sep='')
    }else{
        CI=paste(   '(',   signif( summary$Estimate  - 1.96 * summary$Estimate * summary[,2] ,odsig),' - ',
                signif( summary$Estimate  + 1.96 * summary$Estimate * summary[,2] ,odsig), ')',sep='')
    }
    #effect size: marginal r^2 (Zheng 2000)
    y = mod$y; y_bar=mean(y); y_hat=mod$fitted.values
    mr2=1-sum((y-y_hat)^2)/sum((y-y_bar)^2)
    #degrees of freedom
    df=mod$df.residual
    # anova (dQIC)
    #F=anova(mod, mod0motifs )
    # dfbeta/cooks/influence.measure
    shapiro.wilks=ifelse(logit,NA,shapiro.test(mod$y)$p.value)
    
    if(!is.data.frame(summary)){
    # get matrix
    summary = do.call(rbind,summary);			
    # order and data.frame
    if(glm){
      summary = data.frame(summary[order(rownames(summary),summary[,4]),])
    }else{
      summary = data.frame(summary[order(summary[,4]),])
    }
    }
    
    if(logit){
        summary$'Coef' = signif( exp(summary$Estimate) , odsig )
    }else{
        summary$'Coef' = signif( summary$Estimate , odsig )
    }
    summary$'95 CI' = CI
    summary$'Pr(W)' = summary[,4] # NEJM_pval(summary[,4])

    # gen table 
    xtab <- xtable(summary[,5:7])
    #print(xtable(summary(res)))
    singlevars = data.frame(
        stat=c('Number of observations','Number of Clusters','Marginal R^2','Degrees of Freedom',"Shapiro-Wilks P"),
        value=c(observation_count,signif(group_count,odsig),signif(mr2,odsig),df,shapiro.wilks))
    xtab2<-xtable(singlevars)

    return(list(xtab,xtab2))
}