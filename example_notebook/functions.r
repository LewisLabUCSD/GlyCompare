
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

# marginal entropy, GEE (zhang,2000)
marginal_entropy<-function(mod){
    y = mod$y; 
    alpha1=mean(y); alpha0=1-alpha1; # marginal probablility of response k
    pi1 = mod$fitted.values; pi0 = 1-pi1; #model predicted probability of response k
    n=length(y); T=2

    mH = 1 - (sum(pi1*log(pi1))+sum(pi0*log(pi0)))/(n*T*(sum(alpha1*log(alpha1))+sum(alpha0*log(alpha0))))
    return(mH)
}

# marginal r^2, GEE (zhang,2000)
marginal_rsquared<-function(mod,logit){
    if(logit){
        f=as.formula(paste('y~',as.character(mod$formula)[[3]]))
        data = mod$data
        data$y = mod$fitted.values
        mr2=marginal_rsquared(geeglm(f,id=mod$id,data=data),FALSE)
    }else{
        y = mod$y; 
        y_bar=mean(y);
        y_hat=mod$fitted.values
        mr2=1-sum((y-y_hat)^2)/sum((y-y_bar)^2)
    }
    return(mr2)
}

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
    mr2 = marginal_rsquared(mod,logit)
    mh = ifelse(logit, marginal_entropy(mod),NA)
    #degrees of freedom
    df=mod$df.residual
    # anova (dQIC)
    # dfbeta/cooks/influence.measure
    # test normality
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
        stat=c('Number of observations','Number of Clusters','Marginal R^2','Marginal Entropy','Degrees of Freedom',"Shapiro-Wilks P"),
        value=c(observation_count,signif(group_count,odsig),signif(mr2,odsig),signif(mh,odsig),df,shapiro.wilks))
    xtab2<-xtable(singlevars)

    return(list(xtab,xtab2))
}