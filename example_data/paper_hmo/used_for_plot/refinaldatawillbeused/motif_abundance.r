#### load libraries

#source('functions.R')
#source('unittests.r')
library(reshape2)
library(openxlsx)
library(lme4)
library(heatmap.plus)
library(ggrepel)
# library(ggforce)


#### load data

clean<-function(x,skip=1){
	x = x[,-skip]
	x$Dataset = paste0('d',x$Dataset)
	x$Secretor = ! x$Pateint.ID %in% c('L2','L3')
	x$scretor_col = ifelse(x$Secretor,'black','grey')
	x$DPP_col = heat.colors(length(unique(x$DPP)))[factor(x$DPP)]
	x= droplevels(na.omit(x[ x$Dataset=='d1' ,] ))
	x$subject_col = rainbow(length(levels(x$Pateint.ID)))[x$Pateint.ID]
	x
}

ma = clean( read.csv('/Users/apple/PycharmProjects/GlyCompare/hmo/refinaldatawillbeused/motif_abundance.csv') )
mand = clean( read.csv('/Users/apple/PycharmProjects/GlyCompare/hmo/refinaldatawillbeused/motif_abundance_nodes_dropped.csv') )
mandm = clean( read.csv('/Users/apple/PycharmProjects/GlyCompare/hmo/refinaldatawillbeused/motif_abundance_nodes_dropped_more.csv') )
smandm = clean( read.csv('/Users/apple/PycharmProjects/GlyCompare/hmo/refinaldatawillbeused/smallest_motif_abundance_nodes_dropped_more.csv') )
gly = clean ( tmp<-read.csv('/Users/apple/PycharmProjects/GlyCompare/hmo/Data12.csv'), 21)

vars = c( grep('X',colnames(ma),value=T), colnames(gly)[4:19])

gly_p = gly
gly_p[,colnames(gly)[4:19]] = gly[,colnames(gly)[4:19]]/as.numeric(as.character(gly$SUM))

data=list(motif_abundance=ma,motif_abundance_dropped=mand,motif_abundance_dropped_more=mandm,smallest_motif_abundance_dropped_more=smandm,glycan_concentration=gly,glycan_percent=gly_p)

#### heatmaps

pdf('./secretor.pdf')
lapply(names(data),function(xn){
	print(xn)
	x=data[[xn]]
	data_i = data.matrix(x[,colnames(x)%in%vars])
	cols = as.matrix(x[,rev(c('scretor_col','subject_col','DPP_col'))])
	heatmap.plus( data_i,RowSideColors=cols,main=xn) #, 
#		hclustfun=function(x) hclust(x,method="complete"),distfun=function(x) as.dist((1 - cor(  t(x) ,method='spearman' ))/2))
	heatmap.plus( scale(data_i),RowSideColors=cols,main=paste('z-normalized',xn)) #,
#		hclustfun=function(x) hclust(x,method="complete"),distfun=function(x) as.dist((1 - cor(  t(x) ,method='spearman' ))/2))
	})
dev.off()

pdf('./secretor.pearson_dist.pdf')
lapply(names(data),function(xn){
	print(xn)
	x=data[[xn]]
	data_i = data.matrix(x[,colnames(x)%in%vars])
	cols = as.matrix(x[,rev(c('scretor_col','subject_col','DPP_col'))])
	heatmap.plus( data_i,RowSideColors=cols,main=xn, 
		hclustfun=function(x) hclust(x,method="complete"),distfun=function(x) as.dist((1 - cor(  t(x) ,method='pearson' ))/2))
	heatmap.plus( scale(data_i),RowSideColors=cols,main=paste('z-normalized',xn),
		hclustfun=function(x) hclust(x,method="complete"),distfun=function(x) as.dist((1 - cor(  t(x) ,method='pearson' ))/2))
	})
dev.off()

# changes in correlation
lapply(names(data)[-(1:2)],function(xn){
	print(xn)
	x=data[[xn]]
	cr_sec = cor(data.matrix(x[x$Secretor,colnames(x)%in%vars]),method='spearman')
	cr_nsec = cor(data.matrix(x[!x$Secretor,colnames(x)%in%vars]),method='spearman')

	m1=melt(cr_sec)
	m2=melt(cr_nsec)
	m=cbind(m1,m2[,3])
	m[,1:2] = t(apply( m[,1:2],1,sort))

	colnames(m)[3:4] = c('value_sec','value_nsec')
	m$prod = (m[,3]^2+m[,4]^2)^.5
	m$lab = ifelse(m$prod>.5 & sign(m[,3])!=sign(m[,4]),paste(gsub('X','',m$X1),gsub('X','',m$X2),sep='~'),NA)

	g=ggplot( unique(m) , aes(x=value_sec,y=value_nsec , label=lab)) + geom_label_repel(size=7) + geom_point()+
		geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_abline() +
		theme_classic(base_size=35)
	ggsave(g,filename=paste0('motifAnyl/cor_diff.',xn,'.pdf'))
		#+geom_circle() # +stat_ellipse(level=.5)
})

#### regression
library(lme4)
library(gee)

i=0
out_i=lapply( data, function(x){
	i=i+1
	x$Secretor = as.numeric(x$Secretor)
#	x$Secretor = relevel(factor(ifelse(x$Secretor,'secretor','nonsecretor')),ref='nonsecretor')
	mods_i = lapply( colnames(x)[colnames(x)%in%vars] , function(motif){
		mod=NA
		#try(mod<-gee(as.formula(paste('Secretor ~log(DPP)+scale(',motif,')')),id=Pateint.ID,family='binomial',data=droplevels(x)))
		try(mod<-gee(as.formula(paste('Secretor ~log(DPP)+scale(',motif,')')),id=Pateint.ID,family='binomial',data=droplevels(x),corstr='exchangeable'))
#		if(i>=4){
#			try(mod<-gee(as.formula(paste('Secretor ~log(DPP)+log(',motif,')')),id=Pateint.ID,family='binomial',data=droplevels(x),corstr='exchangeable'))			
#		}else{
#			try(mod<-gee(as.formula(paste('Secretor ~log(DPP)+',motif)),id=Pateint.ID,family='binomial',data=droplevels(x),corstr='exchangeable'))
#		}9632
		#try(mod<-glmer(as.formula(paste('Secretor ~log(DPP)+',motif,'+(1|Pateint.ID)')),family='binomial',data=droplevels(x)))
		#print(coef(summary(mod)))
		if(is.na(mod)){return(NULL)}else{coef(summary(mod))}
	})
	mods=as.data.frame( do.call(rbind,mods_i[!sapply(mods_i,is.null)]) )
	colnames(mods) = c('Coef','naiveSE','naiveZ','robustSE','robustZ')
	mods$vars = gsub('scale|log|\\(|\\)','',rownames(mods))
	mods$response = 'secretor'
	mods$robustPr = 2*pnorm(-abs(mods$robustZ))
	mods1=mods[mods$vars%in%vars,]

	mods_i = lapply( colnames(x)[colnames(x)%in%vars] , function(motif){
		mod=NA
		try(mod<-gee(as.formula(paste('log(DPP) ~Secretor+scale(',motif,')')),id=Pateint.ID,data=droplevels(x),corstr='exchangeable'))
#		if(i>=4){
#			try(mod<-gee(as.formula(paste('log(DPP)~Secretor+log(',motif,')')),id=Pateint.ID,data=droplevels(x),corstr='exchangeable'))			
#		}else{
#			try(mod<-gee(as.formula(paste('log(DPP)~Secretor+',motif)),id=Pateint.ID,data=droplevels(x),corstr='exchangeable'))
#		}		#try(mod<-glmer(as.formula(paste('Secretor ~log(DPP)+',motif,'+(1|Pateint.ID)')),family='binomial',data=droplevels(x)))
		#print(coef(summary(mod)))
		if(is.na(mod)){return(NULL)}else{coef(summary(mod))}
	})
	mods=as.data.frame( do.call(rbind,mods_i[!sapply(mods_i,is.null)]) )
	colnames(mods) = c('Coef','naiveSE','naiveZ','robustSE','robustZ')
	mods$vars = gsub('scale|\\(|\\)','',rownames(mods))
	mods$response = 'DPP'
	mods$robustPr = 2*pnorm(-abs(mods$robustZ))
	mods2=mods[mods$vars%in%vars,]

	rbind(mods1,mods2)
})
out=do.call(rbind,out_i) #[-(4:5)])
out$data = unlist(lapply(strsplit(rownames(out),'\\.'),function(x) x[1]))

#ggplot( out[abs(out$Coef)<15,] , aes(x=Coef,y=-log(robustPr,10),label=vars)) + geom_text(check_overlap=T) + facet_grid(response~data,scale='free') 
g= ggplot( out[abs(out$Coef)<15,] , aes(x=Coef,y=-log(robustPr,10),label=vars,color=vars)) + geom_label_repel() + facet_grid(data~response,scale='free')+geom_hline(yintercept=1.3)+geom_vline(xintercept=0)+theme_classic()
ggsave(plot=g,filename='motifAnyl/motif_assc.pdf',height=15,width=15)
write.csv(out,file='motifAnyl/motif_assc.csv')

library(gridExtra)

g=grid.arrange(grobs=list(
	ggplot(data[[3]],aes(x=log(DPP),y=X80,color=factor(Secretor),shape=Pateint.ID))+geom_point()+stat_smooth(method='lm'),
	ggplot(data[[5]],aes(x=log(DPP),y=DSLNT,color=factor(Secretor),shape=Pateint.ID))+geom_point()+stat_smooth(method='lm'),
	ggplot(data[[5]],aes(x=log(DPP),y=DSLNH,color=factor(Secretor),shape=Pateint.ID))+geom_point()+stat_smooth(method='lm'),
	ggplot(data[[5]],aes(x=log(DPP),y=LSTb,color=factor(Secretor),shape=Pateint.ID))+geom_point()+stat_smooth(method='lm')
	))
ggsave(plot=g,filename='motifAnyl/individual_trajectories.pdf',height=15,width=15)

