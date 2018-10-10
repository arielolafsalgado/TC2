source('tc02_init.R')
## Resumen Tabla 1
datosUtiles <- function(g){
	varNames = c('N','M','z','p','S','p0','nClust','Clustering')
	dt = matrix(NA,ncol=1,nrow=length(varNames))
	rownames(dt)  = varNames
	dt['N',1] = length(V(g))
	dt['M',1] = length(E(g))
	dt['z',1] = length(E(g))/length(V(g))
	dt['p',1] = length(E(g))/(length(V(g))*(length(V(g))-1))
	clt = clusters(g)
	dt['S',1] = max(clt$csize)/length(V(g))
	dt['p0',1] = sum(degree(g)==0)/length(V(g))	
	dt['nClust',1] = length(clt$csize)
	dt['Clustering',1] = transitivity(g)
	return(dt)
}

DaTa = datosUtiles(g_y2h)
DaTa = cbind(DaTa,datosUtiles(g_apms))
DaTa = cbind(DaTa,datosUtiles(g_lit))
DaTa = cbind(DaTa,datosUtiles(g_litr))
colnames(DaTa) = c('y2h','apms','lit','litr')


## Segunda tabla
proporcionDeEges=matrix(NA,nrow=4,ncol=4)
colnames(proporcionDeEges)=c('y2h','apms','lit','litr')
rownames(proporcionDeEges)=c('y2h','apms','lit','litr')
g.list = list(g_y2h, g_apms, g_lit, g_litr)
for (i in 1:4){
	for (j in 1:4){
		proporcionDeEges[i,j]=length(intersect(as_ids(E(g.list[[i]])),as_ids(E(g.list[[j]]))))/length(E(g.list[[i]]))
	}
}



## Grafico zopenko
hub.vs.ess = function(gg, norm=TRUE){
	degs = degree(gg)
	deglist = 0:max(degs)
	ess_count = rep(0, max(degs) + 1)
	for(i in deglist){
		ess_count[i+1] = sum(V(gg)$essential[degs >= i])
	}

	if(norm){
		ess_count = ess_count / ess_count[1]
		deglist = deglist / max(deglist)
	}
	
	data.frame(DEGS = deglist, ESS.CNT = ess_count)
}

g.list = list(g_y2h, g_apms, g_lit, g_litr)
hub.ess.df = lapply(g.list, hub.vs.ess)
pdf('zopenko1.pdf')
plot(hub.ess.df[[1]],col='red',lty='solid',type='l',lwd=3,xlab='HubCutOff',ylab='FracEssCap')
lines(hub.ess.df[[2]],col='green',lwd=3)
lines(hub.ess.df[[3]],col='blue',lwd=3)
lines(hub.ess.df[[4]],col='violet',lwd=3)

plot(hub.ess.df[[1]],col='red',lty='solid',type='l',lwd=3,log='xy',xlab='HubCutOff',ylab='FracEssCap')
lines(hub.ess.df[[2]],col='green',lwd=3)
lines(hub.ess.df[[3]],col='blue',lwd=3)
lines(hub.ess.df[[4]],col='violet',lwd=3)

plot(hub.ess.df[[1]],col='red',lty='solid',type='l',lwd=3,log='x',xlab='HubCutOff',ylab='FracEssCap')
lines(hub.ess.df[[2]],col='green',lwd=3)
lines(hub.ess.df[[3]],col='blue',lwd=3)
lines(hub.ess.df[[4]],col='violet',lwd=3)



dev.off()

