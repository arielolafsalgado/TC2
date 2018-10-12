require(igraph)

fitTheEssential <- function(g){
	degs = degree(g)
	K = NULL
	pE = NULL
	for( k in sort(unique(degs))){
		thisK = degs == k
		pess = sum(V(g)$essential[thisK])/sum(thisK)
		K = c(K,k)
		pE = c(pE,pess)
	}
	lpE = log(1-pE)
	mod = lm(lpE[pE<1 & K<=10]~K[pE<1 & K<=10])
	alfa = summary(mod)$coefficients[2,1]
	beta = summary(mod)$coefficients[1,1]
	alfa = 1-exp(alfa)
	beta = 1-exp(beta)
	return(list('k'=K,'pE'=pE,'lpE'=lpE,'model'=mod,'alfa'=alfa,'beta'=beta))
}

g.list = list('y2h'=g_y2h,'apms'=g_apms,'lit'=g_lit,'litr'=g_litr)

fit.list = lapply(g.list,fitTheEssential)

pdf('puntoDHe.pdf')
for(gn in names(g.list)){
	lpE = fit.list[[gn]]$lpE
	k = fit.list[[gn]]$k
	mod = fit.list[[gn]]$model
	pendiente = round(summary(mod)$coefficients[2,1],3)
	sdpendiente = round(summary(mod)$coefficients[2,2],3)
	intercept = round(summary(mod)$coefficients[1,1],3)
	sdintercept = round(summary(mod)$coefficients[1,2],3)
	legendLine = paste('(',pendiente,'±',sdpendiente,')',sep='')
	legendLine = paste(legendLine,'k')
	legendLine = paste(legendLine,paste('(',intercept,'±',sdintercept,')',sep=''),sep='+')
	plot(k,lpE,pch=18,col='black',xlab='grado',ylab='ln(1-pE)',xlim=c(0,10),main=gn)
	lines(k[is.finite(lpE) & k<=10],predict(mod),col='red')
	legend(1,-.7,legend=legendLine,lty='solid',col='red')
}
dev.off()


### La medición de zotenko

tablaZotenko <- function(g){
	require(Matrix)
	# Total de pares (con vecinos en comun mayores a 3
	X = as_adjacency_matrix(g)
	XX = X%*%X
	diag(XX) = 0
	total.pares = sum(XX>=3)/2
	#Total pares essenciales
	A = V(g)$essential*t(X*V(g)$essential) | (!V(g)$essential*t(X*!V(g)$essential)) 
	total.iguales = sum(A*(XX>=3))/2
}
