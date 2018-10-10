require(igraph)
arranca.nodos <- function(g){

	grados = V(g)$name[order(degree(g),decreasing=T)]
	entre = V(g)$name[order(betweenness(g),decreasing=T)]
	trans = V(g)$name[order(transitivity(g,type='local'),decreasing=T)]
	auto = V(g)$name[order(eigen_centrality(g)$vector,decreasing=T)]
	cualca = V(g)$name[sample(length(V(g)))]

	g0 = g

	S_grados = rep(NA,length(V(g)))
	S_entre = rep(NA,length(V(g)))
	S_trans = rep(NA,length(V(g)))
	S_auto = rep(NA,length(V(g)))
	S_cualca = rep(NA,length(V(g)))

	for(i in 1:length(V(g))){
		nodo = grados[i]
		g = delete_vertices(g,nodo)
		S_grados[i] = max(clusters(g)$csize)/length(V(g))
	}
	g = g0
	for(i in 1:length(V(g))){
		nodo = entre[i]
		g = delete_vertices(g,nodo)
		S_entre[i] = max(clusters(g)$csize)/length(V(g))
	}
	g = g0
	for(i in 1:length(trans)){
		nodo = trans[i]
		g = delete_vertices(g,nodo)
		S_trans[i] = max(clusters(g)$csize)/length(V(g))
	}
	g = g0
	for(i in 1:length(auto)){
		nodo = auto[i]
		g = delete_vertices(g,nodo)
		S_auto[i] = max(clusters(g)$csize)/length(V(g))
	}
	g = g0
	for(i in 1:length(cualca)){
		nodo = cualca[i]
		g = delete_vertices(g,nodo)
		S_cualca[i] = max(clusters(g)$csize)/length(V(g))
	}
	g = g0

	return(data.frame('grados'=S_grados,'entrenes'=S_entre,'transidad'=S_trans,'autovalor'=S_auto,'ramdom'=S_cualca))
}

pdf('puntoC1.pdf')
g.list = list(g_y2h,g_apms,g_lit,g_litr)
g.names = list('Y2H','APMS','LIT','LITR')
for(i in 1:4){
	g = g.list[[i]]
	S = arranca.nodos(g)
	N = length(V(g))
	plot((1:N)/N,S$grados,xlab='Fracción de nodos extraidos',ylab='Tamaño componente gigante',log='x',type='l',col='red',lwd=3,xlim=c(1/N,.7),main=g.names[i])
	lines((1:N)/N,S$autovalor,col='yellow',lwd=3)
	lines((1:N)/N,S$ramdom,col='black',lwd=3)
	lines((1:N)/N,S$entrenes,col='cyan',lwd=3)
	lines((1:N)/N,S$transidad,col='magenta',lwd=3)
	S_ess = max(clusters(delete_vertices(g,V(g)$essential))$csize)/length(V(g))
	points(sum(V(g)$essential)/N,S_ess,col='red',pch=18)
	legend(x=10**-3,y=.4,legend=c('Grado','Betweenness','Eigenvalue','Clustering Coef.','Random'),col=c('red','cyan','yellow','magenta','black'),lty='solid',lwd=3)
}
dev.off()
