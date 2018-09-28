require(igraph)
graphics.off()
data_y2h = as.matrix(read.table('yeast_Y2H.txt'))
data_APMS = as.matrix(read.table('yeast_AP-MS.txt'))
data_LIT = as.matrix(read.table('yeast_LIT.txt'))
data_LITR = as.matrix(read.csv('yeast_LIT_Reguly.txt',sep='\t'))[,c(1,2)]
#data_HE = as.matrix(read.table('Essential_ORFs_paperHE.txt'))

g_y2h = simplify(graph_from_edgelist(data_y2h))	
g_apms = simplify(graph_from_edgelist(data_APMS))
g_lit = simplify(graph_from_edgelist(data_LIT))
g_litr = simplify(graph_from_edgelist(data_LITR))
#g_he = graph_from_edgelist(data_HE)

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
#DaTa = cbind(DaTa,datosUtiles(g_he))
colnames(DaTa) = c('y2h','apms','lit','litr')
