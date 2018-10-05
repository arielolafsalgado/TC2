require(igraph)
require(stringr)
graphics.off()
data_y2h = as.matrix(read.table('yeast_Y2H.txt'))
data_APMS = as.matrix(read.table('yeast_AP-MS.txt'))
data_LIT = as.matrix(read.table('yeast_LIT.txt'))
data_LITR = as.matrix(read.csv('yeast_LIT_Reguly.txt',sep='\t'))[,c(1,2)]
EP = read.csv('Essential_ORFs_paperHe.txt',sep='\t')
EP = as.character(EP$ORF_name)
EP = EP[EP!='']
EP = str_trim(EP)
EP = toupper(EP)

g_y2h = as.undirected(simplify(graph_from_edgelist(data_y2h)))
g_apms = as.undirected(simplify(graph_from_edgelist(data_APMS)))
g_lit = as.undirected(simplify(graph_from_edgelist(data_LIT)))
g_litr = as.undirected(simplify(graph_from_edgelist(data_LITR)))

V(g_y2h)$essential = is.element(V(g_y2h)$name,EP)
V(g_apms)$essential = is.element(V(g_apms)$name,EP)
V(g_lit)$essential = is.element(V(g_lit)$name,EP)
V(g_litr)$essential = is.element(V(g_litr)$name,EP)

V(g_y2h)$color = ifelse(V(g_y2h)$essential,yes='red',no='green')
V(g_apms)$color = ifelse(V(g_apms)$essential,yes='red',no='green')
V(g_lit)$color = ifelse(V(g_lit)$essential,yes='red',no='green')
V(g_litr)$color = ifelse(V(g_litr)$essential,yes='red',no='green')

