library(ppitools)
library(igraph)
library(magrittr)

data<-read.table("~/Documents/Projects/PPI_framework/PPItools/Manual/pvals.txt",header=T,sep = "\t")
pvals<-setNames(data[,2],data[,1])
reference_interactome<-read.table("~/Documents/Projects/PPI_framework/PPI_algorithm/INWEB_IM/InBio_Map_core_2016_09_12/inwebIM_ppi.txt")
net<-igraph::graph_from_edgelist(as.matrix(reference_interactome),
                                 directed=F)
net<-delete_vertices(net,"UBC")
dd<-igraph::degree(net)
net<-igraph::delete_vertices(net,names(dd[which(dd==0)]))

net<-preprocess(pvals=pvals,
                network=net,
                simplify = T,
                scf=bum_score(pvals,
                              threshold_pval=5e-8,
                              type="aggressive",
                              bum_plot=FALSE))

module<-estimate_network(network = net$net)
V(module)$size<-15
plot(tmp,layout=layout_with_fr,vertex.label.cex=1,vertex.label.dist=1.5)
tmp<-module
V(tmp)$name[1:9]<-""
V(tmp)$name[11:25]<-""
V(tmp)$name[27:34]<-""
V(tmp)$color<-c(rep(rgb(169/255,169/255,169/255,0.5),9),rgb(0.8,0,0,0.5),rep(rgb(169/255,169/255,169/255,0.5),15),rgb(0.8,0,0,0.5),rep(rgb(169/255,169/255,169/255,0.5),8))
plot(tmp,layout=layout_with_fr,vertex.label.cex=1.5,vertex.label.dist=1.5,vertex.shape="sphere")
V(tmp)$size<-6
V(tmp)$size[c(10,26)]<-12

ModuleSignificance<-BM_pvalue(pvals = net$pvals, 
                              threads = 4,
                              permutation_method = "igraph",
                              n_permutations = 100,
                              simplify = TRUE,
                              network = net$net,
                              scoring_function = net)

system.time(RndGeneFreq.perm.both<-permutation_freq(pvals = net$pvals,
                              network = net$net,
                              scoring_function = net,
                              simplify = TRUE,
                              threads = 3,
                              permute_interactome = TRUE,
                              permute_pvals = TRUE,
                              n_permutations = 1000,
                              permutation_method = "igraph"))

system.time(RndGeneFreq.perm.both.tmp<-permutation_freq(pvals = net$pvals,solver_time_factor = 0.2,permutation_time_factor = 0.2,
                                        network = net$net,
                                        scoring_function = net,
                                        simplify = TRUE,
                                        threads = 3,
                                        permute_interactome = TRUE,
                                        permute_pvals = TRUE,
                                        n_permutations = 1000,
                                        permutation_method = "ppitools"))

random.connection.prob.two.nodes<-function(TotalGenes,PositiveNodes,Connections){
  A=1
  for (i in 1:Connections){
    A=A*(TotalGenes-PositiveNodes-Connections+i)/(TotalGenes-Connections+i)
  }
  B=A*Connections/(TotalGenes-Connections-PositiveNodes+1)*PositiveNodes
  return(1-A-B)
}
random.connection.prob.single.node<-function(TotalGenes,PositiveNodes,Connections){
  A=1
  for (i in 1:Connections){
    A=A*(TotalGenes-PositiveNodes-Connections+i)/(TotalGenes-Connections+i)
  }
  return(1-A)
}

pr<-0
neg.genes<-names(which(igraph::degree(net$net,v = names(which(net$score(net$pvals)<0)))>1))
for (i in neg.genes){
  pr<-pr+random.connection.prob.two.nodes(TotalGenes = 14925,PositiveNodes = 26,Connections = igraph::degree(net$net,v = i))
}

module1<-estimate_network(network = delete_vertices(net$net,"TMEM38B"))
module2<-estimate_network(network = delete_vertices(net$net,"TAL2"))
module3<-estimate_network(network = delete_vertices(net$net,"RAD23B"))
###TMEM38BvsTAL2vsRAD23B
par(mar=c(6,5,3,2))
plot(setNames(c(sum(net$score(net$pvals[V(module)$name])),sum(net$score(net$pvals[V(module)$name]))-sum(net$score(net$pvals["TAL2"])),sum(net$score(net$pvals[V(module1)$name])),sum(net$score(net$pvals[V(module2)$name])),sum(net$score(net$pvals[V(module3)$name]))),c("Observed","Expected drop","Observed-TMEM38B","Observed-TAL2","Observed-RAD23B")),pch=2,ylab="Cumulative Score",type="b",xaxt="n",xlab="",cex.lab=1.5,cex.axis=1.5)
lines(x = c(1,1),y = c(171.5,174),lty=2)
lines(x = c(2,2),y = c(171.5,174),lty=2)
lines(x = c(3,3),y = c(171.5,174),lty=2)
lines(x = c(4,4),y = c(171.5,174),lty=2)
lines(x = c(4,4),y = c(171.5,174),lty=2)
axis(1,at=1:5,labels=c("Observed","Expected\ndrop","Observed\n-TMEM38B","Observed\n-TAL2","Observed\n-RAD23B"),cex.axis=1.5,las=1,padj = 1)

plot(setNames(c(length(E(module))/length(V(module)),length(E(module1))/length(V(module1)),length(E(module2))/length(V(module2))),c("Observed","Observed-TPCN2","Observed-CCND1")),pch=2,ylab="Edges/Node",type="b",xaxt="n",xlab="",cex.lab=1.5,cex.axis=1.5)
lines(x = c(1,1),y = c(169,174),lty=2)
lines(x = c(2,2),y = c(169,174),lty=2)
lines(x = c(3,3),y = c(169,174),lty=2)
#lines(x = c(4,4),y = c(169,174),lty=2)
axis(1,at=1:3,labels=c("Observed","Observed\n-TPCN2","Observed\n-CCND1"),cex.axis=1.5,las=1,padj = 1)

##CCND1vsTPCN2
module1<-estimate_network(network = delete_vertices(net$net,"TPCN2"))
module2<-estimate_network(network = delete_vertices(net$net,"CCND1"))
module3<-estimate_network(network = delete_vertices(net$net,"TERT"))

par(mar=c(6,5,3,2))
plot(setNames(c(sum(net$score(net$pvals[V(module)$name])),sum(net$score(net$pvals[V(module)$name]))-sum(net$score(net$pvals["TPCN2"])),sum(net$score(net$pvals[V(module1)$name])),sum(net$score(net$pvals[V(module2)$name]))),c("Observed","Expected drop","Observed-TPCN2","Observed-CCND1")),pch=2,ylab="Cumulative Score",type="b",xaxt="n",xlab="",cex.lab=1.5,cex.axis=1.5)
lines(x = c(1,1),y = c(169,174),lty=2)
lines(x = c(2,2),y = c(169,174),lty=2)
lines(x = c(3,3),y = c(169,174),lty=2)
lines(x = c(4,4),y = c(169,174),lty=2)
axis(1,at=1:4,labels=c("Observed","Expected\ndrop","Observed\n-TPCN2","Observed\n-CCND1"),cex.axis=1.5,las=1,padj = 1)

plot(setNames(c(length(E(module))/length(V(module)),length(E(module1))/length(V(module1)),length(E(module2))/length(V(module2))),c("Observed","Observed-TPCN2","Observed-CCND1")),pch=2,ylab="Edges/Node",type="b",xaxt="n",xlab="",cex.lab=1.5,cex.axis=1.5,ylim=c(1.30,1.46))
points(x = c(length(E(module))/length(V(module)),(length(E(module))-1)/(length(V(module))-1),(length(E(module))-4)/(length(V(module))-1)), type="b",pch=17,col="darkgrey")
lines(x = c(1,1),y = c(1.3,1.45),lty=2)
lines(x = c(2,2),y = c(1.3,1.45),lty=2)
lines(x = c(3,3),y = c(1.3,1.45),lty=2)
#lines(x = c(4,4),y = c(169,174),lty=2)
axis(1,at=1:3,labels=c("Observed","Observed\n-TPCN2","Observed\n-CCND1"),cex.axis=1.5,las=1,padj = 1)
legend("bottomleft",col=c("grey","black"),pch=c(17,2),legend = c("Expected","Observed"))
##SETDB1vsARNT
module1<-estimate_network(network = delete_vertices(net$net,"SETDB1"))
module2<-estimate_network(network = delete_vertices(net$net,"ARNT"))

par(mar=c(6,5,3,2))
plot(setNames(c(sum(net$score(net$pvals[V(module)$name])),sum(net$score(net$pvals[V(module)$name]))-sum(net$score(net$pvals["SETDB1"])),sum(net$score(net$pvals[V(module1)$name])),sum(net$score(net$pvals[V(module2)$name]))),c("Observed","Expected drop","Observed-SETDB1","Observed-ARNT")),pch=2,ylab="Cumulative Score",type="b",xaxt="n",xlab="",cex.lab=1.5,cex.axis=1.5)
lines(x = c(1,1),y = c(169,174),lty=2)
lines(x = c(2,2),y = c(169,174),lty=2)
lines(x = c(3,3),y = c(169,174),lty=2)
lines(x = c(4,4),y = c(169,174),lty=2)
axis(1,at=1:4,labels=c("Observed","Expected\ndrop","Observed\n-SETDB1","Observed\n-ARNT"),cex.axis=1.5,las=1,padj = 1)

plot(setNames(c(length(E(module))/length(V(module)),length(E(module1))/length(V(module1)),length(E(module2))/length(V(module2))),c("Observed","Observed-TPCN2","Observed-CCND1")),pch=2,ylab="Edges/Node",type="b",xaxt="n",xlab="",cex.lab=1.5,cex.axis=1.5)
lines(x = c(1,1),y = c(169,174),lty=2)
lines(x = c(2,2),y = c(169,174),lty=2)
lines(x = c(3,3),y = c(169,174),lty=2)
#lines(x = c(4,4),y = c(169,174),lty=2)
axis(1,at=1:3,labels=c("Observed","Observed\n-TPCN2","Observed\n-CCND1"),cex.axis=1.5,las=1,padj = 1)

means<-c(2.3,2.63)
bottom_ci<-c(1.012,1.154216)
upper_ci<-c(4.58191,5.272603)
barCenters<-barplot(height=means,beside=F,las=1,names=c("All GWAS","Mapped GWAS"),ylim=c(0,5),main="",col=c(gray.colors(2)), ylab="", xlab="", border="black",axes=T,legend.text=T,args.legend=list(title="Gene set", x="topright", legend=c("All GWAS","Mapped GWAS")),cex.axis=1.5,cex.names=1.5)
segments(barCenters,bottom_ci,barCenters,upper_ci,lwd=1.5)
arrows(barCenters,bottom_ci,barCenters,upper_ci,lwd=1.5,angle=90,code=3, length=0.05)


pvals_tmp<-pvals
pvals_tmp["ARNT"]<-1e-8
net1<-preprocess(pvals = pvals_tmp,network = ref_net,simplify=T,
                 scf=bum_score(pvals,                                                                              threshold_pval=5e-8,
                               type="aggressive",bum_plot=FALSE))

module_test<-estimate_network(net1$net)
"ARNT"%in%V(module_test)$name