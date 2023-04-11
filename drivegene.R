setwd('./')
getwd()
names(mul_gene_list) <- c('n','genes')  
mul_genenames <- mul_gene_list$genes      


AdjMatrix <- matrix(data=0,nrow=length(mul_genenames),ncol=length(mul_genenames)) 
rownames(AdjMatrix) <- mul_genenames   
colnames(AdjMatrix) <- mul_genenames
for(i in 1:nrow(mul_edge_list))
{
  tem_list <- mul_edge_list[i,]  
  x<-as.numeric(tem_list[1])   
  y<- as.numeric(tem_list[2])
  AdjMatrix[x,y] = 1
  AdjMatrix[y,x] = 1 
}
dele_row <- which(tf_gene[,3]=='Unknown')  
tf_known_r <- nrow(tf_gene) - length(dele_row)
tf_known <- matrix(0,tf_known_r,3)   
j <- 1
for (i in 1:nrow(tf_gene)) {
  if(tf_gene[i,3]=='Unknown')next()   
  tf_known[j,1] <- tf_gene[i,1]
  tf_known[j,2] <- tf_gene[i,2]
  tf_known[j,3] <- tf_gene[i,3]
  j <- j + 1
}
count <- 0
count <- 0     
for (i in 1:nrow(tf_known)) {
  gene1 <- tf_known[i,1]
  gene2 <- tf_known[i,2]
  index1 <- which(mul_gene_list[,2]==gene1)
  index2 <- which(mul_gene_list[,2]==gene2)
  AdjMatrix[index1,index2] <- 1
  AdjMatrix[index2,index1] <- 1
  count <- count + 1
}




degree<-matrix(0,nrow(AdjMatrix),2)
degree[,1] <- mul_genenames   
for (i in 1:nrow(AdjMatrix)) {
  neib_index <- which(AdjMatrix[i,]!=0)
  if(length(neib_index)==0)next()
  degree[i,2] <- length(neib_index)
}
mul_gene_list[,3]<-0
for(i in 1:nrow(mul_gene_list)){   
  gene_name <- mul_gene_list[i,2]   
  gene_index <- which(degree[,1]==gene_name)
  if(length(gene_index)==0)next()
  mul_gene_list[i,3] <- as.numeric(degree[gene_index,2])
}
maxdegree<- max(mul_gene_list[,3])
mul_gene_list[,4]<-0
for(i in 1:nrow(mul_gene_list))
{
  mul_gene_list[i,4] <- mul_gene_list[i,3]/maxdegree   
}



ECC<-matrix(0,nrow(AdjMatrix),ncol(AdjMatrix))
for(i in 1:10389)
{
  for (j in 1:10389){
      number<-0
     if(AdjMatrix[i,j]==1)
     {
       for (k in 1:10389) {
         if(k!=j && k!=i && AdjMatrix[i,k]==1 && AdjMatrix[j,k]==1)
         {
           number<-number+1
          
         }
       }
       if(length(which(AdjMatrix[i,]==1))>1 && length(which(AdjMatrix[j,]==1))>1)
         {
         ECC[i,j]<-number/min(length(which(AdjMatrix[i,]==1)),length(which(AdjMatrix[j,]==1)))
       }
     }
  }
}
SOECC<-matrix(0,nrow(AdjMatrix),1)
SOECC<-rowSums(ECC)
mul_gene_list[,5]<-0
mul_gene_list[,5]<-SOECC
maxECC<- max(mul_gene_list[,5])
mul_gene_list[,6]<-0
for(i in 1:nrow(mul_gene_list))
{
  mul_gene_list[i,6] <- mul_gene_list[i,5]/maxECC     
}



R_mut <- mRNA_EXP
mrna_exp <- R_mut
gene_exp1 <- matrix(data=0,nrow=nrow(mrna_exp),ncol=6)
gene_exp1[,1] <- mrna_exp[,1]
colnames(gene_exp1) <- c("geng","V1","V2","V3","V4","V5")
rownames(mrna_exp) <- mrna_exp[,1]
mrna_exp <- mrna_exp[,-1]
means<-apply(mrna_exp,1,mean)
var<-apply(mrna_exp,1,var)
sd<-apply(mrna_exp,1,sd)
V<-var/(1+var)
G<-means+2.5*sd*V
gene_exp1[,2] <- means
gene_exp1[,3] <- var
gene_exp1[,4] <- sd
gene_exp1[,5] <- V
gene_exp1[,6] <- G
write.table(gene_exp1,'gene_exp1.txt',quote=FALSE, row.names=FALSE)
write.csv(mrna_exp,'mrna_exp.csv')
write.table(means,'means.txt',quote=F,row.names=F,col.names=F)





remrna_exp <- read.table('./remrna_exp.txt',sep=',')
hangming <- rownames(mrna_exp)
lieming <- colnames(mrna_exp)
rownames(remrna_exp)<-hangming
colnames(remrna_exp)<-lieming
gene_exp2<-matrix(0,nrow(remrna_exp),2)
gene_exp2[,1]<-hangming
mrna_col_name <- colnames(remrna_exp)
normal_exp <- grep('.11A|11B',mrna_col_name)
nor <- length(normal_exp)
cer <- nor + 1
total <- length(mrna_col_name)
for(i in 1:nrow(remrna_exp))
{
  mean1 <- mean(as.numeric(remrna_exp[i,1:nor]))      
  mean2 <- mean(as.numeric(remrna_exp[i,cer:total]))   
  var1 <- sd(as.numeric(remrna_exp[i,1:nor]))
  var2 <- sd(as.numeric(remrna_exp[i,cer:total]))
  var_sum <- var1^2+var2^2
  sum1 <- (mean1-mean2)^2/(4*var_sum) 
  sum2 <-  0.5*log(var_sum/(2*var1*var2))   
  gene_exp2[i,2] <- sum1+sum2    
}

mul_gene_list[,7]<-0
for(i in 1:nrow(mul_gene_list)){   
  gene_name1 <- mul_gene_list[i,2]   
  gene_index1 <- which(gene_exp2[,1]==gene_name1)
  if(length(gene_index1)==0)next()
  mul_gene_list[i,7] <- as.numeric(gene_exp2[gene_index1,2])
}
dele_index <- which(mul_gene_list[,7]=='Inf')
mul_gene_list[dele_index,7] <- 0  








gene_name <- matrix(data=0,nrow=nrow(mul_edge_list),ncol=2)
for (i in 1:nrow(mul_edge_list)) {
  index=mul_edge_list[i,1]
  a=mul_gene_list[index,2]
  gene_name[i,1]=a

}
for (i in 1:nrow(mul_edge_list)) {
  index2=mul_edge_list[i,2]
  a=mul_gene_list[index2,2]
  gene_name[i,2]=a

}
tf=tf_gene
tf=tf[,-c(3,4)]
new=rbind(gene_name,tf)

library(igraph,warn.conflicts = F)
graph=graph.data.frame(new)
Mutation <- SM_mution
Mutation=t(Mutation)
colnames(Mutation)=Mutation[1,]
Mutation=Mutation[-1,]
SM_mut<-Mutation

SM_score=matrix(data=0,nrow=nrow(SM_mut),ncol=2)
SM_score[,1]=rownames(SM_mut)

for(i in 1:nrow(SM_mut))
{
  SM_score[i,2] <- sum(as.numeric(SM_mut[i,]))/ncol(SM_mut)  
}

mul_gene_list[,8]<-0
for(i in 1:nrow(mul_gene_list))
{
  gene_find <- mul_gene_list[i,2]
  index_find <- which(SM_score[,1]==gene_find)      
  if(length(index_find)==0)
  {
    mul_gene_list[i,8] <- 0
    next()
  }
  mul_gene_list[i,8] <- as.numeric(SM_score[index_find,2]) 
}
vertex_W=matrix(0,nrow=nrow(Mutation),ncol = ncol(Mutation))
vertexes=rownames(Mutation)
colnames(vertex_W)=colnames(Mutation)
rownames(vertex_W)=rownames(Mutation)

for(i in 1:ncol(Mutation)){

  v_i=intersect(vertexes[which(Mutation[,i]==1)],V(graph)$name)  

  v_diff=setdiff(vertexes[which(Mutation[,i]==1)],V(graph)$name)  


  if(length(v_i)!=0){
    for (j in 1:length(v_i)){
      V=v_i[j]
      A=which(SM_score[,1]==V)
      Fdegree=-(as.numeric(SM_score[A,2]))*log2(as.numeric(SM_score[A,2]))
      vertex_W[V,i]=Fdegree
    }
  }
  vertex_W[v_diff,i]=0.1                         
}
vertex_W[is.na(vertex_W)]<-0
hypergraph=matrix(0,nrow =nrow(vertex_W),ncol=2 )
hypergraph[,1]=rownames(vertex_W)
TT=rowSums(vertex_W)
hypergraph[,2]=TT
mul_gene_list[,9]<-0
for(i in 1:nrow(mul_gene_list))
{
  gene_find1 <- mul_gene_list[i,2]
  index_find1 <- which(hypergraph[,1]==gene_find1) 
  if(length(index_find1)==0)
  {
    mul_gene_list[i,9] <- 0
    next()
  }
  mul_gene_list[i,9] <- as.numeric(hypergraph[index_find1,2])
}
mirna_exp <-mirna_EXP 
gene_exp <- matrix(data=0,nrow=nrow(mirna_exp),ncol=2)
gene_exp[,1] <- mirna_exp[,1]
rownames(mirna_exp) <- mirna_exp[,1]
mirna_exp <- mirna_exp[,-1]
mirna_col_name <- colnames(mirna_exp)
normal_expi <- grep('.11A|11B',mirna_col_name)
nori <- length(normal_expi)
ceri <- nori + 1
totali <- length(mirna_col_name)
for(i in 1:nrow(mirna_exp))
{
  m1 <- mean(as.numeric(mirna_exp[i,1:nori]))      
  m2 <- mean(as.numeric(mirna_exp[i,ceri:totali]))  
  var1 <- sd(as.numeric(mirna_exp[i,1:nori]))
  var2 <- sd(as.numeric(mirna_exp[i,ceri:totali]))
  var_sum <- var1^2+var2^2
  sum1 <- (m1-m2)^2/(4*var_sum)
  sum2 <-  0.5*log(var_sum/(2*var1*var2))
  gene_exp[i,2] <- sum1+sum2
}
dele_index <- which(gene_exp[,2]=='Inf')
gene_exp[dele_index,2] <- 0



mirna_list1<- mirna_List
rna_list <- which(mirna_list1[,3]==1)
mirna_list <- matrix(data=0,nrow = length(rna_list),ncol = 3)
for (i in 1:length(rna_list)) {
  index1 <- rna_list[i]
  for (j in 1:3) {
    mirna_list[i,j] <- mirna_list1[index1,j]
  }
}


mul_gene_list[,10] <- 0
for(i in 1:nrow(mul_gene_list))
{
  have_gene <- which(mirna_list[,2]==mul_gene_list[i,2])
  mul_gene_list[i,10] <- length(have_gene)
}
mul_gene_list[,11] <- 0
for (i in 1:nrow(mul_gene_list)) {
  if(mul_gene_list[i,10]==0)next()
  now_gene <- mul_gene_list[i,2] 
  mirna_index <- which(mirna_list[,2]==now_gene) 
  mirna_value <- 0
  weight <- 0
  cont=0
  for (j in 1:length(mirna_index)) {  
    index1 <- mirna_index[j]
    find_mirna <- mirna_list[index1,1]  
    weight <- mirna_list[index1,3]
    if(weight==0)next()
    index2 <- which(gene_exp[,1]==find_mirna)
    if(length(index2)==0)next()
    if(gene_exp[index2,2]==0)next()
    mirna_value1 <- as.numeric(gene_exp[index2,2])
    cont <- cont+1
    mirna_value <- mirna_value + mirna_value1   
  }
  if(cont==0)cont=1
  mul_gene_list[i,11] <- mirna_value   
}




mul_gene_list[,7][is.nan(mul_gene_list[,7])]<-0    
mul_gene_list[,6][is.nan(mul_gene_list[,6])]<-0    
mul_gene_list[,8][is.nan(mul_gene_list[,8])]<-0    


list<-matrix(0,nrow(mul_gene_list),1)
for(i in 1:nrow(mul_gene_list))
{

  list[i,1] <- mul_gene_list[i,6]*mul_gene_list[i,9]
}
txt<-matrix(0,10389,1)
txt[,1]<-mul_gene_list[,2]
write.table(txt,'genename.txt',quote=FALSE, row.names=FALSE,col.names = FALSE)
write.table(list,'CMCC.txt',quote=FALSE, row.names=FALSE,col.names = FALSE)
write.table(mul_gene_list[,7],'mrna_score.txt',quote=FALSE, row.names=FALSE,col.names = FALSE)
write.table(mul_gene_list[,11],'mirna_score.txt',quote=FALSE, row.names=FALSE,col.names = FALSE)



