library(dplyr)
info <- read.csv("D:/360Downloads/tianchi_mobile_recommend_train_user/tianchi_mobile_recommend_train_user.csv")
matrinfo <- as.matrix(info)
#for (i in c(1:nrow(matrinfo)))
#{
# S <- strsplit(matrinfo[i,6]," ")
# S1 <- as.data.frame(S,stringsAsFactors=F)
#  matrinfo[i,6] <- S1[1]
#}
#days <- unique(matrinfo[,6])
#N <- matrix(0,length(days),1)
#for (i in c(1:length(days)))
#{
#  bydayinfo <- filter(info,behavior_type="4",time==days[i])
#  N[i] <- nrow(bydayinfo)
#}
#plot(c(1:length(days)),N)
N <- 3800
#-----the above is to predict the possible trade number on the new coming day
user <- unique(matrinfo[,1])
item <- unique(matrinfo[,2])
userAttr <- matrix(0,length(user),length(item))
rownames(userAttr) <- unlist(user)
colnames(userAttr) <- unlist(item)
for (i in c(1:length(user)))
{
  indiviInfo <- filter(info,user_id==user[i])
  itemInfo <- unique(as.matrix(indiviInfo[,2]))
  for (j in c(1:length(itemInfo)))
  {
    behaveInfo <- filter(indiviInfo,item_id==itemInfo[j])
    behavetype <- max(as.matrix(behaveInfo[,3]))
    if (behavetype=="1")
    {
      userAttr[i,which(item==itemInfo[j],T)] <- 0.3
    }
    else if (behavetype=="2")
    {
      userAttr[i,which(item==itemInfo[j],T)] <- 0.5
    }
    else if (behavetype=="3")
    {
      userAttr[i,which(item==itemInfo[j],T)] <- 0.8
    }
    else 
    {
      userAttr[i,which(item==itemInfo[j],T)] <- 1.0
    }
  }
}
message(1)
userSimi <- matrix(0,nrow(userAttr),nrow(userAttr))
for (i in c(1:nrow(userAttr)))
{
  for (j in c(1:nrow(userAttr)))
  {
    userSimi[i,j] <- (userAttr[i,]%*%userAttr[j,])/(nrom(userAttr[i,],type="2")*nrom(userAttr[j,],type="2"))
  }
}
message(2)
#following are two functions for clustering
estimateNumberOfClustersGivenGraph <- function(W, NUMC=2:1000) {
  
  #   This function estimates the number of clusters given the two huristics
  #   given in the supplementary materials of our nature method paper
  #   W is the similarity graph
  #   NUMC is a vector which contains the possible choices of number of
  #   clusters.
  #   
  #   
  #   K1 is the estimated best number of clusters according to eigen-gaps
  #   K12 is the estimated SECOND best number of clusters according to eigen-gaps
  #   
  #   K2 is the estimated number of clusters according to rotation cost
  #   K22 is the estimated SECOND number of clusters according to rotation cost
  #   
  #   an example would be [K1, K2, K12,K22] = Estimate_Number_of_Clusters_given_graph(W,
  #                                                                                   [2:5]);
  #   
  #   Note that this function can only give an estimate of the number of
  #   clusters. How to determine the "OPTIMAL" number of clusters, is still an
  #   open question so far. 
  
  if (min(NUMC) == 1) {
    warning('Note that we always assume there are more than one cluster.');
    NUMC = NUMC[NUMC > 1]  
  }
  
  W = (W + t(W))/2
  diag(W) = 0
  
  if (length(NUMC) > 0) {
    degs = rowSums(W)
    
    
    # compute unnormalized Laplacian
    
    degs[degs == 0] = .Machine$double.eps    
    D = diag(degs)    
    L = D - W
    Di = diag(1 / sqrt(degs))
    L = Di %*% L %*% Di
    
    # compute the eigenvectors corresponding to the k smallest
    # eigs$valuess
    eigs = eigen(L)
    eigs_order = sort(eigs$values, index.return=T)$ix
    eigs$values = eigs$values[eigs_order]
    eigs$vectors = eigs$vectors[, eigs_order]
    eigengap = abs(diff(eigs$values))
    eigengap = eigengap * (1 - eigs$values[1:length(eigs$values) - 1] ) / (1 - eigs$values[2:length(eigs$values)])
    
    quality = list()
    for (c_index in 1:length(NUMC)) {
      ck = NUMC[c_index]
      UU = eigs$vectors[, 1:ck]
      EigenvectorsDiscrete <- .discretisation(UU)[[1]]
      EigenVectors = EigenvectorsDiscrete^2
      
      # MATLAB: sort(EigenVectors,2, 'descend');
      temp1 <- EigenVectors[do.call(order, lapply(1:ncol(EigenVectors), function(i) EigenVectors[, i])), ]
      temp1 <- t(apply(temp1, 1, sort, TRUE))  
      
      quality[[c_index]] = (1 - eigs$values[ck + 1]) / (1 - eigs$values[ck]) * 
        sum( sum( diag(1 / (temp1[, 1] + .Machine$double.eps) ) %*% temp1[, 1:max(2, ck-1)] ))
    }
    
    t1 <- sort(eigengap[NUMC], decreasing=TRUE, index.return=T)$ix
    K1 = NUMC[t1[1]]
    K12 = NUMC[t1[2]]
    t2 <- sort(unlist(quality), index.return=TRUE)$ix
    K2 <- NUMC[t2[1]]
    K22 <- NUMC[t2[2]]    
  }
  
  return (list(K1, K12, K2, K22))
}
spectralClustering <- function(affinity, K, type=3) {
  
  ###This function implements the famous spectral clustering algorithms. There are three variants. The default one is the third type. 
  ###THe inputs are as follows:
  
  #affinity: the similarity matrix;
  #K: the number of clusters
  # type: indicators of variants of spectral clustering 
  
  d = rowSums(affinity)
  d[d == 0] = .Machine$double.eps
  D = diag(d)
  L = D - affinity
  if (type == 1) {
    NL = L
  } else if (type == 2) {
    Di = diag(1 / d)
    NL = Di %*% L
  } else if(type == 3) {
    Di = diag(1 / sqrt(d))
    NL = Di %*% L %*% Di
  }
  eig = eigen(NL)
  res = sort(abs(eig$values),index.return = TRUE)
  U = eig$vectors[,res$ix[1:K]]
  normalize <- function(x) x / sqrt(sum(x^2))
  if (type == 3) {
    U = t(apply(U,1,normalize))
  }
  eigDiscrete = .discretisation(U)
  eigDiscrete = eigDiscrete$discrete
  labels = apply(eigDiscrete,1,which.max)
  
  
  
  return(labels)
}
K <- estimateNumberOfClustersGivenGraph(userSimi,NUMC=2:1000)
K <- as.data.frame(K)
K <- as.matrix(K)
C <- K[1]
group <- spectralClustering(userSimi,C,type=3)
names(group) <- unlist(user)
message(3)


#the next is attributes of items for possibility calculation
vectByItem <- matrinfo[,3]
for (i in c(1:length(vectByItem)))
{
  if (vectByItem[i]=="1")
  {
    vectByItem[i] <- 0.3
  }
  else if (vectByItem[i]=="2")
  {
    vectByItem[i] <- 0.5
  }
  else if (vectByItem[i]=="3")
  {
    vectByItem[i] <- 0.7
  }
  else 
  {
    vectByItem[i] <- 1.0
  }
}
message(4)
vectByUser <- vectByItem
possiSum <- sum(vectByItem)
names(vectByItem) <- unlist(matrinfo[,2])
names(vectByUser) <- unlist(matrinfo[,1])
possiItem <- matrix(0,length(item),1)
names(possiItem) <- unlist(item)
possiUser <- matrix(0,length(user),1)
names(possiUser) <- unlist(user)
for (i in c(1:length(item)))
{
  possiItem[i] <- sum(vectByItem[which(rownames(vectByItem)==item[i],T)])/possisum
}
for (i in c(1:length(user)))
{
  possiUser[i] <- sum(vectByUser[which(rownames(vectByUser)==user[i],T)])/possisum
}
message(5)
recom <- matrix(0,(length(user)*3),3)
k <- 0
for (i in c(1:length(user)))
{
  subUser <- user[which(names(group)==user[i],T)]
  info1 <- filter(info,user_id %in% subUser,behavior_type=="4")
  if (nrow(info1)>0)
  {
    item1 <- unique(as.matrix(info1[,2]))
  }
  info2 <- filter(info1,user_id==user[i])
  differ <- setdiff(item1,as.matrix(info2[,2]))
  if (length(differ)>0)
  {
    differ1 <- matrix(0,length(differ),1)
    names(differ1) <- unlist(differ)
    for (j in c(1:nrow(differ1)))
    {
      differ1[j,1] <- possiItem[which(names(possiItem)==differ[i],T)]
    }
    differ1 <- sort(differ1,decreasing=T)
    if (length(differ)>=3)
    {
      recom[(k+1):(k+3),1] <- user[i]
      recom[(k+1):(k+3),2] <- names(differ1)[1:3]
      recom[(k+1):(k+3),3] <- differ1[1:3]
      k <- k+3
    }
    else
    {
      recom[(k+1):(k+length(differ)),1] <- user[i]
      recom[(k+1):(k+length(differ)),2] <- names(differ1)
      recom[(k+1):(k+length(differ)),3] <- differ1
      k <- k+length(differ)
    }
  }
}
recom <- recom[1:k,]
message(6)
for (i in c(1:nrow(recom)))
{
  recom[i,3] <- recom[i,3]*possiUser[which(names(possiUser)==recom[i,1],t)]
}
frequ <- recom[,3]
names(frequ) <- unlist(c(1:nrow(recom)))
frequ <- sort(frequ,decreasing=T)
posi <- names(frequ)[1:N]
finalRecom <- matrix(0,N,2)
for (i in c(1:N))
{
  finalRecom[i,1] <- recom[posi[i],1]
  finalRecom[i,2] <- recom[posi[i],2]
}
message("the end")