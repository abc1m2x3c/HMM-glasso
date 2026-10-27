
# ----------------------------------------
# Random graph generator: use of simulateGraph
# ----------------------------------------


Graph.gen<-function(mod,p){
  if (mod==1){
    # p=4
    #write the neighborhood matrix
    Gr1.NB=matrix(0,p,p)
    Gr2.NB=matrix(0,p,p)
    Gr1.NB[2,3]=1 #makes corresponding elems in upper triangle to be 1
    Gr1.NB[1,4]=1
    Gr1.NB[lower.tri(Gr1.NB)]=t(Gr1.NB)[lower.tri(Gr1.NB)]
    Gr2.NB[3,4]=1
    Gr2.NB[lower.tri(Gr2.NB)]=t(Gr2.NB)[lower.tri(Gr2.NB)]
    Gr1.C=diag(1,p,p)
    Gr2.C=diag(1,p,p)
    Gr1.C[as.logical(Gr1.NB)]=0.7
    Gr2.C[as.logical(Gr2.NB)]=0.8
    Gr1.P=solve(Gr1.C)
    Gr2.P=solve(Gr2.C)
  }
  if (mod==2){
    # p=4
    #write the neighborhood matrix
    Gr1.NB=matrix(0,p,p)
    Gr2.NB=matrix(0,p,p)
    Gr1.NB[2,3]=1 #makes corresponding elems in upper triangle to be 1
    Gr1.NB[lower.tri(Gr1.NB)]=t(Gr1.NB)[lower.tri(Gr1.NB)]
    Gr2.NB[3,4]=1
    Gr2.NB[lower.tri(Gr2.NB)]=t(Gr2.NB)[lower.tri(Gr2.NB)]
    Gr1.C=diag(1,p,p)
    Gr2.C=diag(1,p,p)
    Gr1.C[as.logical(Gr1.NB)]=0.7
    Gr2.C[as.logical(Gr2.NB)]=0.8
    Gr1.P=solve(Gr1.C)
    Gr2.P=solve(Gr2.C)
  }
  if (mod==3){
    # p=2
    #write the neighborhood matrix
    Gr1.NB=matrix(0,p,p)
    Gr2.NB=matrix(0,p,p)
    Gr1.NB[1,2]=1 #makes corresponding elems in upper triangle to be 1
    Gr1.NB[lower.tri(Gr1.NB)]=t(Gr1.NB)[lower.tri(Gr1.NB)]
    Gr2.NB[lower.tri(Gr2.NB)]=t(Gr2.NB)[lower.tri(Gr2.NB)]
    Gr1.C=diag(1,p,p)
    Gr2.C=diag(1,p,p)
    Gr1.C[as.logical(Gr1.NB)]=0.7
    Gr1.P=solve(Gr1.C)
    Gr2.P=solve(Gr2.C)
  }
  if (mod==4){
    # p=1
    #write the neighborhood matrix
    Gr1.NB=matrix(0,p,p)
    Gr2.NB=matrix(0,p,p)
    # Gr1.NB[1,2]=1 #makes corresponding elems in upper triangle to be 1
    # Gr1.NB[lower.tri(Gr1.NB)]=t(Gr1.NB)[lower.tri(Gr1.NB)]
    # Gr2.NB[lower.tri(Gr2.NB)]=t(Gr2.NB)[lower.tri(Gr2.NB)]
    Gr1.C=diag(1,p,p)
    Gr2.C=diag(1,p,p)
    # Gr1.C[as.logical(Gr1.NB)]=0.7
    Gr1.P=solve(Gr1.C)
    Gr2.P=solve(Gr2.C)
  }
  if (mod==5){ # chain vs chain
    set.seed(123) # set seed in network generating for same structure 
    # p<-90
    precision1<- matrix(rep(0,p*p),nrow=p)
    
    cov1<-cov2<-matrix(NA,nrow=p,ncol=p)
    
    a<-0.1
    diag(cov1)<-1
    s<-rep(0,nrow(precision1))
    s[1]<-1
    for (i in 2:length(s))
    {
      s[i]<-s[i-1]+runif(1,0.5,1)
    }
    for (i in 1:(nrow(precision1)-1))
    {
      for (j in (i+1):nrow(precision1))
      {
        cov1[i,j]<-exp(-a*abs(s[j]-s[i]))
        cov1[j,i]<-cov1[i,j]
      }
    }
    precision1<-solve(cov1)
    Gr1.C=cov1
    Gr1.P=precision1
    lower<-precision1*lower.tri(precision1)
    loc<-which(abs(lower)<(10*10^(-4)),arr.ind=TRUE)
    loc1<-which(abs(lower)>(10*10^(-4)),arr.ind=TRUE)
    for (i in 1:nrow(loc))
    {
      lower[loc[i,1],loc[i,2]]<-0
    }
    for (i in 1:nrow(loc1))
    {
      lower[loc1[i,1],loc1[i,2]]<-1
    }
    Gr1.NB<-lower+t(lower)
    # for same networks
    precision2<-precision1
    precision2[2,3]<-precision2[3,2]<-precision2[3,4]<-precision2[4,3]<-0
    cov2<-solve(precision2)
    Gr2.C=cov2
    Gr2.P=precision2
    lower<-precision2*lower.tri(precision2)
    loc<-which(abs(lower)<(10*10^(-4)),arr.ind=TRUE)
    loc1<-which(abs(lower)>(10*10^(-4)),arr.ind=TRUE)
    for (i in 1:nrow(loc))
    {
      lower[loc[i,1],loc[i,2]]<-0
    }
    for (i in 1:nrow(loc1))
    {
      lower[loc1[i,1],loc1[i,2]]<-1
    }
    Gr2.NB<-lower+t(lower)
  }
  if (mod==6){ #nearest neighbor vs nearest neighbor
    set.seed(123) # set seed in network generating for same structure 
    
    # p<-90
    MaxNodes<-4
    #group 1
    xpoints<-runif(p,0,1)
    ypoints<-runif(p,0,1)
    distMat<-matrix(NA,nrow=p,ncol=p)
    G0<-matrix(0,nrow=p,ncol=p)
    num<-rep(0,nrow(distMat))
    
    precision0<-matrix(0,nrow=p,ncol=p)
    
    # this factor is for the calculation of diagonal of precision0.
    #factor<-1
    # the following factor 1.2 was used to generate p=30 sample size 200 using R 3.6.1.
    #factor<-1.2
    # the following factor is used for p=30 for other sample sizes other than 200 using R Studio.
    #factor<-1.3
    # this is for 50 and 70 p
    factor<-1.5
    
    for (i in 1:(p-1))
    {
      
      
      NumCon<-sample.int(min(MaxNodes, (p-i)),1)
      
      for (j in (i+1):p)
      {
        distMat[i,j]<-(xpoints[i]-xpoints[j])^2+(ypoints[i]-ypoints[j])^2
      }
      dist<-distMat[i,(i+1):p]
      toSortDist<-dist
      sorted<-sort(toSortDist,decreasing=FALSE)
      if (num[i]<MaxNodes)
      {
        for (k in 1:NumCon)
        {
          loc<-which(distMat[i,(i+1):p]==sorted[k])
          if (num[i]<MaxNodes && num[loc+i]<MaxNodes)
          {
            num[c(i,(loc+i))]<-num[c(i,(loc+i))]+1
            G0[i,(loc+i)]<-G0[(loc+i),i]<-1
            samp<-sample(c(runif(1,0.5,1),runif(1,-1,-0.5)),1)
            precision0[i,(loc+i)]<-precision0[(loc+i),i]<-samp
          }
        }
      }
    }
    
    for (i in 1:p)
    {
      precision0[i,i]<-factor*sum(abs(precision0[i,]))
    }
    
    for (i in 1:p)
    {
      if (precision0[i,i]!=0){
        precision0[i,i:p]<-precision0[i,i:p]/precision0[i,i]  
      }#else just keep that row all 0's
    }
    
    upper<-precision0*upper.tri(precision0)
    
    precision0<-upper+t(upper)+diag(1,p)
    Gr1.NB=G0
    Gr1.P=precision0
    Gr1.C=solve(precision0)
    
    #group 2
    xpoints<-runif(p,0,1)
    ypoints<-runif(p,0,1)
    distMat<-matrix(NA,nrow=p,ncol=p)
    G0<-matrix(0,nrow=p,ncol=p)
    num<-rep(0,nrow(distMat))
    
    precision0<-matrix(0,nrow=p,ncol=p)
    
    # this factor is for the calculation of diagonal of precision0.
    #factor<-1
    # the following factor 1.2 was used to generate p=30 sample size 200 using R 3.6.1.
    #factor<-1.2
    # the following factor is used for p=30 for other sample sizes other than 200 using R Studio.
    #factor<-1.3
    # this is for 50 and 70 p
    factor<-1.5
    
    for (i in 1:(p-1))
    {
      
      
      NumCon<-sample.int(min(MaxNodes, (p-i)),1)
      
      for (j in (i+1):p)
      {
        distMat[i,j]<-(xpoints[i]-xpoints[j])^2+(ypoints[i]-ypoints[j])^2
      }
      dist<-distMat[i,(i+1):p]
      toSortDist<-dist
      sorted<-sort(toSortDist,decreasing=FALSE)
      if (num[i]<MaxNodes)
      {
        for (k in 1:NumCon)
        {
          loc<-which(distMat[i,(i+1):p]==sorted[k])
          if (num[i]<MaxNodes && num[loc+i]<MaxNodes)
          {
            num[c(i,(loc+i))]<-num[c(i,(loc+i))]+1
            G0[i,(loc+i)]<-G0[(loc+i),i]<-1
            samp<-sample(c(runif(1,0.5,1),runif(1,-1,-0.5)),1)
            precision0[i,(loc+i)]<-precision0[(loc+i),i]<-samp
          }
        }
      }
    }
    
    for (i in 1:p)
    {
      precision0[i,i]<-factor*sum(abs(precision0[i,]))
    }
    
    for (i in 1:p)
    {
      if (precision0[i,i]!=0){
        precision0[i,i:p]<-precision0[i,i:p]/precision0[i,i]  
      }#else just keep that row all 0's
    }
    
    upper<-precision0*upper.tri(precision0)
    
    precision0<-upper+t(upper)+diag(1,p)
    Gr2.NB=G0
    Gr2.P=precision0
    Gr2.C=solve(precision0) 
  }
  if (mod==7){ #chain vs nearest neighbor
    set.seed(123) # set seed in network generating for same structure 
    MaxNodes<-4
    
    # p<-90
    precision1<- matrix(rep(0,p*p),nrow=p)
    
    cov1<-cov2<-matrix(NA,nrow=p,ncol=p)
    
    a<-0.1
    diag(cov1)<-1
    s<-rep(0,nrow(precision1))
    s[1]<-1
    for (i in 2:length(s))
    {
      s[i]<-s[i-1]+runif(1,0.5,1)
    }
    for (i in 1:(nrow(precision1)-1))
    {
      for (j in (i+1):nrow(precision1))
      {
        cov1[i,j]<-exp(-a*abs(s[j]-s[i]))
        cov1[j,i]<-cov1[i,j]
      }
    }
    precision1<-solve(cov1)
    Gr1.C=cov1
    Gr1.P=precision1
    lower<-precision1*lower.tri(precision1)
    loc<-which(abs(lower)<(10*10^(-4)),arr.ind=TRUE)
    loc1<-which(abs(lower)>(10*10^(-4)),arr.ind=TRUE)
    for (i in 1:nrow(loc))
    {
      lower[loc[i,1],loc[i,2]]<-0
    }
    for (i in 1:nrow(loc1))
    {
      lower[loc1[i,1],loc1[i,2]]<-1
    }
    Gr1.NB<-lower+t(lower)
    
    
    
    xpoints<-runif(p,0,1)
    ypoints<-runif(p,0,1)
    distMat<-matrix(NA,nrow=p,ncol=p)
    G0<-matrix(0,nrow=p,ncol=p)
    num<-rep(0,nrow(distMat))
    
    precision0<-matrix(0,nrow=p,ncol=p)
    
    # this factor is for the calculation of diagonal of precision0.
    #factor<-1
    # the following factor 1.2 was used to generate p=30 sample size 200 using R 3.6.1.
    #factor<-1.2
    # the following factor is used for p=30 for other sample sizes other than 200 using R Studio.
    #factor<-1.3
    # this is for 50 and 70 p
    factor<-1.5
    
    for (i in 1:(p-1))
    {
      
      
      NumCon<-sample.int(min(MaxNodes, (p-i)),1)
      
      for (j in (i+1):p)
      {
        distMat[i,j]<-(xpoints[i]-xpoints[j])^2+(ypoints[i]-ypoints[j])^2
      }
      dist<-distMat[i,(i+1):p]
      toSortDist<-dist
      sorted<-sort(toSortDist,decreasing=FALSE)
      if (num[i]<MaxNodes)
      {
        for (k in 1:NumCon)
        {
          loc<-which(distMat[i,(i+1):p]==sorted[k])
          if (num[i]<MaxNodes && num[loc+i]<MaxNodes)
          {
            num[c(i,(loc+i))]<-num[c(i,(loc+i))]+1
            G0[i,(loc+i)]<-G0[(loc+i),i]<-1
            samp<-sample(c(runif(1,0.5,1),runif(1,-1,-0.5)),1)
            precision0[i,(loc+i)]<-precision0[(loc+i),i]<-samp
          }
        }
      }
    }
    
    for (i in 1:p)
    {
      precision0[i,i]<-factor*sum(abs(precision0[i,]))
    }
    
    for (i in 1:p)
    {
      if (precision0[i,i]!=0){
        precision0[i,i:p]<-precision0[i,i:p]/precision0[i,i]  
      }#else just keep that row all 0's
    }
    
    upper<-precision0*upper.tri(precision0)
    
    precision0<-upper+t(upper)+diag(1,p)
    Gr2.NB=G0
    Gr2.P=precision0
    Gr2.C=solve(precision0)
  }
  set.seed(NULL) #stop seed
  return (list(Gr1.C=Gr1.C,Gr2.C=Gr2.C))
}


Graph.gen.multiK<-function(mod,p){
  if (mod==5){ # chain vs chain
    set.seed(123) # set seed in network generating for same structure 
    # p<-90
    precision1<- matrix(rep(0,p*p),nrow=p)
    
    cov1<-cov2<-matrix(NA,nrow=p,ncol=p)
    
    a<-0.1
    diag(cov1)<-1
    s<-rep(0,nrow(precision1))
    s[1]<-1
    for (i in 2:length(s))
    {
      s[i]<-s[i-1]+runif(1,0.5,1)
    }
    for (i in 1:(nrow(precision1)-1))
    {
      for (j in (i+1):nrow(precision1))
      {
        cov1[i,j]<-exp(-a*abs(s[j]-s[i]))
        cov1[j,i]<-cov1[i,j]
      }
    }
    precision1<-solve(cov1)
    Gr1.C=cov1
    Gr1.P=precision1
    lower<-precision1*lower.tri(precision1)
    loc<-which(abs(lower)<(10*10^(-4)),arr.ind=TRUE)
    loc1<-which(abs(lower)>(10*10^(-4)),arr.ind=TRUE)
    for (i in 1:nrow(loc))
    {
      lower[loc[i,1],loc[i,2]]<-0
    }
    for (i in 1:nrow(loc1))
    {
      lower[loc1[i,1],loc1[i,2]]<-1
    }
    Gr1.NB<-lower+t(lower)
    # network 2
    precision2<-precision1
    precision2[2,3]<-precision2[3,2]<-precision2[3,4]<-precision2[4,3]<-0
    cov2<-solve(precision2)
    Gr2.C=cov2
    Gr2.P=precision2
    lower<-precision2*lower.tri(precision2)
    loc<-which(abs(lower)<(10*10^(-4)),arr.ind=TRUE)
    loc1<-which(abs(lower)>(10*10^(-4)),arr.ind=TRUE)
    for (i in 1:nrow(loc))
    {
      lower[loc[i,1],loc[i,2]]<-0
    }
    for (i in 1:nrow(loc1))
    {
      lower[loc1[i,1],loc1[i,2]]<-1
    }
    Gr2.NB<-lower+t(lower)
    # network 3
    precision3<-precision1
    precision3[5,6]<-precision3[6,5]<-precision3[7,8]<-precision3[8,7]<-0
    cov3<-solve(precision3)
    Gr3.C=cov3
    Gr3.P=precision3
    lower<-precision3*lower.tri(precision3)
    loc<-which(abs(lower)<(10*10^(-4)),arr.ind=TRUE)
    loc1<-which(abs(lower)>(10*10^(-4)),arr.ind=TRUE)
    for (i in 1:nrow(loc))
    {
      lower[loc[i,1],loc[i,2]]<-0
    }
    for (i in 1:nrow(loc1))
    {
      lower[loc1[i,1],loc1[i,2]]<-1
    }
    Gr3.NB<-lower+t(lower)
    # # network 4
    # precision4<-precision1
    # precision4[9,10]<-precision4[10,9]<-precision4[11,12]<-precision4[12,11]<-0
    # cov4<-solve(precision4)
    # Gr4.C=cov4
    # Gr4.P=precision4
    # lower<-precision4*lower.tri(precision4)
    # loc<-which(abs(lower)<(10*10^(-4)),arr.ind=TRUE)
    # loc1<-which(abs(lower)>(10*10^(-4)),arr.ind=TRUE)
    # for (i in 1:nrow(loc))
    # {
    #   lower[loc[i,1],loc[i,2]]<-0
    # }
    # for (i in 1:nrow(loc1))
    # {
    #   lower[loc1[i,1],loc1[i,2]]<-1
    # }
    # Gr4.NB<-lower+t(lower)
    # # network 5
    # precision5<-precision1
    # precision5[13,14]<-precision5[14,13]<-precision5[15,16]<-precision5[16,15]<-0
    # cov5<-solve(precision5)
    # Gr5.C=cov5
    # Gr5.P=precision5
    # lower<-precision5*lower.tri(precision5)
    # loc<-which(abs(lower)<(10*10^(-4)),arr.ind=TRUE)
    # loc1<-which(abs(lower)>(10*10^(-4)),arr.ind=TRUE)
    # for (i in 1:nrow(loc))
    # {
    #   lower[loc[i,1],loc[i,2]]<-0
    # }
    # for (i in 1:nrow(loc1))
    # {
    #   lower[loc1[i,1],loc1[i,2]]<-1
    # }
    # Gr5.NB<-lower+t(lower)
  }
  
  if (mod==6){ #nearest neighbor vs nearest neighbor
    set.seed(123) # set seed in network generating for same structure 
    
    # p<-90
    MaxNodes<-4
    #group 1
    xpoints<-runif(p,0,1)
    ypoints<-runif(p,0,1)
    distMat<-matrix(NA,nrow=p,ncol=p)
    G0<-matrix(0,nrow=p,ncol=p)
    num<-rep(0,nrow(distMat))
    
    precision0<-matrix(0,nrow=p,ncol=p)
    
    # this factor is for the calculation of diagonal of precision0.
    #factor<-1
    # the following factor 1.2 was used to generate p=30 sample size 200 using R 3.6.1.
    #factor<-1.2
    # the following factor is used for p=30 for other sample sizes other than 200 using R Studio.
    #factor<-1.3
    # this is for 50 and 70 p
    factor<-1.5
    
    for (i in 1:(p-1))
    {
      
      
      NumCon<-sample.int(min(MaxNodes, (p-i)),1)
      
      for (j in (i+1):p)
      {
        distMat[i,j]<-(xpoints[i]-xpoints[j])^2+(ypoints[i]-ypoints[j])^2
      }
      dist<-distMat[i,(i+1):p]
      toSortDist<-dist
      sorted<-sort(toSortDist,decreasing=FALSE)
      if (num[i]<MaxNodes)
      {
        for (k in 1:NumCon)
        {
          loc<-which(distMat[i,(i+1):p]==sorted[k])
          if (num[i]<MaxNodes && num[loc+i]<MaxNodes)
          {
            num[c(i,(loc+i))]<-num[c(i,(loc+i))]+1
            G0[i,(loc+i)]<-G0[(loc+i),i]<-1
            samp<-sample(c(runif(1,0.5,1),runif(1,-1,-0.5)),1)
            precision0[i,(loc+i)]<-precision0[(loc+i),i]<-samp
          }
        }
      }
    }
    
    for (i in 1:p)
    {
      precision0[i,i]<-factor*sum(abs(precision0[i,]))
    }
    
    for (i in 1:p)
    {
      if (precision0[i,i]!=0){
        precision0[i,i:p]<-precision0[i,i:p]/precision0[i,i]  
      }#else just keep that row all 0's
    }
    
    upper<-precision0*upper.tri(precision0)
    
    precision0<-upper+t(upper)+diag(1,p)
    Gr1.NB=G0
    Gr1.P=precision0
    Gr1.C=solve(precision0)
    
    # Generate Group 2 by modifying two edges in Group 1
    # Function to remove two edges and regenerate if not positive definite
    generate_modified_group <- function(Gr.NB, Gr.P) {
      repeat {
        GrNew.NB <- Gr.NB
        GrNew.P <- Gr.P
        
        removed_edges <- 0
        while (removed_edges < 2) {
          edge <- which(GrNew.NB == 1, arr.ind = TRUE)
          if (nrow(edge) == 0) break
          selected_edge <- edge[sample(1:nrow(edge), 1), ]
          GrNew.NB[selected_edge[1], selected_edge[2]] <- 0
          GrNew.NB[selected_edge[2], selected_edge[1]] <- 0
          GrNew.P[selected_edge[1], selected_edge[2]] <- 0
          GrNew.P[selected_edge[2], selected_edge[1]] <- 0
          removed_edges <- removed_edges + 1
        }
        
        if (all(eigen(GrNew.P, symmetric = TRUE)$values > 0)) {
          return(list(NB = GrNew.NB, P = GrNew.P, C = solve(GrNew.P)))
        }
      }
    }
    # Generate Group 2 by removing two edges from Group 1
    Gr2 <- generate_modified_group(Gr1.NB, Gr1.P)
    
    # Generate Group 3 by removing two edges from Group 2
    Gr3 <- generate_modified_group(Gr2$NB, Gr2$P)
    
    
    Gr2.C = Gr2$C
    Gr3.C = Gr3$C
    Gr2.NB = Gr2$NB
    Gr3.NB = Gr3$NB
  }
  
  
  if (mod==7){ #nearest neighbor vs nearest neighbor
    set.seed(123) # set seed in network generating for same structure 
    
    # p<-90
    MaxNodes<-4
    #group 1
    xpoints<-runif(p,0,1)
    ypoints<-runif(p,0,1)
    distMat<-matrix(NA,nrow=p,ncol=p)
    G0<-matrix(0,nrow=p,ncol=p)
    num<-rep(0,nrow(distMat))
    
    precision0<-matrix(0,nrow=p,ncol=p)
    
    # this factor is for the calculation of diagonal of precision0.
    #factor<-1
    # the following factor 1.2 was used to generate p=30 sample size 200 using R 3.6.1.
    #factor<-1.2
    # the following factor is used for p=30 for other sample sizes other than 200 using R Studio.
    #factor<-1.3
    # this is for 50 and 70 p
    factor<-1.5
    
    for (i in 1:(p-1))
    {
      
      
      NumCon<-sample.int(min(MaxNodes, (p-i)),1)
      
      for (j in (i+1):p)
      {
        distMat[i,j]<-(xpoints[i]-xpoints[j])^2+(ypoints[i]-ypoints[j])^2
      }
      dist<-distMat[i,(i+1):p]
      toSortDist<-dist
      sorted<-sort(toSortDist,decreasing=FALSE)
      if (num[i]<MaxNodes)
      {
        for (k in 1:NumCon)
        {
          loc<-which(distMat[i,(i+1):p]==sorted[k])
          if (num[i]<MaxNodes && num[loc+i]<MaxNodes)
          {
            num[c(i,(loc+i))]<-num[c(i,(loc+i))]+1
            G0[i,(loc+i)]<-G0[(loc+i),i]<-1
            samp<-sample(c(runif(1,0.5,1),runif(1,-1,-0.5)),1)
            precision0[i,(loc+i)]<-precision0[(loc+i),i]<-samp
          }
        }
      }
    }
    
    for (i in 1:p)
    {
      precision0[i,i]<-factor*sum(abs(precision0[i,]))
    }
    
    for (i in 1:p)
    {
      if (precision0[i,i]!=0){
        precision0[i,i:p]<-precision0[i,i:p]/precision0[i,i]  
      }#else just keep that row all 0's
    }
    
    upper<-precision0*upper.tri(precision0)
    
    precision0<-upper+t(upper)+diag(1,p)
    Gr1.NB=G0
    Gr1.P=precision0
    Gr1.C=solve(precision0)
    
    # Generate Group 2 by modifying two edges in Group 1
    # Function to remove two edges and regenerate if not positive definite
    generate_modified_group <- function(Gr.NB, Gr.P) {
      repeat {
        GrNew.NB <- Gr.NB
        GrNew.P <- Gr.P
        
        removed_edges <- 0
        while (removed_edges < 2) {
          edge <- which(GrNew.NB == 1, arr.ind = TRUE)
          if (nrow(edge) == 0) break
          selected_edge <- edge[sample(1:nrow(edge), 1), ]
          GrNew.NB[selected_edge[1], selected_edge[2]] <- 0
          GrNew.NB[selected_edge[2], selected_edge[1]] <- 0
          GrNew.P[selected_edge[1], selected_edge[2]] <- 0
          GrNew.P[selected_edge[2], selected_edge[1]] <- 0
          removed_edges <- removed_edges + 1
        }
        
        if (all(eigen(GrNew.P, symmetric = TRUE)$values > 0)) {
          return(list(NB = GrNew.NB, P = GrNew.P, C = solve(GrNew.P)))
        }
      }
    }
    # Generate Group 2 by removing two edges from Group 1
    Gr2 <- generate_modified_group(Gr1.NB, Gr1.P)

    # Generate Group 3 by removing two edges from Group 2
    Gr3 <- generate_modified_group(Gr2$NB, Gr2$P)
    
    # Gr2.C = 0.7*Gr2$C #the only difference with mod==6
    # Gr3.C = 0.4*Gr3$C
    # modify the above two rows on 202507
    Gr1.C = 0.1*Gr1.C
    Gr2.C = 0.1*Gr2$C 
    Gr3.C = 0.1*Gr3$C
    
    Gr2.NB = Gr2$NB
    Gr3.NB = Gr3$NB
  }
  
  
  if (mod==8){
    set.seed(123)
    
    nedges <- p-1  # desired number of edges
    
    # total possible undirected edges
    total_edges <- choose(p, 2)
    
    # sample exactly p edges (indices in upper-triangle)
    edge_indices <- sample(total_edges, nedges)
    
    # adjacency (binary)
    adj <- matrix(0, p, p)
    adj[upper.tri(adj)][edge_indices] <- 1
    Gr1.NB=adj + t(adj)
    
    
    # assign weights only once (upper-tri)
    weights <- runif(length(edge_indices), 0.3, 0.8)
    adj[upper.tri(adj)][edge_indices] <- weights
    
    # symmetrize
    adj <- adj + t(adj)
    

    # Precision matrix: make it diagonally dominant (so positive-definite)
    # smaller pad -> stronger partial correlations
    tau <- 0.01                    
    W <- adj
    precision <- -W
    diag(precision) <- rowSums(W) + tau
    
    # Covariance matrix
    cov1 <- solve(precision)
    Gr1.P=precision
    Gr1.C=cov1
    
    precision_to_partial <- function(Omega, tol = 1e-12, clamp = TRUE) {
      # Symmetrize defensively (numerical noise)
      if (!isSymmetric(Omega, tol = 1e-8)) {
        Omega <- 0.5 * (Omega + t(Omega))
      }
      
      d <- diag(Omega)
      if (any(d <= tol)) {
        stop("Diagonal entries must be positive (matrix should be positive definite).")
      }
      
      # Denominator sqrt(Omega_ii * Omega_jj)
      denom <- sqrt(outer(d, d))
      R <- -Omega / denom
      diag(R) <- 1
      
      # Optional: clamp tiny numerical overshoots to [-1, 1]
      if (clamp) {
        R <- pmin(pmax(R, -1), 1)
      }
      R
    }
    
    partial_corr = precision_to_partial(precision)
    
    # work only with upper triangle (to avoid duplicates, skip diagonal)
    ut <- upper.tri(partial_corr, diag = FALSE)
    vals <- abs(partial_corr[ut])
    
    # order indices by descending absolute value
    ord <- order(vals, decreasing = TRUE)[1:4]
    
    # convert back to (i,j) positions
    idx_mat <- which(ut, arr.ind = TRUE)
    top4 <- idx_mat[ord, , drop = FALSE]
    
    Gr2.NB=Gr1.NB
    Gr2.P=Gr1.P
    #remove the first and 4th strongest edges
    Gr2.NB[top4[1,1], top4[1,2]] <- 0
    Gr2.NB[top4[1,2], top4[1,1]] <- 0
    Gr2.P[top4[1,1], top4[1,2]] <- 0
    Gr2.P[top4[1,2], top4[1,1]] <- 0
    Gr2.NB[top4[4,1], top4[4,2]] <- 0
    Gr2.NB[top4[4,2], top4[4,1]] <- 0
    Gr2.P[top4[4,1], top4[4,2]] <- 0
    Gr2.P[top4[4,2], top4[4,1]] <- 0
    Gr2.C=solve(Gr2.P)
    
    Gr3.NB=Gr1.NB
    Gr3.P=Gr1.P
    #remove the first and 4th strongest edges
    Gr3.NB[top4[2,1], top4[2,2]] <- 0
    Gr3.NB[top4[2,2], top4[2,1]] <- 0
    Gr3.P[top4[2,1], top4[2,2]] <- 0
    Gr3.P[top4[2,2], top4[2,1]] <- 0
    Gr3.NB[top4[3,1], top4[3,2]] <- 0
    Gr3.NB[top4[3,2], top4[3,1]] <- 0
    Gr3.P[top4[3,1], top4[3,2]] <- 0
    Gr3.P[top4[3,2], top4[3,1]] <- 0
    Gr3.C=solve(Gr3.P)
    
    #check positive after revoming edges
    is.positive.definite <- function(M) {
      eigenvalues <- eigen(M, symmetric = TRUE)$values
      all(eigenvalues > 1e-8)   # allow small tolerance for numerical errors
    }

    if (!(is.positive.definite(Gr1.P) & is.positive.definite(Gr2.P) & is.positive.definite(Gr3.P)))
      stop('Not positive definite')
  }
  
  set.seed(NULL)
  return (list(Gr1.C=Gr1.C,Gr2.C=Gr2.C,Gr3.C=Gr3.C,Gr1.NB=Gr1.NB,Gr2.NB=Gr2.NB,Gr3.NB=Gr3.NB))
}


EM.glasso <- function(X,K,method,rho=NULL,sigma.true=NULL,penalize.diagonal=FALSE,epsilon.EM=1e-8,max.iter.EM=1000){
  if (method==1){
    if (dim(sigma.true)[3]!=K) return(NA)
  }  
  p=ncol(X)
  n=nrow(X)
  # Initial values
  tau <- rep(1/K,K)
  
  # Initial random assignment for mu and sigma
  sample.start <- sample(1:K, n, replace=TRUE)
  mu <- t(sapply(1:K, function(k) colMeans(X[sample.start == k,])))
  
  sigma <- array(dim=c(p,p,K))
  for (k in 1:K) {
    sigma[,,k] <- var(X[sample.start == k, ])
  }
  
  
  # EM Algorithm
  logllk.pre <- -Inf
  count <- 0
  
  
  
  repeat {
    count <- count + 1
    
    # E-step: Calculate responsibilities
    f <- array(dim=c(n,K))
    for (k in 1:K){
      f[,k] <- dmvnorm(X, mean=mu[k,], sigma=sigma[,,k]) 
    }
    
    f_weighted <- sweep(f,2,tau,FUN='*')
    T <- sweep(f_weighted, 1, rowSums(f_weighted), FUN="/")
    
    # M-step: Update parameters
    if (method==1) sigma.true.temp=sigma.true
    for (k in 1:K) {
      tau[k] <- mean(T[,k]) # or sum(T[,k])/sum(T) where sum(T)=1000
      mu[k,] <- colSums(T[,k] * X) / sum(T[,k])
      
      X.centered <- sweep(X, 2, mu[k,])
      sigma.temp <- t(X.centered) %*% (X.centered * T[,k]) / sum(T[,k])
      if (method==2){
        sigma[,,k] <- sigma.temp
      }else if(method==10){
        sigma[,,k] <- glasso(sigma.temp,rho=rho.method10,penalize.diagonal=penalize.diagonal)$w
      }else if(method==11){
        sigma[,,k] <- glasso(sigma.temp,rho=rho.method11,penalize.diagonal=penalize.diagonal)$w
      }else if(method==1){
        distances <- sapply(1:dim(sigma.true.temp)[3], function(j) norm(sigma.temp - sigma.true.temp[,,j], type = "F"))
        sigma[,,k] <- sigma.true.temp[,,which.min(distances)]
        sigma.true.temp <- sigma.true.temp[,,-which.min(distances),drop=FALSE] # remove the layer that has closest distance 
        # in case of repeated use the same sigma
      }
    }
    
    
    # Log-likelihood calculation
    logllk.curr <- sum(log(rowSums(sweep(f,2,tau,FUN='*'))))
    
    # print(paste('Method: ',method,' Iteration: ',count,sep=''))
    # print(paste('Current Q: ',round(logllk.curr,2),'; Current llk ', round(logllk.true.curr,2),sep=''))
    # print(paste('Current Q: ',round(logllk.curr,2),sep=''))
    
    # Convergence check
    if (abs(logllk.curr - logllk.pre) < epsilon.EM || count > max.iter.EM) {
      if (count > max.iter.EM) warning("Maximum iteration reached")
      break
    }
    
    # Prepare for next iteration
    logllk.pre <- logllk.curr
    
  }
  
  if (method==2){
    num_param=K*ncol(X) + #mu
      K*0.5*ncol(X)*(ncol(X)+1) + # sigma
      K-1  #tau
    BIC=-2*logllk.curr+n*log(num_param)
  }
  if (method==10 | method==11){
    num_param=K*ncol(X) + #mu
      K*0.5*ncol(X)*(ncol(X)+1) + # sigma
      K-1 + #tau
      1 #rho
    BIC=-2*logllk.curr+num_param*log(n)
  }
  
  return (list(sigma=sigma,mu=mu,logllk=logllk.curr,tau=tau,T=T,BIC=BIC))
}



cv.cluster_num<-function(X,method,n_folds=5,K_values=1:6){
  # Define the range of K to test
  
  log_likelihood_values <- vector("list", max(K_values))
  names(log_likelihood_values) <- as.character(K_values)
  
  # Number of folds for cross-validation
  
  
  # Create folds
  folds <- sample(rep(1:n_folds, length.out = nrow(X)))
  
  for (k in K_values) {
    log_likelihoods <- numeric(n_folds)
    
    for (i in 1:n_folds) {
      # Split data into training and validation
      train_set <- X[folds != i, , drop = FALSE]
      valid_set <- X[folds == i, , drop = FALSE]
      
      # Fit model on training data
      fit <- EM.glasso(train_set, K = k, method = method)  # Adjust method as needed
      
      # Calculate log-likelihood on validation data
      mu <- fit$mu
      sigma <- fit$sigma
      tau <- fit$tau  # Example; adjust if T is not structured as [samples, clusters]
      
      log_likelihoods[i] <- sum(sapply(1:nrow(valid_set), function(j) {
        log(sum(sapply(1:k, function(m) {
          dmvnorm(valid_set[j, , drop = FALSE], mean = mu[m, ], sigma = sigma[, , m]) * tau[m]
        })))
      }))
    }
    
    log_likelihood_values[[as.character(k)]] <- log_likelihoods
  }
  
  # Calculate average log-likelihood for each K across all folds
  average_log_likelihood <- sapply(log_likelihood_values, mean)
  # optimal_K <- which.max(average_log_likelihood)
  return (average_log_likelihood)
}


# Function to predict new data points membership
##not verified
predict_gmm <- function(X_new, fitted.model) {
  mu=fitted.model$mu
  sigma=fitted.model$sigma
  tau=fitted.model$tau
  n_new <- nrow(X_new)
  K <- length(tau)
  
  # Probability densities and responsibilities (posterior probabilities)
  f_new <- matrix(0, nrow = n_new, ncol = K)
  for (k in 1:K) {
    f_new[, k] <- dmvnorm(X_new, mean = mu[k, ], sigma = sigma[, , k])
  }
  
  # Weight by the mixing probabilities
  weighted_f_new <- sweep(f_new, 2, tau, FUN = "*")
  
  # Normalize to get responsibilities
  T_new <- sweep(weighted_f_new, 1, rowSums(weighted_f_new), FUN = "/")
  
  # Assign to cluster with highest responsibility
  cluster_assignments <- max.col(T_new, ties.method = "first")
  
  return(cluster_assignments)
}


max.accu <- function(mem.mat,true_membership){
  n=nrow(mem.mat)
  predicted_clusters <- max.col(mem.mat, ties.method = "first")
  permutations <- permuteGeneral(length(unique(predicted_clusters)), length(unique(predicted_clusters)))
  contingency_table <- table(predicted_clusters, true_membership)
  if (length(unique(apply(contingency_table,1,which.max)))==1) warning('All clusters belong to the same true cluster')
  
  # Initialize best accuracy
  best_accuracy <- 0
  best_permutation <- NULL
  # Iterate over all permutations to find the best alignment
  for (i in 1:nrow(permutations)) {
    # Map the predicted clusters to true clusters using the current permutation
    perm <- permutations[i, ]
    
    # Permuted contingency table
    permuted_table <- contingency_table[, perm]
    
    # Calculate accuracy for this permutation
    accuracy <- sum(diag(permuted_table)) / n
    
    # Update best accuracy
    if (accuracy > best_accuracy) {
      best_accuracy <- accuracy
      best_permutation <- perm
    }
  }

  # Return the best accuracy and the corresponding permutation
  return(list(
    best_accuracy = best_accuracy,
    best_permutation = best_permutation
  ))
  
  # accuracy <- sum(apply(contingency_table,1,max)) / n
}


# in max.accu, we exhaust all permutations to find the best matches for pred and true membership
# however, if the number of predicted clusters is different from the number of true membership, this won't work
# instead, just assign the predicted cluster to the true cluster which has highest overlap counts
max.accu2 <- function(mem.mat,true_membership){
  n=nrow(mem.mat)
  predicted_clusters <- max.col(mem.mat, ties.method = "first")
  # permutations <- permuteGeneral(length(unique(predicted_clusters)), length(unique(predicted_clusters)))
  contingency_table <- table(predicted_clusters, true_membership)
  if (length(unique(apply(contingency_table,1,which.max)))==1) warning('All clusters belong to the same true cluster')
  best_permutation <- apply(contingency_table,1,which.max)
  # Iterate over all permutations to find the best alignment
  
  
  
  # Permuted contingency table
  permuted_table <- contingency_table[, best_permutation]
  
  # Calculate accuracy for this permutation
  best_accuracy <- sum(diag(permuted_table)) / n
  
  
  # Return the best accuracy and the corresponding permutation
  return(list(
    best_accuracy = best_accuracy,
    best_permutation = best_permutation
  ))
  
  # accuracy <- sum(apply(contingency_table,1,max)) / n
}


EM.glasso.MergeSigma <- function(X,K,method,rho=NULL,sigma.true=NULL,penalize.diagonal=FALSE,epsilon.EM=1e-8,max.iter.EM=1000,sigma_merge_threshold = 1e-3){
  if (method==1){
    if (dim(sigma.true)[3]!=K) return(NA)
  }  
  p=ncol(X)
  n=nrow(X)
  # Initial values
  tau <- rep(1/K,K)
  
  # Initial random assignment for mu and sigma
  sample.start <- sample(1:K, n, replace=TRUE)
  mu <- t(sapply(1:K, function(k) colMeans(X[sample.start == k,])))
  
  sigma <- array(dim=c(p,p,K))
  for (k in 1:K) {
    sigma[,,k] <- var(X[sample.start == k, ])
  }
  
  
  # EM Algorithm
  logllk.pre <- -Inf
  count <- 0
  
  
  
  repeat {
    count <- count + 1
    
    # E-step: Calculate responsibilities
    f <- array(dim=c(n,K))
    for (k in 1:K){
      f[,k] <- dmvnorm(X, mean=mu[k,], sigma=sigma[,,k]) 
    }
    
    f_weighted <- sweep(f,2,tau,FUN='*')
    T <- sweep(f_weighted, 1, rowSums(f_weighted), FUN="/")
    
    # M-step: Update parameters
    if (method==1) sigma.true.temp=sigma.true
    for (k in 1:K) {
      tau[k] <- mean(T[,k]) # or sum(T[,k])/sum(T) where sum(T)=1000
      mu[k,] <- colSums(T[,k] * X) / sum(T[,k])
      
      X.centered <- sweep(X, 2, mu[k,])
      sigma.temp <- t(X.centered) %*% (X.centered * T[,k]) / sum(T[,k])
      if (method==2){
        sigma[,,k] <- sigma.temp
      }else if(method==10){
        sigma[,,k] <- glasso(sigma.temp,rho=rho.method10,penalize.diagonal=penalize.diagonal)$w
      }else if(method==11){
        sigma[,,k] <- glasso(sigma.temp,rho=rho.method11,penalize.diagonal=penalize.diagonal)$w
      }else if(method==1){
        distances <- sapply(1:dim(sigma.true.temp)[3], function(j) norm(sigma.temp - sigma.true.temp[,,j], type = "F"))
        sigma[,,k] <- sigma.true.temp[,,which.min(distances)]
        sigma.true.temp <- sigma.true.temp[,,-which.min(distances),drop=FALSE] # remove the layer that has closest distance 
        # in case of repeated use the same sigma
      }
    }
    
    merged <- rep(FALSE, K)
    for (k in 1:(K-1)) {
      if (!merged[k]) {
        for (l in (k+1):K) {
          if (!merged[l]) {
            distance <- norm(sigma[,,k] - sigma[,,l], type = "F")
            if (distance < sigma_merge_threshold){
              # Merge cluster l into cluster k
              sigma[,,k] <- (sigma[,,k] * tau[k] + sigma[,,l] * tau[l]) / (tau[k] + tau[l])
              merged[l] <- TRUE
            }
          }
        }
      }
    }
    
    
    # Log-likelihood calculation
    logllk.curr <- sum(log(rowSums(sweep(f,2,tau,FUN='*'))))
    
    print(paste('Method: ',method,' Iteration: ',count,sep=''))
    # print(paste('Current Q: ',round(logllk.curr,2),'; Current llk ', round(logllk.true.curr,2),sep=''))
    print(paste('Current Q: ',round(logllk.curr,2),sep=''))
    
    # Convergence check
    if (abs(logllk.curr - logllk.pre) < epsilon.EM || count > max.iter.EM) {
      if (count > max.iter.EM) warning("Maximum iteration reached")
      break
    }
    
    # Prepare for next iteration
    logllk.pre <- logllk.curr
    
  }
  
  if (method==2){
    num_param=K*ncol(X) + #mu
      K*0.5*ncol(X)*(ncol(X)+1) + # sigma
      K-1  #tau
    BIC=-2*logllk.curr+n*log(num_param)
  }
  if (method==10 | method==11){
    num_param=K*ncol(X) + #mu
      K*0.5*ncol(X)*(ncol(X)+1) + # sigma
      K-1 + #tau
      1 #rho
    BIC=-2*logllk.curr+num_param*log(n)
  }
  
  return (list(sigma=sigma,mu=mu,logllk=logllk.curr,tau=tau,T=T,BIC=BIC))
}

EM.glasso.2 <- function(X, K, J, method, rho = NULL, sigma.true = NULL, penalize.diagonal = FALSE, epsilon.EM = 1e-8, max.iter.EM = 1000) {
  if (J > K) {
    stop("J must be less than or equal to K")
  }
  
  if (method == 1) {
    if (!is.null(sigma.true) && dim(sigma.true)[3] != J) {
      return(NA)
    }
  }
  
  p <- ncol(X)
  n <- nrow(X)
  
  # Initial values
  tau <- rep(1 / K, K)
  
  # Initial random assignment for mu
  sample.start <- sample(1:K, n, replace = TRUE)
  mu <- t(sapply(1:K, function(k) colMeans(X[sample.start == k, , drop = FALSE])))
  
  # Initialize sigma with J unique covariance matrices
  sigma <- array(0, dim = c(p, p, J))
  # Assign each of the K clusters to one of the J sigma groups
  sigma_assign <- sample(1:J, K, replace = TRUE)
  
  for (j in 1:J) {
    clusters_in_group <- which(sigma_assign == j)
    if (length(clusters_in_group) == 0) {
      # If no clusters assigned to this sigma group, initialize with identity
      sigma[, , j] <- diag(p)
    } else {
      # Initialize sigma[j] as the average covariance of assigned clusters
      sigma_group <- X[sample.start %in% clusters_in_group, , drop = FALSE]
      if (nrow(sigma_group) > 1) {
        sigma[, , j] <- cov(sigma_group)
      } else {
        sigma[, , j] <- diag(p)
      }
    }
  }
  
  # EM Algorithm
  logllk.pre <- -Inf
  count <- 0
  
  repeat {
    count <- count + 1
    
    # E-step: Calculate responsibilities
    f <- matrix(0, n, K)
    for (k in 1:K) {
      sigma_k <- sigma[, , sigma_assign[k]]
      f[, k] <- dmvnorm(X, mean = mu[k, ], sigma = sigma_k)
    }
    
    f_weighted <- sweep(f, 2, tau, FUN = '*')
    T_resp <- sweep(f_weighted, 1, rowSums(f_weighted), FUN = "/")
    
    # M-step: Update parameters
    tau <- colMeans(T_resp)
    
    for (k in 1:K) {
      mu[k, ] <- colSums(T_resp[, k] * X) / sum(T_resp[, k])
    }
    
    if (method == 1) sigma.true.temp <- sigma.true
    
    for (j in 1:J) {
      assigned_clusters <- which(sigma_assign == j)
      if (length(assigned_clusters) == 0) {
        next
      }
      # Aggregate responsibilities for sigma group j
      weighted_sum <- matrix(0, p, p)
      total_resp <- 0
      for (k in assigned_clusters) {
        X_centered <- sweep(X, 2, mu[k, ])
        weighted_sum <- weighted_sum + t(X_centered) %*% (X_centered * T_resp[, k])
        total_resp <- total_resp + sum(T_resp[, k])
      }
      sigma_temp <- weighted_sum / total_resp
      
      if (method == 2) {
        sigma[, , j] <- sigma_temp
      } else if (method == 10) {
        sigma[, , j] <- glasso(sigma_temp, rho = rho.method10, penalize.diagonal = penalize.diagonal)$w
      } else if (method == 11) {
        sigma[, , j] <- glasso(sigma_temp, rho = rho.method11, penalize.diagonal = penalize.diagonal)$w
      } else if (method == 1) {
        distances <- apply(sigma.true.temp, 3, function(sigma_true_j) norm(sigma_temp - sigma_true_j, type = "F"))
        min_idx <- which.min(distances)
        sigma[, , j] <- sigma.true.temp[, , min_idx]
        if (J > 1) {
          sigma.true.temp <- sigma.true.temp[, , -min_idx, drop = FALSE]
        }
      }
    }
    
    # Log-likelihood calculation
    logllk.curr <- sum(log(rowSums(sweep(f, 2, tau, FUN = '*'))))
    
    print(paste('Method: ', method, ' Iteration: ', count, sep = ''))
    print(paste('Current Q: ', round(logllk.curr, 2), sep = ''))
    
    # Convergence check
    if (abs(logllk.curr - logllk.pre) < epsilon.EM || count > max.iter.EM) {
      if (count > max.iter.EM) warning("Maximum iteration reached")
      break
    }
    
    # Prepare for next iteration
    logllk.pre <- logllk.curr
  }
  
  if (method == 2) {
    num_param <- K * ncol(X) + # mu
      J * 0.5 * ncol(X) * (ncol(X) + 1) + # sigma
      K - 1  # tau
    BIC <- -2 * logllk.curr + n * log(num_param)
  }
  if (method == 10 | method == 11) {
    num_param <- K * ncol(X) + # mu
      J * 0.5 * ncol(X) * (ncol(X) + 1) + # sigma
      K - 1 + # tau
      1 # rho
    BIC <- -2 * logllk.curr + num_param * log(n)
  }
  
  return(list(sigma = sigma, mu = mu, logllk = logllk.curr, tau = tau, T = T_resp, BIC = BIC))
}


hierarchical_gmm <- function(X, K, L, max_iter = 100, tol = 1e-6, rho , penalize.diagonal = FALSE) {
  n <- nrow(X)
  d <- ncol(X)
  
  # Initialize parameters
  pi_k <- rep(1/K, K) # Mixing coefficients for top-level clusters
  pi_kl <- matrix(1/L, nrow = K, ncol = L) # Mixing coefficients for sub-level clusters
  
  # Randomly divide X's rows into K groups and calculate covariance as initial value for Sigma_k
  groups <- split(as.data.frame(X), sample(rep(1:K, length.out = n)))
  Sigma_k <- array(0, dim = c(K, d, d))
  for (k in 1:K) {
    Sigma_k[k, , ] <- cov(as.matrix(groups[[k]]))
  }
  
  # Further divide each group into L subgroups and use the mean as the initial value of mu_kl
  mu_kl <- array(0, dim = c(K, L, d))
  for (k in 1:K) {
    subgroups <- split(as.data.frame(groups[[k]]), rep(1:L, length.out = nrow(groups[[k]])))
    for (l in 1:L) {
      mu_kl[k, l, ] <- colMeans(as.matrix(subgroups[[l]]))
    }
  }
  
  log_likelihood <- NULL
  
  for (iteration in 1:max_iter) {
    # E-Step: Calculate responsibilities
    gamma <- array(0, dim = c(n, K, L))
    
    for (k in 1:K) {
      for (l in 1:L) {
        gamma[, k, l] <- pi_k[k] * pi_kl[k, l] * dmvnorm(X, mean = mu_kl[k, l, ], sigma = Sigma_k[k, , ])
      }
    }
    
    # Normalize to get posterior probabilities for each data point
    gamma_sum <- rowSums(gamma, dims = 1)
    gamma <- gamma / gamma_sum
    
    # M-Step: Update parameters
    # Update top-level mixing coefficients pi_k
    N_k <- apply(gamma, 2, sum)
    pi_k <- N_k / n
    
    # Update sub-level mixing coefficients pi_kl and means mu_kl
    N_kl <- apply(gamma, c(2, 3), sum)
    pi_kl <- N_kl / N_k
    
    for (k in 1:K) {
      for (l in 1:L) {
        mu_kl[k, l, ] <- colSums(gamma[, k, l] * X) / N_kl[k, l]
      }
    }
    
    # Update covariance matrices Sigma_k for each top-level cluster
    for (k in 1:K) {
      Sigma_k[k, , ] <- matrix(0, nrow = d, ncol = d)
      for (l in 1:L) {
        diff <- sweep(X, 2, mu_kl[k, l, ])
        Sigma_k[k, , ] <- Sigma_k[k, , ] + t(diff) %*% (gamma[, k, l] * diff)
      }
      Sigma_k[k, , ] <- Sigma_k[k, , ] / N_k[k]
      Sigma_k[k, , ] <- glasso(Sigma_k[k, , ],rho=rho,penalize.diagonal=penalize.diagonal)$w
      
    }
    
    # Compute log-likelihood for convergence
    # Compute log-likelihood for convergence
    new_log_likelihood <- 0
    for (i in 1:n) {
      likelihood_i <- 0
      for (k in 1:K) {
        for (l in 1:L) {
          likelihood_i <- likelihood_i + pi_k[k] * pi_kl[k, l] * 
            dmvnorm(X[i, ], mean = mu_kl[k, l, ], sigma = Sigma_k[k, , ])
        }
      }
      new_log_likelihood <- new_log_likelihood + log(likelihood_i)
    }
    
    # new_log_likelihood <- sum(log(rowSums(apply(gamma, c(2, 3), function(g, k, l) {
    #   pi_k[k] * pi_kl[k, l] * dmvnorm(X, mean = mu_kl[k, l, ], sigma = Sigma_k[k, , ])
    # }, simplify = TRUE))))
    
    # Check convergence
    if (!is.null(log_likelihood) && abs(new_log_likelihood - log_likelihood) < tol) {
      cat("Converged at iteration", iteration, "\n")
      break
    }
    log_likelihood <- new_log_likelihood
    cat("Iteration", iteration, "Log Likelihood:", log_likelihood, "\n")
  }
  
  # Return parameters and responsibilities
  list(pi_k = pi_k, pi_kl = pi_kl, mu_kl = mu_kl, Sigma_k = Sigma_k, gamma = gamma)
}

# Function to convert precision matrix to neighboring matrix
precision_to_neighboring <- function(precision_matrix, cutoff) {
  # Check if the input is a matrix
  if (!is.matrix(precision_matrix)) {
    stop("Input must be a matrix.")
  }
  
  # Ensure the matrix is square
  if (nrow(precision_matrix) != ncol(precision_matrix)) {
    stop("Precision matrix must be square.")
  }
  
  # Apply the cutoff
  neighboring_matrix <- ifelse(abs(precision_matrix) >= cutoff, 1, 0)
  
  # Ensure diagonal entries are 0 (no self-loops)
  diag(neighboring_matrix) <- 0
  
  # Return the neighboring matrix
  return(neighboring_matrix)
}


calc_sensitivity_specificity <- function(true_matrix, estimated_matrix) {
  # Ensure the matrices are of the same size
  if (!all(dim(true_matrix) == dim(estimated_matrix))) {
    stop("True and estimated matrices must have the same dimensions.")
  }
  
  # Flatten the matrices for element-wise comparison
  true_vec <- as.vector(true_matrix[upper.tri(true_matrix)])
  estimated_vec <- as.vector(estimated_matrix[upper.tri(estimated_matrix)])
  
  # Calculate TP, FP, TN, FN
  TP <- sum(true_vec == 1 & estimated_vec == 1)
  FP <- sum(true_vec == 0 & estimated_vec == 1)
  TN <- sum(true_vec == 0 & estimated_vec == 0)
  FN <- sum(true_vec == 1 & estimated_vec == 0)
  
  # Calculate sensitivity and specificity
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  
  # Return the results
  return(list(
    sensitivity = sensitivity,
    specificity = specificity,
    TP = TP,
    FP = FP,
    TN = TN,
    FN = FN
  ))
}

#this enables each K, the number of sub clsuters may be different, i.e. L_k instead of L
hierarchical_gmm2 <- function(X, K, L_k, max_iter = 100, tol = 1e-6, rho, penalize.diagonal = FALSE) {
  n <- nrow(X)
  d <- ncol(X)
  
  # Initialize parameters
  pi_k <- rep(1/K, K)  # Mixing coefficients for top-level clusters
  pi_kl <- vector("list", K)  # Mixing coefficients for sub-level clusters
  mu_kl <- vector("list", K)  # Means for sub-level clusters
  Sigma_k <- array(0, dim = c(K, d, d))  # Covariance matrices for top-level clusters
  
  # Randomly divide X's rows into K groups and calculate covariance as initial value for Sigma_k
  groups <- split(as.data.frame(X), sample(rep(1:K, length.out = n)))
  for (k in 1:K) {
    Sigma_k[k, , ] <- cov(as.matrix(groups[[k]]))
    
    # Initialize sub-cluster mixing coefficients and means
    pi_kl[[k]] <- rep(1 / L_k[k], L_k[k])  # Uniform initialization
    subgroups <- split(as.data.frame(groups[[k]]), rep(1:L_k[k], length.out = nrow(groups[[k]])))
    mu_kl[[k]] <- matrix(0, nrow = L_k[k], ncol = d)
    for (l in 1:L_k[k]) {
      mu_kl[[k]][l, ] <- colMeans(as.matrix(subgroups[[l]]))
    }
  }
  
  log_likelihood <- NULL
  
  for (iteration in 1:max_iter) {
    # E-Step: Calculate responsibilities
    gamma <- vector("list", K)
    for (k in 1:K) {
      gamma[[k]] <- matrix(0, nrow = n, ncol = L_k[k])
      for (l in 1:L_k[k]) {
        gamma[[k]][, l] <- pi_k[k] * pi_kl[[k]][l] * 
          dmvnorm(X, mean = mu_kl[[k]][l, ], sigma = Sigma_k[k, , ])
      }
    }
    
    # Normalize to get posterior probabilities for each data point
    gamma_sum <- rowSums(do.call(cbind, lapply(gamma, rowSums)))
    for (k in 1:K) {
      gamma[[k]] <- gamma[[k]] / gamma_sum
    }
    
    # M-Step: Update parameters
    # Update top-level mixing coefficients pi_k
    N_k <- sapply(gamma, function(g) sum(rowSums(g)))
    pi_k <- N_k / n
    
    # Update sub-level mixing coefficients pi_kl and means mu_kl
    for (k in 1:K) {
      N_kl <- colSums(gamma[[k]])
      pi_kl[[k]] <- N_kl / N_k[k]
      for (l in 1:L_k[k]) {
        mu_kl[[k]][l, ] <- colSums(gamma[[k]][, l] * X) / N_kl[l]
      }
    }
    
    # Update covariance matrices Sigma_k for each top-level cluster
    for (k in 1:K) {
      Sigma_k[k, , ] <- matrix(0, nrow = d, ncol = d)
      for (l in 1:L_k[k]) {
        diff <- sweep(X, 2, mu_kl[[k]][l, ])
        Sigma_k[k, , ] <- Sigma_k[k, , ] + t(diff) %*% (gamma[[k]][, l] * diff)
      }
      Sigma_k[k, , ] <- Sigma_k[k, , ] / N_k[k]
      Sigma_k[k, , ] <- glasso(Sigma_k[k, , ], rho = rho, penalize.diagonal = penalize.diagonal)$w
    }
    
    # Compute log-likelihood for convergence
    new_log_likelihood <- 0
    for (i in 1:n) {
      likelihood_i <- 0
      for (k in 1:K) {
        for (l in 1:L_k[k]) {
          likelihood_i <- likelihood_i + pi_k[k] * pi_kl[[k]][l] * 
            dmvnorm(X[i, ], mean = mu_kl[[k]][l, ], sigma = Sigma_k[k, , ])
        }
      }
      new_log_likelihood <- new_log_likelihood + log(likelihood_i)
    }
    
    # Check convergence
    if (!is.null(log_likelihood) && abs(new_log_likelihood - log_likelihood) < tol) {
      cat("Converged at iteration", iteration, "\n")
      break
    }
    log_likelihood <- new_log_likelihood
    cat("Iteration", iteration, "Log Likelihood:", log_likelihood, "\n")
  }
  
  # Return parameters and responsibilities
  list(pi_k = pi_k, pi_kl = pi_kl, mu_kl = mu_kl, Sigma_k = Sigma_k, gamma = gamma)
}

#based on gmm2, can initialize based on the prior information of assignments
hierarchical_gmm3 <- function(X, K, max_iter = 100, tol = 1e-6, rho, 
                              penalize.diagonal = FALSE, 
                              subcluster_prior = NULL, 
                              layer_assignments = NULL) {
  n <- nrow(X)
  d <- ncol(X)
  

  
  if (is.null(layer_assignments)) {
    stop("layer_assignments is required for initializing the first layer clusters.")
  }
  

  # Determine which subclusters belong to each top-level cluster
  subcluster_indices <- split(1:length(layer_assignments), layer_assignments)
  L_k <- sapply(subcluster_indices, length)
  K <- length(subcluster_indices)
  
  # Initialize parameters
  pi_k <- rep(1 / K, K)  # Mixing coefficients for top-level clusters
  pi_kl <- vector("list", K)  # Mixing coefficients for sub-level clusters
  mu_kl <- vector("list", K)  # Means for sub-level clusters
  Sigma_k <- array(0, dim = c(K, d, d))  # Covariance matrices for top-level clusters
  
  if (!is.null(subcluster_prior)) {
    # Initialize using prior probabilities
    for (k in 1:K) {
      subclusters_k <- subcluster_indices[[k]]
      pi_kl[[k]] <- if (is.matrix(subcluster_prior[, subclusters_k])) {
        colSums(subcluster_prior[, subclusters_k]) / n
      } else {
        subcluster_prior[, subclusters_k] / n
      }
      mu_kl[[k]] <- t(sapply(subclusters_k, function(l) {
        weight <- if (is.vector(subcluster_prior[, l])) {
          subcluster_prior[, l]
        } else {
          subcluster_prior[, l, drop = FALSE]
        }
        colSums(weight * X) / sum(weight)
      }))
      
      # Initialize Sigma_k using weighted covariance
      for (j in seq_along(subclusters_k)) {
        l <- subclusters_k[j]
        weight <- if (is.vector(subcluster_prior[, l])) {
          subcluster_prior[, l]
        } else {
          subcluster_prior[, l, drop = FALSE]
        }
        diff <- sweep(X, 2, mu_kl[[k]][j, ])
        Sigma_k[k, , ] <- Sigma_k[k, , ] + t(diff) %*% (weight * diff)
      }
      Sigma_k[k, , ] <- Sigma_k[k, , ] / sum(subcluster_prior[, subclusters_k])
    }
  } else {
    # Uniform initialization if no prior provided
    groups <- split(as.data.frame(X), sample(rep(1:K, length.out = n)))
    for (k in 1:K) {
      Sigma_k[k, , ] <- cov(as.matrix(groups[[k]]))
      pi_kl[[k]] <- rep(1 / L_k[k], L_k[k])
      subgroups <- split(as.data.frame(groups[[k]]), rep(1:L_k[k], length.out = nrow(groups[[k]])))
      mu_kl[[k]] <- matrix(0, nrow = L_k[k], ncol = d)
      for (l in 1:L_k[k]) {
        mu_kl[[k]][l, ] <- colMeans(as.matrix(subgroups[[l]]))
      }
    }
  }
  log_likelihood <- NULL
  
  for (iteration in 1:max_iter) {
    # E-Step: Calculate responsibilities
    gamma <- vector("list", K)
    for (k in 1:K) {
      gamma[[k]] <- matrix(0, nrow = n, ncol = L_k[k])
      for (l in 1:L_k[k]) {
        gamma[[k]][, l] <- pi_k[k] * pi_kl[[k]][l] * 
          dmvnorm(X, mean = mu_kl[[k]][l, ], sigma = Sigma_k[k, , ])
      }
    }
    
    # Normalize to get posterior probabilities for each data point
    gamma_sum <- rowSums(do.call(cbind, lapply(gamma, rowSums)))
    for (k in 1:K) {
      gamma[[k]] <- gamma[[k]] / gamma_sum
    }
    
    # M-Step: Update parameters
    # Update top-level mixing coefficients pi_k
    N_k <- sapply(gamma, function(g) sum(rowSums(g)))
    pi_k <- N_k / n
    
    # Update sub-level mixing coefficients pi_kl and means mu_kl
    for (k in 1:K) {
      N_kl <- colSums(gamma[[k]])
      pi_kl[[k]] <- N_kl / N_k[k]
      for (l in 1:L_k[k]) {
        mu_kl[[k]][l, ] <- colSums(gamma[[k]][, l] * X) / N_kl[l]
      }
    }
    
    # Update covariance matrices Sigma_k for each top-level cluster
    for (k in 1:K) {
      Sigma_k[k, , ] <- matrix(0, nrow = d, ncol = d)
      for (l in 1:L_k[k]) {
        diff <- sweep(X, 2, mu_kl[[k]][l, ])
        Sigma_k[k, , ] <- Sigma_k[k, , ] + t(diff) %*% (gamma[[k]][, l] * diff)
      }
      Sigma_k[k, , ] <- Sigma_k[k, , ] / N_k[k]
      Sigma_k[k, , ] <- glasso(Sigma_k[k, , ], rho = rho, penalize.diagonal = penalize.diagonal)$w
    }
    
    # Compute log-likelihood for convergence
    new_log_likelihood <- 0
    for (i in 1:n) {
      likelihood_i <- 0
      for (k in 1:K) {
        for (l in 1:L_k[k]) {
          likelihood_i <- likelihood_i + pi_k[k] * pi_kl[[k]][l] * 
            dmvnorm(X[i, ], mean = mu_kl[[k]][l, ], sigma = Sigma_k[k, , ])
        }
      }
      new_log_likelihood <- new_log_likelihood + log(likelihood_i)
    }
    
    # Check convergence
    if (!is.null(log_likelihood) && abs(new_log_likelihood - log_likelihood) < tol) {
      cat("Converged at iteration", iteration, "\n")
      break
    }
    log_likelihood <- new_log_likelihood
    # cat("Iteration", iteration, "Log Likelihood:", log_likelihood, "\n")
  }
  
  # Return parameters and responsibilities
  list(pi_k = pi_k, pi_kl = pi_kl, mu_kl = mu_kl, Sigma_k = Sigma_k, gamma = gamma, llk=new_log_likelihood) 
  
 }
