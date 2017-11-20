# compute  covariance or correlation matrix
cov.compute = function(X, covType = "pearson"){
  
  # number of samples
  n = dim(X)[1]
  
  # compute covariance matrix if covType=="pearson"
  if(covType=="pearson"){
    S = (n-1)*cov(X, method = "pearson")/n
  }
  
  # compute correlation matrix if covType=="kendall"
  if (covType == "kendall"){
    S = cor(X, method = "kendall")
    # connect Kendall's correlation to Pearson's correlation using bridge function sin(S*pi/2)
    S = sin(S*pi/2)
    S[is.na(S)] = 0
    diag(S) = 1
    # project S into the cone of positive semidefinite matrices using Matrix::nearPD()
    S = as.matrix(nearPD(S, corr = TRUE)$mat)
  }
  
  # compute correlation matrix if covType=="spearman"
  if (covType == "spearman"){
    S = cor(X, method = "spearman")
    # connect Spearman's correlation to Pearson's correlation using bridge function 2*sin(S*pi/6)
    S = 2*sin(S*pi/6)
    S[is.na(S)] = 0
    diag(S) = 1
    # project S into the cone of positive semidefinite matrices using Matrix::nearPD()
    S = as.matrix(nearPD(S, corr = TRUE)$mat)
  }
  
  S
}





Dtrace = function(X, lambda, covType = "pearson", tol = 1e-5, maxiter = 500, rho = 0.1, rho.incr = 1.05, rho.max = 1e5){
  
  # number of states
  G = length(X)
  # number of genes
  p = dim(X[[1]])[2]
  
  if (is.data.frame(X[[1]])){
    for (g in seq(G)){
      X[[g]] = as.matrix(X[[g]])
    }
  }
  
  # assign gene names if none exist
  if(length(dimnames(X[[1]])[[2]])==0){
    for(g in seq(G)){
      dimnames(X[[g]])[[2]]=paste("V",1:p,sep="")
    }
  }
  
  # obtain the cov/cor matrices
  S = list()
  if (isSymmetric(X[[1]])){
    S = X
  }else{
    try(if (covType %in% c("pearson","kendall", "spearman") == FALSE) stop("The cov/cor type you provide is not include in this package. Please use your own function to obtain the list of cov/cor and use them as the input of FGL()"))
    for (g in seq(G)){
      S[[g]] = cov.compute(X[[g]], covType)
    }  
  }
  
  # solve optimization model using Dtrace.solve()
  result = Dtrace.solve(S, lambda, tol, maxiter, rho, rho.incr, rho.max)
  
  # create differential network from the estimated difference between two precision matrices, Delta = Theta_2^{-1} - Theta_1^{-1}
  Delta.graph = (abs(result$Delta)>1e-5)*1
  diag(Delta.graph) = 0
  Degree = apply(Delta.graph, 1, sum)
  Delta.graph.connected = Delta.graph[Degree>0, Degree>0]
  result$Delta.graph.full =  graph_from_adjacency_matrix(Delta.graph, mode = "undirected", weighted = TRUE, diag = FALSE)
  result$Delta.graph.connected =  graph_from_adjacency_matrix(Delta.graph.connected, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  result   
}





### ADMM for Dtrace
Dtrace.solve = function(S, lambda, tol = 1e-5, maxiter = 500, rho = 0.1, rho.incr = 1.05, rho.max = 1e10){
  
  # number of genes
  p = dim(S[[1]])[1]
  G = length(S)
  # assign gene names if none exist
  if(length(dimnames(S[[1]])[[2]])==0){
    for(g in seq(G)){
      dimnames(S[[g]])[[1]]=paste("V",1:p,sep="")
      dimnames(S[[g]])[[2]]=paste("V",1:p,sep="")
    }
  }
  
  
  # initialize 
  Delta = diag(p)
  Delta_p = diag(p)
  A = matrix(0, p, p)
  
  # compute eigenvalue decomposition of 0.5*(S1'*S2 + S2'*S1)
  Q = 0.5*(t(S[[1]])%*%S[[2]] + t(S[[2]])%*%S[[1]])
  W = S[[1]] - S[[2]]
  edecomp = svd(Q)
  SigmaQ = edecomp$d
  UQ = edecomp$u
  
  # ADMM iterations
  for(i in seq(maxiter)){
    
    # update Delta
    Delta.prev = Delta
    Delta = solve_G(UQ, SigmaQ, W + rho*Delta_p - A, rho)
    
    # update Delta_p
    Delta_p = soft(A/rho + Delta,lambda/rho, penalize.diagonal=TRUE)
    
    # updata A
    A = A + rho*(Delta - Delta_p)
    
    # check convergence condition
    diff1 = sum(abs(Delta - Delta.prev))
    diff2 = sum(abs(Delta - Delta_p))
    norm_value = sum(abs(Delta))
    if (max(diff1, diff2) < tol*max(norm_value, 1)){
      break
    }
    #update rho
    rho = min(rho*rho.incr,rho.max)
    
  }
  
  #assign gene names
  Delta = Delta_p
  diag(Delta) = 0
  row.names(Delta) = row.names(S[[1]])
  colnames(Delta) = colnames(S[[1]])
  
  result = list(Delta = Delta)
  result
}






# -------------------------- The Expand Operator ----------------------------- #
# --- Expand_n(A, rho, n) = argmin_{Theta} -n\log\det(Theta) + (rho/2)*||Theta - A||_F^2    						       #	
# --- Let A = U*D*U' 		       #	
# Set D2_{ii} = 0.5*(D_{ii} + sqrt(D_{ii}^2 + 4n/rho));            #  			 
# Expand_n(A, rho,n) = U*D2*U'					       #	
#									       

expand = function(A, rho, n){
  edecomp = eigen(A)
  D = edecomp$values
  U = edecomp$vectors
  D2 = 0.5*(D + sqrt(D^2 + 4*n/rho ))
  Theta = U %*% diag(D2) %*% t(U)
  Theta
}




### infer differential network using FGL
FGL = function(X, lambda1, lambda2, covType = "pearson", weights="equal",  penalize.diagonal = FALSE, tol = 1e-5, maxiter = 500, rho = 0.1, rho.incr = 1.05, rho.max = 1e5){
  
  # number of states
  G = length(X)
  # number of genes
  p = dim(X[[1]])[2]
  # sample sizes
  n = c()
  if (is.data.frame(X[[1]])){
    for (g in seq(G)){
      X[[g]] = as.matrix(X[[g]])
      n[g] = dim(X[[g]])[1]
    }
  }
  
  # assign gene names if none exist
  if(length(dimnames(X[[1]])[[2]])==0){
    for(g in seq(G)){
      dimnames(X[[g]])[[2]]=paste("V",1:p,sep="")
    }
  }
  
  
  # obtain the cov/cor matrices
  S = list()
  if (isSymmetric(X[[1]])){
    S = X
  }else{
    try(if (covType %in% c("pearson","kendall", "spearman") == FALSE) stop("The cov/cor type you provide is not include in this package. Please use your own function to obtain the list of cov/cor and use them as the input of FGL()"))
    for (g in seq(G)){
      S[[g]] = cov.compute(X[[g]], covType)
    }  
  }
  
  # set weights
  if(length(weights)==1){
    if(weights == "equal"){
      weights = rep(1,G)
    }
  }
  if(length(weights)==1){
    if(weights == "sample.size"){
      weights = n/sum(n)
    }
  }
  
  # solve optimization model using FGL.solve()
  result = FGL.solve(S, lambda1, lambda2, weights, penalize.diagonal, tol, maxiter, rho, rho.incr, rho.max)
  
  # create differential network from the estimated difference between two precision matrices, Delta = Theta_2^{-1} - Theta_1^{-1} 
  Delta.graph = (abs(result$Delta)>1e-5)*1
  diag(Delta.graph) = 0
  Degree = apply(Delta.graph, 1, sum)
  Delta.graph.connected = Delta.graph[Degree>0, Degree>0]
  result$Delta.graph.full =  graph_from_adjacency_matrix(Delta.graph, mode = "undirected", weighted = TRUE, diag = FALSE)
  result$Delta.graph.connected =  graph_from_adjacency_matrix(Delta.graph.connected, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  # create state-specific gene networks from the estimated  precision matrices, Theta
  result$Theta.graph.full = list()
  result$Theta.graph.connected = list()
  for(g in 1:G){
    Theta.graph = (abs(result$Theta[[g]])>1e-5)*1
    diag(Theta.graph) = 0
    Degree = apply(Theta.graph, 1, sum)
    Theta.graph.connected = Theta.graph[Degree>0, Degree>0]
    result$Theta.graph.full[[g]] =  graph_from_adjacency_matrix(Theta.graph, mode = "undirected", weighted = TRUE, diag = FALSE)
    result$Theta.graph.connected[[g]] =  graph_from_adjacency_matrix(Theta.graph.connected, mode = "undirected", weighted = TRUE, diag = FALSE)
    
  }
  
  result   
}



### ADMM for FGL
FGL.solve = function(S, lambda1, lambda2, n = c(1,1), penalize.diagonal = FALSE, tol = 1e-5, maxiter = 500, rho = 0.1, rho.incr = 1.05, rho.max = 1e10){
  
  # number of genes
  p = dim(S[[1]])[1]
  G = length(S)
  # assign gene names if none exist:
  if(length(dimnames(S[[1]])[[2]])==0){
    for(g in seq(G)){
      dimnames(S[[g]])[[1]]=paste("V",1:p,sep="")
      dimnames(S[[g]])[[2]]=paste("V",1:p,sep="")
    }
  }
  
  # initialize 
  Theta = list()
  Z = list()
  U = list()
  Theta[[1]] = diag(p)
  Theta[[2]] = diag(p)
  Z[[1]] = diag(p)
  Z[[2]] = diag(p)
  U[[1]] = matrix(0, p, p)
  U[[2]] = matrix(0, p, p)
  
  # ADMM iterations
  for (i in seq(maxiter)){
    
    # update theta
    Theta.prev = Theta
    Theta[[1]] = expand(Z[[1]] - (n[1]*S[[1]]+U[[1]])/rho, rho, n[1])
    Theta[[2]] = expand(Z[[2]] - (n[2]*S[[2]]+U[[2]])/rho, rho, n[2])
    
    # update Z
    A = list()
    A[[1]] = Theta[[1]] + U[[1]]/rho
    A[[2]] = Theta[[2]] + U[[2]]/rho
    Z = flsa2(A, rho, lambda1, lambda2, penalize.diagonal)
    
    # update U
    U[[1]] = U[[1]] + rho*(Theta[[1]] - Z[[1]])
    U[[2]] = U[[2]] + rho*(Theta[[2]] - Z[[2]])
    
    # check the convergence condition
    diff1 = sum(abs(Theta[[1]] - Theta.prev[[1]])) + sum(abs(Theta[[2]] - Theta.prev[[2]]))
    diff2 = sum(abs(Theta[[1]] - Z[[1]])) + sum(abs(Theta[[2]] - Z[[2]]))
    norm_value = sum(abs(Theta[[1]])) + sum(abs(Theta[[2]]))
    if (max(diff1, diff2) < tol*max(norm_value, 1)){
      break
    }
    
    # update rho
    rho = min(rho*rho.incr,rho.max)
    
  }
  
  # compute precision matrix difference
  Theta = Z
  Delta = Theta[[2]] - Theta[[1]]
  diag(Delta) = 0
  
  # assign gene names
  row.names(Theta[[1]]) = row.names(S[[1]])
  colnames(Theta[[1]]) = colnames(S[[1]])
  row.names(Theta[[2]]) = row.names(S[[1]])
  colnames(Theta[[2]]) = colnames(S[[1]])
  row.names(Delta) = row.names(S[[1]])
  colnames(Delta) = colnames(S[[1]])
  
  result = list(Delta = Delta, Theta = Theta)
  result
}




## -------------------- fused lasso signal approximator ------------------- ##
##
##---------------- minimize rho/2*(||Z_1 - A_1||_2^2 + ||Z_2 - A_2||_2^2)+ lambda1*(||Z_1||_1 + ||Z_2||_1) + lambda2*||Z_1 - Z_2||_1 -------- ##

flsa2 <-
  function(A, rho, lam1,lam2,penalize.diagonal)  #A is a list of 2 matrices from which we apply an L2 penalty to departures
  {
    
    S1 = abs(A[[1]]-A[[2]])<=2*lam2/rho
    X1 = (A[[1]]+A[[2]])/2
    Y1 = X1
    
    S2 = (A[[1]] > A[[2]]+2*lam2/rho)
    X2 = A[[1]] - lam2/rho
    Y2 = A[[2]] + lam2/rho
    
    S3 = (A[[2]] > A[[1]]+2*lam2/rho)
    X3 = A[[1]] + lam2/rho
    Y3 = A[[2]] - lam2/rho
    
    Z = list()
    Z[[1]] = soft(S1*X1 + S2*X2 + S3*X3, lam1/rho, penalize.diagonal)
    Z[[2]] = soft(S1*Y1 + S2*Y2 + S3*Y3, lam1/rho, penalize.diagonal)
    
    Z
  }





### ADMM for the subproblem of pDNA 
pDNA.admm.iters = function(S, lambda, Sigma.svd = NULL, Delta.init = NULL, V.init = NULL, tol = 1e-5, maxiter = 500, rho = 0.1, rho.incr = 1.05, rho.max = 1e6){
  
  # number of genes
  p = dim(S[[1]])[1]
  
  # compute eigenvalue decomposition of 0.5*(S_{k1}'*S_{k2} + S_{k2}'*S_{k1}) if not exist
  if (is.null(Sigma.svd)){
    Sigma.svd = list()
    A = 0.5*(t(S[[1]])%*%S[[2]] + t(S[[2]])%*%S[[1]])
    edecomp = svd(Q)
    Sigma.svd$D = edecomp$d
    Sigma.svd$U = edecomp$u
  }
  
  # initialize Delta
  if (is.null(Delta.init)){
    Delta = diag(p)
  }else{
    Delta = (Delta.init + t(Delta.init))/2
  }
  
  # initialize V
  if (is.null(V.init)){
    V = diag(p)
  }else{
    V = V.init
  }
  W = t(V)
  
  # initialize P and Q
  P = matrix(0, p, p)
  Q = matrix(0, p, p)
  
  
  # ADMM iterations
  for(i in seq(maxiter)){
    
    # update Delta
    Delta.prev = Delta;
    Delta = solve_G(Sigma.svd$U, Sigma.svd$D, S[[1]]-S[[2]]+rho*(V+W)-P, rho)
    # update H
    H = 0.5*(Delta - W + t(W)) + 0.5*(P - Q)/rho
    # update V
    V = soft(H, lambda/(2*rho), penalize.diagonal=TRUE)
    diag(V) = 0
    # update W
    W = 0.5*(Delta - V + t(V)) + 0.5*(P + t(Q))/rho
    
    # update P
    P = P + rho*(Delta - (V + W))
    # update Q
    Q = Q + rho*(V - t(W))
    
    # check convergence condition
    diff1 = sum(abs(Delta - Delta.prev))
    diff2 = sum(abs(Delta - (V + W)))
    diff3 = sum(abs(V - t(W)))
    norm_value = sum(abs(Delta))   
    if (max(c(diff1, diff2, diff2)) < tol*max(norm_value, 1)){
      break
    }
    
    # update rho
    rho = min(rho*rho.incr,rho.max)
  }
  
  Delta = V + t(V)
  result = list(Delta = Delta, V = V)
  result
}




### infer differential network using pDNA
pDNA = function(X, lambda, covType = "pearson", tol = 1e-5, maxiter = 500, rho = 0.1, rho.incr = 1.05, rho.max = 1e10){
  
  # number of data types
  K = dim(X)[1]
  # number of states
  G = dim(X)[2]
  # number of genes
  p = dim(X[[1,1]])[2]
  
  if (is.data.frame(X[[1,1]])){
    for (k in seq(G)){
      for (g in seq(G)){
        X[[k,g]] = as.matrix(X[[k,g]])
      }
    }
  }
  
  # assign gene names if none exist
  if (is.data.frame(X[[1,1]])){
    for (k in seq(G)){
      for (g in seq(G)){
        dimnames(X[[k,g]])[[2]]=paste("V",1:p,sep="")
        
      }
    }
  }
  
  # obtain the cov/cor matrices
  S = matrix(list(), K, G)
  if (isSymmetric(X[[1,1]])){
    S = X
  }else{
    try(if (covType %in% c("pearson","kendall", "spearman") == FALSE) stop("The cov/cor type you provide is not include in this package. Please use your own function to obtain the list of cov/cor and use them as the input of FGL()"))
    for (k in seq(K)){
      for (g in seq(G)){
        S[[k,g]] = cov.compute(X[[k,g]], covType)        
      }
    }
  }
  
  # solve optimization model using pDNA.solve()
  result = pDNA.solve(S, lambda, tol, maxiter, rho, rho.incr, rho.max)
  
  # create differential network from the estimated difference between two precision matrices, Delta_k = Theta_{k2}^{-1} - Theta_{k1}^{-1} 
  result$Delta.graph.full = list()
  result$Delta.graph.connected = list()
  Delta.weight = matrix(0,p,p)  
  for(k in seq(K)){
    Delta.graph = (abs(result$Delta[[k]])>1e-5)*1
    diag(Delta.graph) = 0
    rownames(Delta.graph) = colnames(X[[1,1]])
    colnames(Delta.graph) = colnames(X[[1,1]])
    Degree = apply(Delta.graph, 1, sum)
    Delta.graph.connected = Delta.graph[Degree>0, Degree>0]
    result$Delta.graph.full[[k]] =  graph_from_adjacency_matrix(Delta.graph, mode = "undirected", weighted = TRUE, diag = FALSE)
    result$Delta.graph.connected[[k]] =  graph_from_adjacency_matrix(Delta.graph.connected, mode = "undirected", weighted = TRUE, diag = FALSE)
    Delta.weight = Delta.weight + Delta.graph
  }  
  
  Delta.weight = Delta.weight/K
  result$Delta.weight = Delta.weight
  rownames(Delta.weight) = colnames(X[[1,1]])
  colnames(Delta.weight) = colnames(X[[1,1]])
  Degree = apply((Delta.weight!=0)*1, 1, sum)
  Delta.weight.connected = Delta.weight[Degree>0, Degree>0]
  result$Delta.graph.weight.full = graph_from_adjacency_matrix(Delta.weight, mode = "undirected", weighted = TRUE, diag = FALSE)
  result$Delta.graph.weight.connected = graph_from_adjacency_matrix(Delta.weight.connected, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  result   
}



### local linear approximation for pDNA
pDNA.solve = function(S, lambda, tol = 1e-5, maxiter = 500, rho = 0.1, rho.incr = 1.05, rho.max = 1e6){
  
  # number of data types
  K = dim(S)[1]
  G = dim(S)[2]
  
  # number of genes
  p = dim(S[[1,1]])[1]
 
  
  # assign gene names if none exist
  if(length(dimnames(S[[1]])[[2]])==0){
    for(k in seq(K)){
      for(g in seq(G)){
        dimnames(S[[k,g]])[[1]]=paste("V",1:p,sep="")
        dimnames(S[[k,g]])[[2]]=paste("V",1:p,sep="")
      }
    }
  }
  
  # compute eigenvalue decomposition of 0.5*(S_{k1}'*S_{k2} + S_{k2}'*S_{k1})
  Sigma.svd = list()
  for (k in seq(K)){
    Sigma.svd[[k]] = list()
    A = 0.5*(t(S[[k,1]])%*%S[[k,2]] + t(S[[k,2]])%*%S[[k,1]])
    edecomp = svd(A)
    Sigma.svd[[k]]$D = edecomp$d
    Sigma.svd[[k]]$U = edecomp$u
  }
  
  # initialize paramters using the solution to a lasso-type model
  Delta = list()
  V = list()
  for(k in seq(K)){
    temp.result = pDNA.admm.iters(S[k,], lambda, Sigma.svd = Sigma.svd[[k]], Delta.init = NULL, V.init = NULL, tol, maxiter, rho, rho.incr, rho.max)
    Delta[[k]] = temp.result$Delta
    V[[k]] = temp.result$V
  }
  
  
  # LLA iterations
  for(i in seq(maxiter)){
    
    # compute penalty weights, omega and psi 
    temp = matrix(0,p,p)
    for(k in seq(K)){
      temp =  temp +  abs(V[[k]])
    }
    temp.omega = sqrt(temp)
    temp.psi = sqrt(apply(temp.omega, 2, sum))
    omega = 1/(temp.omega+(1e-4))
    psi = 1/(temp.psi+(1e-4))
    
    # compute weighted penalties
    lambda.mat = lambda*omega%*%diag(psi)/4
    
    # update parameters using pDNA.admm.iters()
    Delta.prev = Delta  
    for(k in seq(K)){
      temp.result = pDNA.admm.iters(S[k,], lambda.mat, Sigma.svd = Sigma.svd[[k]], Delta.init = Delta[[k]], V.init = V[[k]], tol, maxiter, rho, rho.incr, rho.max)
      Delta[[k]] = temp.result$Delta
      V[[k]] = temp.result$V
    }
    
    # check convergence condition
    diff = 0
    norm_value = 0
    for (k in 1:K){
      diff = diff+sum(abs(Delta[[k]]-Delta.prev[[k]]))
      norm_value = norm_value + sum(abs(Delta[[k]]))
    }
    if(diff<tol*max(norm_value,1)){
      break
    }  
    
  }
  
  # compute precision matrix differences
  Delta.weight = matrix(0,p,p)
  for(k in 1:K){
    Delta[[k]] = V[[k]] + t(V[[k]])
    diag(Delta[[k]]) = 0
    row.names(V[[k]]) = row.names(S[[1]])
    colnames(V[[k]]) = colnames(S[[1]])
    row.names(Delta[[k]]) = row.names(S[[1]])
    colnames(Delta[[k]]) = colnames(S[[1]])
  }
  
  result = list(Delta = Delta)
  result
}






### infer differential network using PNJGL
PNJGL = function(X, lambda1, lambda2, covType = "pearson", weights="equal",  penalize.diagonal = FALSE, tol = 1e-5, maxiter = 500, rho = 0.1, rho.incr = 1.05, rho.max = 1e5){
  
  # number of states
  G = length(X)
  # number of genes
  p = dim(X[[1]])[2]
  # sample sizes
  n = c()
  if (is.data.frame(X[[1]])){
    for (g in seq(G)){
      X[[g]] = as.matrix(X[[g]])
      n[g] = dim(X[[g]])[1]
    }
  }
  
  # assign gene names if none exist
  if(length(dimnames(X[[1]])[[2]])==0){
    for(g in seq(G)){
      dimnames(X[[g]])[[2]]=paste("V",1:p,sep="")
    }
  }
  
  
  # obtain the cov/cor matrices
  S = list()
  if (isSymmetric(X[[1]])){
    S = X
  }else{
    try(if (covType %in% c("pearson","kendall", "spearman") == FALSE) stop("The cov/cor type you provide is not include in this package. Please use your own function to obtain the list of cov/cor and use them as the input of FGL()"))
    for (g in seq(G)){
      S[[g]] = cov.compute(X[[g]], covType)
    }  
  }
  
  # set weights
  if(length(weights)==1){
    if(weights == "equal"){
      weights = rep(1,G)
    }
  }
  if(length(weights)==1){
    if(weights == "sample.size"){
      weights = n/sum(n)
    }
  }
  
  # solve optimization model using PNJGL.solve()
  result = PNJGL.solve(S, lambda1, lambda2, weights, penalize.diagonal, tol, maxiter, rho, rho.incr, rho.max)
  
  # create differential network from the estimated difference between two precision matrices, Delta = Theta_2^{-1} - Theta_1^{-1}
  Delta.graph = (abs(result$Delta)>1e-5)*1
  diag(Delta.graph) = 0
  Degree = apply(Delta.graph, 1, sum)
  Delta.graph.connected = Delta.graph[Degree>0, Degree>0]
  result$Delta.graph.full =  graph_from_adjacency_matrix(Delta.graph, mode = "undirected", weighted = TRUE, diag = FALSE)
  result$Delta.graph.connected =  graph_from_adjacency_matrix(Delta.graph.connected, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  # create state-specific gene networks from the estimated  precision matrices, Theta
  result$Theta.graph.full = list()
  result$Theta.graph.connected = list()
  for(g in 1:G){
    Theta.graph = (abs(result$Theta[[g]])>1e-5)*1
    diag(Theta.graph) = 0
    Degree = apply(Theta.graph, 1, sum)
    Theta.graph.connected = Theta.graph[Degree>0, Degree>0]
    result$Theta.graph.full[[g]] =  graph_from_adjacency_matrix(Theta.graph, mode = "undirected", weighted = TRUE, diag = FALSE)
    result$Theta.graph.connected[[g]] =  graph_from_adjacency_matrix(Theta.graph.connected, mode = "undirected", weighted = TRUE, diag = FALSE)  
  }
  
  
  result   
}





### ADMM for PNJGL
PNJGL.solve = function(S, lambda1, lambda2, n = c(1,1), penalize.diagonal = FALSE, tol = 1e-5, maxiter = 500, rho = 0.1, rho.incr = 1.05, rho.max = 1e5){
  
  # number of genes
  p = dim(S[[1]])[1]
  
  #initialize primal variables
  Theta = list()
  Z = list()
  Theta[[1]] = diag(p)
  Theta[[2]] = diag(p)
  Z[[1]] = diag(p)
  Z[[2]] = diag(p)
  
  # initialize dual Variables
  V = matrix(0, p, p)
  W = matrix(0, p, p)
  F = matrix(0, p, p)
  H = matrix(0, p, p)
  Q = list()
  Q[[1]] = matrix(0, p, p)
  Q[[2]] = matrix(0, p, p)
  
  # ADMM iterations
  for (i in seq(maxiter)){
    
    # update theta
    Theta.prev = Theta
    Theta[[1]] = expand(1/(2*rho)*(rho*(Theta[[2]] + V + W+ Z[[1]]) - (Q[[1]] + (S[[1]]*n[1]) + F)), 2*rho, n[1])
    Theta[[2]] = expand(1/(2*rho)*(rho*(Theta[[1]] - (V + W) + Z[[2]]) - (Q[[2]] + (S[[2]]*n[2]) - F)), 2*rho, n[2])
    
    # update Z
    Z[[1]] = soft(Theta[[1]] + Q[[1]]/rho, lambda1/rho, penalize.diagonal)
    Z[[2]] = soft(Theta[[2]] + Q[[2]]/rho, lambda1/rho, penalize.diagonal)
    
    # update V
    V = soft_scal(0.5*(t(W) - W + Theta[[1]]- Theta[[2]]) + 1/(2*rho)*(F - H), lambda2/(2*rho))
    
    # update W
    W = 0.5*(Theta[[1]] - Theta[[2]] + t(V) - V) + 1/(2*rho)*(F + t(H))
    
    # update F
    F = F + rho*(Theta[[1]] - Theta[[2]] - (V + W))
    
    # update H
    H = H + rho*(V - t(W))
    
    # update Q
    Q[[1]] = Q[[1]] + rho*(Theta[[1]] - Z[[1]])
    Q[[2]] = Q[[2]] + rho*(Theta[[2]] - Z[[2]])
    
    # check the convergence condition
    diff1 = sum(abs(Theta[[1]] - Theta.prev[[1]])) + sum(abs(Theta[[2]] - Theta.prev[[2]]))
    diff2 = sum(abs(Theta[[1]] - Z[[1]])) + sum(abs(Theta[[2]] - Z[[2]]))
    norm_value = sum(abs(Theta[[1]])) + sum(abs(Theta[[2]]))
    if (max(diff1, diff2) < tol*max(norm_value, 1)){
      break
    }
    
    # update rho
    rho = min(rho*rho.incr,rho.max)
    
  }
  
  # compute precision matrix difference
  Theta = Z
  Delta = Theta[[2]] - Theta[[1]]
  diag(Delta) = 0
  
  # assign gene names
  row.names(Theta[[1]]) = row.names(S[[1]])
  colnames(Theta[[1]]) = colnames(S[[1]])
  row.names(Theta[[2]]) = row.names(S[[1]])
  colnames(Theta[[2]]) = colnames(S[[1]]) 
  row.names(Delta) = row.names(S[[1]])
  colnames(Delta) = colnames(S[[1]])
  
  result = list(Delta = Delta, Theta = Theta)
  result
  
}




## -------------------- SOFT-THRESHOLDING OPERATOR ------------------- ##
##
##---------------- minimize 1/2*||X - Y||_2^2 + lambda*||X||_1 -------- ##

soft = function(Y,lambda, penalize.diagonal = FALSE){ # if last argument is FALSE, soft-threshold Y matrix but don't penalize the diagonal
  X = sign(Y)*pmax(0, abs(Y)-lambda)
  if(!penalize.diagonal) diag(X) = diag(Y)
  X
}




## ------------------- SOFT-SCALING OPERATOR ----------------------- ##

# ---------------- minimize_X \frac{1}{2}||X - Y||_2^2 + \lambda*sum_i||X_i||_2

# ---------------- solution is  X_i = max(||Y_i||_2 - \lambda, 0)*Y_i/||Y_i||_2

soft_scal = function(Y, lambda){
  m = dim(Y)[1]
  n = dim(Y)[2]
  
  X = matrix(0, m, n)
  
  for(i in seq(n)){
    norm_Y = sqrt(sum((Y[,i])^2))
    if (norm_Y>lambda){
      X[,i] = (1 - lambda/norm_Y)*Y[,i]
    }  
  }
  X
}




## ------------------- solve_G OPERATOR ----------------------- ##

# ---------------- minimize_{X =X'} \frac{1}{2} trace(Delta'*(A+rho*I)*Delta) - trace(Delta*B)

# ---------------- solution is  UA*((UA'*B*UA)*C)*UA', where A = UA*diag(SigmaA)*UA'  and C_{ij} = 2/(SigmaA_ii+rho+SigmaA_jj+rho)


solve_G = function(UA, SigmaA, B, rho){ 
  D = length(SigmaA);
  Sigma = SigmaA + rho;
  Sigma.mat = matrix(rep(Sigma,D),D)
  C = 2/(Sigma.mat+t(Sigma.mat))
  X = UA%*%((t(UA)%*%B%*%UA)*C)%*%t(UA)          
}


