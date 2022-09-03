library(HyperG)
library(combinat)
library(parallel)
library(Rcsdp)

############################
# Model build
############################
inci_m <-
  function(x,v)
  {
    s <- length(x)
    if(missing(v))
      v <- unique(unlist(x))
    n <- length(v)
    M <- Matrix::Matrix(0,nrow=s,ncol=n,dimnames=list(NULL,as.character(v)))
    for(i in 1:s){
      M[i,match(x[[i]],v)] <- 1
    }
    return(M)
  }

to_adja = function(M) {
  A <- Matrix::t(M) %*% M
  diag(A) <- 0
  return(A)
}

# hypergraph setting
# numCore = detectCores() / 2
numCore=1
# n=20;m=3;p=0.45;q=0.05;alpha=0.9
Gen_hyperG = function(n, m, p, q, alpha, numCore) {
  n1 <- sample(n, n/2) # random select n/2 nodes as community 1 nodes
  nf <- sequence(c(0, n), by = 1) # generate a full sequenceset.seed(2022)
  n2 <- nf[!nf %in% n1] # get n2 as community 2
  comb1 <- combn(n1, m) # hyperedges formed by vertices from community 1
  comb2 <- combn(n2, m) # hyperedges formed by vertices from community 2
  
  comb3_base1 <- combn(n1, m-1) # sample 1, randomly select cross-cluster nodes to form hyperedges
  comb3_base2 <- combn(n2, 1)
  
  upper <- replicate(ncol(comb3_base2), comb3_base1) # combine two group of nodes toghther
  comb3_upper <- mclapply(1:ncol(comb3_base2), function(d) rbind(upper[,,d], comb3_base2[d]), mc.cores = numCore)
  
  comb3_base1s <- combn(n2, m-1) # sample 2, randomly select cross-cluster nodes to form hyperedges
  comb3_base2s <- combn(n1, 1)
  
  lower <- replicate(ncol(comb3_base2s), comb3_base1s) # combine two group of nodes toghther
  comb3_lower <- mclapply(1:ncol(comb3_base2s), function(d) rbind(lower[,,d], comb3_base2s[d]), mc.cores = numCore)
  
  comb3_uppers_binded <- do.call(cbind, mclapply(comb3_upper, as.matrix, mc.cores = numCore))
  comb3_lowers_bined <- do.call(cbind, mclapply(comb3_lower, as.matrix, mc.cores = numCore))
  
  # comb3_uppers_binded <- do.call(cbind, mclapply(comb3_upper, as.data.frame, mc.cores = numCore))
  # comb3_lowers_bined <- do.call(cbind, mclapply(comb3_lower, as.data.frame, mc.cores = numCore))
  comb3 <- cbind(comb3_uppers_binded, comb3_lowers_bined)
  
  # # cencored hyperedges from different comnunities
  # comb3_cencored <- comb3[,rbinom(ncol(comb3), 1, alpha) * sequence(c(0, ncol(comb3)), by = 1) == 0]
  # # cencored hyperedges from same comnunity
  # comb1_cencored <- comb1[,rbinom(ncol(comb1), 1, alpha) * sequence(c(0, ncol(comb1)), by = 1) == 0]
  # comb2_cencored <- comb2[,rbinom(ncol(comb2), 1, alpha) * sequence(c(0, ncol(comb2)), by = 1) == 0]
  
  # generate valid hyperedge index with probability q times censored rate
  validhyeq <- rbinom(ncol(comb3), 1, alpha*q) * sequence(c(0, ncol(comb3)), by = 1)
  # generate valid hyperedge index for in community nodes with probability p times censored rate
  validhye <- rbinom(ncol(comb1), 1, alpha*p) * sequence(c(0, ncol(comb1)), by = 1)  
  validhye2 <- rbinom(ncol(comb2), 1, alpha*p) * sequence(c(0, ncol(comb2)), by = 1)  
  
  comb1s <- comb1[,validhye] # screened hyperedges
  comb2s <- comb2[,validhye2] # screened hyperedges
  comb3s <- comb3[,validhyeq] # screened hyperedges
  
  comb_all <- cbind(comb1s, comb2s)
  comb_all <- cbind(comb_all, comb3s) # hyperedges from in-community and cross-community
  
  comb_list <- mclapply(1:ncol(comb_all), function(c) list(comb_all[,c]), mc.cores = numCore)
  comb_list = unlist(comb_list, recursive = FALSE)
  
  # censored_all <- cbind(comb1_cencored, comb2_cencored)
  # censored_all <- cbind(censored_all, comb3_cencored) # all censored hyperedges
  # censored_list <- mclapply(1:ncol(censored_all), function(c) list(censored_all[,c]), mc.cores = numCore)
  # censored_list = unlist(censored_list, recursive = FALSE)
  
  # h <- hypergraph_from_edgelist(comb_list)
  hm <- inci_m(comb_list)
  out = list()
  # out$h = h
  # out$censored_list = censored_list
  out$h = hm
  out$n1 = n1
  return(out)
}

############################
# Algorithm applying
############################
# semi-definite algorithm
sdp = function(n,m,s,alpha,lambda,numCore) {
  p = s+lambda; q=s-lambda
  out <- Gen_hyperG(n,m,p,q,alpha,numCore)
  A_tilde <- out$h
  adj_A_tilde <- to_adja(A_tilde)
  ns <- seq(from=1,to=n,by=1)
  len_A <- nrow(adj_A_tilde)
  if (len_A != n) { # compensate the missing rows and columns
    miss_set <- ns %!in% as.integer(colnames(adj_A_tilde)) * ns
    miss_n <-  miss_set[!miss_set %in% 0]
    for (i in 1:length(miss_n)) {
      if (i == 1) {
        ferry_ini <- cbind(adj_A_tilde, 0)
        ferry <- rbind(ferry_ini, 0)
      } else {
        ferry = cbind(ferry, 0)
        ferry = rbind(ferry, 0)
      }
    }
    col_names <- as.integer(colnames(adj_A_tilde))
    col_names = append(col_names, miss_n)
    colnames(ferry) <- col_names
    rownames(ferry) <- col_names
    G <- ferry[, order(as.integer(colnames(ferry)))]
    G <- G[order(as.integer(rownames(G))),]
  } else {
    G <- adj_A_tilde[, order(as.integer(colnames(adj_A_tilde)))]
    G <- G[order(as.integer(rownames(G))),]
  }
  
  # input for semi-definite algorithm
  J <- matrix(1,n,n)
  C <- list(G)
  S_list = list()
  for (i in 1:n) {
    S <- matrix(0,n,n)
    S[i,i] <- 1
    S_list = append(S_list, list(list(S)))
  }
  A <- append(list(list(J)),S_list)
  b <- c(0,rep(1,n))
  K <- list(type=c("s"),size=c(n))
  
  # checking result
  n1 <- out$n1
  `%!in%` <- Negate(`%in%`)
  ns <- seq(from=1,to=n,by=1)
  sigma <- ns %!in% sort(n1) * rep(1, n) #true labels
  new_sig <- sigma*2-1
  Y <- t(matrix(new_sig,1,n))%*%matrix(new_sig,1,n) #true label matrix
  
  # result from semi-definite
  Y_hat <- csdp(C,A,b,K)$X[[1]] # Y_hat
  Y_hat = apply(Y_hat, 2, round)
  res <- Y_hat == Y
  true_num <- sum(res, na.rm = TRUE)
  error_rate <- ((n*n) - true_num) / (n*n)
  return(error_rate)
}

# simulation
n_list <- as.numeric(rep(50,4))
alpha_list <- rep(c(0.013,0.012,0.011,0.01), 1)
lambda1_list <- rep(0.15, 4)
lambda2_list <- rep(0.1, 4)
s = rep(0.25, 4)
sim_list <- data.frame(n_list,s,alpha_list,lambda1_list,lambda2_list)
nax <- c("n", "s", "alpha", "l1", "l2")
colnames(sim_list) <- nax

err_list1 = as.numeric()
sd_list1 = as.numeric()
for (i in 1:nrow(sim_list)) {
  rt <- replicate(25,sdp(sim_list[i,1], 3, sim_list[i,2], sim_list[i,3], sim_list[i,4], numCore))
  err_list1 = append(err_list1, mean(rt))
  sd_list1 = append(sd_list1, sd(rt))
}

err_list2 = as.numeric()
sd_list2 = as.numeric()
for (i in 1:nrow(sim_list)) {
  rt <- replicate(25,sdp(sim_list[i,1], 3, sim_list[i,2], sim_list[i,3], sim_list[i,5], numCore))
  err_list2 = append(err_list2, mean(rt))
  sd_list2 = append(sd_list2, sd(rt))
}

final_df <- data.frame(err_list1, sd_list1, err_list2, sd_list2)
write.table(final_df, 'sdp_results_add_50n_p0.425.csv', col.names = c('err_pq1','sd_pq1', 'err_pq2','sd_pq2'), 
            row.names = sim_list$n, sep = ',')
