# library(HyperG)
library(combinat)
library(parallel)
library(magrittr)
library(dplyr)
library(caret)

setwd('E:/Dissertation/Simulation')
############################
# equal-sized k-means
############################
kMeanAdj = function(kdat, iter = t, k) {
  data1 <- kMeansInit(kdat, k)
  initc <- data1$centers
  for (i in 2:iter) {
    if (i == 2) {
      data2 <- kMeansIter(kdat, initc, k)
      update_c <- data2$centers
    } else {
      data2 <- kMeansIter(kdat, update_c, k)
      update_c <- data2$centers
    }
  }
  return(data2)
}

kMeansInit = function(kdat, k) {
  kdat = as.data.frame(kdat)
  kdat %>% kmeans(k) -> kclust
  kdist = function(x1, y1, x2, y2){
    sqrt((x1-x2)^2 + (y1-y2)^2)
  }
  centers = kclust$centers
  
  kdat %<>% 
    mutate(D1 = kdist(kdat[,1], kdat[,2], centers[1,1], centers[1,2]))
  kdat %<>% 
    mutate(D2 = kdist(kdat[,1], kdat[,2], centers[2,1], centers[2,2]))
  
  kdat$assigned = 0
  kdat$index = 1:nrow(kdat)
  working = kdat
  nassign = trunc(nrow(kdat)/k)
  
  for(i in 1:nassign){
    ind1 = working$index[which(working$D1 == min(working$D1))[1]]
    kdat$assigned[kdat$index == ind1] = 1
    working %<>% filter(!index == ind1)
    
    ind2 = working$index[which(working$D2 == min(working$D2))[1]]
    kdat$assigned[kdat$index == ind2] = 2
    working %<>% filter(!index == ind2)
  }
  
  NewCenters <- kdat %>% filter(assigned == 1) %>% 
    select(V1, V2) %>%
    kmeans(1) %$% centers
  
  NewCenters %<>% rbind(kdat %>% 
                          filter(assigned == 2) %>%
                          select(V1, V2) %>%
                          kmeans(1) %$% centers)
  
  NewCenters %<>% as.data.frame()
  
  out <- list()
  out$centers = NewCenters
  out$Data = kdat
  return(out)
}

kMeansIter = function(kdat, center, k) {
  kdat = as.data.frame(kdat)
  kdat %>% kmeans(center) -> kclust
  kdist = function(x1, y1, x2, y2){
    sqrt((x1-x2)^2 + (y1-y2)^2)
  }
  centers = kclust$centers
  
  kdat %<>% 
    mutate(D1 = kdist(kdat[,1], kdat[,2], centers[1,1], centers[1,2]))
  kdat %<>% 
    mutate(D2 = kdist(kdat[,1], kdat[,2], centers[2,1], centers[2,2]))
  
  kdat$assigned = 0
  kdat$index = 1:nrow(kdat)
  working = kdat
  nassign = trunc(nrow(kdat)/k)
  
  for(i in 1:nassign){
    ind1 = working$index[which(working$D1 == min(working$D1))[1]]
    kdat$assigned[kdat$index == ind1] = 1
    working %<>% filter(!index == ind1)
    
    ind2 = working$index[which(working$D2 == min(working$D2))[1]]
    kdat$assigned[kdat$index == ind2] = 2
    working %<>% filter(!index == ind2)
  }
  
  NewCenters <- kdat %>% filter(assigned == 1) %>% 
    select(V1, V2) %>%
    kmeans(1) %$% centers
  
  NewCenters %<>% rbind(kdat %>% 
                          filter(assigned == 2) %>%
                          select(V1, V2) %>%
                          kmeans(1) %$% centers)
  
  NewCenters %<>% as.data.frame()
  
  out <- list()
  out$centers = NewCenters
  out$Data = kdat
  return(out)
}

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

adja_m = function(x, v) {
  s <- length(x)
  if(missing(v))
    v <- unique(unlist(x))
  n <- length(v)
  M <- Matrix::Matrix(0,nrow=s,ncol=n,dimnames=list(NULL,as.character(v)))
  for(i in 1:s){
    M[i,match(x[[i]],v)] <- 1
  }
  A <- Matrix::t(M) %*% M
  diag(A) <- 0
  return(A)
}

############################
# Model build
############################
# hypergraph setting
numCore = detectCores() / 2
# numCore = 1
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
ES_algo = function(n, m, p, q, h, n1) {
  ns <- sequence(c(0, n), by = 1)
  M <- combn(ns, m)
  S1_index <- rbinom(ncol(M), 1, 0.2) # randomly form S1 with selecting hyperedges with probability loglogn/logn
  S1 <- M[, S1_index * sequence(c(0, ncol(M)), by = 1)]
  S <- rbinom(ncol(M), 1, 1)
  S2_index <- S - S1_index
  S2 <- M[, S2_index * sequence(c(0, ncol(M)), by = 1)]
  
  # subgraph A_tilde and A_bar
  valid_egde <- unlist(mclapply(1:nrow(h), function(row) list(sort(as.integer(colnames(h)[h[row,]==1]))),
                                mc.cores = numCore), recursive = FALSE) # extract formed hyperedges
  
  comb_S1 <- unlist(mclapply(1:ncol(S1), function(c) list(S1[,c]), mc.cores = numCore), recursive = FALSE)
  A_tilde_edge <- intersect(valid_egde, comb_S1) # hyperedges formed by vertices in set S1
  # A1 <- hypergraph_from_edgelist(A_tilde_edge) # subgraph A_tilde
  
  comb_S2 <- unlist(mclapply(1:ncol(S2), function(c) list(S2[,c]), mc.cores = numCore), recursive = FALSE)
  A_bar_edge <- intersect(valid_egde, comb_S2) # hyperedges formed by vertices in set S2
  
  # A2 <- hypergraph_from_edgelist(A_bar_edge) # subgraph A_bar
  # A_bar_complete <- append(A_bar_edge, censored_list) # add censored edges onto sub-hypergraph A_bar
  
  ############################
  # apply HSC on A_tilde
  ha_a1 <- adja_m(A_tilde_edge)
  # ha_a1 <- as.matrix(hypergraph_as_adjacency_matrix(A1)) # Adjacency matrix of A_tilde
  # ns <- seq(from=1,to=n,by=1)
  `%!in%` <- Negate(`%in%`)
  
  # obtain the similarity matrix, we use the following codes to compensate missing rows and cols
  len_A <- nrow(ha_a1)
  if (len_A != n) {
    miss_set <- ns %!in% as.integer(colnames(ha_a1)) * ns
    miss_n <-  miss_set[!miss_set %in% 0]
    for (i in 1:length(miss_n)) {
      if (i == 1) {
        ferry_ini <- cbind(ha_a1, 0)
        ferry <- rbind(ferry_ini, 0)
      } else {
        ferry = cbind(ferry, 0)
        ferry = rbind(ferry, 0)
      }
    }
    col_names <- as.integer(colnames(ha_a1))
    col_names = append(col_names, miss_n)
    colnames(ferry) <- col_names
    rownames(ferry) <- col_names
    A <- ferry[, order(as.integer(colnames(ferry)))]
    A <- A[order(as.integer(rownames(A))),]
  } else {
    A <- ha_a1[, order(as.integer(colnames(ha_a1)))]
    A <- A[order(as.integer(rownames(A))),]
  }
  
  # k-means clustering
  U_zero <- t(eigen(A)$vectors[1:2,])
  # U_zero <- t(eigen(A)$vectors[1:2,]) # k largest eigenvector
  data2 = kMeanAdj(U_zero, iter = 10, 2)
  cluster_res <- data2$Data$assigned # cluster result for vertices
  I_posi <- seq(from=1, to=n, by=1) * (cluster_res==1)
  I_posi = I_posi[I_posi!=0]
  I_nega <- seq(from=1, to=n, by=1) * (cluster_res==2)
  I_nega = I_nega[I_nega!=0]
  
  # refinement
  calc_e_posi = function(x) {
    nc = rbind(combn(setdiff(I_posi, x), 2), as.numeric(x))
    reconstruct <- unlist(lapply(1:ncol(nc), function(c) list(sort(nc[,c]))), 
                          recursive = FALSE)
    if (length(intersect(reconstruct, A_bar_edge)) > 0) {
      e_i_res = log(p/q) * length(intersect(reconstruct, A_bar_edge)) +
        log((1-p)/(1-q)) * (length(A_bar_edge) - length(intersect(reconstruct, A_bar_edge)))
    } else {
      e_i_res = log((1-p)/(1-q)) * length(A_bar_edge)
    }
    return(e_i_res)
  }
  
  calc_e_nega = function(x) {
    nc = rbind(combn(setdiff(I_nega, x), 2), as.numeric(x))
    reconstruct <- unlist(lapply(1:ncol(nc), function(c) list(sort(nc[,c]))), 
                          recursive = FALSE)
    if (length(intersect(reconstruct, A_bar_edge)) > 0) {
      e_i_res = log(p/q) * length(intersect(reconstruct, A_bar_edge)) +
        log((1-p)/(1-q)) * (length(A_bar_edge) - length(intersect(reconstruct, A_bar_edge)))
    } else {
      e_i_res = log((1-p)/(1-q)) * length(A_bar_edge)
    }
    return(e_i_res)
  }
  
  e_i_pos = unlist(mclapply(I_posi, calc_e_posi, mc.cores = numCore), recursive = FALSE)
  e_i_neg = unlist(mclapply(I_posi, calc_e_nega, mc.cores = numCore), recursive = FALSE)
  e_j_pos = unlist(mclapply(I_nega, calc_e_posi, mc.cores = numCore), recursive = FALSE)
  e_j_neg = unlist(mclapply(I_nega, calc_e_nega, mc.cores = numCore), recursive = FALSE)
  
  # e_i_pos <- as.numeric()
  # for (i in 1:length(I_posi)) {
  #   nc = rbind(combn(setdiff(I_posi, I_posi[i]), 2), I_posi[i]) # i and other two nodes form a hyperedge
  #   reconstruct <- unlist(mclapply(1:ncol(nc), function(c) list(sort(nc[,c])), mc.cores = numCore), 
  #                         recursive = FALSE) # to separate lists
  #   if (length(intersect(reconstruct, A_bar_edge)) > 0) {
  #     e_i_pos = append(e_i_pos, log(p/q) * length(intersect(reconstruct, A_bar_edge)) +
  #                        log((1-p)/(1-q)) * (length(A_bar_edge) - length(intersect(reconstruct, A_bar_edge))))
  #   } else {
  #     e_i_pos = append(e_i_pos, log((1-p)/(1-q)) * length(A_bar_edge))
  #   }
  # }
  
  # e_i_neg = as.numeric()
  # for (i in 1:length(I_posi)) {
  #   nc = rbind(combn(setdiff(I_nega, I_posi[i]), 2), I_posi[i])
  #   reconstruct <- unlist(mclapply(1:ncol(nc), function(c) list(sort(nc[,c])), mc.cores = numCore), 
  #                         recursive = FALSE)
  #   if (length(intersect(reconstruct, A_bar_edge)) > 0) {
  #     e_i_neg = append(e_i_neg, log(p/q) * length(intersect(reconstruct, A_bar_edge)) +
  #                        log((1-p)/(1-q)) * (length(A_bar_edge) - length(intersect(reconstruct, A_bar_edge))))
  #   } else {
  #     e_i_neg = append(e_i_neg, log((1-p)/(1-q)) * length(A_bar_edge))
  #   }
  # }
  
  # e_j_neg = as.numeric()
  # for (i in 1:length(I_nega)) {
  #   nc = rbind(combn(setdiff(I_nega, I_nega[i]), 2), I_nega[i])
  #   reconstruct <- unlist(mclapply(1:ncol(nc), function(c) list(sort(nc[,c])), mc.cores = numCore), 
  #                         recursive = FALSE)
  #   if (length(intersect(reconstruct, A_bar_edge)) > 0) {
  #     e_j_neg = append(e_j_neg, log(p/q) * length(intersect(reconstruct, A_bar_edge)) +
  #                        log((1-p)/(1-q)) * (length(A_bar_edge) - length(intersect(reconstruct, A_bar_edge))))
  #   } else {
  #     e_j_neg = append(e_j_neg, log((1-p)/(1-q)) * length(A_bar_edge))
  #   }
  # }
  
  # e_j_pos = as.numeric()
  # for (i in 1:length(I_nega)) {
  #   nc = rbind(combn(setdiff(I_posi, I_nega[i]), 2), I_nega[i])
  #   reconstruct <- unlist(mclapply(1:ncol(nc), function(c) list(sort(nc[,c])), mc.cores = numCore), 
  #                         recursive = FALSE)
  #   if (length(intersect(reconstruct, A_bar_edge)) > 0) {
  #     e_j_pos = append(e_j_pos, log(p/q) * length(intersect(reconstruct, A_bar_edge)) +
  #                        log((1-p)/(1-q)) * (length(A_bar_edge) - length(intersect(reconstruct, A_bar_edge))))
  #   } else {
  #     e_j_pos = append(e_j_pos, log((1-p)/(1-q)) * length(A_bar_edge))
  #   }
  # }
  
  i_m <- matrix(c(e_i_pos, e_i_neg), ncol = 2) # < then flip the label
  j_m <- matrix(c(e_j_neg, e_j_pos), ncol = 2) # < then flip the label
  
  # flip the labels
  I_posi_new <- unlist(mcmapply(function(l1, l2, vec, c) if (l1 < l2) {vec = vec[!vec %in% vec[c]]}, 
                                e_i_pos, e_i_neg, I_posi, 1:length(e_i_pos), mc.cores = numCore), recursive = FALSE)
  I_posi_new = I_posi[!I_posi %in% I_posi_new]
  
  I_nega_new <- unlist(mcmapply(function(l1, l2, vec, c) if (l1 < l2) {vec = vec[!vec %in% vec[c]]}, 
                                e_j_neg, e_j_pos, I_nega, 1:length(e_j_neg), mc.cores = numCore), recursive = FALSE)
  I_nega_new = I_nega[!I_nega %in% I_nega_new]
  
  I_nega_fliped <- sort(append(I_nega_new, setdiff(I_posi, I_posi_new)))
  I_posi_fliped <- sort(append(I_posi_new, setdiff(I_nega, I_nega_new)))
  
  # group the extra vertices to the other community
  len1 <- length(I_posi_fliped)
  len2 <- length(I_nega_fliped)
  
  if (len1 > len2) {
    lendiff <- (len1-len2)/2
    dif <- as.numeric(i_m[,1]) - as.numeric(i_m[,2])
    if (0 %in% dif) {
      # nl = as.numeric()
      # vt = as.numeric()
      # for (i in 1:length(dif)) {
      #   if (dif[i] == 0) {
      #     nl <- append(nl, as.numeric(i_m[i,1]))
      #     vt <- append(vt, as.numeric(I_posi[i]))
      #   }
      # }
      nl <- unlist(mcmapply(function(x, y) if (x==0) {as.numeric(y)}, dif, i_m[,1], mc.cores = numCore), recursive = FALSE)
      vt <- unlist(mcmapply(function(x, y) if (x==0) {as.numeric(y)}, dif, I_posi, mc.cores = numCore), recursive = FALSE)
      
      dt <- data.frame(nl,vt)
      if (length(vt) < lendiff) {
        res_node <- dt$vt
        I_nega_fliped = sort(append(I_nega_fliped, res_node))
        I_posi_fliped = I_posi_fliped[! I_posi_fliped %in% res_node]
        needmore <- lendiff - length(vt)
        # ds = as.numeric()
        # nd = as.numeric()
        # for (i in 1:length(dif)) {
        #   if (dif[i] > 0) {
        #     ds <- append(ds, as.numeric(i_m[i,1]))
        #     nd <- append(nd, as.numeric(I_posi[i]))
        #   }
        # }
        ds <- unlist(mcmapply(function(x, y) if (x>0) {as.numeric(y)}, dif, i_m[,1], mc.cores = numCore), recursive = FALSE)
        nd <- unlist(mcmapply(function(x, y) if (x>0) {as.numeric(y)}, dif, I_posi, mc.cores = numCore), recursive = FALSE)
        
        dt2 <- data.frame(ds,nd)
        more_node <- dt2$nd[1:needmore] # if multiple found, return the first several ones
        I_nega_fliped = sort(append(I_nega_fliped, more_node))
        I_posi_fliped = I_posi_fliped[! I_posi_fliped %in% more_node]
        
      } else if (length(vt) > lendiff) {
        dt_sort = dt[order(dt$nl, decreasing = TRUE),]
        res_node <- dt_sort$vt[1:lendiff] # if multiple found, return the first several one
        I_nega_fliped = sort(append(I_nega_fliped, res_node))
        I_posi_fliped = I_posi_fliped[! I_posi_fliped %in% res_node]
      } else {
        res_node <- dt$vt
        I_nega_fliped = sort(append(I_nega_fliped, res_node))
        I_posi_fliped = I_posi_fliped[! I_posi_fliped %in% res_node]
      }
    } else {
      # nl = as.numeric()
      # vt = as.numeric()
      # for (i in 1:length(dif)) {
      #   if (dif[i] > 0) {
      #     nl <- append(nl, as.numeric(i_m[i,1]))
      #     vt <- append(vt, as.numeric(I_posi[i]))
      #   }
      # }
      nl <- unlist(mcmapply(function(x, y) if (x>0) {as.numeric(y)}, dif, i_m[,1], mc.cores = numCore), recursive = FALSE)
      vt <- unlist(mcmapply(function(x, y) if (x>0) {as.numeric(y)}, dif, I_posi, mc.cores = numCore), recursive = FALSE)
      
      dt <- data.frame(nl,vt)
      res_node <- dt$vt[1:lendiff]
      I_nega_fliped = sort(append(I_nega_fliped, res_node))
      I_posi_fliped = I_posi_fliped[! I_posi_fliped %in% res_node]
    }
  } else if (len1 < len2) {
    lendiff <- (len2-len1)/2
    dif <- as.numeric(j_m[,1]) - as.numeric(j_m[,2])
    if (0 %in% dif) {
      # nl = as.numeric()
      # vt = as.numeric()
      # for (i in 1:length(dif)) {
      #   if (dif[i] == 0) {
      #     nl <- append(nl, as.numeric(j_m[i,1]))
      #     vt <- append(vt, as.numeric(I_nega[i]))
      #   }
      # }
      
      nl <- unlist(mcmapply(function(x, y) if (x==0) {as.numeric(y)}, dif, j_m[,1], mc.cores = numCore), recursive = FALSE)
      vt <- unlist(mcmapply(function(x, y) if (x==0) {as.numeric(y)}, dif, I_nega, mc.cores = numCore), recursive = FALSE)
      
      dt <- data.frame(nl,vt)
      if (length(vt) < lendiff) {
        res_node <- dt$vt
        I_posi_fliped = sort(append(I_posi_fliped, res_node))
        I_nega_fliped = I_nega_fliped[! I_nega_fliped %in% res_node]
        needmore <- lendiff - length(vt)
        # ds = as.numeric()
        # nd = as.numeric()
        # for (i in 1:length(dif)) {
        #   if (dif[i] > 0) {
        #     ds <- append(ds, as.numeric(j_m[i,1]))
        #     nd <- append(nd, as.numeric(I_nega[i]))
        #   }
        # }
        ds <- unlist(mcmapply(function(x, y) if (x>0) {as.numeric(y)}, dif, j_m[,1], mc.cores = numCore), recursive = FALSE)
        nd <- unlist(mcmapply(function(x, y) if (x>0) {as.numeric(y)}, dif, I_nega, mc.cores = numCore), recursive = FALSE)
        
        dt2 <- data.frame(ds,nd)
        more_node <- dt2$nd[1:needmore] # if multiple found, return the first several ones
        I_posi_fliped = sort(append(I_posi_fliped, more_node))
        I_nega_fliped = I_nega_fliped[! I_nega_fliped %in% more_node]
        
      } else if (length(vt) > lendiff) {
        dt_sort = dt[order(dt$nl, decreasing = TRUE),]
        res_node <- dt_sort$vt[1:lendiff] # if multiple found, return the first one
        I_posi_fliped = sort(append(I_posi_fliped, res_node))
        I_nega_fliped = I_nega_fliped[! I_nega_fliped %in% res_node]
      } else {
        res_node <- dt$vt
        I_posi_fliped = sort(append(I_posi_fliped, res_node))
        I_nega_fliped = I_nega_fliped[! I_nega_fliped %in% res_node]
      }
    } else {
      # nl = as.numeric()
      # vt = as.numeric()
      # for (i in 1:length(dif)) {
      #   if (dif[i] > 0) {
      #     nl <- append(nl, as.numeric(j_m[i,1]))
      #     vt <- append(vt, as.numeric(I_nega[i]))
      #   }
      # }
      nl <- unlist(mcmapply(function(x, y) if (x>0) {as.numeric(y)}, dif, j_m[,1], mc.cores = numCore), recursive = FALSE)
      vt <- unlist(mcmapply(function(x, y) if (x>0) {as.numeric(y)}, dif, I_nega, mc.cores = numCore), recursive = FALSE)
      
      dt <- data.frame(nl,vt)
      res_node <- dt$vt[1:lendiff]
      I_posi_fliped = sort(append(I_posi_fliped, res_node))
      I_nega_fliped = I_nega_fliped[! I_nega_fliped %in% res_node]
    }
  }
  
  #Creates vectors having data points
  expected <- ns %!in% sort(n1) * rep(1, n)
  predicted <- ns %!in% sort(I_posi_fliped) * rep(1, n)
  predicted_flip <- (predicted -1) * (-1) # another possible assignment
  
  expected_value <- factor(expected)
  predicted_value <- factor(predicted)
  predicted_flip_value <- factor(predicted_flip)
  
  # Creating confusion matrix
  confusion_res1 = confusionMatrix(data=predicted_value, reference = expected_value)
  confusion_res2 = confusionMatrix(data=predicted_flip_value, reference = expected_value)
  accu1 = confusion_res1$overall[1] # we are in equal community setting, so obtain accuracy
  accu2 = confusion_res2$overall[1]
  
  # Display results 
  return(1-max(accu1, accu2))
}

# create a list for simulation task
n_list <- as.numeric(rep(50,4))
alpha_list <- rep(c(0.9,0.8,0.7,0.6), 1)
p1_list <- rep(0.4, 4)
q1_list <- rep(0.1, 4)
p2_list <- rep(0.35, 4)
q2_list <- rep(0.15, 4)
sim_list <- data.frame(n_list,alpha_list,p1_list,q1_list,p2_list,q2_list)
nax <- c("n", "alpha", "p1", "q1", "p2", "q2")
colnames(sim_list) <- nax

# RepParallel <- function(n, expr, simplify = "array",...) {
#   answer <-
#     mclapply(integer(n), eval.parent(substitute(function(...) expr)), mc.cores = numCore, ...)
#   if (!identical(simplify, FALSE) && length(answer))
#     return(simplify2array(answer, higher = (simplify == "array")))
#   else return(answer)
# }

err_list1 = as.numeric()
sd_list1 = as.numeric()
for (i in 1:nrow(sim_list)) {
  mF = function(x, r1, r3) {
    ES_algo(sim_list[i,1], 3, sim_list[i,3], sim_list[i,4], x[[r1]], x[[r3]])
  }
  out<-replicate(50, Gen_hyperG(sim_list[i,1],3,sim_list[i,3],sim_list[i,4],sim_list[i,2],numCore))
  rt<-apply(out, 2, mF, r1='h', r3='n1')
  err_list1 = append(err_list1, mean(rt))
  sd_list1 = append(sd_list1, sd(rt))
}

err_list2 = as.numeric()
sd_list2 = as.numeric()
for (i in 1:nrow(sim_list)) {
  mF = function(x, r1, r3) {
    ES_algo(sim_list[i,1], 3, sim_list[i,5], sim_list[i,6], x[[r1]], x[[r3]])
  }
  out<-replicate(50, Gen_hyperG(sim_list[i,1],3,sim_list[i,5],sim_list[i,6],sim_list[i,2],numCore))
  rt<-apply(out, 2, mF, r1='h', r3='n1')
  err_list2 = append(err_list2, mean(rt))
  sd_list2 = append(sd_list2, sd(rt))
}

final_df <- data.frame(err_list1, sd_list1, err_list2, sd_list2)
write.table(final_df, 'results_50n_additional_p0.425.csv', col.names = c('err_pq1','sd_pq1', 'err_pq2','sd_pq2'),
            row.names = sim_list$n, sep = ',')
