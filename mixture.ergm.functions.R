########################################################################
########################################################################
########################################################################
################### This script is written for Bayesian mixtures of ERGMs
#########################################################################
#########################################################################
#########################################################################
#wants <- c("network", "ergm","coda", "Matrix", "mvtnorm",
#           "foreach", "doParallel",  "dplyr", "mclust", "mcclust", "gtools",
#           "flexclust", "ClusterR","Bergm", "igraph", "intergraph","brainGraph") 
wants <- c("network", "ergm", "coda", "Matrix", "mvtnorm",
           "foreach", "doParallel",  "dplyr", "mclust", "mcclust", "gtools",
           "flexclust", "ClusterR", "Bergm") 
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
lapply(wants, library, character.only=TRUE)
#########################################################################
expit <- function(x) 1/(1 + exp(-x))
#########################################################################
#########################################################################
logPL <- function(theta, y, X, weights, size, size.offset=TRUE){
  # theta_transf <- c(calibr.info$W %*% (theta - calibr.info$Theta_MLE) + 
  #                     calibr.info$Theta_PL)
  if(size.offset){
    xtheta <- c(X %*% theta)
    log.pl <- sum(dbinom(weights * y, weights, expit(-log(size) + xtheta), 
                         log = TRUE))
  }
  else{
    xtheta <- c(X %*% theta)
    log.pl <- sum(dbinom(weights * y, weights, expit(xtheta), 
                         log = TRUE)) 
  }
  return(log.pl)
}
##########################################################################
logL <- function(theta,
                 mplesetup=NULL, 
                 size.offset=TRUE,...){
  # true likelihood
  #browser()
  temp.dat <- mplesetup$dat
  direct.flag <- network::is.directed(temp.dat)
  size = mplesetup$size
  if(!size.offset){
    # if no krivitsky offset # LHS of formula is temp.dat
    log.l = ergm.bridge.llr(mplesetup$form,
                            from = rep(0,length(theta)), to = theta, 
                            control = control.ergm.bridge(MCMC.burnin = 4*size^2),
                            llronly = TRUE) + ifelse(direct.flag, -2*choose(size,2)*log(2), 
                                                     -choose(size,2)*log(2) )
    # MCMC.burnin = 4*size^2,
    # MCMC.interval = 2*size^2,
    
  }
  else{
    # update the formula to include an offset at the end
    offset.form <- statnet.common::nonsimp_update.formula(mplesetup$form, 
                                                          .~.+offset(edges), 
                                                          from.new = TRUE)
    # offset
    log.l = ergm.bridge.llr(offset.form,
                            from = rep(0, length(theta)+1), to = c(theta,-log(size)), 
                            control = control.ergm.bridge(MCMC.burnin = 4*size^2),
                            llronly = TRUE) + ifelse(direct.flag, -2*choose(size,2)*log(2), 
                                                     -choose(size,2)*log(2) )
  }
  return(log.l)
}
#################################################
logAPL <- function(theta, 
                   y, X, 
                   weights, size, 
                   size.offset=TRUE, 
                   theta_MLE,
                   theta_MPLE, W, 
                   logC, ...){
  # browser()
  # theta_transf <- c(calibr.info$W %*% (theta - calibr.info$Theta_MLE) + 
  #                     calibr.info$Theta_PL)
  if(size.offset){
    # browser()
    # xtheta <- c(X %*% theta)
    theta_transf <- theta_MPLE + W %*% (theta - theta_MLE)
    xtheta_transf <- c(X %*% theta_transf)
    log.apl <- logC + sum(dbinom(weights * y, weights, 
                                 expit(-log(size) + xtheta_transf), 
                                 log = TRUE))
  }
  else{
    # xtheta <- c(X %*% theta)
    theta_transf <- theta_MPLE + W %*% (theta - theta_MLE)
    xtheta_transf <- c(X %*% theta_transf)
    log.apl <- logC + sum(dbinom(weights * y, weights, 
                                 expit(xtheta_transf), 
                                 log = TRUE)) 
  }
  return(log.apl)
}
#################################################
# logL <- function(theta, y, X, weights, 
#                  size, 
#                  mplesetup_list=NULL, 
#                  size.offset=TRUE, pl=TRUE,...){
#  
#   if(pl){
#    if(size.offset){
#      xtheta <- c(X %*% theta)
#      log.pl <- sum(dbinom(weights * y, weights, expit( -log(size) + xtheta), 
#                          log = TRUE))
#    }
#    else{
#      xtheta <- c(X %*% theta)
#      log.pl <- sum(dbinom(weights * y, weights, expit(xtheta), 
#                          log = TRUE)) 
#    }
#     log.l=log.pl #use pseudo-likelihood approximation
#   }
#   else{
#     # true likelihood
#     #browser()
#     temp.dat <- mplesetup_list$dat
#     direct.flag <- network::is.directed(temp.dat)
#     if(!size.offset){
#       # if no krivitsky offset # LHS of formula is temp.dat
#       log.l = ergm.bridge.llr(mplesetup_list$form,
#                               from = rep(0,length(theta)), to = theta, 
#                               control = control.ergm.bridge(MCMC.burnin = 4*size^2,
#                                                             MCMC.interval = 2*size^2,
#                                                             MCMC.samplesize = 10000),
#                               llronly = TRUE) + ifelse(direct.flag, -2*choose(size,2)*log(2), 
#                                                       -choose(size,2)*log(2) )
#       
#     }
#     else{
#       # update the formula to include an offset at the end
#       offset.form <- nonsimp_update.formula(mplesetup_list$form, 
#                                             .~.+offset(edges), 
#                                             from.new = TRUE)
#       log.l = ergm.bridge.llr(offset.form,
#                               from = rep(0, length(theta)+1), to = c(theta,-log(size)), 
#                               control = control.ergm.bridge(MCMC.burnin = 4*size^2,
#                                                             MCMC.interval = 2*size^2,
#                                                             MCMC.samplesize = 10000),
#                               llronly = TRUE) + ifelse(direct.flag, -2*choose(size,2)*log(2), 
#                                                       -choose(size,2)*log(2) )
#     }
#   }
#   return(log.l)
# }
########################################################################
### Inverse Fisher Z transformation
########################################################################
inv.FZT <- function(z){
  r <- (exp(2*z)-1) / (exp(2*z)+1)
  return(r)
}
#########################################################################
################# Thresholding weighted matrix to binarize it
#########################################################################
w2b_thr <- function(w, thr=0.2,equal=TRUE,...){
  # input : a weighted matrix, a threshold
  # output : a binary matrix calculated based on 
  # thresholding the weighted matrix
  b <- w
  if(equal){
    b[which(w>=thr, arr.ind = T)] <- 1
  } else {
    b[which(w>thr, arr.ind = T)] <- 1
  }
  b[which(w<thr, arr.ind = T)] <- 0
  
  return(b)
}
######################################################
w2b_mdeg <- function(w, mdeg=NULL,equal=FALSE,...){
  # convert a weighted matrix into binary matrix given threshold or mean deg constraint
  # using binary search to identify the threshold according to mean degree
  # constraint
  nsize <- nrow(w)
  # assuming it is undirected
  # stop when the lb and ub lead to identical binary adjacency matrix
  lb <- min(w, na.rm=TRUE)
  ub <- max(w, na.rm = TRUE)
  temp_thr <- (lb+ub)/2
  b2 <- w2b_thr(w, thr = lb,equal=equal)
  b1 <- w2b_thr(w, thr = ub,equal=equal) 
  
  # 
  while(sum(b2 - b1)>0){
    # 
    temp_b <- w2b_thr(w, thr = temp_thr,equal=equal)
    #
    cat("Current threshold : ", temp_thr, "\n")
    # calculate mean degree
    temp_b_mdeg <- mean(sna::degree(temp_b, gmode = "graph", cmode="indegree")) 
    # default : adjacency matrix undirected
    if(abs(temp_b_mdeg - mdeg)<=(2/nsize)){
      break
    }
    # check if the constraint is larger than temp_b_mdeg?
    if(mdeg >= temp_b_mdeg){
      # need more edges, decrease upper bound 
      ub <- temp_thr
    } else{
      # need less edges, increase lower bound
      lb <- temp_thr
    }
    temp_thr <- (lb+ub)/2
    b2 <- w2b_thr(w, thr = lb, equal = equal)
    b1 <- w2b_thr(w, thr = ub, equal = equal) 
  }
  
  return(list(b=temp_b, thr = temp_thr, mdeg = mdeg) )
}
#########################################################################
######## clustering networks based on mple estimates, return ordered labels and estimates
#########################################################################
####### mple_kmeans
#########################################################################
mple_kmeans <- function(dat, form, p,K,num.cores=1,size.offset = TRUE,center.flag=TRUE,scale.flag=TRUE,...){
  require(doParallel)
  # sample size
  n <- length(dat)
  # network size
  netsize <- do.call("c", lapply(dat, function(x) network::network.size(x)))
  # mple
  if(num.cores>1){
    registerDoParallel(cores = num.cores)
    mple_coef_list <- foreach(i = 1:n) %dopar% {
      temp_dat <- dat[[i]]
      temp_form <- statnet.common::nonsimp_update.formula(form, temp_dat~., from.new = TRUE)
      temp_mple <- ergm(temp_form, estimate = "MPLE")
      temp_mple$coef
    }
  } else{
    mple_coef_list <- lapply(dat, function(x) ergm((statnet.common::nonsimp_update.formula(form, x~., from.new = TRUE)), 
                                                   estimate = "MPLE")$coef )
    
  }
  mple_coef <- do.call(rbind, mple_coef_list)
  
  if(size.offset){
    mple_coef[,1] <- mple_coef[,1] + log(netsize)
  }
  
  cl <- kmeans(scale(mple_coef,center = center.flag, scale=scale.flag), centers = K)
  # labels, cl$cluster
  # the cluster with smallest value on first element of the parameter vector 
  # got the smallest label
  # 
  mple_coef_cluster <- data.frame(mple_coef, clabel = cl$cluster)
  # mean of edges terms, assume there is edges term
  
  mple_coef_cluster_summary <- mple_coef_cluster %>% group_by(clabel) %>% summarise(medges = mean(edges)) %>% arrange(desc(medges))
  mple_coef_cluster_summary$olabel <- K:1
  # browser()
  # ord <- mple_coef_cluster_summary$clabel[order(mple_coef_cluster_summary$medges, decreasing = FALSE)] 
  
  colabel <- data.frame(clabel = mple_coef_cluster_summary$clabel, olabel = mple_coef_cluster_summary$olabel)
  # data.frame(clabel = mple_coef_cluster_summary$clabel, ord = ord, olabel=1:3)
  
  # merge the two tables by clabel and then keep the olabel
  
  mple_coef_colabel <- left_join(mple_coef_cluster, colabel, by = c("clabel"= "clabel"))
  # ord is the true order
  return(mple_coef_colabel)
}
#############################################################################################
#############################################################################################
#############################################################################################
#########################################################################
########################## Posterior classification probability
########################## input : posterior samples, mple info, samp size
########################## output : Monte Carlo approximation to the posterior class probabilites
########################################################################
class_prob <- function(mplesetup_list, 
                       tau, 
                       theta, 
                       samp.size, 
                       size.offset=TRUE,
                       lik=c("PL","APL","full"),
                       num.cores=1,
                       theta_MLE_list=NULL,
                       theta_MPLE_list=NULL,
                       W_list=NULL,
                       logC_list=NULL,...){
  # 
  n <- length(mplesetup_list)
  K <- ncol(tau)
  #########################################
  calc.id <- 1:nrow(tau)
  # if(nrow(tau) >= samp.size){
  #   calc.id <- sample(1:nrow(tau), size=samp.size,
  #                         replace = FALSE)
  # } else{
  #   calc.id <- 1:nrow(tau)
  # }
  # 
  log_prob_array <- array(NA, dim=c(length(calc.id), n, K)) # each entry equals to log(eta_j^m * f(yi, theta_j^m))
  if(num.cores==1){
    for(m in 1:length(calc.id)){
      for(i in 1:n){
        for(k in 1:K){
          if(match.arg(lik) == "PL"){
            temp_log_lik <- logPL(theta = theta[calc.id[m],,k],
                                  y= mplesetup_list[[i]]$response,
                                  X = mplesetup_list[[i]]$predictor,
                                  weights = mplesetup_list[[i]]$weights,
                                  size = mplesetup_list[[i]]$size,
                                  size.offset = size.offset)
            # logPL <- function(theta, y, X, weights, size, size.offset=TRUE)
          } else if(match.arg(lik) == "APL"){
            temp_log_lik <- logAPL(theta = theta[calc.id[m],,k],
                                   y= mplesetup_list[[i]]$response,
                                   X = mplesetup_list[[i]]$predictor,
                                   weights = mplesetup_list[[i]]$weights,
                                   size = mplesetup_list[[i]]$size,
                                   size.offset = size.offset,
                                   theta_MLE = theta_MLE_list[[i]],
                                   theta_MPLE = theta_MPLE_list[[i]],
                                   W = W_list[[i]],
                                   logC = logC_list[[i]])
            # logAPL <- function(
            #                    theta, 
            #                    y, X, 
            #                    weights, size, 
            #                    size.offset=TRUE, 
            #                    theta_MLE,
            #                    theta_MPLE, W, 
            #                    logC, ...)
          } else if(match.arg(lik) == "full"){
            temp_log_lik <- logL(theta = theta[calc.id[m],,k],
                                 mplesetup = mplesetup_list[[i]],
                                 size.offset = size.offset)
            # logL <- function(theta,
            #                  mplesetup=NULL, 
            #                  size.offset=TRUE,...)
          }
          
          log_prob_array[m,i,k] <- log(tau[calc.id[m],k]) + temp_log_lik
          
          # logL(theta = theta[calc.id[m],,k],
          #                                                     y= mplesetup_list[[i]]$response,
          #                                                     X = mplesetup_list[[i]]$predictor,
          #                                                     weights = mplesetup_list[[i]]$weights,
          #                                                     mplesetup_list = mplesetup_list[[i]],
          #                                                     size = mplesetup_list[[i]]$size,
          #                                                     size.offset = size.offset, pl = pl)
        }
      }
    }
  } else if(num.cores>1){
    for(m in 1:length(calc.id)){
      log_prob_array[m,,] <- foreach(i = 1:n, .combine = rbind)%dopar%{
        temp_vec <- rep(NA, K)
        for(k in 1:K){
          
          if(match.arg(lik) == "PL"){
            temp_log_lik <- logPL(theta = theta[calc.id[m],,k],
                                  y= mplesetup_list[[i]]$response,
                                  X = mplesetup_list[[i]]$predictor,
                                  weights = mplesetup_list[[i]]$weights,
                                  size = mplesetup_list[[i]]$size,
                                  size.offset = size.offset)
            # logPL <- function(theta, y, X, weights, size, size.offset=TRUE)
          } else if(match.arg(lik) == "APL"){
            temp_log_lik <- logAPL(theta = theta[calc.id[m],,k],
                                   y= mplesetup_list[[i]]$response,
                                   X = mplesetup_list[[i]]$predictor,
                                   weights = mplesetup_list[[i]]$weights,
                                   size = mplesetup_list[[i]]$size,
                                   size.offset = size.offset,
                                   theta_MLE = theta_MLE_list[[i]],
                                   theta_MPLE = theta_MPLE_list[[i]],
                                   W = W_list[[i]],
                                   logC = logC_list[[i]])
            # logAPL <- function(
            #                    theta, 
            #                    y, X, 
            #                    weights, size, 
            #                    size.offset=TRUE, 
            #                    theta_MLE,
            #                    theta_MPLE, W, 
            #                    logC, ...)
          } else if(match.arg(lik) == "full"){
            temp_log_lik <- logL(theta = theta[calc.id[m],,k],
                                 mplesetup = mplesetup_list[[i]],
                                 size.offset = size.offset)
            # logL <- function(theta,
            #                  mplesetup=NULL, 
            #                  size.offset=TRUE,...)
          }
          
          # log_prob_array[m,i,k] <- log(eta[calc.id[m],k]) + temp_log_lik
          # log_prob_array[m,i,k]
          temp_vec[k] <- log(tau[calc.id[m],k]) + temp_log_lik
        }
        temp_vec
      }
    }
  }
  # convert log_prob_array to prob_matrix
  prob_matrix <- matrix(NA, nrow = n, ncol = K)
  for(i in 1:n){
    for(k in 1:K){
      prob_matrix[i,k] <- mean(exp(log_prob_array[,i,k] - apply(matrix(log_prob_array[,i,],ncol=K),
                                                                1,
                                                                function(x) sna::logSum(x[which(x != -Inf)]) ))) 
    }
  }
  return(prob_matrix)
}
########################################################################
########################################################################
########################################################################
DIC_mixture_ergm <- function(mplesetup_list, 
                             tau, 
                             theta, 
                             size.offset=TRUE, 
                             lik = c("PL", "APL", "full"),
                             theta_MLE_list=NULL,
                             theta_MPLE_list=NULL,
                             W_list=NULL,
                             logC_list=NULL,...){
  #####
  n <- length(mplesetup_list)
  K <- ncol(tau)
  calc.id <- 1:nrow(tau)
  # browser()
  # the theta and eta are all thinned
  # figure out the index of the samples used to calculate DIC
  # if(nrow(tau) >= DIC.samp.size){
  #   calc.id <- sample(1:nrow(tau), size=DIC.samp.size,
  #                         replace = FALSE)
  # } else{
  #   calc.id <- 1:nrow(tau)
  # }
  # the matrix for stroing intermediat values for the calculation of posterior expectation of loglikelihood 
  DIC.matrix <- matrix(NA, 
                       nrow = length(calc.id),
                       ncol = n)
  #
  for(l in 1:length(calc.id)){
    for(i in 1:n){
      temp_val <- rep(NA, length=K)
      for(k in 1:K){
        if(match.arg(lik) == "PL"){
          temp_log_lik <- logPL(theta = theta[calc.id[l],,k],
                                y= mplesetup_list[[i]]$response,
                                X = mplesetup_list[[i]]$predictor,
                                weights = mplesetup_list[[i]]$weights,
                                size = mplesetup_list[[i]]$size,
                                size.offset = size.offset)
          # logPL <- function(theta, y, X, weights, size, size.offset=TRUE)
        } else if(match.arg(lik) == "APL"){
          temp_log_lik <- logAPL(theta = theta[calc.id[l],,k],
                                 y= mplesetup_list[[i]]$response,
                                 X = mplesetup_list[[i]]$predictor,
                                 weights = mplesetup_list[[i]]$weights,
                                 size = mplesetup_list[[i]]$size,
                                 size.offset = size.offset,
                                 theta_MLE = theta_MLE_list[[i]],
                                 theta_MPLE = theta_MPLE_list[[i]],
                                 W = W_list[[i]],
                                 logC = logC_list[[i]])
          # logAPL <- function(
          #                    theta, 
          #                    y, X, 
          #                    weights, size, 
          #                    size.offset=TRUE, 
          #                    theta_MLE,
          #                    theta_MPLE, W, 
          #                    logC, ...)
        } else if(match.arg(lik) == "full"){
          temp_log_lik <- logL(theta = theta[calc.id[l],,k],
                               mplesetup = mplesetup_list[[i]],
                               size.offset = size.offset)
          # logL <- function(theta,
          #                  mplesetup=NULL, 
          #                  size.offset=TRUE,...)
        }
        temp_val[k] <- log(tau[calc.id[l],k]) + temp_log_lik
        
        
        # logL(theta = theta[DIC.calc.id[l],,k], 
        #                                                     y= mplesetup_list[[i]]$response,
        #                                                     X = mplesetup_list[[i]]$predictor, 
        #                                                     weights = mplesetup_list[[i]]$weights,
        #                                                     size = mplesetup_list[[i]]$size, 
        #                                                     size.offset = size.offset,
        #                                                     mplesetup_list = mplesetup_list[[i]],
        #                                                     pl=pl)
        
      }
      DIC.matrix[l,i] <- sna::logSum(temp_val[which(temp_val != -Inf)]) 
      # which(unnorm_logp!=-Inf)
    }
  }
  Dbar <- (-2)*sum(DIC.matrix)/length(calc.id)
  cat("Dbar equals to", Dbar, "\n")
  # calculating marginal density, using logSum function
  log.marginal.den.vec <- rep(NA, length = n)
  for(i in 1:n){
    temp_mat <- matrix(NA, nrow = length(calc.id),
                       ncol = K)
    for(l in 1:length(calc.id)){
      for(k in 1:K){
        if(match.arg(lik) == "PL"){
          temp_log_lik <- logPL(theta = theta[calc.id[l],,k],
                                y= mplesetup_list[[i]]$response,
                                X = mplesetup_list[[i]]$predictor,
                                weights = mplesetup_list[[i]]$weights,
                                size = mplesetup_list[[i]]$size,
                                size.offset = size.offset)
          # logPL <- function(theta, y, X, weights, size, size.offset=TRUE)
        } else if(match.arg(lik) == "APL"){
          temp_log_lik <- logAPL(theta = theta[calc.id[l],,k],
                                 y= mplesetup_list[[i]]$response,
                                 X = mplesetup_list[[i]]$predictor,
                                 weights = mplesetup_list[[i]]$weights,
                                 size = mplesetup_list[[i]]$size,
                                 size.offset = size.offset,
                                 theta_MLE = theta_MLE_list[[i]],
                                 theta_MPLE = theta_MPLE_list[[i]],
                                 W = W_list[[i]],
                                 logC = logC_list[[i]])
          # logAPL <- function(
          #                    theta, 
          #                    y, X, 
          #                    weights, size, 
          #                    size.offset=TRUE, 
          #                    theta_MLE,
          #                    theta_MPLE, W, 
          #                    logC, ...)
        } else if(match.arg(lik) == "full"){
          temp_log_lik <- logL(theta = theta[calc.id[l],,k],
                               mplesetup = mplesetup_list[[i]],
                               size.offset = size.offset)
          # logL <- function(theta,
          #                  mplesetup=NULL, 
          #                  size.offset=TRUE,...)
        }
        
        temp_mat[l,k] <- log(tau[calc.id[l],k]/length(calc.id)) + temp_log_lik
        
        
        
        # logL(theta = theta[DIC.calc.id[l],,k], 
        #                                                                       y= mplesetup_list[[i]]$response,
        #                                                                       X = mplesetup_list[[i]]$predictor, 
        #                                                                       weights = mplesetup_list[[i]]$weights,
        #                                                                      size = mplesetup_list[[i]]$size, 
        #                                                                      size.offset = size.offset,
        #                                                                     mplesetup_list = mplesetup_list[[i]],
        #                                                                     pl = pl) 
      }
    }
    log.marginal.den.vec[i] <- sna::logSum(temp_mat[which(temp_mat != -Inf, arr.ind = T)])
  }
  DIC <- 2*Dbar + 2*sum(log.marginal.den.vec)
  cat("DIC3 equals to", DIC, "\n")
  return(DIC)
}
############################################################################
###### Bayesian information criterion
#############################################################################
BIC_mixture_ergm <- function(mplesetup_list,
                             p,
                             tau, 
                             theta, 
                             size.offset=TRUE, 
                             lik = c("PL", "APL", "full"),
                             theta_MLE_list=NULL,
                             theta_MPLE_list=NULL,
                             W_list=NULL,
                             logC_list=NULL,...){
  #####
  n <- length(mplesetup_list)
  K <- ncol(tau)
  
  ##### 2*lhat - (G*M + G -1)*log(N)
  ##### lhat calculated from posterior mean
  ##### N: number of networks
  ##### G: number of groups
  ##### M: number of summary statistics
  
  ##### figure out the posterior mean estimate
  theta_post_mean <- apply(theta, c(2,3), mean)
  tau_post_mean <- apply(tau, 2, mean)
  ##### calculate the loglikelihood evaluated at the theta_post_mean
  #####
  BIC.vec <- rep(NA, n)
  for(i in 1:n){
    temp_val <- rep(NA, length=K)
    for(k in 1:K){
      if(match.arg(lik) == "PL"){
        temp_log_lik <- logPL(theta = theta_post_mean[,k],
                              y= mplesetup_list[[i]]$response,
                              X = mplesetup_list[[i]]$predictor,
                              weights = mplesetup_list[[i]]$weights,
                              size = mplesetup_list[[i]]$size,
                              size.offset = size.offset)
        # logPL <- function(theta, y, X, weights, size, size.offset=TRUE)
      } else if(match.arg(lik) == "APL"){
        temp_log_lik <- logAPL(theta = theta_post_mean[,k],
                               y= mplesetup_list[[i]]$response,
                               X = mplesetup_list[[i]]$predictor,
                               weights = mplesetup_list[[i]]$weights,
                               size = mplesetup_list[[i]]$size,
                               size.offset = size.offset,
                               theta_MLE = theta_MLE_list[[i]],
                               theta_MPLE = theta_MPLE_list[[i]],
                               W = W_list[[i]],
                               logC = logC_list[[i]])
        # logAPL <- function(
        #                    theta, 
        #                    y, X, 
        #                    weights, size, 
        #                    size.offset=TRUE, 
        #                    theta_MLE,
        #                    theta_MPLE, W, 
        #                    logC, ...)
      } else if(match.arg(lik) == "full"){
        temp_log_lik <- logL(theta = theta_post_mean[,k],
                             mplesetup = mplesetup_list[[i]],
                             size.offset = size.offset)
        # logL <- function(theta,
        #                  mplesetup=NULL, 
        #                  size.offset=TRUE,...)
      }
      temp_val[k] <- log(tau[k]) + temp_log_lik
      # logL(theta = theta[DIC.calc.id[l],,k], 
      #                                                     y= mplesetup_list[[i]]$response,
      #                                                     X = mplesetup_list[[i]]$predictor, 
      #                                                     weights = mplesetup_list[[i]]$weights,
      #                                                     size = mplesetup_list[[i]]$size, 
      #                                                     size.offset = size.offset,
      #                                                     mplesetup_list = mplesetup_list[[i]],
      #                                                     pl=pl)
      
    }
    BIC.vec[i] <- sna::logSum(temp_val[which(temp_val != -Inf)]) 
  }
  ###
  return(-(2*sum(BIC.vec) - (K*p + K - 1)*log(n) ))
}
############################################################################
############################################################################
##### calculate the likelihood of a network under the mixture of ergms model
##### given parameter values
##### 
############################################################################
############################################################################
loglik_mixture_ergm <- function(dat,
                                mplesetup,
                                p,
                                tau, 
                                theta, 
                                size.offset=TRUE, 
                                lik = c("PL", "APL", "full"),
                                theta_MLE=NULL,
                                theta_MPLE=NULL,
                                W=NULL,
                                logC=NULL,...){
  
  # number of clusters
  K <- length(tau)
  # tau is a vector of length K
  # theta is a matrix of size p * K, with each row representing a specific component
  temp_val <- rep(NA, length=K)
  for(k in 1:K){
    if(match.arg(lik) == "PL"){
      temp_log_lik <- logPL(theta = theta[,k],
                            y= mplesetup$response,
                            X = mplesetup$predictor,
                            weights = mplesetup$weights,
                            size = mplesetup$size,
                            size.offset = size.offset)
      # logPL <- function(theta, y, X, weights, size, size.offset=TRUE)
    } else if(match.arg(lik) == "APL"){
      temp_log_lik <- logAPL(theta = theta[,k],
                             y= mplesetup$response,
                             X = mplesetup$predictor,
                             weights = mplesetup$weights,
                             size = mplesetup$size,
                             size.offset = size.offset,
                             theta_MLE = theta_MLE,
                             theta_MPLE = theta_MPLE,
                             W = W,
                             logC = logC)
      # logAPL <- function(
      #                    theta, 
      #                    y, X, 
      #                    weights, size, 
      #                    size.offset=TRUE, 
      #                    theta_MLE,
      #                    theta_MPLE, W, 
      #                    logC, ...)
    } else if(match.arg(lik) == "full"){
      temp_log_lik <- logL(theta = theta[,k],
                           mplesetup = mplesetup,
                           size.offset = size.offset)
      # logL <- function(theta,
      #                  mplesetup=NULL, 
      #                  size.offset=TRUE,...)
    }
    temp_val[k] <- log(tau[k]) + temp_log_lik
    # logL(theta = theta[DIC.calc.id[l],,k], 
    #                                                     y= mplesetup_list[[i]]$response,
    #                                                     X = mplesetup_list[[i]]$predictor, 
    #                                                     weights = mplesetup_list[[i]]$weights,
    #                                                     size = mplesetup_list[[i]]$size, 
    #                                                     size.offset = size.offset,
    #                                                     mplesetup_list = mplesetup_list[[i]],
    #                                                     pl=pl)
    
  }
  return(sna::logSum(temp_val[which(temp_val != -Inf)]))
}

############################################################################
#### bayesian mixture of ERGMs using pseudolikelihood
# use pseudolikelihood when updating membership indicators and calculate DICs
# use exchange step when updating cluster-specific parameters
###########################################################################
bayes_mixture_ergm <- function(dat, form, p, K, prior.mean=c(-1,rep(0,p-1)), 
                               prior.sigma=diag(25,p), alpha=rep(3,K),
                               main.iters=1000,burn.in=(main.iters)/2,
                               thin=1, exchange.step = FALSE, aux.iters = 1000,
                               sigma.epsilon=diag(0.001,p), pre.assign=NULL,
                               num.cores=1,
                               DIC.calc=TRUE,DIC.samp.size=500, post.class.calc=TRUE,
                               random.init=TRUE, init.method = "mple", order.constraints = TRUE,
                               edge.init.value=-2, size.offset = TRUE, pl=TRUE,
                               ...){
  init.start.time <- Sys.time()
  # sample size
  n <- length(dat)
  # netsize <- network::network.size(dat[[1]])
  netsize <- do.call(c,lapply(dat, function(x) network::network.size(x) ) )
  # sufficient stats for each observed data
  sy <- matrix(NA, nrow = n, ncol = p)
  mplesetup_list <- vector(mode="list", length = n)
  data.glm.initial_list <- vector(mode="list", length=n)
  # mstats_list <- vector(mode="list", length=n)
  for(i in 1:n){
    # sufficient statistics
    temp.dat <- dat[[i]]
    temp.form <- statnet.common::nonsimp_update.formula(form, temp.dat~., from.new = TRUE) # LHS of formula is temp.dat
    sy[i,] <- summary_formula(temp.form)
    mplesetup <- ergmMPLE(temp.form)
    mplesetup$size <- netsize[i]
    mplesetup$form <- temp.form
    mplesetup$dat <- dat[[i]]
    mplesetup_list[[i]] <- mplesetup
    data.glm.initial_list[[i]] <- cbind(mplesetup$response, mplesetup$weights, 
                                        mplesetup$predictor)
    colnames(data.glm.initial_list[[i]]) <- c("responses", "weights", colnames(mplesetup$predictor))
  }
  # we have got sy as the 
  # Initialize the matrix for storing sampled parameters
  eta <- matrix(NA, nrow = 1+main.iters, ncol = K)
  theta <- array(NA, dim = c(1+main.iters, p, K))
  Z <- matrix(NA, nrow = 1+main.iters, ncol = n)
  acc.counts <- matrix(FALSE, nrow=1+main.iters, ncol = K)
  #browser()
  mple_coef_colabel <- NA
  # merge the two tables by clabel and then keep the olabel
  if(!random.init){
    if(init.method=="mple"){
      mple_coef_colabel <- mple_kmeans(dat, form = form, p=p, K=K, num.cores = num.cores, size.offset = size.offset)
      Z[1,] <- mple_coef_colabel$olabel
      eta[1,] <- table(factor(Z[1,], levels = 1:K))/n
      for(k in 1:K){
        # k is the index for subgroups
        theta[1,,k] <- apply(filter(mple_coef_colabel, olabel == k),2,mean)[1:p]
      }
    }
    else if(init.method=="stats_scaled") {
      Z[1,] <- kmeans( scale(sy, center = TRUE, scale = TRUE), centers = K)$cluster
      
      eta[1,] <- table(factor(Z[1,], levels = 1:K))/n
      for(k in 1:K){
        theta[1,,k] <- runif(n=p, min=-0.1,max=0.1)
        # apply(filter(mple_coef_colabel, olabel == k),2,mean)[1:p]
      }
    } else if(init.method == "stats_mdist"){
      # Z[1,] <- kmeans( scale(sy, center = TRUE, scale = TRUE), centers = K)$cluster
      # use hierarchical clustering
      mean_sy <- apply(sy, 2, mean)
      cov_sy <- cov(sy)
      mdist <- mahalanobis(sy, center = mean_sy, cov = cov_sy)
      # hierarchical clustering
      hclust_sy <- hclust(mdist, method = "ward.D2")
      Z[1,] <- cutree(hclust_sy, k=K)
      
      eta[1,] <- table(factor(Z[1,], levels = 1:K))/n
      for(k in 1:K){
        theta[1,,k] <- runif(n=p, min=-0.1,max=0.1)
        # apply(filter(mple_coef_colabel, olabel == k),2,mean)[1:p]
      }
    }
    ############ old method
    # mple_coef_colabel <- mple_kmeans(dat, form = form, p=p, K=K, num.cores = num.cores, size.offset = size.offset)
    # # left_join(mple_coef_cluster, colabel, by = c("clabel"= "clabel"))
    # # Initialize the Z^0 vector and parameters for \theta_j^0, eta0
    # Z[1,] <- mple_coef_colabel$olabel
    # eta[1,] <- table(factor(mple_coef_colabel$olabel, levels = 1:K))/n
    # Z[1,] <- kmeans( scale(sy, center = TRUE, scale = TRUE), centers = K)$cluster
    # for(k in 1:K){
    #   # k is the index for subgroups
    #   theta[1,,k] <- apply(filter(mple_coef_colabel, olabel == k),2,mean)[1:p]
    # }
    ############ old method
    
    ## new method
    # eta[1,] <- table(factor(Z[1,], levels = 1:K))/n
    # for(k in 1:K){
    #   theta[1,,k] <- runif(n=p, min=-0.1,max=0.1)
    #   # apply(filter(mple_coef_colabel, olabel == k),2,mean)[1:p]
    # }
    ## new method
  }
  else {
    # random initialization
    Z[1,] <- sample(x=1:K,size=n,replace = TRUE,prob=alpha/sum(alpha))
    eta[1,] <- alpha/sum(alpha)
    # 
    for(k in 1:K){
      theta[1,,k] <- runif(n=p, min=-0.1,max=0.1)
      # apply(filter(mple_coef_colabel, olabel == k),2,mean)[1:p]
    }
    theta[1,1,] <- theta[1,1,]+edge.init.value
  }
  
  # if pre.assign not null, use pre-assignments
  if(!is.null(pre.assign)){
    Z[1,] <- pre.assign
  }
  
  # let edge parameter start at somewhere that is negative and order them to make sure the first cluster has the
  # smallest initial value on edge parameter
  ########
  theta_order <- order(theta[1,1,], decreasing = FALSE)
  theta[1,,] <- theta[1,,theta_order]
  #########
  # metropolis within gibbs sampler
  # online algorithm to deal with label switching
  # start updating the algorithm from j = 2
  j = 2
  # main iterations
  samp.start.time <- Sys.time()
  if(exchange.step){
    # setup the list for simulating from ergms
    # y <- ergm.getnetwork(mod2$formula)
    control <- control.ergm(MCMC.burnin = aux.iters, MCMC.interval = 1, 
                            MCMC.samplesize = 1)
    # model_list <- vector(mode="list", length=n)
    Clist_list <- vector(mode = "list",length=n)
    proposal_list <- vector(mode = "list", length=n)
    # sy is a matrix
    for(i in 1:n){
      temp_model <- ergm::ergm_model(form, dat[[i]])
      Clist_list[[i]] <- ergm::ergm.Cprepare(dat[[i]], temp_model)
      proposal_list[[i]] <- ergm::ergm_proposal(object = ~., constraints = ~., 
                                                arguments = control$MCMC.prop.args,
                                                nw = dat[[i]])
    }
  }
  while(j <= 1+main.iters){
    # update the assignment label (total of n labels)
    for(i in 1:n){
      unnorm_p <- rep(NA, length = K)
      unnorm_logp <- rep(NA, length = K)
      for(k in 1:K){
        # ignore the constant C
        cat("Previous eta", eta[j-1,k], "\n")
        cat("Previous theta", theta[j-1,,k], "\n")
        temp_logPL <- logPL(theta = c(theta[j-1,,k]), y=mplesetup_list[[i]]$response,
                            X = mplesetup_list[[i]]$predictor, weights = mplesetup_list[[i]]$weights,
                            size = mplesetup_list[[i]]$size, size.offset = size.offset)
        cat("Corrected pseudo-like", temp_logPL, "\n")
        unnorm_logp[k] <- log(eta[j-1,k])  + temp_logPL
        cat("unnormalized log prob", unnorm_logp[k], "\n")
      }
      # A = max(unnorm_logp)
      lsum = sna::logSum(unnorm_logp)
      unnorm_p <- exp(unnorm_logp-lsum)
      cat("The unnormalized probability", unnorm_p, "\n")
      # Update Z...
      Z[j,i] <- sample(x=1:K, size=1, prob = unnorm_p)
    }
    
    # update eta (from dirichlet), now we have Z[j,]
    eta[j,] <- rdirichlet(n=1, alpha = (alpha+table(factor(Z[j,], levels = 1:K))) ) # K-dimensional vector
    cat("eta", eta[j,], "\n")
    
    # update theta, using random walk metropolis update or exchange algorithm
    for(k in 1:K){
      # proposing a new theta for k-th component
      theta1 <- rmvnorm(1,mean=theta[j-1,,k], sigma = sigma.epsilon)
      # log prior vector
      pr <- dmvnorm( rbind(theta1, theta[j-1,,k]), mean= prior.mean, sigma = prior.sigma, log=TRUE)
      # adjusted pseudolikelihood vector
      # lr <- c( logPL.corr(theta = theta[j-1,,k]) )
      k_id <- which(Z[j,] == k) # index of observed data whose corresponding assigned cluster is k 
      cat("The",k,"-th cluster:",k_id, "\n")
      if(length(k_id) >= 1){
        ll_mat <- matrix(NA, nrow=2,ncol = length(k_id))
        for(l in 1:length(k_id)){
          # network index k_id[l]
          if(!exchange.step){
            ll_mat[2,l] <- logPL(theta = theta[j-1,,k], y=mplesetup_list[[k_id[l]]]$response,
                                 X = mplesetup_list[[k_id[l]]]$predictor, weights = mplesetup_list[[k_id[l]]]$weights,
                                 size = mplesetup_list[[k_id[l]]]$size, size.offset = size.offset)
            ll_mat[1,l] <- logPL(theta = c(theta1), y=mplesetup_list[[k_id[l]]]$response,
                                 X = mplesetup_list[[k_id[l]]]$predictor, weights = mplesetup_list[[k_id[l]]]$weights,
                                 size = mplesetup_list[[k_id[l]]]$size, size.offset = size.offset)
          }
          else{
            # exchange step for updating cluster-specific parameters, the proposed theta for cluster k is theta1
            # simulation step, simulate network from theta1
            ## figure out network size
            temp_coef <- theta1
            
            if(size.offset){
              temp_size <- network::network.size(dat[[ k_id[l] ]])
              temp_coef[1] <- theta1[1] - log(temp_size)
            }
            
            temp_stats <- ergm_MCMC_slave(Clist = Clist_list[[ k_id[l] ]], 
                                          proposal = proposal_list[[ k_id[l] ]], 
                                          eta = temp_coef, 
                                          control = control, 
                                          verbose = FALSE)$s + sy[k_id[l],] # set up them first
            # obs stats is named sy, which is of size n by p
            # likelihood ratio calculation step
            ll_mat[1,l] <- sum(theta[j-1,,k] * temp_stats) + sum( theta1 * sy[k_id[l],] )
            ll_mat[2,l] <- sum(theta[j-1,,k] * sy[k_id[l],]) + sum( theta1 * temp_stats )
          }
        }
        lr <- apply(ll_mat,1,sum) # row-sum
      }
      else {
        lr <- c(0,0)
      }
      # log acceptance ratio, proposal is symmetric
      beta <- (lr[1] - lr[2]) + (pr[1] - pr[2]) 
      # accept if greater than acceptance ratio 
      if(beta >= log(runif(1))){
        theta[j,,k] <- theta1
        acc.counts[j,k] <- TRUE
      }
      else {
        theta[j,,k] <- theta[j-1,,k]
      }
    }
    # re-order parameters to deal with label switching issue, from small to large, using the order constraints on 
    # edge parameter
    if(order.constraints){
      theta_order <- order(theta[j,1,], decreasing = FALSE)
      theta[j,,] <- theta[j,,theta_order]
      eta[j,] <- eta[j,theta_order]
      # also need to permute the latent cluster assignment... 
      temp_Z <- data.frame(old_label = Z[j,])
      reference_Z <- data.frame(old_label=1:K, correct_label = theta_order)
      correct_Z <- left_join(temp_Z, reference_Z, by = c("old_label" = "old_label"))
      Z[j,] <- correct_Z$correct_label
    }
    # 
    cat("current theta", j, "-th iteration", theta[j,,], "\n")
    # update j
    j = j+1
  }
  samp.end.time <- Sys.time()
  #### DIC calculation
  ####
  DIC.calc.start.time <- Sys.time()
  DIC3 = NA
  DIC4 = NA
  # Thinning the MCMC samples
  id.thin <- seq(from=burn.in+1, to = main.iters+1, by = thin)
  theta.thin <- array(theta[id.thin,,], dim=c(length(id.thin),p,K))
  eta.thin <- matrix(eta[id.thin,],ncol=K)
  Z.thin <- Z[id.thin,]
  # browser()
  # if DIC.calc is true, we are going to calculate the DIC value
  if(DIC.calc){
    if(is.null(DIC.samp.size)){
      DIC.samp.size <- nrow(eta.thin)
    }
    DIC3 = DIC3(mplesetup_list = mplesetup_list,
                eta = eta.thin,
                theta = theta.thin,
                DIC.samp.size = DIC.samp.size, 
                size.offset = size.offset, 
                pl = pl, num.cores = num.cores)  
  }
  DIC.calc.end.time <- Sys.time()
  ############################# posterior class probability
  # class_prob_PL <- function(mplesetup_list, eta, theta, samp.size,...)
  post_class_prob <- NA
  post_class_label <- rep(NA, length = n)
  #browser()
  if(post.class.calc){
    post_class_prob <- class_prob(mplesetup_list = mplesetup_list,
                                  eta = eta.thin,
                                  theta = theta.thin, 
                                  samp.size = DIC.samp.size, 
                                  size.offset = size.offset, pl = pl, num.cores=num.cores)
    post_class_prob <- matrix(post_class_prob, ncol = K)
    for(i in 1:n){
      post_class_label[i] <- which(post_class_prob[i,] == max(post_class_prob[i,]))
    }
  }
  #############################
  # browser()
  out = list(samp.Time = difftime(samp.end.time, samp.start.time),
             init.Time = difftime(samp.start.time, init.start.time),
             DIC.calc.Time = difftime(DIC.calc.end.time,DIC.calc.start.time),
             form = form, theta = theta, eta = eta, Z = Z, 
             acc.counts = acc.counts, DIC3 = DIC3, DIC4 = DIC4,
             burn.in = burn.in, main.iters=main.iters, thin = thin,
             mple.coef.colabel = mple_coef_colabel, id.thin = id.thin, 
             eta.thin = eta.thin, exchange.step = exchange.step, aux.iters=aux.iters,
             theta.thin = theta.thin,
             Z.thin = Z.thin, post.class.label = post_class_label, 
             post.class.prob = post_class_prob, size.offset = size.offset)
  return(out)
}
####################################################################
############# functions for simulating from mixture of ERGMs
############# input : network size, size of clusters, parameters
############# covariates_list, length : p indicating the number of nodal covariates
############# output : a list of networks
#####################################################################
sim_mixture_ergm <- function(form, network_size, cluster_size, 
                             params, covariates_list = NULL,
                             seed, directed=FALSE, tprob=0.1,...){
  # initialize a network
  y <- network::network.initialize(n = network_size, directed = directed)
  if(!is.null(covariates_list)){
    p = length(covariates_list)
    for(j in 1:p){
      network::set.vertex.attribute(y,attrname = covariates_list[[j]]$attrname,
                                    value = covariates_list[[j]]$value)
      
    }
  }
  form <- statnet.common::nonsimp_update.formula(form, y~., from.new = TRUE)
  # check if nrow(params) matches the length of cluster size vector
  if(nrow(params) != length(cluster_size)){
    stop("The number of clusters does not match the number of vectors of parameters!")
  }
  # 
  sim_networks_sub <- vector(mode="list", length = length(cluster_size))
  for(i in 1:length(cluster_size)){
    set.seed(seed*i)
    temp_graph <- sna::rgraph(n = network_size, tprob = tprob, mode = "graph")
    y[,] <- temp_graph[,]
    
    sim_networks_sub[[i]] <- simulate_formula(form,
                                              nsim=cluster_size[i],
                                              seed = seed*i,
                                              coef=params[i,],
                                              control=control.simulate.formula(MCMC.burnin = 20*network_size^2,
                                                                               MCMC.interval = 4*network_size^2))
    
  }
  sim_networks <- do.call("c", sim_networks_sub)
  out <- list(form = form, network_size = network_size, cluster_size = cluster_size,
              params = params, seed = seed, directed = directed, sim_networks = sim_networks)
  return(out)
}
####################
triangles.proportion <- function(dat, 
                                 x,...){
  # triangles_proportion <- rep(NA, 3)
  # figure out n0 and n1
  n <- table(dat%v%x)
  # rearrange the original adjacency matrix
  ordered_dat <- dat[order(dat%v%x),
                     order(dat%v%x)]
  # 
  D <- dat
  D[,] <- 0
  D[1:n[1],1:n[1]] <- ordered_dat[1:n[1], 1:n[1]]
  # 
  R <- dat
  R[,] <- 0
  R[(n[1]+1):(n[1]+n[2]), (n[1]+1):(n[1]+n[2])] <- ordered_dat[(n[1]+1):(n[1]+n[2]), (n[1]+1):(n[1]+n[2])]
  # 
  D_triangles <- sna::triad.census(D[,], mode="graph")[4]
  R_triangles <- sna::triad.census(R[,], mode="graph")[4]
  total_triangles <- sna::triad.census(dat[,], mode="graph")[4]
  DR_triangles <- total_triangles-D_triangles-R_triangles
  #
  return(data.frame(D = D_triangles/choose(n[1],3), 
                    R = R_triangles/choose(n[2],3), 
                    DR = DR_triangles/(choose(n[1],2)*choose(n[2],1) + choose(n[2],2)*choose(n[1],1)) ))
}
###############################################
# post check mixture ergm
# given a list of networks, return the summary statistics
post_check_mixture_ergm <- function(dat, ...){
  # figure out the network size 
  n <- length(dat)
  # netsize <- network::network.size(dat[[1]])
  # netsize is a vector giving the size of each network in the dat
  netsize <- do.call(c,lapply(dat, function(x) network::network.size(x) ) )
  
  # and check if network is directed
  directed.flag <- network::is.directed(dat[[1]])
  
  # initialize a list for results
  # post_check_rlt_list <- list()
  # return the summary stats
  
  # observed data
  obs_networks <- dat
  # post_check_obs_rlt_list[[i]]$
  # 
  if(directed.flag){
    tot_triangles <- factorial(netsize)/factorial(netsize-3)
  } else{
    tot_triangles <- choose(netsize,3)
  }
  # 
  obs_triangles_freq <- sna::triad.census(obs_networks,
                                          mode=ifelse(directed.flag,"digraph","graph"))[,4]/tot_triangles
  
  # 
  obs_gtrans <- do.call("c", lapply(obs_networks, function(x) sna::gtrans(x, mode=ifelse(directed.flag,"digraph","graph")) ))
  # 
  obs_assort <- do.call("c", lapply(obs_networks, function(x) igraph::assortativity_degree(intergraph::asIgraph(x),
                                                                                           directed = directed.flag)))
  # 
  obs_mean_ec <- do.call("c", lapply(obs_networks, function(x) igraph::eigen_centrality(intergraph::asIgraph(x),
                                                                                        directed = directed.flag)$vector %>% mean() ) ) 
  #
  obs_deg_sd <- do.call("c", lapply(obs_networks, function(x) sna::degree(x, 
                                                                          gmode = ifelse(directed.flag,"digraph","graph"),
                                                                          cmode = "indegree") %>% sd() )) 
  
  # 
  obs_deg_mean <- do.call("c", lapply(obs_networks, function(x) sna::degree(x, 
                                                                            gmode = ifelse(directed.flag,"digraph","graph"),
                                                                            cmode = "indegree") %>% mean() ))
  
  # 
  obs_mean_i_path_length <- rep(NA, length(obs_networks))
  for(j in 1:length(obs_networks)){
    temp_geodist <- sna::geodist(obs_networks[[j]])$gdist
    obs_mean_i_path_length[j] <- mean((1/temp_geodist)[lower.tri(1/temp_geodist)])
  }
  # global efficiency
  obs_ge <- do.call("c", lapply(obs_networks, function(x) brainGraph::efficiency(intergraph::asIgraph(x),
                                                                                 type = "global")))
  #
  obs_le <- do.call("c", lapply(obs_networks, function(x) brainGraph::efficiency(intergraph::asIgraph(x),
                                                                                 type = "local") %>% mean() ))
  # post_check_rlt_list$triangles_freq <- obs_triangles_freq
  # post_check_rlt_list$gtrans <- obs_gtrans
  # post_check_rlt_list$assort <- obs_assort
  # post_check_rlt_list$mean_ec <- obs_mean_ec
  # post_check_rlt_list$deg_sd <- obs_deg_sd
  # post_check_rlt_list$mean_i_path_length <- obs_mean_i_path_length
  
  return(data.frame(rbind(data.frame(value = obs_triangles_freq, metric = "Triangles Frequency"),
                          data.frame(value = obs_gtrans, metric = "Transitivity"),
                          data.frame(value = obs_assort, metric = "Assortativity"),
                          data.frame(value = obs_mean_ec, metric = "Mean Eigencentrality"),
                          data.frame(value = obs_deg_sd, metric = "Degree Standard Deviation"),
                          data.frame(value = obs_deg_mean, metric = "Degree Mean"),
                          data.frame(value = obs_ge, metric = "Global Efficiency"),
                          data.frame(value = obs_le, metric = "Local Efficiency"),
                          data.frame(value = obs_mean_i_path_length, metric = "Average Inverse Path Length")) ) )
  
  # 
  # return(list(triangles_freq = obs_triangles_freq,
  #             gtrans = obs_gtrans,
  #             assort = obs_assort,
  #             mean_ec = obs_mean_ec,
  #             deg_sd = obs_deg_sd,
  #             mean_i_path_length = obs_mean_i_path_length ) )
}
###############################################
#### bayesian mixture of ERGMs using pseudolikelihood
# use pseudolikelihood when updating membership indicators and calculate DICs
# use exchange step when updating cluster-specific parameters
# use overclustering to determine the number of clusters
# K-centroid clustering
# https://rdrr.io/cran/flexclust/man/kcca.html
####
#### input: 
### dat : a list of network data, each element is a network
### form: formula
### p: model size (number of parameters)
### K: number of clusters to fit
### prior.mean: mean of the Gaussian prior 
### prior.sigma: variance of the Gaussian prior
### alpha: concentration parameter of the dirichlet distribution
### a: alpha ~ Gamma(a, a*K)
### exchange.step = perform exchange algorithm when updating ERGM parameters?
### pseudo.step = use pseudolikelihood when updating cluster memberships
### aux.iters = number of iterations used for simulating from ERGMs
### num.cores: number of cores when simulating networks from ERGMs
### random.init: do you want to randomly assign parameter values to ERGM parameters?
### init.method: method used for initialization
### size.offset: do you want to introduce Krivitsky offset?
###########################################################################
bayes_mixture_ergm_overclustering <- function(dat, form, p, K, 
                                              prior.mean=c(-1,rep(0,p-1)), 
                                              prior.sigma=diag(25,p), a=10,
                                              main.iters=1000, burn.in=(main.iters)/2,
                                              thin=1, 
                                              Z.step = c("PL", "APL", "full"), 
                                              theta.step = c("PL", "APL", "exchange"),
                                              aux.iters = 1000,
                                              sigma.epsilon=diag(0.001,p), 
                                              alpha.epsilon=0.01^2,
                                              pre.assign=NULL,
                                              num.cores=1,
                                              DIC.calc=TRUE,
                                              BIC.calc=TRUE,
                                              DIC.samp.size=500,
                                              post.class.calc=TRUE,
                                              post.check = TRUE,
                                              random.init=TRUE, init.method = "mple", 
                                              order.constraints = TRUE,
                                              edge.init.value=-2, size.offset = TRUE,...){
  init.start.time <- Sys.time()
  # sample size
  n <- length(dat)
  # netsize <- network::network.size(dat[[1]])
  # netsize is a vector giving the size of each network in the dat
  netsize <- do.call(c,lapply(dat, function(x) network::network.size(x) ) )
  # sufficient stats for each observed data
  sy <- matrix(NA, nrow = n, ncol = p)
  mplesetup_list <- vector(mode="list", length = n)
  data.glm.initial_list <- vector(mode="list", length=n)
  # mstats_list <- vector(mode="list", length=n)
  # get the information needed for calculating the pseudolikelihood
  for(i in 1:n){
    # sufficient statistics
    temp.dat <- dat[[i]]
    temp.form <- statnet.common::nonsimp_update.formula(form, temp.dat~., from.new = TRUE) # LHS of formula is temp.dat
    sy[i,] <- summary_formula(temp.form)
    mplesetup <- ergmMPLE(temp.form)
    mplesetup$size <- netsize[i]
    mplesetup$form <- temp.form
    mplesetup$dat <- dat[[i]]
    # MPLE calculation
    temp.mple.fit <- ergm(temp.form, estimate = "MPLE")
    mplesetup$theta_MPLE <- temp.mple.fit$coef
    # 
    mplesetup_list[[i]] <- mplesetup
    data.glm.initial_list[[i]] <- cbind(mplesetup$response, mplesetup$weights, 
                                        mplesetup$predictor)
    colnames(data.glm.initial_list[[i]]) <- c("responses", "weights", colnames(mplesetup$predictor))
  }
  # we have got sy as the 
  # Initialize the matrix for storing sampled parameters
  alpha <- rep(NA, 1+main.iters)
  tau <- matrix(NA, nrow = 1+main.iters, ncol = K)
  theta <- array(NA, dim = c(1+main.iters, p, K))
  Z <- matrix(NA, nrow = 1+main.iters, ncol = n)
  acc.counts <- matrix(FALSE, nrow=1+main.iters, ncol = K)
  K_nonzero <- rep(NA, 1+main.iters) # keep track of the number of non-zero components
  #browser()
  # mple_coef_colabel <- NA
  
  # randomly initialize the parameters
  # according to the prior
  alpha[1] <- rgamma(n=1, shape = a, rate = a*K)
  Z[1,] <- sample(x=1:K,size=n,replace=TRUE,prob=rep(alpha[1]/(alpha[1] * K), K))
  tau[1,] <- rdirichlet(n=1, alpha = rep(alpha[1], K)) 
  # rep(alpha[1]/(alpha[1] * K), K)
  for(k in 1:K){
    theta[1,,k] <- runif(n=p, min=-0.1,max=0.1)
    # apply(filter(mple_coef_colabel, olabel == k),2,mean)[1:p]
  }
  # theta[1,1,] <- theta[1,1,]+edge.init.value
  K_nonzero[1] <- length(table(Z[1,]))
  ######
  # if pre.assign not null, use pre-assignments
  if(!is.null(pre.assign)){
    Z[1,] <- pre.assign
  }
  
  # let edge parameter start at somewhere that is negative and order them to make sure the first cluster has the
  # smallest initial value on edge parameter
  ########
  # theta_order <- order(theta[1,1,], decreasing = FALSE)
  # theta[1,,] <- theta[1,,theta_order]
  #########
  # metropolis within gibbs sampler
  # online algorithm to deal with label switching
  # start updating the algorithm from j = 2
  j = 2
  # main iterations
  samp.start.time <- Sys.time()
  # prepare for the exchange algorithm
  #########################
  if(match.arg(theta.step) == "exchange"){
    # setup the list for simulating from ergms
    # y <- ergm.getnetwork(mod2$formula)
    control <- control.ergm(MCMC.burnin = aux.iters, 
                            MCMC.interval = 1, 
                            MCMC.samplesize = 1)
    # model_list <- vector(mode="list", length=n)
    Clist_list <- vector(mode = "list",length=n)
    proposal_list <- vector(mode = "list", length=n)
    # sy is a matrix
    for(i in 1:n){
      temp_model <- ergm::ergm_model(form, dat[[i]])
      Clist_list[[i]] <- ergm::ergm.Cprepare(dat[[i]], temp_model)
      proposal_list[[i]] <- ergm::ergm_proposal(object = ~., constraints = ~., 
                                                arguments = control$MCMC.prop.args,
                                                nw = dat[[i]])
    }
  }
  
  # prepare for the adjusted pseudolikelihood
  APL.calc.start.time <- Sys.time()
  # APL.calc.Time = difftime(APL.calc.end.time, APL.calc.start.time),
  W_list <- vector(mode="list", length=n)
  logC_list <- vector(mode="list", length=n)
  Theta_MLE_list <- vector(mode="list", length=n)
  Theta_MPLE_list <- vector(mode="list", length=n)
  #########################
  if(match.arg(Z.step) == "APL" | match.arg(theta.step) == "APL"){
    # use Bergm::ergmAPL to get the parameters
    for(i in 1:n){
      # 
      temp.dat <- dat[[i]]
      temp.form <- statnet.common::nonsimp_update.formula(form, temp.dat~., from.new = TRUE) # LHS of formula is temp.dat
      
      ### try MCMLE first
      temp_ergmAPL_rlt <- try(Bergm::ergmAPL(formula = temp.form,
                                             ladder = 500,
                                             aux.iters = 2*netsize[i]^2,
                                             n.aux.draws = 1500,
                                             aux.thin = 50, estimate = "MLE",
                                             control = control.ergm(main.method = c("MCMLE"),
                                                                    MCMC.samplesize = 1024) ))
      
      ### then CD 
      temp_ergmAPL_rlt <- try(Bergm::ergmAPL(formula = temp.form,
                                             ladder = 500,
                                             aux.iters = 2*netsize[i]^2,
                                             n.aux.draws = 1500,
                                             aux.thin = 50, estimate = "CD" ))
      
      ### then Robbins-Monro
      if(class(temp_ergmAPL_rlt) == "try-error"){
        temp_ergmAPL_rlt <- try(Bergm::ergmAPL(formula = temp.form,
                                               ladder = 500,
                                               aux.iters = 2*netsize[i]^2,
                                               n.aux.draws = 1500,
                                               aux.thin = 50, estimate = "MLE",
                                               control = control.ergm(main.method = c("Robbins-Monro"),
                                                                      MCMC.samplesize = 1024) ))
      }
      
      ### then stochastic approximation
      if(class(temp_ergmAPL_rlt) == "try-error"){
        temp_ergmAPL_rlt <- try(Bergm::ergmAPL(formula = temp.form,
                                               ladder = 500,
                                               aux.iters = 2*netsize[i]^2,
                                               n.aux.draws = 1500,
                                               aux.thin = 50, estimate = "MLE",
                                               control = control.ergm(main.method = c("Stochastic-Approximation"),
                                                                      MCMC.samplesize = 1024) ))
      }
      
      ### then stepping 
      if(class(temp_ergmAPL_rlt) == "try-error"){
        temp_ergmAPL_rlt <- try(Bergm::ergmAPL(formula = temp.form,
                                               ladder = 500,
                                               aux.iters = 2*netsize[i]^2,
                                               n.aux.draws = 1500,
                                               aux.thin = 50, estimate = "MLE",
                                               control = control.ergm(main.method = c("Stepping"),
                                                                      MCMC.samplesize = 1024) ))
      }
      
      # 4*size^2
      W_list[[i]] <- temp_ergmAPL_rlt$W
      logC_list[[i]] <- temp_ergmAPL_rlt$logC
      Theta_MLE_list[[i]] <- temp_ergmAPL_rlt$Theta_MLE
      Theta_MPLE_list[[i]] <- temp_ergmAPL_rlt$Theta_PL
    }
  }
  ####################
  APL.calc.end.time <- Sys.time()
  ####################
  # main MCMC sampling 
  while(j <= 1+main.iters){
    # -- update tau --#
    # update tau (from dirichlet), now we have Z[j,]
    tau[j,] <- gtools::rdirichlet(n=1, 
                                  alpha = (rep(alpha[j-1],K)+table(factor(Z[j-1,], levels = 1:K))) ) # K-dimensional vector
    cat("tau", tau[j,], "\n")
    # -- tau updated --# 
    
    # update theta, using random walk metropolis update or exchange algorithm
    # update theta_k in random order
    nK <- sample(1:K, K)
    for(k in nK){
      # proposing a new theta for k-th component
      theta1 <- rmvnorm(1,mean=theta[j-1,,k], sigma = sigma.epsilon)
      # log prior vector
      pr <- dmvnorm( rbind(theta1, theta[j-1,,k]), mean= prior.mean, sigma = prior.sigma, log=TRUE)
      # adjusted pseudolikelihood vector
      # lr <- c( logPL.corr(theta = theta[j-1,,k]) )
      k_id <- which(Z[j-1,] == k) # index of observed data whose corresponding assigned cluster is k 
      # 
      cat("The",k,"-th cluster:",k_id, "\n")
      if(length(k_id) >= 1){
        ll_mat <- matrix(NA, nrow=2,ncol = length(k_id))
        for(l in 1:length(k_id)){
          # network index k_id[l]
          if(match.arg(theta.step)=="PL"){
            ll_mat[2,l] <- logPL(theta = theta[j-1,,k], y=mplesetup_list[[k_id[l]]]$response,
                                 X = mplesetup_list[[k_id[l]]]$predictor, weights = mplesetup_list[[k_id[l]]]$weights,
                                 size = mplesetup_list[[k_id[l]]]$size, size.offset = size.offset)
            ll_mat[1,l] <- logPL(theta = c(theta1), y=mplesetup_list[[k_id[l]]]$response,
                                 X = mplesetup_list[[k_id[l]]]$predictor, weights = mplesetup_list[[k_id[l]]]$weights,
                                 size = mplesetup_list[[k_id[l]]]$size, size.offset = size.offset)
          }
          else if(match.arg(theta.step)=="exchange"){
            # exchange step for updating cluster-specific parameters, the proposed theta for cluster k is theta1
            # simulation step, simulate network from theta1
            ## figure out network size
            temp_coef <- theta1
            if(size.offset){
              temp_size <- network::network.size(dat[[ k_id[l] ]])
              temp_coef[1] <- theta1[1] - log(temp_size)
            }
            
            temp_stats <- ergm_MCMC_slave(Clist = Clist_list[[ k_id[l] ]], 
                                          proposal = proposal_list[[ k_id[l] ]], 
                                          eta = temp_coef, 
                                          control = control, 
                                          verbose = FALSE)$s + sy[k_id[l],] # set up them first
            # obs stats is named sy, which is of size n by p
            # likelihood ratio calculation step
            ll_mat[1,l] <- sum(theta[j-1,,k] * temp_stats) + sum( theta1 * sy[k_id[l],] )
            ll_mat[2,l] <- sum(theta[j-1,,k] * sy[k_id[l],]) + sum( theta1 * temp_stats )
          } else if(match.arg(theta.step) == "APL"){
            # adjusted pseudo-likelihood
            ll_mat[2,l] <- logAPL(theta = theta[j-1,,k], 
                                  y=mplesetup_list[[k_id[l]]]$response,
                                  X = mplesetup_list[[k_id[l]]]$predictor, 
                                  weights = mplesetup_list[[k_id[l]]]$weights,
                                  size = mplesetup_list[[k_id[l]]]$size, 
                                  size.offset = size.offset,
                                  theta_MLE = Theta_MLE_list[[k_id[l]]],
                                  theta_MPLE = Theta_MPLE_list[[k_id[l]]],
                                  W = W_list[[k_id[l]]],
                                  logC = logC_list[[k_id[l]]])
            #####
            ll_mat[1,l] <- logAPL(theta = c(theta1), y=mplesetup_list[[k_id[l]]]$response,
                                  X = mplesetup_list[[k_id[l]]]$predictor, 
                                  weights = mplesetup_list[[k_id[l]]]$weights,
                                  size = mplesetup_list[[k_id[l]]]$size, 
                                  size.offset = size.offset,
                                  theta_MLE = Theta_MLE_list[[k_id[l]]],
                                  theta_MPLE = Theta_MPLE_list[[k_id[l]]],
                                  W = W_list[[k_id[l]]],
                                  logC = logC_list[[k_id[l]]])
          }
        }
        lr <- apply(ll_mat,1,sum) # row-sum
      }
      else {
        lr <- c(0,0)
      }
      # log acceptance ratio, proposal is symmetric
      beta.theta <- (lr[1] - lr[2]) + (pr[1] - pr[2]) 
      # accept if greater than acceptance ratio 
      if(beta.theta >= log(runif(1))){
        theta[j,,k] <- theta1
        acc.counts[j,k] <- TRUE
      }
      else {
        theta[j,,k] <- theta[j-1,,k]
      }
    }
    ### theta updated
    
    # update the assignment label (total of n labels)
    for(i in 1:n){
      unnorm_p <- rep(NA, length = K)
      unnorm_logp <- rep(NA, length = K)
      if(match.arg(Z.step)=="PL"){
        # pseudo.step = TRUE indicates that we want to use pseudo likelihood to update
        # the cluster membership
        for(k in 1:K){
          # ignore the constant C
          cat("Previous tau", tau[j,k], "\n")
          cat("Previous theta", theta[j,,k], "\n")
          # calculating log pseudolikelihood
          temp_logPL <- logPL(theta = c(theta[j,,k]), y=mplesetup_list[[i]]$response,
                              X = mplesetup_list[[i]]$predictor, 
                              weights = mplesetup_list[[i]]$weights,
                              size = mplesetup_list[[i]]$size, size.offset = size.offset)
          cat("Corrected pseudo-like", temp_logPL, "\n")
          unnorm_logp[k] <- log(tau[j,k]) + temp_logPL
          cat("unnormalized log prob", unnorm_logp[k], "\n")
        }
        # A = max(unnorm_logp)
        cat("The unnormalized log-probability", unnorm_logp, "\n")
        # -Inf -Inf -884.7563 -1165.712 -422.0812 -517.0914 
        # lead to NaN
        # 
        lsum = sna::logSum(unnorm_logp[which(unnorm_logp!=-Inf)])
        unnorm_p <- exp(unnorm_logp-lsum)
        cat("The unnormalized probability", unnorm_p, "\n")
        cat("lsum=", lsum, "\n")
        print(paste("j=",j, "i=",i, "k=",k))
        cat("\n")
        # Update Z...
        Z[j,i] <- sample(x=1:K, size=1, prob = unnorm_p)
      } else if(match.arg(Z.step) == "full"){
        # if using full likelihood
        for(k in 1:K){
          cat("Previous tau", tau[j,k], "\n")
          cat("Previous theta", theta[j,,k], "\n")
          # calculating log pseudolikelihood
          temp_logL <- logL(theta = c(theta[j,,k]),
                            mplesetup = mplesetup_list[[i]],
                            size.offset = size.offset)
          # 
          cat("full likelihood", temp_logL, "\n")
          unnorm_logp[k] <- log(tau[j,k]) + temp_logL
          cat("unnormalized log prob", unnorm_logp[k], "\n")
        }
        # A = max(unnorm_logp)
        # lsum = sna::logSum(unnorm_logp)
        lsum = sna::logSum(unnorm_logp[which(unnorm_logp!=-Inf)])
        unnorm_p <- exp(unnorm_logp-lsum)
        cat("The unnormalized probability", unnorm_p, "\n")
        # Update Z...
        Z[j,i] <- sample(x=1:K, size=1, prob = unnorm_p)
      } else if(match.arg(Z.step) == "APL"){
        # adjusted pseudo likelihood
        for(k in 1:K){
          cat("Previous tau", tau[j,k], "\n")
          cat("Previous theta", theta[j,,k], "\n")
          # calculating log pseudolikelihood
          temp_logAPL <- logAPL(theta = c(theta[j,,k]),
                                y=mplesetup_list[[i]]$response,
                                X = mplesetup_list[[i]]$predictor,
                                weights = mplesetup_list[[i]]$weights,
                                size = mplesetup_list[[i]]$size,
                                size.offset = size.offset,
                                theta_MLE = Theta_MLE_list[[i]],
                                theta_MPLE = Theta_MPLE_list[[i]],
                                W = W_list[[i]],
                                logC = logC_list[[i]])
          
          # logPL(theta = c(theta[j-1,,k]), y=mplesetup_list[[i]]$response,
          #                   X = mplesetup_list[[i]]$predictor, 
          #                   weights = mplesetup_list[[i]]$weights,
          #                   size = mplesetup_list[[i]]$size, size.offset = size.offset)
          cat("Corrected pseudo-like", temp_logAPL, "\n")
          unnorm_logp[k] <- log(tau[j,k]) + temp_logAPL
          cat("unnormalized log prob", unnorm_logp[k], "\n")
        }
        # A = max(unnorm_logp)
        cat("The unnormalized log-probability", unnorm_logp, "\n")
        # -Inf -Inf -884.7563 -1165.712 -422.0812 -517.0914 
        # lead to NaN
        # 
        lsum = sna::logSum(unnorm_logp[which(unnorm_logp!=-Inf)])
        unnorm_p <- exp(unnorm_logp-lsum)
        cat("The unnormalized probability", unnorm_p, "\n")
        cat("lsum=", lsum, "\n")
        print(paste("j=",j, "i=",i, "k=",k))
        cat("\n")
        # Update Z...
        Z[j,i] <- sample(x=1:K, size=1, prob = unnorm_p)
      }
    }
    
    #-- update alpha using random-walk metropolis --#
    # propose a new value
    alpha1 <- exp(rnorm(1, mean=log(alpha[j-1]), sd = sqrt(alpha.epsilon)))
    # calculate the logarithm of acceptance probability
    # browser()
    beta.alpha <- (sum( (alpha1-1)*log(tau[j,]) ) + log(gamma(K * alpha1)) + 
                     log(dgamma(alpha1,shape = a,rate=a*K)) - K*log(gamma(alpha1))) - (sum((alpha[j-1])*log(tau[j,]))) + log(gamma(K * alpha[j-1])) + log(dgamma(alpha[j-1], shape=a, rate=a*K)) - K*log(gamma(alpha[j-1])) + log(1/alpha[j-1]) - log(1/alpha1)
    # accept the proposal based on the value of beta.alpha
    if(beta.alpha >= log(runif(1))){
      # accept the proposal
      alpha[j] <- alpha1
      # acc.counts[j,k] <- TRUE
    }
    else{
      # reject the proposal
      # theta[j,,k] <- theta[j-1,,k]
      alpha[j] <- alpha[j-1]
    }
    # -- alpha updated --#
    
    # re-order parameters to deal with label switching issue, from small to large, using the order constraints on 
    # edge parameter
    if(order.constraints){
      theta_order <- order(theta[j,1,], decreasing = FALSE)
      theta[j,,] <- theta[j,,theta_order]
      tau[j,] <- tau[j,theta_order]
      # also need to permute the latent cluster assignment... 
      temp_Z <- data.frame(old_label = Z[j,])
      reference_Z <- data.frame(old_label=1:K, correct_label = theta_order)
      correct_Z <- dplyr::left_join(temp_Z, reference_Z, by = c("old_label" = "old_label"))
      Z[j,] <- correct_Z$correct_label
    }
    else{
      # randomly permute their labels
      theta_order <- sample(1:K,size=K,replace = FALSE)
      theta[j,,theta_order] <- theta[j,,]
      tau[j,theta_order] <- tau[j,]
      # also need to permute the latent cluster assignment... 
      temp_Z <- data.frame(old_label = Z[j,])
      reference_Z <- data.frame(old_label=1:K, correct_label = theta_order)
      correct_Z <- dplyr::left_join(temp_Z, reference_Z, by = c("old_label" = "old_label"))
      Z[j,] <- correct_Z$correct_label
    }
    # 
    cat("current theta", j, "-th iteration", theta[j,,], "\n")
    # update K_nonzero: number of components with at least one observation
    K_nonzero[j] <- length(table(Z[j,])) # only clusters with observations are counted
    # update j
    j = j+1
  }
  samp.end.time <- Sys.time()
  #### figure out the number of clusters using posterior mode of K
  # Khat_table = table(factor(K_nonzero[(burn.in+1):(main.iters+1)], 
  #                           levels=1:K))
  # thinning the MCMC draws 
  # browser()
  K_nonzero.thin = K_nonzero[seq(burn.in+1, main.iters+1, by=thin)]
  tau.thin = tau[seq(burn.in+1, main.iters+1,by=thin),]
  Z.thin = Z[seq(burn.in+1, main.iters+1,by=thin),]
  theta.thin = theta[seq(burn.in+1, main.iters+1,by=thin),,]
  alpha.thin = alpha[seq(burn.in+1, main.iters+1, by=thin)]
  # 
  Khat_table = table(factor(K_nonzero.thin, 
                            levels=1:K))
  #### 
  Khat = as.numeric(which(Khat_table == max(Khat_table))[1])
  ####
  #### K-centroid clustering 
  # first identify those iterations with the number of non-zero components = Khat
  sub_id <- which(K_nonzero.thin == Khat)
  # keep those after burn_in
  # sub_id <- sub_id[sub_id >= (burn.in+1)]
  M0 <- length(sub_id)
  # pull out the theta.thin's to form a matrix : p columns, Khat * M0 rows
  theta_sub <- matrix(NA, nrow = M0*Khat, ncol = p)
  tau_sub <- rep(NA, length=M0*Khat)
  alpha_sub <- rep(NA, length=M0)
  #browser()
  for(jj in 1:M0){
    theta_sub[(1+(jj-1)*Khat):(jj*Khat),] <- t(theta.thin[sub_id[jj],,unique(Z.thin[sub_id[jj],])])
    tau_sub[(1+(jj-1)*Khat):(jj*Khat)] <- tau.thin[sub_id[jj],unique(Z.thin[sub_id[jj],])]
    alpha_sub[jj] <- alpha.thin[jj]
  }
  # k-centroid clustering
  # kcca_rlt <- kcca(Nclus, k=Khat, family=kccaFamily("kmedians"),
  #             control=list(initcent="kmeanspp"))
  #
  kcca_rlt <- Cluster_Medoids(theta_sub, 
                              clusters = Khat, 
                              distance_metric = "mahalanobis")
  # 
  # kcca_rlt$clusters is a sequence of cluster membership indicators
  # check if the clustering result is valid for each iteration of parameters
  theta_sub_sub <-  array(NA, dim = c(M0, p, Khat))
  tau_sub_sub <- matrix(NA, nrow = M0, ncol = Khat)
  alpha_sub_sub <- rep(NA, length = M0)
  M0_sub = 0 # number of iterations retained 
  #
  for(jj in 1:M0){
    # if the cluster assignment is a permutation
    if(all(sort(kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)]) == c(1:Khat)) ){
      theta_sub_sub[M0_sub+1,,kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)]] <- t(theta_sub[(1+(jj-1)*Khat):(jj*Khat),])
      alpha_sub_sub[M0_sub+1] <- alpha_sub[jj]
      tau_sub_sub[M0_sub+1,kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)]] <- tau_sub[(1+(jj-1)*Khat):(jj*Khat)]
      # tau_sub_sub[jj,kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)]] <- tau.thin[jj, ] 
      # theta[sub_id[jj], ,unique(Z[sub_id[jj],])]
      M0_sub = M0_sub + 1
    }
  }
  # only keep the complete.cases
  theta_sub_sub = theta_sub_sub[1:M0_sub,,]
  alpha_sub_sub = alpha_sub_sub[1:M0_sub]
  tau_sub_sub = tau_sub_sub[1:M0_sub,]
  # 
  tau_sub_sub <- sweep(tau_sub_sub, 1, apply(tau_sub_sub,1,sum), "/")
  # 
  #### the cluster with smaller value on edge parameter should have smaller cluster id number
  ####
  cluster_rank <- rank(kcca_rlt$medoids[,1])
  ####
  for(jj in 1:nrow(tau_sub_sub)){
    tau_sub_sub[jj,cluster_rank] <- tau_sub_sub[jj,]
    theta_sub_sub[jj,,cluster_rank] <- theta_sub_sub[jj,,]
  }
  # 
  # browser()
  #### DIC calculation
  ####
  DIC.calc.start.time <- Sys.time()
  DIC3 = NA
  DIC4 = NA
  # Thinning the MCMC samples
  # id.thin <- seq(from=burn.in+1, to = main.iters+1, by = thin)
  # theta.thin <- array(theta[id.thin,,], dim=c(length(id.thin),p,K))
  # tau.thin <- matrix(tau[id.thin,],ncol=K)
  # Z.thin <- Z[id.thin,]
  # browser()
  # if DIC.calc is true, we are going to calculate the DIC value
  # DIC calculation does not require Z here
  if(DIC.calc){
    # if(is.null(DIC.samp.size)){
    #   DIC.samp.size <- nrow(tau_sub_sub)
    # }
    # DIC3 = DIC3(mplesetup_list = mplesetup_list,
    #             eta = tau_sub_sub,
    #             theta = theta_sub_sub,
    #             DIC.samp.size = DIC.samp.size, 
    #             size.offset = size.offset)  
    
    DIC3 = DIC_mixture_ergm(mplesetup_list = mplesetup_list,
                            tau = tau_sub_sub,
                            theta = theta_sub_sub,
                            DIC.samp.size = nrow(tau_sub_sub),
                            size.offset = size.offset,
                            lik = Z.step,
                            theta_MLE_list = Theta_MLE_list,
                            theta_MPLE_list = Theta_MPLE_list,
                            W_list = W_list,
                            logC_list = logC_list)
    
    # DIC3 <- function(mplesetup_list, 
    #                  tau, 
    #                  theta, 
    #                  DIC.samp.size=500, 
    #                  size.offset=TRUE, 
    #                  lik = c("PL", "APL", "full"),
    #                  theta_MLE_list=NULL,
    #                  theta_MPLE_list=NULL,
    #                  W_list=NULL,
    #                  logC_list=NULL,...)
  }
  DIC.calc.end.time <- Sys.time()
  ############################# BIC calculation
  #############################
  BIC.calc.start.time <- Sys.time()
  #############################
  BIC1 <- NA
  if(BIC.calc){
    BIC1 <- BIC_mixture_ergm(mplesetup_list = mplesetup_list,
                             p=p,
                             tau = tau_sub_sub,
                             theta = theta_sub_sub,
                             size.offset = size.offset,
                             lik = Z.step,
                             theta_MLE_list = Theta_MLE_list,
                             theta_MPLE_list = Theta_MPLE_list,
                             W_list = W_list,
                             logC_list = logC_list)
    # BIC_mixture_ergm <- function(mplesetup_list,
    #                              p,
    #                              tau, 
    #                              theta, 
    #                              size.offset=TRUE, 
    #                              lik = c("PL", "APL", "full"),
    #                              theta_MLE_list=NULL,
    #                              theta_MPLE_list=NULL,
    #                              W_list=NULL,
    #                              logC_list=NULL,...)
  }
  #############################
  BIC.calc.end.time <- Sys.time()
  ############################# posterior class probability
  # class_prob_PL <- function(mplesetup_list, eta, theta, samp.size,...)
  post_class_prob <- NA
  post_class_label <- rep(NA, length = n)
  #browser()
  if(post.class.calc){
    post_class_prob <- class_prob(mplesetup_list = mplesetup_list,
                                  tau = tau_sub_sub,
                                  theta = theta_sub_sub,
                                  samp.size = nrow(tau_sub_sub),
                                  size.offset = size.offset,
                                  lik = Z.step,
                                  num.cores = num.cores,
                                  theta_MLE_list = Theta_MLE_list,
                                  theta_MPLE_list = Theta_MPLE_list,
                                  W_list = W_list,
                                  logC_list = logC_list)
    ###########################
    post_class_prob <- matrix(post_class_prob, ncol = Khat)
    for(i in 1:n){
      post_class_label[i] <- which(post_class_prob[i,] == max(post_class_prob[i,]))
    }
  }
  #############################
  #############################
  # posterior predictive checks
  # summary stats calculated on the observed data
  # summary stats calculated on simulated data
  obs_post_check_stats <- NA
  sim_post_check_stats <- NA
  if(post.check){
    # observed data summary stats
    obs_post_check_stats <- post_check_mixture_ergm(dat)
    
    # 
    sim_dat <- vector(mode="list", length = n)
    # 
    for(i in 1:length(sim_dat)){
      temp_id <- i
      y <- dat[[i]]
      # y <- temp_sim_mixture_ergm$sim_networks[[temp_id]]
      temp_membership <- as.numeric(factor(y%v%"x") )
      temp_Z <- sample(1:Khat, size=1, prob = tau_sub_sub[(i-1)%%M0_sub+1,]) 
      # which(rmultinom(n=1, size = 1, prob = tau_sub_sub[i,])==1)
      temp_network_size <- network::network.size(y)
      temp_coef <- theta_sub_sub[(i-1)%%M0_sub+1,,temp_Z]
      # rlt_list[[ii]][[ Khat[ii] ]]$theta.thin[i,,temp_Z]
      # temp_coef[1] <- temp_coef[1] - log(temp_network_size)
      temp_sim_network <- simulate(form,
                                   nsim=1,
                                   seed = 1234*i,
                                   coef = temp_coef,
                                   basis = y,
                                   control=control.simulate.formula(MCMC.burnin = 50*temp_network_size^2,
                                                                    MCMC.interval = 1))
      # sim_modularity_1[i,j] <- modularity(intergraph::asIgraph(temp_sim_network),
      #                                      membership = factor(y%v%"party") )
      sim_dat[[i]] <- temp_sim_network 
    }
    
    # simulated data, stats
    sim_post_check_stats <- post_check_mixture_ergm(sim_dat)
  }
  # browser()
  out = list(samp.Time = difftime(samp.end.time, samp.start.time),
             init.Time = difftime(samp.start.time, init.start.time),
             DIC.calc.Time = difftime(DIC.calc.end.time,DIC.calc.start.time),
             BIC.calc.Time = difftime(BIC.calc.end.time, BIC.calc.start.time),
             APL.calc.Time = difftime(APL.calc.end.time, APL.calc.start.time),
             form = form, 
             acc.counts = acc.counts, DIC3 = DIC3, DIC4 = DIC4, BIC1 = BIC1,
             burn.in = burn.in, main.iters=main.iters, thin = thin,
             tau.out = tau_sub_sub, 
             theta.out = theta_sub_sub, 
             aux.iters=aux.iters,
             theta.thin = theta.thin,
             post.class.label = post_class_label, 
             post.class.prob = post_class_prob, 
             size.offset = size.offset,
             Khat = Khat,
             sy = sy,
             mplesetup.list = mplesetup_list,
             theta.MLE.list = Theta_MLE_list,
             W.list = W_list,
             logC.list = logC_list,
             theta.MPLE.list = Theta_MPLE_list,
             sim.post.check.stats = sim_post_check_stats,
             obs.post.check.stats = obs_post_check_stats,
             M0 = M0,
             M0_sub = M0_sub)
  #              mple.coef.colabel = mple_coef_colabel, 
  # theta = theta, tau = tau, Z = Z
  return(out)
}
#############################
#############################
#############################
bayes_mixture_ergm_overclustering_fixedalpha <- function(dat, form, p, K, 
                                                         prior.mean=c(-1,rep(0,p-1)), 
                                                         prior.sigma=diag(25,p), alpha=0.005,
                                                         main.iters=1000, burn.in=(main.iters)/2,
                                                         thin=1, 
                                                         Z.step = c("PL", "APL", "full"), 
                                                         theta.step = c("PL", "APL", "exchange"),
                                                         aux.iters = 1000,
                                                         sigma.epsilon=diag(0.001,p), 
                                                         pre.assign=NULL,
                                                         num.cores=1,
                                                         DIC.calc=TRUE,
                                                         BIC.calc=TRUE,
                                                         DIC.samp.size=500,
                                                         post.class.calc=TRUE,
                                                         random.init=TRUE,
                                                         post.check=TRUE,
                                                         init.method = "mple", 
                                                         order.constraints = TRUE,
                                                         edge.init.value=-2, 
                                                         size.offset = TRUE,
                                                         APL.precomuted = FALSE,
                                                         W.list=NULL,
                                                         logC.list=NULL,
                                                         Theta.MLE.list=NULL,
                                                         Theta.MPLE.list=NULL,...){
  init.start.time <- Sys.time()
  # sample size
  n <- length(dat)
  # netsize <- network::network.size(dat[[1]])
  # netsize is a vector giving the size of each network in the dat
  netsize <- do.call(c,lapply(dat, function(x) network::network.size(x) ) )
  # sufficient stats for each observed data
  sy <- matrix(NA, nrow = n, ncol = p)
  mplesetup_list <- vector(mode="list", length = n)
  data.glm.initial_list <- vector(mode="list", length=n)
  # mple <- matrix(NA, nrow=n, ncol=p)
  # mstats_list <- vector(mode="list", length=n)
  # get the information needed for calculating the pseudolikelihood
  for(i in 1:n){
    # sufficient statistics
    temp.dat <- dat[[i]]
    temp.form <- statnet.common::nonsimp_update.formula(form, temp.dat~., from.new = TRUE) # LHS of formula is temp.dat
    sy[i,] <- summary_formula(temp.form)
    mplesetup <- ergmMPLE(temp.form)
    mplesetup$size <- netsize[i]
    mplesetup$form <- temp.form
    mplesetup$dat <- dat[[i]]
    # browser()
    # mple calculation
    temp.mple.fit <- ergm(temp.form, estimate = "MPLE")
    mplesetup$theta_MPLE <- temp.mple.fit$coef
    #
    mplesetup_list[[i]] <- mplesetup
    data.glm.initial_list[[i]] <- cbind(mplesetup$response, 
                                        mplesetup$weights, 
                                        mplesetup$predictor)
    colnames(data.glm.initial_list[[i]]) <- c("responses", "weights", colnames(mplesetup$predictor))
  }
  # we have got sy as the 
  # Initialize the matrix for storing sampled parameters
  # alpha <- rep(NA, 1+main.iters)
  tau <- matrix(NA, nrow = 1+main.iters, ncol = K)
  theta <- array(NA, dim = c(1+main.iters, p, K))
  Z <- matrix(NA, nrow = 1+main.iters, ncol = n)
  acc.counts <- matrix(FALSE, nrow=1+main.iters, ncol = K)
  K_nonzero <- rep(NA, 1+main.iters) # keep track of the number of non-zero components
  # browser()
  # mple_coef_colabel <- NA
  # randomly initialize the parameters
  # according to the prior
  # alpha[1] <- rgamma(n=1, shape = a, rate = a*K)
  Z[1,] <- sample(x=1:K,size=n,replace=TRUE,prob=rep(alpha[1]/(alpha[1] * K), K))
  tau[1,] <- rdirichlet(n=1, alpha = rep(alpha[1], K)) 
  # rep(alpha[1]/(alpha[1] * K), K)
  for(k in 1:K){
    theta[1,,k] <- runif(n=p, min=-0.1,max=0.1)
    # apply(filter(mple_coef_colabel, olabel == k),2,mean)[1:p]
  }
  # theta[1,1,] <- theta[1,1,]+edge.init.value
  K_nonzero[1] <- length(table(Z[1,]))
  ######
  
  # if pre.assign not null, use pre-assignments
  if(!is.null(pre.assign)){
    Z[1,] <- pre.assign
  }
  
  # let edge parameter start at somewhere that is negative and order them to make sure the first cluster has the
  # smallest initial value on edge parameter
  #########
  # metropolis within gibbs sampler
  # online algorithm to deal with label switching
  # start updating the algorithm from j = 2
  # main iterations
  samp.start.time <- Sys.time()
  # prepare for the exchange algorithm
  #########################
  if(match.arg(theta.step) == "exchange"){
    # setup the list for simulating from ergms
    # y <- ergm.getnetwork(mod2$formula)
    control <- control.ergm(MCMC.burnin = aux.iters, 
                            MCMC.interval = 1, 
                            MCMC.samplesize = 1)
    # model_list <- vector(mode="list", length=n)
    Clist_list <- vector(mode = "list",length=n)
    proposal_list <- vector(mode = "list", length=n)
    # sy is a matrix
    for(i in 1:n){
      temp_model <- ergm::ergm_model(form, dat[[i]])
      Clist_list[[i]] <- ergm:::ergm.Cprepare(dat[[i]], temp_model)
      proposal_list[[i]] <- ergm:::ergm_proposal(object = ~., constraints = ~., 
                                                 arguments = control$MCMC.prop.args,
                                                 nw = dat[[i]])
    }
  }
  
  # prepare for the adjusted pseudolikelihood
  #########################
  APL.calc.start.time <- Sys.time()
  #
  APL.calc.time.vec <- rep(NA, n)
  if(!APL.precomuted){
    #########################
    W_list <- vector(mode="list", length=n)
    logC_list <- vector(mode="list", length=n)
    Theta_MLE_list <- vector(mode="list", length=n)
    Theta_MPLE_list <- vector(mode="list", length=n)
    #########################
    if(match.arg(Z.step) == "APL" | match.arg(theta.step) == "APL"){
      #
      # temp.APL.calc.time1 <- Sys.time()
      # use Bergm::ergmAPL to get the parameters
      #
      for(i in 1:n){
        cat("\n")
        cat("calculating the adjusted pseudolikelihood info for network", i)
        cat("\n")
        #
        temp.APL.calc.time1 <- Sys.time()
        # 
        temp.dat <- dat[[i]]
        temp.form <- statnet.common::nonsimp_update.formula(form, temp.dat~., from.new = TRUE) # LHS of formula is temp.dat
        
        ### try MCMLE first
        temp_ergmAPL_rlt <- try(Bergm::ergmAPL(formula = temp.form,
                                               ladder = 500,
                                               aux.iters = 2*netsize[i]^2,
                                               n.aux.draws = 1500,
                                               aux.thin = 50, 
                                               estimate = "MLE",
                                               control = control.ergm(main.method = "MCMLE",
                                                                      MCMC.samplesize = 4096,
                                                                      MCMC.burnin = 4*network::network.size(temp.dat)^2,
                                                                      SAN.maxit = 16) ))
        
        ### then CD
        if(class(temp_ergmAPL_rlt) == "try-error"){
          temp_ergmAPL_rlt <- try(Bergm::ergmAPL(formula = temp.form,
                                                 ladder = 500,
                                                 aux.iters = 2*netsize[i]^2,
                                                 n.aux.draws = 1500,
                                                 aux.thin = 50, estimate = "CD",
                                                 control = control.ergm(CD.nsteps.obs = 1024*4,
                                                                        CD.nsteps = 32) ))
        }
        
        ### then Robbins-Monro
        if(class(temp_ergmAPL_rlt) == "try-error"){
          temp_ergmAPL_rlt <- try(Bergm::ergmAPL(formula = temp.form,
                                                 ladder = 500,
                                                 aux.iters = 2*netsize[i]^2,
                                                 n.aux.draws = 1500,
                                                 aux.thin = 50, estimate = "MLE",
                                                 control = control.ergm(main.method = c("Robbins-Monro"),
                                                                        MCMC.samplesize = 1024*4,
                                                                        MCMC.burnin = 4*network::network.size(temp.dat)^2) ))
        }
        
        ### then stochastic approximation
        if(class(temp_ergmAPL_rlt) == "try-error"){
          temp_ergmAPL_rlt <- try(Bergm::ergmAPL(formula = temp.form,
                                                 ladder = 500,
                                                 aux.iters = 2*netsize[i]^2,
                                                 n.aux.draws = 1500,
                                                 aux.thin = 50, estimate = "MLE",
                                                 control = control.ergm(main.method = c("Stochastic-Approximation"),
                                                                        MCMC.samplesize = 1024*4,
                                                                        MCMC.burnin = 4*network::network.size(temp.dat)^2) ))
        }
        
        ### then stepping 
        if(class(temp_ergmAPL_rlt) == "try-error"){
          temp_ergmAPL_rlt <- try(Bergm::ergmAPL(formula = temp.form,
                                                 ladder = 500,
                                                 aux.iters = 2*netsize[i]^2,
                                                 n.aux.draws = 1500,
                                                 aux.thin = 50, estimate = "MLE",
                                                 control = control.ergm(main.method = c("Stepping"),
                                                                        MCMC.samplesize = 1024*4,
                                                                        MCMC.burnin = 4*network::network.size(temp.dat)^2) ))
        }
        #
        temp.APL.calc.time2 <- Sys.time()
        #
        APL.calc.time.vec[i] <- difftime(temp.APL.calc.time2,
                                         temp.APL.calc.time1,
                                         units = "mins")
        
        W_list[[i]] <- temp_ergmAPL_rlt$W
        logC_list[[i]] <- temp_ergmAPL_rlt$logC
        Theta_MLE_list[[i]] <- temp_ergmAPL_rlt$Theta_MLE
        Theta_MPLE_list[[i]] <- temp_ergmAPL_rlt$Theta_PL
        
        # 
        print(temp_ergmAPL_rlt$Theta_MLE)
        print(temp_ergmAPL_rlt$Theta_PL)
        print(temp_ergmAPL_rlt$W)
        print(temp_ergmAPL_rlt$logC)
      }
    }
  } else{
    W_list <- W.list
    logC_list <- logC.list
    Theta_MLE_list <- Theta.MLE.list
    Theta_MPLE_list <- Theta.MPLE.list
  }
  ####################
  if(!(match.arg(Z.step) == "APL" | match.arg(theta.step) == "APL")){
    Theta_MPLE_list <- lapply(mplesetup_list, function(x) x$theta_MPLE)
  }
  ####################
  APL.calc.end.time <- Sys.time()
  ####################
  ####################
  # main MCMC sampling 
  j = 2
  while(j <= 1+main.iters){
    
    # -- update tau --#
    # update tau (from dirichlet), now we have Z[j,]
    tau[j,] <- gtools::rdirichlet(n=1, 
                                  alpha = (rep(alpha[1],K) + table(factor(Z[j-1,], levels = 1:K))) ) # K-dimensional vector
    cat("tau", tau[j,], "\n")
    # -- tau updated --# 
    
    # update theta, using random walk metropolis update or exchange algorithm
    # update theta_k in random order
    nK <- sample(1:K, K)
    for(k in nK){
      # proposing a new theta for k-th component
      # browser()
      theta1 <- rmvnorm(1,mean=theta[j-1,,k], sigma = sigma.epsilon)
      # log prior vector
      pr <- dmvnorm( rbind(theta1, theta[j-1,,k]), mean= prior.mean, sigma = prior.sigma, log=TRUE)
      # adjusted pseudolikelihood vector
      # lr <- c( logPL.corr(theta = theta[j-1,,k]) )
      k_id <- which(Z[j-1,] == k) # index of observed data whose corresponding assigned cluster is k 
      # 
      cat("The",k,"-th cluster:",k_id, "\n")
      if(length(k_id) >= 1){
        ll_mat <- matrix(NA, nrow=2,ncol = length(k_id))
        for(l in 1:length(k_id)){
          # network index k_id[l]
          if(match.arg(theta.step)=="PL"){
            ll_mat[2,l] <- logPL(theta = theta[j-1,,k], y=mplesetup_list[[k_id[l]]]$response,
                                 X = mplesetup_list[[k_id[l]]]$predictor, weights = mplesetup_list[[k_id[l]]]$weights,
                                 size = mplesetup_list[[k_id[l]]]$size, size.offset = size.offset)
            ll_mat[1,l] <- logPL(theta = c(theta1), y=mplesetup_list[[k_id[l]]]$response,
                                 X = mplesetup_list[[k_id[l]]]$predictor, weights = mplesetup_list[[k_id[l]]]$weights,
                                 size = mplesetup_list[[k_id[l]]]$size, size.offset = size.offset)
          }
          else if(match.arg(theta.step)=="exchange"){
            # exchange step for updating cluster-specific parameters, the proposed theta for cluster k is theta1
            # simulation step, simulate network from theta1
            ## figure out network size
            temp_coef <- theta1
            if(size.offset){
              temp_size <- network::network.size(dat[[ k_id[l] ]])
              temp_coef[1] <- theta1[1] - log(temp_size)
            }
            
            temp_stats <- ergm_MCMC_slave(Clist = Clist_list[[ k_id[l] ]], 
                                          proposal = proposal_list[[ k_id[l] ]], 
                                          eta = temp_coef, 
                                          control = control, 
                                          verbose = FALSE)$s + sy[k_id[l],] # set up them first
            # obs stats is named sy, which is of size n by p
            # likelihood ratio calculation step
            ll_mat[1,l] <- sum(theta[j-1,,k] * temp_stats) + sum( theta1 * sy[k_id[l],] )
            ll_mat[2,l] <- sum(theta[j-1,,k] * sy[k_id[l],]) + sum( theta1 * temp_stats )
          } else if(match.arg(theta.step) == "APL"){
            # adjusted pseudo-likelihood
            ll_mat[2,l] <- logAPL(theta = theta[j-1,,k], 
                                  y=mplesetup_list[[k_id[l]]]$response,
                                  X = mplesetup_list[[k_id[l]]]$predictor, 
                                  weights = mplesetup_list[[k_id[l]]]$weights,
                                  size = mplesetup_list[[k_id[l]]]$size, 
                                  size.offset = size.offset,
                                  theta_MLE = Theta_MLE_list[[k_id[l]]],
                                  theta_MPLE = Theta_MPLE_list[[k_id[l]]],
                                  W = W_list[[k_id[l]]],
                                  logC = logC_list[[k_id[l]]])
            #####
            ll_mat[1,l] <- logAPL(theta = c(theta1), y=mplesetup_list[[k_id[l]]]$response,
                                  X = mplesetup_list[[k_id[l]]]$predictor, 
                                  weights = mplesetup_list[[k_id[l]]]$weights,
                                  size = mplesetup_list[[k_id[l]]]$size, 
                                  size.offset = size.offset,
                                  theta_MLE = Theta_MLE_list[[k_id[l]]],
                                  theta_MPLE = Theta_MPLE_list[[k_id[l]]],
                                  W = W_list[[k_id[l]]],
                                  logC = logC_list[[k_id[l]]])
          }
        }
        lr <- apply(ll_mat,1,sum) # row-sum
      }
      else {
        lr <- c(0,0)
      }
      # log acceptance ratio, proposal is symmetric
      beta.theta <- (lr[1] - lr[2]) + (pr[1] - pr[2]) 
      print(beta.theta)
      cat("\n")
      # 
      if(is.nan(beta.theta)){
        # browser()
        beta.theta <- (pr[1] - pr[2]) 
      }
      # browser()
      # accept if greater than acceptance ratio 
      if(beta.theta >= log(runif(1))){
        theta[j,,k] <- theta1
        acc.counts[j,k] <- TRUE
      }
      else {
        theta[j,,k] <- theta[j-1,,k]
      }
    }
    
    # update the assignment label (total of n labels) Z
    for(i in 1:n){
      unnorm_p <- rep(NA, length = K)
      unnorm_logp <- rep(NA, length = K)
      if(match.arg(Z.step)=="PL"){
        # pseudo.step = TRUE indicates that we want to use pseudo likelihood to update
        # the cluster membership
        for(k in 1:K){
          # ignore the constant C
          cat("Previous tau", tau[j,k], "\n")
          cat("Previous theta", theta[j,,k], "\n")
          # calculating log pseudolikelihood
          temp_logPL <- logPL(theta = c(theta[j,,k]), 
                              y=mplesetup_list[[i]]$response,
                              X = mplesetup_list[[i]]$predictor, 
                              weights = mplesetup_list[[i]]$weights,
                              size = mplesetup_list[[i]]$size, size.offset = size.offset)
          cat("Corrected pseudo-like", temp_logPL, "\n")
          unnorm_logp[k] <- log(tau[j,k]) + temp_logPL
          cat("unnormalized log prob", unnorm_logp[k], "\n")
        }
        # A = max(unnorm_logp)
        cat("The unnormalized log-probability", unnorm_logp, "\n")
        # -Inf -Inf -884.7563 -1165.712 -422.0812 -517.0914 
        # lead to NaN
        # 
        lsum = sna::logSum(unnorm_logp[which(unnorm_logp!=-Inf)])
        unnorm_p <- exp(unnorm_logp-lsum)
        cat("The unnormalized probability", unnorm_p, "\n")
        cat("lsum=", lsum, "\n")
        print(paste("j=",j, "i=",i, "k=",k))
        cat("\n")
        # Update Z...
        Z[j,i] <- sample(x=1:K, size=1, prob = unnorm_p)
      } else if(match.arg(Z.step) == "full"){
        # if using full likelihood
        for(k in 1:K){
          cat("Previous tau", tau[j,k], "\n")
          cat("Previous theta", theta[j,,k], "\n")
          # calculating log pseudolikelihood
          temp_logL <- logL(theta = c(theta[j,,k]),
                            mplesetup = mplesetup_list[[i]],
                            size.offset = size.offset)
          # 
          cat("full likelihood", temp_logL, "\n")
          unnorm_logp[k] <- log(tau[j,k]) + temp_logL
          cat("unnormalized log prob", unnorm_logp[k], "\n")
        }
        # A = max(unnorm_logp)
        # lsum = sna::logSum(unnorm_logp)
        lsum = sna::logSum(unnorm_logp[which(unnorm_logp!=-Inf)])
        unnorm_p <- exp(unnorm_logp-lsum)
        cat("The unnormalized probability", unnorm_p, "\n")
        # Update Z...
        Z[j,i] <- sample(x=1:K, size=1, prob = unnorm_p)
      } else if(match.arg(Z.step) == "APL"){
        # adjusted pseudo likelihood
        for(k in 1:K){
          cat("Previous tau", tau[j,k], "\n")
          cat("Previous theta", theta[j,,k], "\n")
          # calculating log pseudolikelihood
          temp_logAPL <- logAPL(theta = c(theta[j,,k]),
                                y=mplesetup_list[[i]]$response,
                                X = mplesetup_list[[i]]$predictor,
                                weights = mplesetup_list[[i]]$weights,
                                size = mplesetup_list[[i]]$size,
                                size.offset = size.offset,
                                theta_MLE = Theta_MLE_list[[i]],
                                theta_MPLE = Theta_MPLE_list[[i]],
                                W = W_list[[i]],
                                logC = logC_list[[i]])
          
          # logPL(theta = c(theta[j-1,,k]), y=mplesetup_list[[i]]$response,
          #                   X = mplesetup_list[[i]]$predictor, 
          #                   weights = mplesetup_list[[i]]$weights,
          #                   size = mplesetup_list[[i]]$size, size.offset = size.offset)
          cat("Corrected pseudo-like", temp_logAPL, "\n")
          unnorm_logp[k] <- log(tau[j,k]) + temp_logAPL
          cat("unnormalized log prob", unnorm_logp[k], "\n")
        }
        # A = max(unnorm_logp)
        cat("The unnormalized log-probability", unnorm_logp, "\n")
        # -Inf -Inf -884.7563 -1165.712 -422.0812 -517.0914 
        # lead to NaN
        # 
        lsum = sna::logSum(unnorm_logp[which(unnorm_logp!=-Inf)])
        unnorm_p <- exp(unnorm_logp-lsum)
        cat("The unnormalized probability", unnorm_p, "\n")
        cat("lsum=", lsum, "\n")
        print(paste("j=",j, "i=",i, "k=",k))
        cat("\n")
        
        # 
        if(all(unnorm_logp == -Inf)){
          Z[j,i] <- sample(x=1:K, size=1)
        } else{
          # Update Z...
          Z[j,i] <- sample(x=1:K, size=1, prob = unnorm_p)
        }
      }
    }
    
    # #-- update alpha using random-walk metropolis --#
    # # propose a new value
    # alpha1 <- exp(rnorm(1, mean=log(alpha[1]), sd = sqrt(alpha.epsilon)))
    # # calculate the logarithm of acceptance probability
    # # browser()
    # beta.alpha <- (sum( (alpha1-1)*log(tau[j,]) ) + log(gamma(K * alpha1)) + log(dgamma(alpha1,shape = a,rate=a*K)) - K*log(gamma(alpha1))) - (sum((alpha[j-1])*log(tau[j,]))) + log(gamma(K * alpha[j-1])) + log(dgamma(alpha[j-1], shape=a, rate=a*K)) - K*log(gamma(alpha[j-1])) + log(1/alpha[j-1]) - log(1/alpha1)
    # # accept the proposal based on the value of beta.alpha
    # if(beta.alpha >= log(runif(1))){
    #   # accept the proposal
    #   alpha[j] <- alpha1
    #   # acc.counts[j,k] <- TRUE
    # }
    # else{
    #   # reject the proposal
    #   # theta[j,,k] <- theta[j-1,,k]
    #   alpha[j] <- alpha[j-1]
    # }
    # # -- alpha updated --#
    # re-order parameters to deal with label switching issue, from small to large, using the order constraints on 
    # edge parameter
    if(order.constraints){
      theta_order <- order(theta[j,1,], decreasing = FALSE)
      theta[j,,] <- theta[j,,theta_order]
      tau[j,] <- tau[j,theta_order]
      # also need to permute the latent cluster assignment... 
      temp_Z <- data.frame(old_label = Z[j,])
      reference_Z <- data.frame(old_label=1:K, correct_label = theta_order)
      correct_Z <- dplyr::left_join(temp_Z, reference_Z, by = c("old_label" = "old_label"))
      Z[j,] <- correct_Z$correct_label
    }
    else{
      # randomly permute their labels
      theta_order <- sample(1:K,size=K,replace = FALSE)
      theta[j,,theta_order] <- theta[j,,]
      tau[j,theta_order] <- tau[j,]
      # also need to permute the latent cluster assignment... 
      temp_Z <- data.frame(old_label = Z[j,])
      reference_Z <- data.frame(old_label=1:K, correct_label = theta_order)
      correct_Z <- dplyr::left_join(temp_Z, reference_Z, by = c("old_label" = "old_label"))
      Z[j,] <- correct_Z$correct_label
    }
    # 
    cat("current theta", j, "-th iteration", theta[j,,], "\n")
    # update K_nonzero: number of components with at least one observation
    K_nonzero[j] <- length(table(Z[j,])) # only clusters with observations are counted
    # update j
    j = j+1
  }
  samp.end.time <- Sys.time()
  #### figure out the number of clusters using posterior mode of K
  # Khat_table = table(factor(K_nonzero[(burn.in+1):(main.iters+1)], 
  #                           levels=1:K))
  # thinning the MCMC draws 
  # browser()
  K_nonzero.thin = K_nonzero[seq(burn.in+1, main.iters+1, by=thin)]
  tau.thin = tau[seq(burn.in+1, main.iters+1,by=thin),]
  Z.thin = Z[seq(burn.in+1, main.iters+1,by=thin),]
  theta.thin = theta[seq(burn.in+1, main.iters+1,by=thin),,]
  
  #
  if(p == 1){
    theta.thin <- array(theta.thin, dim = c(length(seq(burn.in+1, main.iters+1,by=thin)),
                                            p,
                                            K))
  }
  
  # alpha.thin = alpha[seq(burn.in+1, main.iters+1, by=thin)]
  # 
  Khat_table = table(factor(K_nonzero.thin, 
                            levels=1:K))
  #### 
  Khat = as.numeric(which(Khat_table == max(Khat_table))[1])
  ####
  ####
  #### K-centroid clustering 
  # first identify those iterations with the number of non-zero components = Khat
  sub_id <- which(K_nonzero.thin == Khat)
  # keep those after burn_in
  # sub_id <- sub_id[sub_id >= (burn.in+1)]
  M0 <- length(sub_id)
  # pull out the theta.thin's to form a matrix : p columns, Khat * M0 rows
  theta_sub <- matrix(NA, nrow = M0*Khat, ncol = p)
  tau_sub <- rep(NA, length=M0*Khat)
  # browser()
  # alpha_sub <- rep(NA, length=M0)
  for(jj in 1:M0){
    if(p == 1){
      theta_sub[(1+(jj-1)*Khat):(jj*Khat),] <- theta.thin[sub_id[jj], ,unique(Z.thin[sub_id[jj],])]
    } else{
      theta_sub[(1+(jj-1)*Khat):(jj*Khat),] <- t(theta.thin[sub_id[jj], ,unique(Z.thin[sub_id[jj],])])
    }
    tau_sub[(1+(jj-1)*Khat):(jj*Khat)] <- tau.thin[sub_id[jj],unique(Z.thin[sub_id[jj],])]
    # alpha_sub[jj] <- alpha.thin[jj]
  }
  # k-centroid clustering
  # kcca_rlt <- kcca(Nclus, k=Khat, family=kccaFamily("kmedians"),
  #             control=list(initcent="kmeanspp"))
  # do clustering if and only if Khat > 1
  # browser()
  if(Khat > 1){
    #
    kcca_rlt <- Cluster_Medoids(theta_sub, 
                                clusters = Khat, 
                                distance_metric = "mahalanobis")
    # 
    # kcca_rlt$clusters is a sequence of cluster membership indicators
    # check if the clustering result is valid for each iteration of parameters
    theta_sub_sub <-  array(NA, dim = c(M0, p, Khat))
    tau_sub_sub <- matrix(NA, nrow = M0, ncol = Khat)
    alpha_sub_sub <- rep(NA, length = M0)
    M0_sub = 0 # number of iterations retained 
    ###########
    for(jj in 1:M0){
      # if the cluster assignment is a permutation
      if(all(sort(kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)]) == c(1:Khat)) ){
        if(p == 1){
          theta_sub_sub[M0_sub+1,,kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)]] <- (theta_sub[(1+(jj-1)*Khat):(jj*Khat),])
        } else{
          theta_sub_sub[M0_sub+1,,kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)]] <- t(theta_sub[(1+(jj-1)*Khat):(jj*Khat),])
        }
        
        # alpha_sub_sub[M0_sub+1] <- alpha_sub[jj]
        tau_sub_sub[M0_sub+1,kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)]] <- tau_sub[(1+(jj-1)*Khat):(jj*Khat)]
        # tau_sub_sub[jj,kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)]] <- tau.thin[jj, ] 
        # theta[sub_id[jj], ,unique(Z[sub_id[jj],])]
        M0_sub = M0_sub + 1
      }
    }
    
    # only keep the complete.cases
    theta_sub_sub = theta_sub_sub[1:M0_sub,,]
    # alpha_sub_sub = alpha_sub_sub[1:M0_sub]
    tau_sub_sub = tau_sub_sub[1:M0_sub,]
    # 
    tau_sub_sub <- sweep(tau_sub_sub, 1, apply(tau_sub_sub,1,sum), "/")
    # normalize the tau's to make sure rowsum is 1
    #### the cluster with smaller value on edge parameter should have smaller cluster id number
    ####
    cluster_rank <- rank(kcca_rlt$medoids[,1])
    ####
    for(jj in 1:nrow(tau_sub_sub)){
      tau_sub_sub[jj,cluster_rank] <- tau_sub_sub[jj,]
      theta_sub_sub[jj,,cluster_rank] <- theta_sub_sub[jj,,]
    }
    #
  } else if(Khat == 1){
    ###########
    theta_sub_sub <-  array(NA, dim = c(M0, p, Khat))
    tau_sub_sub <- matrix(NA, nrow = M0, ncol = Khat)
    alpha_sub_sub <- rep(NA, length = M0)
    M0_sub = 0 # number of iterations retained 
    ###########
    for(jj in 1:M0){
      # when Khat = 1, no need to check permutation
      # if the cluster assignment is a permutation
      # if(all(sort(kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)]) == c(1:Khat)) ){
      # theta_sub_sub[M0_sub+1,, kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)] ] <- t(theta_sub[(1+(jj-1)*Khat):(jj*Khat),])
      # # alpha_sub_sub[M0_sub+1] <- alpha_sub[jj]
      # tau_sub_sub[M0_sub+1, kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)] ] <- tau_sub[(1+(jj-1)*Khat):(jj*Khat)]
      # tau_sub_sub[jj,kcca_rlt$clusters[(1+(jj-1)*Khat):(jj*Khat)]] <- tau.thin[jj, ] 
      # theta[sub_id[jj], ,unique(Z[sub_id[jj],])]
      if(p == 1){
        theta_sub_sub[M0_sub+1,, 1:Khat ] <- theta_sub[(1+(jj-1)*Khat):(jj*Khat),]
      } else{
        theta_sub_sub[M0_sub+1,, 1:Khat ] <- t(theta_sub[(1+(jj-1)*Khat):(jj*Khat),])
      }
      
      # alpha_sub_sub[M0_sub+1] <- alpha_sub[jj]
      tau_sub_sub[M0_sub+1, 1:Khat ] <- tau_sub[(1+(jj-1)*Khat):(jj*Khat)]
      M0_sub = M0_sub + 1
      # }
    }
    tau_sub_sub <- sweep(tau_sub_sub, 1, apply(tau_sub_sub,1,sum), "/")
    # normalize the tau's to make sure rowsum is 1
    #### the cluster with smaller value on edge parameter should have smaller cluster id number
    ####
    # cluster_rank <- rank(kcca_rlt$medoids[,1])
    # ####
    # for(jj in 1:nrow(tau_sub_sub)){
    #   tau_sub_sub[jj,cluster_rank] <- tau_sub_sub[jj,]
    #   theta_sub_sub[jj,,cluster_rank] <- theta_sub_sub[jj,,]
    # }
  }
  # 
  # # convert the dimension in case Khat = 1
  # if(Khat == 1){
  #   theta_sub_sub <- array(theta_sub_sub, dim = c(dim(theta_sub_sub), 1))
  #   # tau_sub_sub <- array(tau_sub_sub, dim = c(length(tau_sub_sub), 1) )
  #   tau_sub_sub <- matrix(tau_sub_sub, nrow = length(tau_sub_sub), ncol = 1)
  #   kcca_rlt$medoids <- matrix(kcca_rlt$medoids, ncol = 1)
  # }
  # 
  # tau_sub_sub <- sweep(tau_sub_sub, 1, apply(tau_sub_sub,1,sum), "/")
  # # normalize the tau's to make sure rowsum is 1
  # #### the cluster with smaller value on edge parameter should have smaller cluster id number
  # ####
  # cluster_rank <- rank(kcca_rlt$medoids[,1])
  # ####
  # for(jj in 1:nrow(tau_sub_sub)){
  #   tau_sub_sub[jj,cluster_rank] <- tau_sub_sub[jj,]
  #   theta_sub_sub[jj,,cluster_rank] <- theta_sub_sub[jj,,]
  # }
  ####
  # browser()
  #### 
  #### DIC calculation
  ####
  DIC.calc.start.time <- Sys.time()
  DIC3 = NA
  DIC4 = NA
  # Thinning the MCMC samples
  # id.thin <- seq(from=burn.in+1, to = main.iters+1, by = thin)
  # theta.thin <- array(theta[id.thin,,], dim=c(length(id.thin),p,K))
  # tau.thin <- matrix(tau[id.thin,],ncol=K)
  # Z.thin <- Z[id.thin,]
  # browser()
  # if DIC.calc is true, we are going to calculate the DIC value
  # DIC calculation does not require Z here
  if(DIC.calc){
    # if(is.null(DIC.samp.size)){
    #   DIC.samp.size <- nrow(tau_sub_sub)
    # }
    # DIC3 = DIC3(mplesetup_list = mplesetup_list,
    #             eta = tau_sub_sub,
    #             theta = theta_sub_sub,
    #             DIC.samp.size = DIC.samp.size, 
    #             size.offset = size.offset)  
    
    DIC3 = DIC_mixture_ergm(mplesetup_list = mplesetup_list,
                            tau = tau_sub_sub,
                            theta = theta_sub_sub,
                            DIC.samp.size = nrow(tau_sub_sub),
                            size.offset = size.offset,
                            lik = Z.step,
                            theta_MLE_list = Theta_MLE_list,
                            theta_MPLE_list = Theta_MPLE_list,
                            W_list = W_list,
                            logC_list = logC_list)
    
    # DIC3 <- function(mplesetup_list, 
    #                  tau, 
    #                  theta, 
    #                  DIC.samp.size=500, 
    #                  size.offset=TRUE, 
    #                  lik = c("PL", "APL", "full"),
    #                  theta_MLE_list=NULL,
    #                  theta_MPLE_list=NULL,
    #                  W_list=NULL,
    #                  logC_list=NULL,...)
  }
  DIC.calc.end.time <- Sys.time()
  ############################# BIC calculation
  BIC.calc.start.time <- Sys.time()
  #############################
  BIC1 <- NA
  if(BIC.calc){
    BIC1 <- BIC_mixture_ergm(mplesetup_list = mplesetup_list,
                             p=p,
                             tau = tau_sub_sub,
                             theta = theta_sub_sub,
                             size.offset = size.offset,
                             lik = Z.step,
                             theta_MLE_list = Theta_MLE_list,
                             theta_MPLE_list = Theta_MPLE_list,
                             W_list = W_list,
                             logC_list = logC_list)
    cat("BIC value:", BIC1)
    cat("\n")
    # BIC_mixture_ergm <- function(mplesetup_list,
    #                              p,
    #                              tau, 
    #                              theta, 
    #                              size.offset=TRUE, 
    #                              lik = c("PL", "APL", "full"),
    #                              theta_MLE_list=NULL,
    #                              theta_MPLE_list=NULL,
    #                              W_list=NULL,
    #                              logC_list=NULL,...)
  }
  #############################
  BIC.calc.end.time <- Sys.time()
  ############################# posterior class probability
  # class_prob_PL <- function(mplesetup_list, eta, theta, samp.size,...)
  post_class_prob <- NA
  post_class_label <- rep(NA, length = n)
  # browser()
  if(post.class.calc){
    post_class_prob <- class_prob(mplesetup_list = mplesetup_list,
                                  tau = tau_sub_sub,
                                  theta = theta_sub_sub,
                                  samp.size = nrow(tau_sub_sub),
                                  size.offset = size.offset,
                                  lik = Z.step,
                                  num.cores = num.cores,
                                  theta_MLE_list = Theta_MLE_list,
                                  theta_MPLE_list = Theta_MPLE_list,
                                  W_list = W_list,
                                  logC_list = logC_list)
    ###########################
    post_class_prob <- matrix(post_class_prob, ncol = Khat)
    for(i in 1:n){
      post_class_label[i] <- which(post_class_prob[i,] == max(post_class_prob[i,]))
    }
    ###########################
  }
  #############################
  cat("post_class_label", post_class_label)
  #
  cat("\n")
  #############################
  rm(Z)
  # rm(theta)
  rm(tau)
  # rm(acc.counts)
  gc()
  # posterior predictive checks
  # summary stats calculated on the observed data
  # summary stats calculated on simulated data
  obs_post_check_stats <- NA
  sim_post_check_stats <- NA
  cat("obs_post_check_stats", obs_post_check_stats)
  cat("\n")
  ##############################
  obs_post_check_stats <- post_check_mixture_ergm(dat)
  ##############################
  if(post.check){
    # obs_post_check_stats <- post_check_mixture_ergm(dat)
    # browser()
    # observed data summary stats
    # obs_post_check_stats <- post_check_mixture_ergm(dat)
    # browser()
    # data simulated from the posterior samples 
    sim_dat <- vector(mode="list", length = n)
    # 
    for(i in 1:length(sim_dat)){
      temp_id <- i
      y <- dat[[i]]
      # y <- temp_sim_mixture_ergm$sim_networks[[temp_id]]
      temp_membership <- as.numeric(factor(y%v%"x"))
      temp_Z <- sample(1:Khat, size=1, prob = tau_sub_sub[(i-1)%%M0_sub+1,]) 
      # which(rmultinom(n=1, size = 1, prob = tau_sub_sub[i,])==1)
      temp_network_size <- network::network.size(y)
      temp_coef <- theta_sub_sub[(i-1)%%M0_sub+1,,temp_Z]
      # rlt_list[[ii]][[ Khat[ii] ]]$theta.thin[i,,temp_Z]
      if(size.offset){
        temp_coef[1] <- temp_coef[1] - log(temp_network_size)
      }
      # temp_coef[1] <- temp_coef[1] - log(temp_network_size)
      temp_sim_network <- simulate(form,
                                   nsim=1,
                                   seed = 1234*i,
                                   coef = temp_coef,
                                   basis = y,
                                   control=control.simulate.formula(MCMC.burnin = 20*temp_network_size^2,
                                                                    MCMC.interval = 1))
      # sim_modularity_1[i,j] <- modularity(intergraph::asIgraph(temp_sim_network),
      #                                      membership = factor(y%v%"party") )
      sim_dat[[i]] <- temp_sim_network 
    }
    # simulated data, stats
    sim_post_check_stats <- post_check_mixture_ergm(sim_dat)
    rm(sim_dat) # sim_dat is no longer needed
    gc()
  }
  # remove unneccessary objects and garbage collection to reduce memory usage
  #############################
  # browser()
  out = list(samp.Time = difftime(samp.end.time, samp.start.time),
             init.Time = difftime(samp.start.time, init.start.time),
             DIC.calc.Time = difftime(DIC.calc.end.time,DIC.calc.start.time),
             APL.calc.Time = difftime(APL.calc.end.time, APL.calc.start.time),
             form = form, 
             DIC3 = DIC3, 
             DIC4 = DIC4, 
             BIC1 = BIC1,
             burn.in = burn.in, 
             main.iters=main.iters, 
             thin = thin,
             tau.out = tau_sub_sub, 
             theta.out = theta_sub_sub, 
             aux.iters=aux.iters,
             post.class.label = post_class_label, 
             post.class.prob = post_class_prob, 
             size.offset = size.offset,
             Khat = Khat,
             sy = sy,
             # mplesetup.list = mplesetup_list,
             APL.calc.time.vec = APL.calc.time.vec,
             theta.MLE.list = Theta_MLE_list,
             W.list = W_list,
             logC.list = logC_list,
             theta.MPLE.list = Theta_MPLE_list,
             sim.post.check.stats = sim_post_check_stats,
             obs.post.check.stats = obs_post_check_stats,
             M0 = M0,
             M0_sub = M0_sub,
             acc.counts = acc.counts,
             theta = theta)
  # theta = theta, tau = tau, Z = Z,
  # W_list <- vector(mode="list", length=n)
  # logC_list <- vector(mode="list", length=n)
  # Theta_MLE_list <- vector(mode="list", length=n)
  # Theta_MPLE_list <- vector(mode="list", length=n)
  return(out)
  #              mple.coef.colabel = mple_coef_colabel, 
}
########################################
# mple setup list?
