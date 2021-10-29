##################
##################
### This script is written to demonstrate how to run our code
### toy example
##################
##################
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("scales",
              "foreach", "doParallel", 
              "dplyr", "xtable",
              "pwr", "ergm", 
              "coda", "flexclust", "ClusterR", 
              "gtools", "sna")
ipak(packages)
###################
# params
source("./mixture.ergm.functions.R")
###
form_sim <- y ~ edges + gwesp(0.25, fixed=TRUE) + nodematch("x", diff=F)
###
thin <- 50
main_iters <- 60000
sigma2 <- 0.0025
###
init_method <- "mple"
burn_in <- ceiling(main_iters*0.6)
####################
params_sim_40 <- matrix(c(-0.85, -0.10, -0.10,
                          -3.45, 0.75, 2,
                          -5.10, 2.5, 0.5),
                        nrow=3, byrow=T)
####################
num_replicates <- 1
seed_sim <- 12345
p=ncol(params_sim_40)
### iii is the indicator variable for the setting
iii = 1
num_cores = 1
alpha = 0.001
###
temp_time1 <- Sys.time()
###
cluster_size <- 10
network_size <- 40
cluster_num <- 3
###
tau_sim <- rep(1/cluster_num,
               cluster_num)
###
#####################################################
registerDoParallel(num_cores)
#####################################################
sim_networks <- readRDS("./demo_sim_networks.rds") 

# sim_networks is a list of networks simulated from mixtures of ERGMs 
# total sample size: 30; 10 from each component defined by the model specification form_sim
# and model parameters given in params_sim_40

#########################################
# ---- run the model ---- #
# Pseudo likelihood
time1 <- Sys.time()
rlt_PL <- try(bayes_mixture_ergm_overclustering_fixedalpha(sim_networks,
                                                             form=form_sim,
                                                             p=p,
                                                             K=cluster_num*2,
                                                             prior.mean=c(-1,0,0),
                                                             prior.sigma=diag(25,p),
                                                             alpha = alpha,
                                                             main.iters = main_iters,
                                                             num_cores=1,
                                                             sigma.epsilon = diag(rep(0.0025,p)),
                                                             thin = thin,
                                                             burn.in = burn_in,
                                                             Z.step = "PL",
                                                             theta.step = "PL",
                                                             order.constraints = FALSE,
                                                             size.offset = FALSE,
                                                             post.check = FALSE))
time2 <- Sys.time()
#################################

## posterior mean (each column corresponds to a cluster)
rlt_PL$theta.out %>% apply(c(2,3), mean)




