source("3.metric_functions.R")

############################
###MPLRJ function#######
############################

#this function runs the MJ-MCMC algorithm

#It takes as input:
#- data from BDgraph. This "data" type contains data, a true graph. 
#- iter_vec_thin. The iterations at which it needs to produce output.
#- The jump_speed (i.e. the value of epsilon, between 0 and 1).
#- The max_jump (i.e. the percentage of the total edges that are allowed to flip in a single iteration, between 0 and 1).
#- Diminishing (0 for homogeneous, 1 for slow decay, 2 for fast decay).


#It outputs:
#- the AUC-PR at every iteration number indicated in iter_vec_thin
#- the AUC-ROC at every iteration number indicated in iter_vec_thin
#- the p+ at every iteration number indicated in iter_vec_thin
#- the p- at every iteration number indicated in iter_vec_thin
#- the running time at every iteration number indicated in iter_vec_thin
#- other metrics

MPLRJ_jump_solve = function(data,iter_vec_thin,jump_speed,max_jump,diminishing){
  
  #obtain true graph
  G_true = data$G
  response = G_true[upper.tri(G_true)]
  
  #obtain amount of observations and variables
  n = nrow(data$data)
  p = ncol(data$data)
  density = data$density
  graph = data$graph
  rep = data$rep
  algorithm_name = "MPLRJ_jump"
  
  #parameters
  burnin = 0
  g.start = "empty"
  g.prior = 0.005
  jump = "variable"
  save = FALSE
  cores = 8
  MCMCstart = g.start
  plinks_old = 0
  weights_old = 0
  verbose = FALSE
  
  #run the algorithm, produce output at every iteration in the iter_vec_thin vector
  len = length(iter_vec_thin)
  AUC_vec = rep(0,len)
  pplus_vec = rep(0,len)
  pmin_vec = rep(0,len)
  AUC_PR_vec = rep(0,len)
  time_vec = rep(0,len)
  total_edges_vec = rep(0,len)
  jump_speed_vec = rep(0,len)
  size_selected_edges_vec = c()
  
  for (j in 1:len){
    
    #print progress
    print(iter_vec_thin[j])
    
    #determine amount of iterations
    olditer = 0
    if (j > 1){
      olditer = iter_vec_thin[j-1]
    }
    newiter = iter_vec_thin[j]
    iter = newiter-olditer
    
    #run bdgraph for maxiter iterations
    sample  = bdgraph.mpl( data = data, algorithm = "rjmcmc", diminishing = diminishing, iter_dim_start = olditer, max_jump=max_jump, jump_speed = jump_speed, iter = iter, burnin = burnin, jump = jump, save = save,cores=cores,g.start=MCMCstart,g.prior=g.prior,verbose=verbose)  
      
    #calculate new p matrix
    weights_new = iter
    weights_old = olditer
    plinks_new = (weights_old*plinks_old + weights_new*sample$p_links)/(weights_old+weights_new)
    predictor = plinks_new[upper.tri(plinks_new)]
    
    #calculate the AUC and the pplus and pmin
    AUC = calc_AUC_ROC(predictor=predictor,response=response)
    obj = calc_pplus_pmin(predictor=predictor,response=response)
    p_plus = obj$p_plus
    p_min = obj$p_min
    AUC_PR = calc_AUC_PR(predictor=predictor,response=response)
    
    #calculate the runtime
    time_init = 0
    time_end = 0
    if (j==1){time_init = sample$time_init}
    if (j==len){time_end = sample$time_end}
    runtime = time_init + sample$time_all_iterations + time_end
    
    #save metrics
    AUC_vec[j] = AUC
    pplus_vec[j] = p_plus
    pmin_vec[j] = p_min
    time_vec[j] = runtime
    AUC_PR_vec[j] = AUC_PR
    total_edges_vec[j] = sum(sample$last_graph)/2
    jump_speed_vec[j] = sample$final_jump_speed
    size_selected_edges_vec = c(size_selected_edges_vec,sample$size_selected_edges_vec)
    
    #store data for next run
    plinks_old = plinks_new
    MCMCstart = sample$last_graph
    
    #save output per iteration in Rdata file (only for large scale problems)
    if (p > 250){
      output_per_iter = list()
      output_per_iter$AUC = AUC
      output_per_iter$pplus = p_plus
      output_per_iter$pmin = p_min
      output_per_iter$time = sum(time_vec)/10^6
      output_per_iter$iter = newiter
      output_per_iter$AUC_PR = AUC_PR
      output_per_iter$total_edges = sum(sample$last_graph)/2
      output_per_iter$jump_speed_vec = jump_speed_vec
      output_per_iter$size_selected_edges_vec = sample$size_selected_edges_vec
      filename = paste0("yyz_results_p",p,"_n",n,"_",graph,"_",density,"_",algorithm_name,"_speed",jump_speed,"_maxjump",max_jump,"_dimin",diminishing,"_rep",rep,"_iter",j,".Rdata")
      save(output_per_iter,file = filename)
    }
    
  }
  time_vec = cumsum(time_vec)
  time_vec= time_vec/10^6 #bdgraph outputs the time in microseconds. Here we convert it to seconds
  
#add the values at time t = 0
  time_vec = c(0,time_vec)
  iter_vec_thin = c(0,iter_vec_thin)
  
  #recover the graph at time t = 0
  if (is.character(g.start)){
    if (g.start == "empty"){
      predictor_0 = rep(0,p*(p-1)/2)
    }
    if (g.start == "full"){
      predictor_0 = rep(1,p*(p-1)/2)
    }
  } else  {
    predictor_0 = g.start[upper.tri(g.start)]
  }
  
  #calculate the AUC and the pplus and pmin at time t = 0
  AUC_0 = calc_AUC_ROC(predictor=predictor_0,response=response)
  obj = calc_pplus_pmin(predictor=predictor_0,response=response)
  p_plus_0 = obj$p_plus
  p_min_0 = obj$p_min
  AUC_PR_0 = calc_AUC_PR(predictor=predictor_0,response=response)
  
  #add the initial values to the vector
  pplus_vec=c(p_plus_0,pplus_vec)
  pmin_vec=c(p_min_0,pmin_vec)
  AUC_vec = c(AUC_0,AUC_vec)
  AUC_PR_vec = c(AUC_PR_0,AUC_PR_vec)
  total_edges_vec = c(sum(predictor_0),total_edges_vec)
  
  return(list(AUC_vec=AUC_vec,
              pplus_vec=pplus_vec,
              pmin_vec=pmin_vec,
              AUC_PR_vec=AUC_PR_vec,
              time_vec=time_vec,
              iter_vec_thin=iter_vec_thin,
              predictor=predictor,
              total_edges_vec = total_edges_vec,
              jump_speed_vec = jump_speed_vec,
              size_selected_edges_vec = size_selected_edges_vec))
}  