############################
###MPLRJ function#######
############################

#this function runs the MPL_BD algorithm

#It takes as input:
#- data (an nxp matrix)
#- iter_vec_thin. The iterations at which it needs to produce output.

#It outputs:
# the running time (at every iteration of iter_vec_thin)
# iter_vec_thin itself
# the edge inclusion probabilities after all iterations 
# the amount of edges in the graph (at every iteration of iter_vec_thin)
# the final graph of the Markov chain


MPLBD_solve = function(data,iter_vec_thin){
  
  #parameters
  p = ncol(data)
  algorithm_name = "MPLBD"
  burnin = 0
  MCMCstart = "empty"
  jump = 1
  save = FALSE
  cores = 8
  g.prior = 0.01
  verbose = FALSE
  
  #intialize outputs
  len_iter = length(iter_vec_thin)
  edge_vec = rep(0,len_iter)
  time_vec = rep(0,len_iter)
  
  #run the MCMC iterations 
  weights_old = 0
  plinks_old = 0
  for (j in 1:len_iter){
    
    #determine amount of iterations
    olditer = 0
    if (j > 1){
      olditer = iter_vec_thin[j-1]
    }
    newiter = iter_vec_thin[j]
    iter = newiter-olditer
    
    #run bdgraph for maxiter iterations
    sample_mpl_bd  = bdgraph.mpl( data = data, algorithm = "bdmcmc", iter = iter, burnin = burnin, jump = jump, save = save,cores=cores,g.start=MCMCstart,g.prior=g.prior,verbose=verbose)  
    
    #calculate new p matrix
    weights_new = sample_mpl_bd$sum_weights
    plinks_new = (weights_old*plinks_old + weights_new*sample_mpl_bd$p_links)/(weights_old+weights_new)
    predictor = plinks_new[upper.tri(plinks_new)]
    
    #calculate the runtime
    time_init = 0
    time_end = 0
    if (j==1){time_init = sample_mpl_bd$time_init}
    if (j==len_iter){time_end = sample_mpl_bd$time_end}
    runtime = time_init + sample_mpl_bd$time_all_iterations + time_end
    
    #save metrics
    edge_vec[j] = sum(sample_mpl_bd$last_graph)/2
    time_vec[j] = runtime
    
    #save data for next run
    plinks_old = plinks_new
    weights_old = weights_old+weights_new
    MCMCstart = sample_mpl_bd$last_graph
    
    #save data in .Rdata file, only for large problems
    if (p > 250){
      output_per_iter = list()
      output_per_iter$time = sum(time_vec)/10^6
      output_per_iter$iter = newiter
      output_per_iter$total_edges = sum(sample_mpl_bd$last_graph)/2
      output_per_iter$weight = weights_old
      output_per_iter$predictor = predictor
      filename = paste0("yyz_results_",algorithm_name,"_iter",j,".Rdata")
      save(output_per_iter,file = filename)
    }
    
  }
  
  time_vec = cumsum(time_vec)
  time_vec = time_vec/(10^6)
  return(list(sum_weights=weights_old,edge_vec=edge_vec,time_vec=time_vec,iter_vec_thin=iter_vec_thin,predictor=predictor,last_graph=sample_mpl_bd$last_graph))
}