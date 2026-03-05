
#######################################
#sanity tests on the true graphs#######
#######################################

#Goal: checking if the created graphs have the density that I want them to have

p = 1000
n_vec = c(400,1000)
graph = "random"
density_vec = c("sparse","dense")
rep_vec = 1:8

for (n in n_vec){
  for (density in density_vec){
    print(n)
    print(density)
    for (rep in rep_vec){
      filename = paste0("yy_truth_p",p,"_n",n,"_",graph,"_",density,"_rep",rep,".Rdata")
      load(file=filename)
      G = data_output$G
      total_edges = sum(G)/2
      edge_density = total_edges/(p*(p-1)/2)
      max_degree = max(colSums(G))
      text_print = paste0("Replication ",rep," has edge density " ,edge_density," and max degree ", max_degree)
      print(text_print)
    }
  }
}


############################################################
############retrieve iter_vec_thin##########################
############################################################

#Goal: checking if the algorithm ran for the amount MCMC iterations that I told them too

p = 1000
density = "dense"
n = 1000
graph = "random"
jump_speed_vec = c(0.001,0.01,0.05,0.25,0.75,0.25,0.25,1)
max_jump_vec = c(1,1,0.0025,0.01,0.05,0.0025,0.0025,1)
dimin_vec = c(0,0,0,0,0,1,2,0)
algorithm_vec = c(rep("MPLRJ_jump",7),"MPLBD")
rep = 1

for (i in 1:length(jump_speed_vec)){
  jump_speed = jump_speed_vec[i]
  max_jump = max_jump_vec[i]
  dimin = dimin_vec[i]
  algorithm = algorithm_vec[i]
  filename = paste0("zz_results_p",p,"_n",n,"_",graph,"_",density,"_",
                  algorithm,"_speed",jump_speed,"_maxjump",max_jump,"_dimin",dimin,"_rep",rep,".Rdata")
  load(file=filename)
  print(paste(output$iter_vec_thin,collapse=","))
}


############################################################
###did the max_jump cap the maximum available edges#########
############################################################

#Goal: it is important that the maximum jump size (i.e. "max_jump") only caps the amount of flips in the first couple
#of MCMC iterations. We check if this is indeed the case.

p = 1000
density = "dense"
n = 1000
graph = "random"

jump_speed = 0.75
max_jump = 0.01
dimin = 0
algorithm = "MPLRJ_jump"
jump_max = p*(p-1)/2*max_jump
jump_max

rep_vec = 1:8
for (rep in rep_vec){
  filename = paste0("zz_results_p",p,"_n",n,"_",graph,"_",density,"_",
                  algorithm,"_speed",jump_speed,"_maxjump",max_jump,"_dimin",dimin,"_rep",rep,".Rdata")
  load(file=filename)
  max_jump_obs = max(output$size_selected_edges_vec)
  indices = which(output$size_selected_edges_vec==max_jump_obs)
  
  print(max(indices))
}


############################################################
###is there a lot of variance between replications#########
############################################################

p = 1000
n = 1000            
graph = "random"
density = "sparse"         #("sparse" or "dense")
algorithm = "MPLRJ_jump"  #(MPLRJ_jump or MPL_BD)
jump_speed =  0.001
max_jump = 1
diminishing = 0
rep_vec = 1:8

#y-axis
AUC_on = FALSE
AUC_PR_on = TRUE
total_edges_on = FALSE

#x-axis
time_on = TRUE
iter_on = FALSE

#show each replication
rep_on = TRUE

#create empty plot
if (time_on){
  if (AUC_on){plot(NA,xlab="time (seconds)",ylab="AUC",xlim=c(0,max_time*1.2),ylim=c(0.5,1))}
  if (AUC_PR_on){plot(NA,xlab="time (seconds)",ylab="AUC PR",xlim=c(0,4000),ylim=c(0,1))}
  if (total_edges_on){plot(NA,xlab="time (seconds)",ylab="#edges",xlim=c(0,100),ylim=c(0,4000))}
}
if (iter_on){
  if (AUC_on){plot(NA,xlab="iterations",ylab="AUC",xlim=c(0,100),ylim=c(0.5,1))}
  if (AUC_PR_on){plot(NA,xlab="iterations",ylab="AUC PR",xlim=c(0,200),ylim=c(0,1))}
  if (total_edges_on){plot(NA,xlab="iterations",ylab="#edges",xlim=c(0,50),ylim=c(0,25000))}
}

#determine length vectors
rep = 1
filename = paste0("zz_results_p",p,"_n",n,"_",graph,"_",density,"_",algorithm,"_speed",jump_speed,"_maxjump",max_jump,"_dimin",diminishing,"_rep",rep,".Rdata")
load(file=filename)
len = length(output$iter_vec_thin)
sum_time = rep(0,len)
sum_AUC = rep(0,len)
sum_AUC_PR = rep(0,len)
sum_total_edges = rep(0,len)

#plot every single replication
col = 1
count = 0
for (rep in rep_vec){
  print(rep)
  filename = paste0("zz_results_p",p,"_n",n,"_",graph,"_",density,"_",algorithm,"_speed",jump_speed,"_maxjump",max_jump,"_dimin",diminishing,"_rep",rep,".Rdata")
  if (file.exists(file=filename)){
    count = count +1
    load(file=filename)
    sum_time = sum_time + output$time_vec
    sum_AUC = sum_AUC + output$AUC_vec
    sum_AUC_PR = sum_AUC_PR + output$AUC_PR_vec
    sum_total_edges = sum_total_edges + output$total_edges_vec
    
    if (rep_on){
      if (time_on){
        if (AUC_on){points(x=output$time_vec,y=output$AUC_vec,type="l",col=col)}
        if (AUC_PR_on){points(x=output$time_vec,y=output$AUC_PR_vec,type="l",col=col)}
        if (total_edges_on){points(x=output$time_vec,y=output$total_edges_vec,type="l",col=col)}
      }
      if (iter_on){
        if (AUC_on){points(x=output$iter_vec_thin,y=output$AUC_vec,type="l",col=col)}
        if (AUC_PR_on){points(x=output$iter_vec_thin,y=output$AUC_PR_vec,type="l",col=col)}
        if (total_edges_on){points(x=output$iter_vec_thin,y=output$total_edges_vec,type="l",col=col)}
      }
    }
  }
  else {
    print_text = paste0("Rep",rep,"does not exist")
    print(print_text)
  }
}

#plot the average of all the replications
if (time_on){
  if (AUC_on){points(x=sum_time/count,y=sum_AUC/count,type="l",col=col,lw=2)}
  if (AUC_PR_on){points(x=sum_time/count,y=sum_AUC_PR/count,type="l",col=col,lw=2)}
  if (total_edges_on){points(x=sum_time/count,y=sum_total_edges/count,type="l",col=col,lw=2)}
}
if (iter_on){
  if (AUC_on){points(x=output$iter_vec_thin,y=sum_AUC/count,type="l",col=col,lw=2)}
  if (AUC_PR_on){points(x=output$iter_vec_thin,y=sum_AUC_PR/count,type="l",col=col,lw=2)}
  if (total_edges_on){points(x=output$iter_vec_thin,y=sum_total_edges/count,type="l",col=col,lw=2)}
}
   
#######################################################################
###scatterplots to test the similarity in predictor vecs###############
#######################################################################

#Goal: for small epsilon, the edge inclusion probabilities should be similar to the BD-MPL algorithm. 
#we verify if this is actually the case by making a scatterplot

p = 1000
n = 400
graph = "random"
density = "dense"   #("sparse" or "dense")

#download MJ-MCMC
jump_speed = 0.001
max_jump = 1
dimin = 0
algorithm = "MPLRJ_jump"
filename = paste0("zz_results_p",p,"_n",n,"_",graph,"_",density,"_",algorithm,"_speed",jump_speed,"_maxjump",max_jump,"_dimin",diminishing,"_rep",rep,".Rdata")
load(file=filename)
mplrj_jump_pred = output$predictor

#download MPL-BD
jump_speed = 1
max_jump = 1
dimin = 0
algorithm = "MPLBD"
filename = paste0("zz_results_p",p,"_n",n,"_",graph,"_",density,"_",algorithm,"_speed",jump_speed,"_maxjump",max_jump,"_dimin",diminishing,"_rep",rep,".Rdata")
load(file=filename)
mplbd_pred = output$predictor

plot(x=mplrj_jump_pred,y=mplbd_pred,col = rgb(0, 0, 0, alpha = 0.03),cex = 0.2,asp=1,pch=16,xlim=c(0,1),ylim=c(0,1))

#######################################################################
###create figures of supplementary material############################
#######################################################################

#set parameter instance
p = 1000
n = 400
graph = "random"
density = "dense"  #("sparse" or "dense")

#make list of algorithms to compare
jump_speed_vec = c(1,0.001,0.01,0.05,0.25,0.75,0.25,0.25)
max_jump_vec = c(1,1,1,0.01,0.05,0.05,0.0025,0.0025)
dimin_vec = c(0,0,0,0,0,0,1,2)
algorithm_vec = c("MPLBD",rep("MPLRJ_jump",7))
df =data.frame(
  algorithm = algorithm_vec,
  dimin = dimin_vec,
  jump_speed = jump_speed_vec,
  max_jump = max_jump_vec
)
df
rep_vec = 1:8

#set final plot parameters
col = 0
par(mar = c(5, 4, 2, 2)) #these are the default margins, alternatively "par(mar = c(5, 4, 2, 10))" to make margins to create space for legend

#select what metric we want on the y-axis, only one can be set to TRUE
AUC_on = FALSE
AUC_PR_on = FALSE
p_plus_on = FALSE
p_min_on = TRUE
total_edges_on = FALSE

#make empty plot
if (AUC_on){plot(NA,xlab="time (seconds)",ylab="AUC ROC",xlim=c(0,1000),ylim=c(0.5,1))}
if (AUC_PR_on){plot(NA,xlab="time (seconds)",ylab="AUC PR",xlim=c(0,1000),ylim=c(0,1))}
if (p_plus_on){plot(NA,xlab="time (seconds)",ylab="p+",xlim=c(0,1000),ylim=c(0,1))}
if (p_min_on){plot(NA,xlab="time (seconds)",ylab="p-",xlim=c(0,15000),ylim=c(0,0.02))}
if (total_edges_on){plot(NA,xlab="time (seconds)",ylab="#edges",xlim=c(0,750),ylim=c(0,10000))}

#fill plot
for (i in 1:length(jump_speed_vec)){
  
  #set counters
  count = 0
  col = col + 1
  
  #set max_jump
  max_jump = max_jump_vec[i]
  jump_speed = jump_speed_vec[i]
  diminishing = dimin_vec[i]
  algorithm = algorithm_vec[i]
  
  #determine length vectors
  rep_test = 1
  filename = paste0("zz_results_p",p,"_n",n,"_",graph,"_",density,"_",algorithm,"_speed",jump_speed,"_maxjump",max_jump,"_dimin",diminishing,"_rep",rep_test,".Rdata")
  load(file=filename)
  len = length(output$iter_vec_thin)
  sum_time = rep(0,len)
  sum_AUC = rep(0,len)
  sum_AUC_PR = rep(0,len)
  sum_pplus = rep(0,len)
  sum_pmin = rep(0,len)
  sum_total_edges = rep(0,len)
  
  rep_vec = 1:8
  for (rep in rep_vec){
    print(rep)
    filename = paste0("zz_results_p",p,"_n",n,"_",graph,"_",density,"_",algorithm,"_speed",jump_speed,"_maxjump",max_jump,"_dimin",diminishing,"_rep",rep,".Rdata")
    if (file.exists(file=filename)){
      count = count +1
      load(file=filename)
      sum_time = sum_time + output$time_vec
      sum_AUC = sum_AUC + output$AUC_vec
      sum_AUC_PR = sum_AUC_PR + output$AUC_PR_vec
      sum_total_edges = sum_total_edges + output$total_edges_vec
      sum_pplus = sum_pplus + output$pplus_vec
      sum_pmin = sum_pmin + output$pmin_vec
    }
    else {
      print_text = paste0("Rep",rep,"does not exist")
      print(print_text)
    }
  }
  if (AUC_on){points(x=sum_time/count,y=sum_AUC/count,type="l",col=col,lw=2)}
  if (AUC_PR_on){points(x=sum_time/count,y=sum_AUC_PR/count,type="l",col=col,lw=2)}
  if (total_edges_on){points(x=sum_time/count,y=sum_total_edges/count,type="l",col=col,lw=2)}
  if (p_plus_on){points(x=sum_time/count,y=sum_pplus/count,type="l",col=col,lw=2)}
  if (p_min_on){points(x=sum_time/count,y=sum_pmin/count,type="l",col=col,lw=2)}
}

legend_vec = c("BD-MPL",expression(epsilon[s]==0.001),
                 expression(epsilon[s]==0.01),
                 expression(epsilon[s]==0.05),
                 expression(epsilon[s]==0.25),
                 expression(epsilon[s]==0.75),
                 expression(epsilon[s]== "slow decay"),
                 expression(epsilon[s]== "quick decay"))
len = length(legend_vec)
legend( "bottomright", inset=c(-0.31,0),xpd=TRUE, legend_vec, col = 1:len,lw = rep(2,len),lty = rep(1,len), cex = 1, box.lty = 0 )



