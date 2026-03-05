

############################
###set parameters ##########
############################

#--read the arguments --- #
p = 1000
n = 400
graph = "random" 
density = "dense" #either "sparse" or "dense"
algorithm = "MPLRJ_jump" #either "MPLRJ_jump" (our MJ-MCMC algorithm) or "MPLBD" (for the MPL-BD algorithm)
jump_speed = 0.25 #value between 0 and 1
max_jump = 0.05 #value between 0 and 1
diminishing = 2 #0 for homogeneous chain, 1 for slow decay, 2 for quick decay
rep = 5

#--download the necessary libraries and load the run_experiments file --
library(PRROC) #to calculate the AUC PR
library(BDgraph) #for MJ-MCMC and MPL-BD 

# download the supporting functions
source("1.MPLRJ_jump_functions.R") #contains the code to run the MJ-MCMC algorithm
source("2.MPLBD_functions.R") #contains the code to run the MPL-BD algorithm
source("3.metric_functions.R") #contains the calculate the AUC-PR, AUC-ROC, p+, p-, and other metrics

#--set value of "size", determining the amount of edges in the true graph
if (density=="sparse"){size = 1000}
if (density=="dense"){size = 5000}

#set remaining parameters for intances creation
type = "Gaussian"
prob = 0.1 #this value does not matter but needs to be specified
vis = FALSE


###########################################################################################
#set the iter_vec_thin. In other words, at which iterations do our algorithms need to output a metric
###########################################################################################

#p=1000, n=400, sparse
if ((p==1000)*(n==400)*(density=="sparse")){
    if (algorithm=="MPLRJ_jump"){
        if (diminishing==0){
            if ((jump_speed==0.001)*(max_jump==1)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,75,100,150,200,250,300,400,500,600,750,1000,2500,5000,10000,15000,20000,25000,30000,35000,40000)}
            if ((jump_speed==0.01)*(max_jump==1)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,75,100,125,150,175,200,250,500,750,1000,1250,1500,1750,2000,2500,3000,3500,4000)}
            if ((jump_speed==0.05)*(max_jump==1)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,75,100,150,200,250,300,400,500,600,750)}
            if ((jump_speed==0.3)*(max_jump==1)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,75,100,125,150,175,200)}
            if ((jump_speed==0.7)*(max_jump==0.0025)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,60,70,80,90,100)}
            if ((jump_speed==0.95)*(max_jump==0.0025)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,60,70,80,90,100)}
        }
        if ((diminishing==1)*(jump_speed==0.3)*(max_jump==1)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,75,100,125,150,175,200,250,300,350,400,450,500)}
        if ((diminishing==2)*(jump_speed==0.3)*(max_jump==1)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,75,100,150,200,250,300,400,500,600,750,1000,2500,5000,10000,15000,20000,25000,30000,35000,40000)}
    }
    if (algorithm=="MPLBD"){iter_vec_thin = c(10,50,100,500,1000,5000,25000,50000,1e+05,150000,250000,300000)}
}

#p=1000, n=1000, sparse
if ((p==1000)*(n==1000)*(density=="sparse")){
    if (algorithm=="MPLRJ_jump"){
        if (diminishing==0){
            if ((jump_speed==0.001)*(max_jump==1)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,75,100,150,200,250,300,400,500,600,750,1000,2500,5000,10000,15000,20000,25000,30000,35000,40000)}
            if ((jump_speed==0.01)*(max_jump==1)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,75,100,125,150,175,200,250,500,750,1000,1250,1500,1750,2000,2500,3000,3500,4000)}
            if ((jump_speed==0.05)*(max_jump==1)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,75,100,150,200,250,300,400,500,600,750)}
            if ((jump_speed==0.3)*(max_jump==0.0025)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,75,100,125,150,175,200,300,400,500,600,700,800,900,1000)}
            if ((jump_speed==0.7)*(max_jump==0.0025)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,60,70,80,90,100)}
            if ((jump_speed==0.95)*(max_jump==0.0025)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,60,70,80,90,100)}
        }
        if ((diminishing==1)*(jump_speed==0.3)*(max_jump==0.0025)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,75,100,125,150,175,200,250,300,350,400,450,500)}
        if ((diminishing==2)*(jump_speed==0.3)*(max_jump==0.0025)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,75,100,150,200,250,300,400,500,600,750,1000,2500,5000,10000,15000,20000,25000,30000,35000,40000)}
    }
    if (algorithm=="MPLBD"){iter_vec_thin = c(10,50,100,500,1000,5000,25000,50000,1e+05,150000,250000,300000)}
}



#p=1000, n=400, dense
if ((p==1000)*(n==400)*(density=="dense")){
    if (algorithm=="MPLRJ_jump"){
        if (diminishing==0){
            if ((jump_speed==0.001)*(max_jump==1)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,75,100,150,200,250,500,750,1000,2000,3000,4000,5000,7500,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000)}
            if ((jump_speed==0.01)*(max_jump==1)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,75,100,150,200,250,500,750,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)}
            if ((jump_speed==0.05)*(max_jump==0.01)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,75,100,150,200,250,300,400,500)}
            if ((jump_speed==0.25)*(max_jump==0.05)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,60,70,80,90,100)}
            if ((jump_speed==0.75)*(max_jump==0.05)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,75,100,125,150,175,200)}
        }
        if ((diminishing==1)*(jump_speed==0.25)*(max_jump==0.0025)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,25,50,75,100,150,200,300,400,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000)}
        if ((diminishing==2)*(jump_speed==0.25)*(max_jump==0.0025)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,100,200,500,1000,2500,5000,7500,10000,15000,20000,25000,30000,35000,40000,50000,60000,70000,80000,90000,100000)}
    }
    if (algorithm=="MPLBD"){iter_vec_thin = c(10,50,100,500,1000,5000,10000,15000,20000,25000,35000,50000,1e+05,150000,250000,300000,400000,500000,600000,700000,800000,900000,1000000)}
}

#p=1000, n=1000, dense
if ((p==1000)*(n==1000)*(density=="dense")){
    if (algorithm=="MPLRJ_jump"){
        if (diminishing==0){
            if ((jump_speed==0.001)*(max_jump==1)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,75,100,150,200,250,500,750,1000,2000,3000,4000,5000,7500,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000)}
            if ((jump_speed==0.01)*(max_jump==1)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,75,100,150,200,250,500,750,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)}
            if ((jump_speed==0.05)*(max_jump==0.0025)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,75,100,150,200,250,500,750,1000,2000,3000,4000,5000)}
            if ((jump_speed==0.25)*(max_jump==0.01)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,60,70,80,90,100)}
            if ((jump_speed==0.75)*(max_jump==0.05)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,75,100,125,150,175,200)}
        }
        if ((diminishing==1)*(jump_speed==0.25)*(max_jump==0.0025)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,25,50,75,100,150,200,300,400,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000)}
        if ((diminishing==2)*(jump_speed==0.25)*(max_jump==0.0025)){iter_vec_thin = c(1,2,3,4,6,8,10,12,15,20,30,40,50,100,200,500,1000,2500,5000,7500,10000,15000,20000,25000,30000,35000,40000,50000,60000,70000,80000,90000,100000)}
    }
    if (algorithm=="MPLBD"){iter_vec_thin = c(10,50,100,500,1000,5000,10000,15000,20000,25000,35000,50000,1e+05,150000,250000,300000,400000,500000)}
}

####################################  
#create and save the "true" instance######
####################################

#create
set.seed(rep)
data.sim = bdgraph.sim( p = p, n = n, graph = graph, prob = prob, size = size, type = type, vis = vis )
data.sim$density = density
data.sim$rep = rep
data_output = list(G = data.sim$G,p=p,n=n,graph = graph,density=density,rep=rep)

#save
if (algorithm=="MPLBD"){ #all algorithms use the same instance, so we only save it for MPL_BD
    filename = paste0("yy_truth_p",p,"_n",n,"_",graph,"_",density,"_rep",rep,".Rdata")
    save(data_output,file = filename)
}

######################################
#run algorithms and save output ######
######################################

if (algorithm=="MPLRJ_jump"){output = MPLRJ_jump_solve(data=data.sim,iter_vec_thin=iter_vec_thin,jump_speed=jump_speed,max_jump=max_jump,diminishing=diminishing)}
if (algorithm=="MPLBD"){output = MPLBD_solve(data=data.sim,iter_vec_thin=iter_vec_thin)}
filename = paste0("zz_results_p",p,"_n",n,"_",graph,"_",density,"_",algorithm,"_speed",jump_speed,"_maxjump",max_jump,"_dimin",diminishing,"_rep",rep,".Rdata")
save(output,file = filename)













