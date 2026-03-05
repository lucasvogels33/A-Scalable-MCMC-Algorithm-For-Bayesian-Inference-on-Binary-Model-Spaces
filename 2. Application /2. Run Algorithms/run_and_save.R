
#################################################
#####download libraries and supporting files#####
#################################################

#set parameters
args = commandArgs(trailingOnly = TRUE)
algorithm = args[1]
jump_speed = as.numeric(args[2])
diminishing = as.numeric(args[3])
max_jump = as.numeric(args[4])

#download libraries and functions
library(BDgraph) #used by both MJ-MCMC and BD-MPL
source("1.MPLRJ_jump_functions.R") #contains MJ-MCMC function
source("2.MPL_BD_functions.R") #contains BD-MPL function


#######################################
#####load data  ########################
######################################
filename=paste0("cleaned_data.Rdata") #the file "cleaned_data.Rdata" is created in the previous step. Make sure this .Rdata file is saved in the same folder as this R file.
load(file=filename)
data = output$data 


#######################################
#####set iter_vec_thin  ##################
######################################


if ((jump_speed == 0.0001)*(diminishing==0)){iter_vec_thin = c(1,5,10,15,20,25,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500,750,1000,1500,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,150000,200000,250000,300000,400000,500000,600000,700000,800000,900000,1000000,1500000,2000000,2500000,3000000,3500000,4000000)}
if ((jump_speed == 0.001)*(diminishing==0)){iter_vec_thin = c(1,5,10,15,20,25,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500,750,1000,1500,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,150000,200000,250000,300000,400000,500000,600000,700000,800000,900000,1000000)}
if ((jump_speed == 0.01)*(diminishing==0)){iter_vec_thin = c(1,5,10,15,20,25,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500,750,1000,1500,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000)}
if ((jump_speed == 0.3)*(diminishing==0)){iter_vec_thin = c(1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,60,70,80,90,100,110,120,130,140,150,200,250,300,350,400,450,500,550,600,650,700,750,800,900,1000,1100,1200,1300,1400,1500,1750,2000)}
if ((jump_speed == 0.3)*(diminishing==1)){iter_vec_thin = c(1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500,650,700,750,800,900,1000,1250,1500,1750,2000,2500,3000,3500,4000)}
if ((jump_speed == 0.3)*(diminishing==2)){iter_vec_thin = c(1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,100,150,200,250,300,350,400,450,500,750,1000,2000,5000,10000,15000,20000,25000,30000,35000,40000,50000,60000,70000,80000,90000,100000,125000,150000,175000,200000)}
if ((jump_speed == 0.6)*(diminishing==0)){iter_vec_thin = c(1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500,550,600,650,700,800,900,1000,1250,1500,1750,2000)}
if ((jump_speed == 0.9)*(diminishing==0)){iter_vec_thin = c(1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500)}
if (algorithm == "MPLBD"){iter_vec_thin = c(1,5,10,15,20,25,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500,750,1000,1500,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,150000,200000,250000,300000,400000,500000,600000,700000,800000,900000,1000000,1500000,2000000,2500000,3000000,3500000,4000000)}


#######################################
#####run algorithm  ##################
######################################
set.seed(1)
if (algorithm == "MPLRJ_jump"){
  for (G_start in c(0)){
    output = MPLRJ_jump_solve(data=data,iter_vec_thin=iter_vec_thin,jump_speed=jump_speed,max_jump=max_jump,diminishing=diminishing,G_start=G_start)
    filename = paste0("zz_results_MPLRJ_jump_speed",jump_speed,"_maxjump",max_jump,"_dimin",diminishing,"_Gstart",G_start,".Rdata")
    save(output,file = filename)
  }
}
set.seed(1)
if (algorithm == "MPLBD"){
    output = MPLBD_solve(data=data,iter_vec_thin=iter_vec_thin)
    filename = paste0("zz_results_MPLBD.Rdata")
    save(output,file = filename)
}

















