

####
#this Rscript reads the .Rdata files created in the step "2. Run algorithms" 
#and creates the figures in the section 5.2 of the article


##########################################
####test if iterations are as I set them
##########################################

jump_speed_vec = c("1e-04",0.001,0.01,0.3,0.3,0.3,0.6,0.9,1)
dimin_vec = c(0,0,0,0,1,2,0,0,0)
algorithm_vec = c(rep("MPLRJ_jump",8),"MPLBD")
max_jump_vec = c(rep(0.025,7),0.1,1)
G_start = 0

for (i in 1:length(jump_speed_vec)){
  jump_speed = jump_speed_vec[i]
  dimin = dimin_vec[i]
  algorithm = algorithm_vec[i]
  max_jump = max_jump_vec[i]
  if (algorithm == "MPLRJ_jump"){filename = paste0("zz_results_",algorithm,"_speed",jump_speed,"_maxjump",max_jump,"_dimin",dimin,"_Gstart",G_start,".Rdata")}
  if (algorithm == "MPLBD"){filename = paste0("zz_results_MPLBD.Rdata")}
  load(file=filename)
  
  iter_vec_thin = output$iter_vec_thin
  last_iter = iter_vec_thin[length(iter_vec_thin)]
  print(algorithm)
  print(jump_speed)
  print(dimin)
  print(last_iter)
}

#########################################
############edges flipped###############
#########################################
#this code verifies that algorithms with a higher jump_speed flip more edges at every single iteration

jump_speed_vec = c("1e-04",0.001,0.01,0.3,0.6,0.9,0.3,0.3)
dimin_vec = c(0,0,0,0,0,0,1,2)
algorithm_vec = rep("MPLRJ_jump",8)
max_jump_vec = c(0.025,0.025,0.025,0.025,0.025,0.1,0.025,0.025)
G_start = 0

plot_until = 500
col = 0
#par(mar = c(5, 4, 4, 2)) #default settings
par(mar = c(5, 4, 4, 10)) #set margins to create space for legend
plot(NA,xlim=c(0,plot_until),ylim=c(0,2000),xlab="MCMC iterations",ylab="edges flipped")

for (i in 1:length(jump_speed_vec)){
  col = col + 1
  jump_speed = jump_speed_vec[i]
  dimin = dimin_vec[i]
  algorithm = algorithm_vec[i]
  max_jump = max_jump_vec[i]
  filename = paste0("zz_results_",algorithm,"_speed",jump_speed,"_maxjump",max_jump,"_dimin",dimin,"_Gstart",G_start,".Rdata")
  load(file=filename)
  
  edges_flipped_vec = output$size_selected_edges_vec[1:plot_until]
  points(x=1:plot_until,y=edges_flipped_vec,col=col,type="l")
}

legend_vec = c(expression(epsilon[s]==0.0001),
               expression(epsilon[s]==0.001),
               expression(epsilon[s]==0.01),
               expression(epsilon[s]==0.3),
               expression(epsilon[s]==0.6),
               expression(epsilon[s]==0.9),
               expression(epsilon[s]== "slow decay"),
               expression(epsilon[s]== "quick decay"))
legend( "bottomright", inset=c(-0.6,0),xpd=TRUE, legend_vec, col = 1:col,lw = rep(2,col),lty = rep(1,col), cex = 1, box.lty = 0 )

##########################################
####plot correlation plots###############
##########################################
#this code verifies that the MJ-MCMC algorithm becomes 'similar' to the BD-MPL algorithm.
#it does so by plotting a measure of distance (MSE, correlation, absolute difference) over time. 
#these distance measures should all grow close to zero as the running time increases
#Note that these plots do not appear in the paper

#download MPL-BD final prediction
filename ="zz_results_MPLBD_seed1.Rdata"
load(file=filename)
MPLBD_pred = output$predictor
cor_upper = 1.1
cor_lower = -0.01
pred_index = which(MPLBD_pred>cor_lower& MPLBD_pred<cor_upper) #used if we want to exclude certain edge incl. prob.
MPLBD_pred = MPLBD_pred[pred_index]
length(MPLBD_pred)/(623*622/2)

#what algorithms do we include
jump_speed_vec = c(1,0.001,0.01,0.3,0.6,0.9,0.3,0.3)
dimin_vec = c(rep(0,6),1,2)
algorithm_vec = c("MPLBD",rep("MPLRJ_jump",7))
max_jump_vec = c(1,rep(0.025,4),0.1,0.025,0.025)
G_start = 0

#set colour
col = 0

#set parameters x-axis
time_on = TRUE  #set to TRUE for time on the x-axis
iter_on = FALSE #set to TRUE for iterations on the x-axis
if (time_on){
  xlab = "time (seconds)"
  xlim = c(0,100)
}
if (iter_on){
  xlab = "MCMC iterations"
  xlim = c(0,1000)
}

#set parameters y-axis
MSE_on = FALSE
cor_on = TRUE
diff_on = FALSE
cutoff = 0.02 #convergence is defined as the moment the distance measure (MSE, correlation, or abs difference) is smaller than "cutoff"
max_on = FALSE
if (MSE_on){
  ylab = "MSE with MPL-BD"
  ylim=c(0,0.05)
}
if (cor_on){
  ylab = "Corr distance with MPL-BD"
  ylim = c(0,1)
}
if (max_on){
  ylab = "Max diff with MPL-BD"
  ylim = c(0,1)
}
if (diff_on){
  ylab = "Average abs dist with MPL-BD"
  ylim = c(0,0.05)
}

#plot 
#par(mar = c(5, 4, 4, 2)) #default settings
par(mar = c(5, 4, 4, 10)) #set margins to create space for legend
plot(NA,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim)
runtime_vec =c()
iter_needed_vec = c()
final_dist_vec = c()


for (i in 1:length(jump_speed_vec)){
  
  col = col + 1
  
  #find parameters
  jump_speed = jump_speed_vec[i]
  dimin = dimin_vec[i]
  algorithm = algorithm_vec[i]
  max_jump = max_jump_vec[i]
  
  #initialize output
  dist_vec = c()
  
  #download iter_vec_thin length
  if (algorithm == "MPLRJ_jump"){
    filename = paste0("zz_results_",algorithm,"_speed",jump_speed,"_maxjump",max_jump,"_dimin",dimin,"_Gstart",G_start,".Rdata")
    load(file=filename)
    iter_vec_thin = output$iter_vec_thin[2:length(output$iter_vec_thin)]  
    time_vec = output$time_vec[2:length(output$time_vec)] 
  }
  if (algorithm == "MPLBD"){
    filename = paste0("zz_results_MPLBD_seed2.Rdata")
    load(file=filename)
    iter_vec_thin = output$iter_vec_thin
    time_vec = output$time_vec
  }
  len = length(iter_vec_thin)
  
  #fill the correlation vector
  for (iter in 1:len){
    if (algorithm == "MPLRJ_jump"){filename = paste0("yyz_results_",algorithm,"_speed",jump_speed,"_maxjump",max_jump,"_dimin",dimin,"_Gstart",G_start,"_iter",iter,".Rdata")}
    if (algorithm == "MPLBD"){filename = paste0("yyz_results_",algorithm,"_iter",iter,"_seed2.Rdata")}
    load(file=filename)
    
    predictor = output_per_iter$predictor[pred_index]
    if (cor_on){dist_vec = c(dist_vec,1 - cor(MPLBD_pred,predictor))}
    if (MSE_on){dist_vec = c(dist_vec,mean((MPLBD_pred-predictor)^2))}
    if (diff_on){dist_vec = c(dist_vec,mean(abs(MPLBD_pred-predictor)))}
    if (max_on){dist_vec = c(dist_vec,max(abs(MPLBD_pred-predictor)))}
    
  }
  
  #save runtime and iterations needed for "convergence"
  indexes = which(dist_vec<cutoff)
  runtime_vec = c(runtime_vec,min(time_vec[indexes]))
  iter_needed_vec = c(iter_needed_vec,min(iter_vec_thin[indexes]))
  final_dist_vec = c(final_dist_vec,dist_vec[length(dist_vec)])
  
  #plot 
  if (time_on){x_vec = time_vec}
  if (iter_on){x_vec = iter_vec_thin}
  points(x=x_vec,y=dist_vec,type="l",col=col,lw=2)
}

#create legend for convergence plot
legend_vec = c("BD-MPL",
               expression(epsilon[s]==0.001),
               expression(epsilon[s]==0.01),
               expression(epsilon[s]==0.3),
               expression(epsilon[s]==0.6),
               expression(epsilon[s]==0.9),
               expression(epsilon[s]== "slow decay"),
               expression(epsilon[s]== "quick decay"))
legend( "bottomright", inset=c(-0.4,0),xpd=TRUE, legend_vec, col = 1:col,lw = rep(2,col),lty = rep(1,col), cex = 1, box.lty = 0 )

##################################
#plot running time barchart#######
##################################
#this creates the barplot for running time until "convergence" in the article

runtime_vec = round(runtime_vec/60,1) #the vector runtime_vec comes from the code right above 
runtime_vec[6] = 0 #for jump_speed 0.9 we we have NA
values = runtime_vec
par(mar = c(5, 8, 1, 4))  # bottom, left, top, right

#make the barchart
bp= barplot(
  rev(values),
  horiz = TRUE,
  names.arg = rev(expression("BD-MPL", epsilon[s]*" = 0.001", epsilon[s]*" = 0.01", epsilon[s]*" = 0.3",epsilon[s]*" = 0.6",epsilon[s]*" = 0.9",
                             epsilon[s]*"= slow decay",epsilon[s]*" = quick decay")),
  las = 1,
  xlab = "running time (minutes)",
  col = "lightblue",
  border = "lightblue"
)

#write the text values in the plots
text_values = runtime_vec
text_values[6] = "NA" #for jump_speed 0.9 we we have NA
text(
  x = rev(values),
  y = bp,
  labels = rev(text_values),
  pos = 4,      # to the right of the bar
  xpd = TRUE,
  font = 3
)


##############################################
#plot MCMC iterations until convergence#######
##############################################
#same as the above code, but for MCMC iterations

iter_needed_vec[6] = 0 #for jump_speed 0.9 we we have NA
values = iter_needed_vec
par(mar = c(5, 8, 1, 4))  # bottom, left, top, right

bp= barplot(
  rev(values),
  horiz = TRUE,
  names.arg = rev(expression("BD-MPL", epsilon[s]*" = 0.001", epsilon[s]*" = 0.01", epsilon[s]*" = 0.3",epsilon[s]*" = 0.6",epsilon[s]*" = 0.9",
                             epsilon[s]*"= slow decay",epsilon[s]*" = quick decay")),
  las = 1,
  xlab = "MCMC iterations",
  col = "lightblue",
  border = "lightblue"
)

text_values = iter_needed_vec
text_values = format(text_values, big.mark = ",", scientific = FALSE)
text_values[6] = "NA" #for jump_speed 0.9 we we have NA
text(
  x = rev(values),
  y = bp,
  labels = rev(text_values),
  pos = 4,      # to the right of the bar
  xpd = TRUE,
  font = 3
)

###############################################################
##############distances between vectors#######################
############################################################

#set the parameters of all algorithm that ran
jump_speed_vec = c(0.001,0.01,0.3,0.3,0.3,0.6,0.9)
dimin_vec = c(0,0,0,1,2,0,0)
algorithm_vec = rep("MPLRJ_jump",7)
max_jump_vec = c(rep(0.025,6),0.1)

#load the edge incl prob of the MPL-BD vector
filename = "zz_results_MPLBD_seed1.Rdata"
load(file=filename)
mplbd_pred = output$predictor

#1. calculate the distance (correlation, mean difference, max_difference) at each iteration
#2. Make the scatterplot
cor_vec = c()
mean_diff_vec = c()
max_diff_vec= c()
plot_this = TRUE
for (i in c(3)){ #selected c(1),c(2),...,c(7) depending on the MJ-MCMC variant you want to compare to MPL-BD
  
  #load data
  jump_speed = jump_speed_vec[i]
  dimin = dimin_vec[i]
  algorithm = algorithm_vec[i]
  max_jump = max_jump_vec[i]
  filename = paste0("zz_results_",algorithm,"_speed",jump_speed,"_maxjump",max_jump,"_dimin",dimin,"_Gstart0.Rdata")
  load(file=filename)
  mplrj_pred = output$predictor
  
  #calculate distances
  correlation = round(cor(mplrj_pred,mplbd_pred),2)
  cor_vec = c(cor_vec,correlation)
  mean_diff_vec = c(mean_diff_vec,mean(abs(mplrj_pred - mplbd_pred)))
  max_diff_vec = c(max_diff_vec,max(abs(mplrj_pred - mplbd_pred)))
  
  #make scatterplot of predictors
  if (plot_this==TRUE){
    par(mar = c(5, 5, 1, 4))  # bottom, left, top, right
    #title = paste0("jump = ",jump_speed,"; dimin = ",dimin,"; cor = ",correlation)
    plot(x=mplbd_pred,y=mplrj_pred,ylab = expression("edge incl. prob — "*epsilon[s]*" = 0.3"),xlab = "edge incl. prob — BD-MPL",col = rgb(0, 0, 0, alpha = 0.06),cex = 0.2,asp=1,pch=16,xlim=c(0,1),ylim=c(0,1))
    segments(
      x0 = 0, y0 = 0,
      x1 = 1, y1 = 1,
      col = "red",
      lwd = 1
    )
  }
}







