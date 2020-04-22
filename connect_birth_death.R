graphics.off()

library("TDA")

# Set directory to work in
setwd("C:/Users/veron/Dropbox/TDA Rings paper/Revision/TDA sampling/standard_sim1_density") # Windows

samp_bead_density <- 0.1 # dataset corresponds to locations of 10% of the beads along filaments
L = 2000       # range of the coordinate of the points in the data set (in nm)
maxBetti <- 1;
rbw_plot <- 1; # index for plotting all paths 

filenamex<-paste("samplings/bead_dataxv_samp",as.character(samp_bead_density),".csv",sep='')
filenamey<-paste("samplings/bead_datayv_samp",as.character(samp_bead_density),".csv",sep='')
filenamez<-paste("samplings/bead_datazv_samp",as.character(samp_bead_density),".csv",sep='')

# Load time data: matrix with each column corresponding to coordinates at a specific time
allx <- read.csv(filenamex,
                 header = FALSE, 
                 quote="\"", 
                 stringsAsFactors= TRUE, 
                 strip.white = TRUE)

ally <- read.csv(filenamey, 
                 header = FALSE, 
                 quote="\"", 
                 stringsAsFactors= TRUE, 
                 strip.white = TRUE)

allz <- read.csv(filenamez, 
                 header = FALSE, 
                 quote="\"", 
                 stringsAsFactors= TRUE, 
                 strip.white = TRUE)

# Determine number of time frames in the dataset
nr_times = dim(allx)[2]

# Function for computing the persistent homology of each point cloud
m_births_deaths <- function(tb) {
  print(tb)
  matxB <- allx[tb] 
  matyB <- ally[tb]
  matzB <- allz[tb]
  # Combine coordinates into a matrix and calculate the persistence diagram of the Rips filtration built on top of the point cloud.
  dataB = cbind(matxB,matyB,matzB)
  DiagB = ripsDiag(X = dataB, maxdimension = maxBetti, maxscale = L, library = "GUDHI", printProgress = FALSE)
  # Extract matrix of birth and death values from diagram object
  mB<-matrix(DiagB[[1]])
  m_birth <- DiagB[["diagram"]][,2];
  m_death <- DiagB[["diagram"]][,3];
  if (maxBetti == 1) {
    feat0 = sum(DiagB[["diagram"]][,1]==0);
    feat1 = sum(DiagB[["diagram"]][,1]==1);
    
    m_birthB <- m_birth[(feat0+1):(feat0+feat1)];
    m_deathB <- m_death[(feat0+1):(feat0+feat1)];
  } else {
    feat0 = sum(DiagB[["diagram"]][,1]==0);
    feat1 = sum(DiagB[["diagram"]][,1]==1);
    feat2 = sum(DiagB[["diagram"]][,1]==2);
    
    m_birthB <- m_birth[(feat0+feat1+1):(feat0+feat1+feat2)];
    m_deathB <- m_death[(feat0+feat1+1):(feat0+feat1+feat2)];
  }
  
  list(births = m_birthB, deaths = m_deathB)
}

# Function for connecting consecutve times
connection_btwn_two <- function(tsA,tsB,ptsA,ptsB){
  # Making matrix with distance between all the points in each time step
  num_ptsA = length(ptsA)/2
  num_ptsB = length(ptsB)/2
  distance_mat = matrix(, nrow = num_ptsA, ncol = num_ptsB)  # Sets up empty matrix for loops
  
  for(i in 1:num_ptsA){
    for(j in 1:num_ptsB){
      xA = ptsA[i,1]
      yA = ptsA[i,2]
      
      xB = ptsB[j,1]
      yB = ptsB[j,2]
      
      dist_temp = ((xB-xA)^2 + (yB-yA)^2)^(1/2)
      distance_mat[i,j] = dist_temp
    }
  }
  
  ##Finding the closest points (minimum distance)
  pair_list = c() 
  cont_diff = 370 #Continuity difference, maximum distance at which we should still attach two triangles
  
  if(num_ptsA < num_ptsB){
    for(i in 1:num_ptsA){
      min_dist = min(distance_mat[i,], na.rm = TRUE)
      if(min_dist < cont_diff){
        closest_B = which.min(distance_mat[i,])
        distance_mat[i,] <- NA
        distance_mat[,closest_B] <- NA
        pair = cbind(i,closest_B)
        pair_list <- c(pair_list, pair)
      }
      else{
        distance_mat[i,] <- NA
        pair = cbind(i,NA)
        pair_list <- c(pair_list,pair)
      }
    }
    pairs_mat = matrix(pair_list, nrow = (length(pair_list))/2, ncol = 2, byrow = TRUE)
    
    for(j in 1:num_ptsB){
      if(!(j %in% pairs_mat[,2])){
        pair = cbind(NA,j)
        pair_list <- c(pair_list, pair)
      }
    }
    pairs_mat = matrix(pair_list, nrow = (length(pair_list))/2, ncol = 2, byrow = TRUE)  
  }
  
  if(num_ptsB < num_ptsA){
    for(j in 1:num_ptsB){
      min_dist = min(distance_mat[,j], na.rm = TRUE)
      if(min_dist < cont_diff){
        closest_A = which.min(distance_mat[,j])
        distance_mat[,j] <- NA
        distance_mat[closest_A,] <- NA
        pair = cbind(closest_A, j)
        pair_list <- c(pair_list, pair)
      }
      else{
        distance_mat[,j] <- NA
        pair = cbind(NA,j)
        pair_list <- c(pair_list,pair)
      }
    }
    pairs_mat = matrix(pair_list, nrow = (length(pair_list))/2, ncol = 2, byrow = TRUE)
    
    for(i in 1:num_ptsA){
      if(!(i %in% pairs_mat[,1])){
        pair = cbind(i,NA)
        pair_list <- c(pair_list, pair)
      }
    }
    pairs_mat = matrix(pair_list, nrow = (length(pair_list))/2, ncol = 2, byrow = TRUE)
  }
  
  if(num_ptsA == num_ptsB){
    for(j in 1:num_ptsB){
      min_dist = min(distance_mat[,j], na.rm = TRUE)
      if(min_dist < cont_diff){
        closest_A = which.min(distance_mat[,j])
        distance_mat[,j] <- NA
        distance_mat[closest_A,] <- NA
        pair = cbind(closest_A, j)
        pair_list <- c(pair_list, pair)
      }
      else{
        distance_mat[,j] <- NA
        pair = cbind(NA,j)
        pair_list <- c(pair_list,pair)
      }
    }
    pairs_mat = matrix(pair_list, nrow = (length(pair_list))/2, ncol = 2, byrow = TRUE)
  }
  return(pairs_mat)
}


#### Creating path matrices ###########################################

# Set up list structure for all birth-death pairs
desired_length <- 202 # or whatever length you want
empty_list <- vector(mode = "list", length = desired_length);
myList <- str(empty_list);

## Get birth-death pair points for first time step 
myList[[1]] <- m_births_deaths(1)
triangle_pointsA = cbind(myList[[1]]$births,myList[[1]]$deaths)

## Get birth-death pair points for second time step 
myList[[2]] <- m_births_deaths(2)
triangle_pointsB = cbind(myList[[2]]$births,myList[[2]]$deaths)

## Set up connection matrix
ab = connection_btwn_two(1,2,triangle_pointsA,triangle_pointsB)
l_mat = ab

x_pairs = c()
y_pairs = c()
for(k in 1:length(ab[,1])){ 
  cA = ab[k,1]
  cB = ab[k,2]
  ptA = triangle_pointsA[cA,]
  ptB = triangle_pointsB[cB,]
  x_pair = cbind(ptA[1],ptB[1])
  y_pair = cbind(ptA[2],ptB[2])
  x_pairs = cbind(x_pairs,x_pair)
  y_pairs = cbind(y_pairs,y_pair)
}
all_x_conn = matrix(x_pairs,nrow = length(ab[,1]), ncol = 2, byrow=TRUE)
all_y_conn = matrix(y_pairs,nrow = length(ab[,1]), ncol = 2, byrow = TRUE)



### All subsequent times - use lapply to find all birth and death pairs
trials <- 2:nr_times
myList[trials] <- lapply(trials, m_births_deaths)
    
for(tb in 2:(nr_times-1)){ # loop through time  
    ## Extract birth-death pair points for first time step 
    triangle_pointsB = cbind(myList[[tb]]$births,myList[[tb]]$deaths)
    
    ## Extract birth-death pair points for second time step 
    triangle_pointsC = cbind(myList[[tb+1]]$births,myList[[tb+1]]$deaths)
    
    ## Connections to be added
    bc = connection_btwn_two(tb,tb+1,triangle_pointsB,triangle_pointsC)
    
    ## L matrix with number labels for each triangle
    # Sets up matrix for connections, indirectly handles death of orphans
    n_emp_row = cbind(rep(NA,nrow(l_mat)))
    l_mat = cbind(l_mat,n_emp_row)
    all_x_conn = cbind(all_x_conn,n_emp_row)
    all_y_conn = cbind(all_y_conn,n_emp_row)
    
    for(p in 1:length(bc[,1])){
      cB = bc[p,1]
      cC = bc[p,2]
      ptB = triangle_pointsB[cB,]
      ptC = triangle_pointsC[cC,]
      
      if(!is.na(bc[p,1])){  # Matches up all pairs
        if(bc[p,1] %in% l_mat[,ncol(l_mat)-1]){
          pos = which(l_mat[,ncol(l_mat)-1]==bc[p,1])
          l_mat[pos,ncol(l_mat)] = bc[p,2]
          all_x_conn[pos,ncol(all_x_conn)] = ptC[1]
          all_y_conn[pos,ncol(all_y_conn)] = ptC[2]
          
        }
      }
      else{  # Handles birth of orphans
        orphan_birth = c(rep(NA,ncol(l_mat)-1),bc[p,2])
        orphan_birth = matrix(orphan_birth,ncol = length(orphan_birth))
        l_mat = rbind(l_mat,orphan_birth)
        
        orphan_birth_x = c(rep(NA,ncol(all_x_conn)-1),ptC[1])
        orphan_birth_x = matrix(orphan_birth_x,ncol = length(orphan_birth_x))
        all_x_conn = rbind(all_x_conn,orphan_birth_x)
        
        orphan_birth_y = c(rep(NA,ncol(all_y_conn)-1),ptC[2])
        orphan_birth_y = matrix(orphan_birth_y,ncol = length(orphan_birth_y))
        all_y_conn = rbind(all_y_conn,orphan_birth_y)
      }
    }
  }

#### Plot Rainbow Paths #######################################
if (rbw_plot == 1){
  # Make a blank plot
  dev.new()
  dev.new()
  print("Plotting paths in rainbow colors")
  test_mat = matrix(0,nrow=1, ncol=1) 
  testDiag = ripsDiag(X = test_mat, maxdimension = maxBetti, maxscale = L, library = "GUDHI", printProgress = FALSE)
  plot(testDiag[["diagram"]], main = "Time-connected paths")
  
  #Plot the paths
  for(path_num in 1:(length(all_x_conn[,1]))){
    #   print(path_num)
    for(pts in 1:(length(all_x_conn[path_num,])-1)){
      segments(all_x_conn[path_num,pts],all_y_conn[path_num,pts],all_x_conn[path_num,pts+1],all_y_conn[path_num,pts+1],col = rainbow(length(all_x_conn[,1]))[path_num] )
    }
  }
}

#### Isolate Significant Path #############################
# Determine the significant path
sig_index_vec = c()
max_dist_vec = c()

for (col in 1:ncol(l_mat)){
  
  curr_mat  = c()
  empty_col = c()
  
  curr_mat  = cbind(l_mat[,col],all_x_conn[,col],all_y_conn[,col])
  empty_col = cbind(rep(NA,nrow(curr_mat)))
  curr_mat <- cbind(curr_mat, empty_col)
  
  for(i in 1:nrow(curr_mat)){
    if(!is.na(curr_mat[i,1])){
      curr_x = curr_mat[i,2]
      curr_y = curr_mat[i,3]
      curr_dist = curr_y - curr_x
      curr_mat[i,4] = curr_dist
    }
  }
  
  # Finding max distance
  max_dist_vec[col] = max(curr_mat[,4], na.rm = TRUE)
  sig_index_vec[col] = which.max(curr_mat[,4])
}

max_dist = max(max_dist_vec, na.rm = TRUE)
max_ind = which.max(max_dist_vec)
sig_path_index = sig_index_vec[max_ind]


## Plotting the significant path in persistence diagram space
dev.new()
# Make a blank plot
test_mat = matrix(0,nrow=1, ncol=1) 
testDiag = ripsDiag(X = test_mat, maxdimension = maxBetti, maxscale = L, library = "GUDHI", printProgress = FALSE)
plot(testDiag[["diagram"]], main = "Significant path")

# Plot the significant path
for(path_num in sig_path_index){
  print(path_num)
  for(pts in 1:(length(all_x_conn[path_num,])-1)){
    segments(all_x_conn[path_num,pts],all_y_conn[path_num,pts],all_x_conn[path_num,pts+1],all_y_conn[path_num,pts+1],col = "red" )
  }
}


### Plotting Distance vs Time
all_sig_pts_dist  = c()
all_sig_pts_dist1 = c()
all_sig_times = c()
all_diffs = c()

# Significance threshold above the diagonal
gt = 301.1

for(path_num in sig_path_index){
  print(path_num)
  for(pts in 1:(length(all_x_conn[path_num,])-1)){
    curr_sig_x = all_x_conn[path_num,pts]
    curr_sig_y = all_y_conn[path_num,pts]
    curr_sig_dist = curr_sig_y-curr_sig_x
    all_sig_pts_dist1 = c(all_sig_pts_dist1,curr_sig_y-curr_sig_x)
    all_sig_times = c(all_sig_times, pts)
    diff = curr_sig_dist - gt
    all_diffs = c(all_diffs, diff)
  }
}


# Find timing of ring formation
last_neg_diff = max(which(all_diffs<0))
first_sig_time_index = last_neg_diff + 1

dev.new()
plot(all_sig_times*10,all_sig_pts_dist1,pch=20,ylim = c(0,1500),xlab = "Time",ylab = "Distance from Baseline",col="blue", main = "Significant path through time")
lines(c(0,2000),c(gt,gt),lty = "dashed")
if(first_sig_time_index*10-10<2000){
  lines(c(first_sig_time_index*10-10,first_sig_time_index*10-10),c(0,2000),lty="dashed")
}

print(paste0("Time of ring formation: ", first_sig_time_index*10-10))



