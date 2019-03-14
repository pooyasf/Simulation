# Simplex Generator ----
s <- 10  ##length of random unit simplex (sum=1)
s_v <- c(0)
for (si in 1:(s-1)){
  s_v <- cbind(s_v,runif(1,0.0,1.0))
}
s_v <- cbind(s_v,1)
s_v <- sort(s_v,decreasing = FALSE)
z_v <- c()
for (si in 1:s){
  z_v <- cbind(z_v,(s_v[si+1]-s_v[si]))
}


# Calculate Q -----

P <- z_v#z_v  ## probability vector
n_P <- length(P)
Q <- matrix(0,nrow = n_P-1 , ncol = n_P)
P_temp <- matrix(c(P,integer(n_P*n_P-n_P)),nrow = n_P , ncol = n_P , byrow = T)

n <- n_P
for (ln in 1:(n_P-1)){

  n <- n_P-ln+1
  i <- which( replace(P_temp[ln,] , P_temp[ln,] == 0, 1) == min(replace(P_temp[ln,] , P_temp[ln,] == 0, 1)) , arr.ind = T)
  j <- which( P_temp[ln,] == max(P_temp[ln,]) , arr.ind = T)
  if (length(i)>1){
    i <- i[1]
  }
  if (length(j)>1){
    j <- j[2]
  }
  Q[ln,i] <- (n-1)*P_temp[ln,i]
  Q[ln,j] <- 1 - (n-1)*P_temp[ln,i]
  P_temp[ln+1,i] <- 0
  P_temp[ln+1,j] <- ((n-1)/(n-2))*(P_temp[ln,i]+P_temp[ln,j] - (1/(n-1)))
  e_loop <- setdiff((1:n_P), c(i,j))
  for(ii in e_loop){
    P_temp[ln+1,ii] <- ((n-1)/(n-2))*P_temp[ln,ii]
  }
  
}


# Functions for finding i and j position and its values -----

q_i <- function(m){
  i <- which( ( (m < max(m)) & (m > 0)) , arr.ind = TRUE)
  if(length(i) == 0){
    i <- which( m == max(m) , arr.ind = TRUE)
    i <- i[1]
  }
  return(c(m[i],i))
}

q_j <- function(m){
  i <- which( m == max(m) , arr.ind = TRUE)
  if(length(i) > 1){
    i <- i[2]
  }
  return(c(m[i],i))
}

## sim loop ----
h <- c()
nsim <- 10000
for (i in 1:nsim){
  u1 <- runif(1,0.0,1.0)
  
  N <- 1 + floor((n_P-1)*u1)
  u2 <- runif(1,0.0,1.0)
  
  q_i_v <- q_i(Q[N,])
  q_j_v <- q_j(Q[N,])
  
  if ( u2 <  q_i_v[1]){
    x <- q_i_v[2]
  } else {
    x <- q_j_v[2]
  }
  h <- cbind(h,x)
}
# Plot side by side
par(mfrow = c(1, 2))
hist(h,breaks = c(0:n_P),freq=FALSE,main = "Histogram")
plot(1:s,z_v,type = 'h',main = "Orginal Probability Vector")
#mtext("My 'Title' in a strange place", side = 3, line = -2, outer = TRUE)
