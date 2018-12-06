mider.plot <- function(dist, con_array, T, labels){
  z_dist <- apply(dist, 2, function(col) (col - mean(col))/sd(col))
  D <- dist(z_dist, diag = TRUE, upper = TRUE)
  W <- diag(ntotal) - matrix(1/ntotal, ntotal, ntotal)
  M <- W %*% as.matrix(-0.5 * D * D) %*% W
  
  eig <- eigen((M + t(M))/2)
  fe <- eig$values
  V <- eig$vectors
  Y <- V[,fe>0] %*% diag(sqrt(fe[fe>0]))
  Y[, 1] <- -Y[, 1] 
  
  plot(Y[, 1], Y[, 2])
  
  for (i in 1:ntotal){
    for (j in 1:ntotal){
      if (con_array[i, j] > 0){
        if (T[i,j] > 0){
          arrows(Y[i,1], Y[i,2], Y[j,1], Y[j,2], length = 0.1)
        }
        if (T[i,j] < 0){
          arrows(Y[j,1], Y[j,2], Y[i,1], Y[i,2], length = 0.1)
        }
        if (T[i,j] == 0){
          segments(Y[j,1], Y[j,2], Y[i,1], Y[i,2])
        }
      }
    }
  }
  
  text(Y[, 1] - 0.1, Y[, 2]-0.1, labels=procname,cex=0.5, font=2)
}

pids <- names(logs.mem.time)[!is.na(colSums(logs.mem.time))]
procname = c("WindowServer",
             "mds_stores",
             "Finder",
             "Spotlight",
             "DropBox",
             "CloudApp",
             "reversetemplated",
             "Google Chrome",
             "Google Chrome Helper",
             "RStudio",
             "WebKit",
             "WebKit",
             "Messages",
             "WebKit")

text(Y[, 1] - 0.1, Y[, 2]-0.1, labels=procname,cex=0.5, font=2)


