estimateH4 <- function(w, x, y, z, pb, q){
  
  leng = length(x)
  binsize = ceiling(sqrt(pb*leng))   #% Estimated bin size for ~pb points/bin
  b       = floor(leng/binsize)   #% Number of bins
  binsize = floor(leng/b)         #% Integer number of elements per bin
  
  #% Initialize matrices to zero
  his = matrix(0, b, b)
  his3 <- array(0, dim = c(b, b, b))
  his4 <- array(0, dim = c(b, b, b, b))
  
  x0  = matrix(0, b, 1)
  x1  = matrix(0, b, 1)
  x2  = matrix(0, b, 1)
  x3  = matrix(0, b, 1)
  
  #% Sort values of X and Y (descending order)
  #[sx,ix] = sort(x,2,'descend');
  tmp = sort(w, index.return=TRUE, decreasing = TRUE); sw = tmp$x; iw=tmp$ix
  tmp = sort(x, index.return=TRUE, decreasing = TRUE); sx = tmp$x; ix=tmp$ix 
  tmp = sort(y, index.return=TRUE, decreasing = TRUE); sy = tmp$x; iy=tmp$ix 
  tmp = sort(z, index.return=TRUE, decreasing = TRUE); sz = tmp$x; iz=tmp$ix 
  
  #% Replace elements of X and Y by their rank order
  linsp  = 1:leng
  sw[iw] = linsp;
  sx[ix] = linsp;
  sy[iy] = linsp;
  sz[iz] = linsp;
  
  #% Put sorted elements into bins
  sw = ceiling(sw/binsize);
  sx = ceiling(sx/binsize);
  sy = ceiling(sy/binsize);
  sz = ceiling(sz/binsize);
  
  # % Put largest values in highest bin
  for (p in 1:leng){
    if (sw[p] > b) sw[p] = b
    if (sx[p] > b) sx[p] = b
    if (sy[p] > b) sy[p] = b
    if (sz[p] > b) sz[p] = b
  }
  
  #% Populate histograms
  for (i in 1:leng){
    w = sw[i]; x = sx[i]; y = sy[i]; z = sz[i]; #print(c(x,y,z))
    his[y,x] = his[y,x]+1; his3[z,y,x] = his3[z, y,x]+1; his4[z,y,x,w] = his4[z,y,x,w]+1
    
    x0[w] = x0[w]+1
    x1[x] = x1[x]+1
    x2[y] = x2[y]+1
    x3[z] = x3[z]+1
  }
  
  
  #% Continue with calculation of Mutual Information and H2
  his = his/leng;
  his4 <- his4/leng;
  H4 <- 0
  
  if (q > 1){
    for (j in 1:b){
      for (k in 1:b){
        for (l in 1:b){
          for (i in 1:b){
            if (his4[j,k,l,i] > 0){
                an = -his4[j,k,l,i] * (his4[j,k,l,i]^(1-q)-1)/(1-q);
                H4 = H4 + an
            }
          }
        }
      }
    }
  }
  
  if (q <= 1){
    for (j in 1:b){
      for (k in 1:b){
        for (l in 1:b){
          for (i in 1:b){
            if (his4[j,k,l,i] > 0){
              an = -his4[j,k,l,i] * log2(his4[j,k,l,i])
              H4 = H4 + an
            }
          }
        }
      }
    }
  }
  
  
  #% Calculate fracn
  a = dim(his)[1]; b = dim(his)[2];
  his = his*leng;
  m = 1
  hnz <- c()
  
  for (k in 1:a){
    for (j in 1:b){
      if (his[k,j]>0){
        hnz[m]=his[k,j]
        m=m+1
      }
    }
  }
  lenz = length(hnz); 
  num5 = 0; 
  for (k  in 1:lenz){
    if (hnz[k] >=5) num5 = num5 + 1
  }
  fracn = num5/lenz
  return(list("H4" = H4, "fracn" = fracn))
}

# estimateH4(X[, 1], X[, 2], X[, 3], X[, 4], 5, 1)
