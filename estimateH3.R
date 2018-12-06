estimateH3 <- function(x, y, z, pb, q){
  
  leng = length(x)
  binsize = ceiling(sqrt(pb*leng))   #% Estimated bin size for ~pb points/bin
  b       = floor(leng/binsize)   #% Number of bins
  binsize = floor(leng/b)         #% Integer number of elements per bin
  
  #% Initialize matrices to zero
  his = matrix(0, b, b)
  his3 <- lapply(1:b, function(bb) matrix(0, b, b))
  
  x1  = matrix(0, b, 1)
  x2  = matrix(0, b, 1)
  x3  = matrix(0, b, 1)
  
  #% Sort values of X and Y (descending order)
  #[sx,ix] = sort(x,2,'descend');
  tmp = sort(x, index.return=TRUE, decreasing = TRUE); sx = tmp$x; ix=tmp$ix
  tmp = sort(y, index.return=TRUE, decreasing = TRUE); sy = tmp$x; iy=tmp$ix 
  tmp = sort(z, index.return=TRUE, decreasing = TRUE); sz = tmp$x; iz=tmp$ix 

  #% Replace elements of X and Y by their rank order
  linsp  = 1:leng
  sx[ix] = linsp;
  sy[iy] = linsp;
  sz[iz] = linsp;
  
  #% Put sorted elements into bins
  sx = ceiling(sx/binsize);
  sy = ceiling(sy/binsize);
  sz = ceiling(sz/binsize);
  
  # % Put largest values in highest bin
  for (p in 1:leng){
    if (sx[p] > b) sx[p] = b
    if (sy[p] > b) sy[p] = b
    if (sz[p] > b) sz[p] = b
    }
  
  #% Populate histograms
  for (i in 1:leng){
    x = sx[i]; y = sy[i]; z = sz[i]; #print(c(x,y,z))
    his[y,x] = his[y,x]+1;
    #print(his3)
    his3[[x]][z,y] = his3[[x]][z,y]+1;
    x1[x] = x1[x]+1
    x2[y] = x2[y]+1
    x3[z] = x3[z]+1
  }
  
  
  #% Continue with calculation of Mutual Information and H2
  his = his/leng;
  his3 <- lapply(his3, function(his) his/leng)
  H3 <- 0

  if (q > 1){
    for (j in 1:b){
      for (k in 1:b){
        for (l in 1:b){
          if (his3[[j]][k,l] > 0){
            an = -his3[[j]][k,l] * (his3[[j]][k,l]^(1-q)-1)/(1-q);
            H3 = H3 + an
          }
        }
      }
    }
  }
  
  if (q <= 1){
    for (j in 1:b){
      for (k in 1:b){
        for (l in 1:b){
          if (his3[[j]][k,l] > 0){
            an = -his3[[j]][k,l] * log2(his3[[j]][k,l])
            H3 = H3 + an
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
  return(list('H3' = H3, 'fracn' = fracn))
}

# estimateH3(X[, 1], X[, 2], X[, 3], 5, 1)
