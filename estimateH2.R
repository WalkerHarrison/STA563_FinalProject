estimateH2 <- function(x, y, pb, q){
  
  #x <- as.matrix(x); y <- as.matrix(y)
  #dim.x <- dim(x); if (dim.x[1] > dim.x[2]) {x = t(x)}
  #dim.y <- dim(y); if (dim.y[1] > dim.y[2]) {y = t(y)}

  #leng    =  nrow(x)
  leng = length(x)
  #print(leng)
  binsize = ceiling(sqrt(pb*leng))   #% Estimated bin size for ~pb points/bin
  #print(binsize)
  b       = floor(leng/binsize)   #% Number of bins
  #print(b)
  binsize = floor(leng/b)         #% Integer number of elements per bin
  
  #% Initialize matrices to zero
  his = matrix(0, b, b)
  x1  = matrix(0, b, 1)
  x2  = matrix(0, b, 1)
  
  #% Sort values of X and Y (descending order)
  #[sx,ix] = sort(x,2,'descend');
  tmp = sort(x, index.return=TRUE, decreasing = TRUE); sx = tmp$x; ix=tmp$ix
  tmp = sort(y, index.return=TRUE, decreasing = TRUE); sy = tmp$x; iy=tmp$ix 
  
  #% Replace elements of X and Y by their rank order
  linsp  = 1:leng
  sx[ix] = linsp;
  sy[iy] = linsp;
  
  #% Put sorted elements into bins
  sx = ceiling(sx/binsize);
  sy = ceiling(sy/binsize);
  
 # % Put largest values in highest bin
  for (p in 1:leng){
    if (sx[p] > b) sx[p] = b
    if (sy[p] > b) sy[p] = b
  }
  
  #% Populate histograms
  for (i in 1:leng){
    x = sx[i]; y = sy[i];
    his[y,x] = his[y,x]+1;
    x1[x] = x1[x]+1
    x2[y] = x2[y]+1
  }
  
  
  #% Continue with calculation of Mutual Information and H2
  his = his/leng;
  x1  = x1/leng;
  x2  = x2/leng;
  res = 0; 
  res2 = 0;
  if (q > 1){
    for (j in 1:b){
      for (k in 1:b){
        if (his[j,k] > 0){
          an2  = his[j,k] * (his[j,k]^(1-q)-1)/(1-q)               #% H2
          res2 = res2 + an2;                                      #% H2
          an   = his[j,k]*(his[j,k]/(x1[k]*x2[j])^(1-q)-1)/(1-q) ;#% MI
          res  = res + an;                                        #% MI 
        }
      }
    }
  }
  
  if (q <= 1){
    for (j in 1:b){
      for (k in 1:b){
        if (his[j,k] > 0){
          an2  = his[j,k]*log2(his[j,k]);                        #% H2
          res2 = res2 + an2;                                      #% H2
          an   = his[j,k]*log2(his[j,k]/(x1[k]*x2[j]));          #% MI
          res  = res + an;                                        #% MI 
        }
      }
    }
  }
  
  #print("here")
  mutinfo = res;
  H2      = -res2;
  
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
  return(list("mutinfo" = mutinfo, "fracn" = fracn, "H2" = H2))
}
