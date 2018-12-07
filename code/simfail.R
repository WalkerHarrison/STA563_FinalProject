set.seed(4)

n <- 100
X1 <- rnorm(n, mean = 0, sd = 1)
X2 <- rnorm(n, mean = 0, sd = 1)
X3 <- rnorm(n, mean = 0, sd = 1)
X4 <- rnorm(n, mean = 0, sd = 1)
X5 <- rnorm(n, mean = 0, sd = 1)
X6 <- rnorm(n, mean = 0, sd = 1)
X7 <- rnorm(n, mean = 0, sd = 1)
Y <- X1 + X2 + X3 + X4 + rnorm(n, mean = 0, sd = 1)

1/2*log2(2*pi*exp(1)*5)
1/2*log2(2*pi*exp(1)*4)
1/2*log2(2*pi*exp(1)*3)
1/2*log2(2*pi*exp(1)*2)
1/2*log2(2*pi*exp(1)*1)

X.real <- sapply(1:8, function(x) ifelse(x<6, 1/2*log2(2*pi*exp(1)*((6-x))), 
                                    1/2*log2(2*pi*exp(1)*1)))

spurr.cor <- sqrt((1 - cor(Y - X1 - X2 - X3 - X4, X5)^2))
X.sim <- c()
X.sim <- c(X.sim, 1/2*log2(2*pi*exp(1)*((var(Y)))))
X.sim <- c(X.sim, 1/2*log2(2*pi*exp(1)*((var(Y - X1)))))
X.sim <- c(X.sim, 1/2*log2(2*pi*exp(1)*((var(Y - X1 - X2)))))
X.sim <- c(X.sim, 1/2*log2(2*pi*exp(1)*((var(Y - X1 - X2 - X3)))))
X.sim <- c(X.sim, 1/2*log2(2*pi*exp(1)*((var(Y - X1 - X2 - X3 - X4)))))
X.sim <- c(X.sim, 1/2*log2(2*pi*exp(1)*(spurr.cor*(var(Y - X1 - X2 - X3 - X4)))))
X.sim <- c(X.sim, 1/2*log2(2*pi*exp(1)*((spurr.cor - 4*(1-spurr.cor))*(var(Y - X1 - X2 - X3 - X4)))))
X.sim <- c(X.sim, 1/2*log2(2*pi*exp(1)*((spurr.cor - 8*(1-spurr.cor))*(var(Y - X1 - X2 - X3 - X4)))))

plot(X.real, ylim = c(0, 4)); lines(X.real)
points(X.sim, ylim = c(0, 4), col = "blue"); lines(X.sim, col = "blue")
