# compile nodag.f with
# R CMD SHLIB nodag.f -llapack -lblas 
dyn.load("nodag.so")

############################
##   generate data from 
##     X1       X2
##      \       /
##       v     v 
##          X3
##
############################

B <- matrix(nrow = 3, ncol = 3, 
	    data = c(0, 0, 0, 0, 0, 0, -1, 2, 0))

n <- 100
p <- 3

X <- matrix(nrow = n, ncol = p, rnorm(n*p))
for (j in 1:p){
	for (i in 1:n){
		X[i, j] = X[i, j] + sum(X[i,]*B[,j]) 
	}
}

### penalization coefficient
lambda <- 0.3
out <- .Fortran("NODAG", as.integer(p), as.double(cor(X)),
                as.double(diag(p)), as.double(lambda), 
                as.double(1e-5), as.double(0.5), 
                as.integer(1000))

A <- matrix(ncol = p, nrow = p, data = out[[3]])
message("algorithm terminates after ", out[[7]], " iterations")
message("True adjacency matrix:")
print(sign(abs(B)))
message("Estimated A matrix ")
print(A)

message("Estimated adjaceny matrix")
print(0 + (A - diag(diag(A)) !=0)) 

