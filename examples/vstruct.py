## example of recovery of a v-structure
import numpy as np
## compile nodag.f with 
## f2py -llapack -c -m nodag nodag.f
import nodag as nd

############################
##   generate data from 
##     X1       X2
##      \       /
##       v     v 
##          X3
##
############################
B = np.asarray([[0, 0, -1],
                [0, 0, 2],
                [0, 0, 0]]) 
n, p = 100, 3
X = np.random.randn(n, p)
for j in range(1, p):
    for i in range(1, n):
        X[i, j] = X[i, j] + X[i,:].dot(B[:,j])
Sigma = np.corrcoef(X.T)
Aest, diff, value, itr = nd.nodag(Sigma, lambd = 0.5)

print('terminated with function value {0} (diff {1})  after {2} iterations'.format(value, diff, itr)) 
print('Estimated A matrix')
print(Aest)

print('Estimated adjacency matrix')

print( Aest - np.diag(Aest.diagonal()) != 0) 

print('True adjacency matrix')
print(B != 0)

