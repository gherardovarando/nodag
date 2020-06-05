# nodag
## structure learning without imposing acyclicity 


This repository contains fortran implementation of a 
proximal gradient method to estimate the structure of 
a strucural equation model (SEM)
 without imposing the acyclicity constraint. 
See the preprint on [arXiv:2006.03005](https://arxiv.org/abs/2006.03005).

The subrutine `NODAG` solve the following l1-penalized 
minus log-likelihood minimization:

```
argmin  -2log(det(A)) + trace( Sigma AA^t) + lambda ||A||_1 
```

### use 

The fortran subrutine `NODAG` can be easily used both from python and R 

* python: compile `nodag.f` with `f2py` using `f2py -llapack -c -m nodag nodag.f` 
* R: compile `nodag.f` with `R CMD SHLIB nodag.f -llapack -lblas`

Check the provided examples to see how to load and call the subroutine.  
 



