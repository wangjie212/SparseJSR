# SparseJSR
SparseJSR is an efficient tool to compute joint spetral radii based on sparse SOS decompositions. To use SparseJSR in Julia, run
```Julia
pkg> add https://github.com/wangjie212/SparseJSR
 ```

## Dependencies
- Julia
- MOSEK
- JuMP

SparseJSR has been tested on WINDOW 10, Julia 1.2, JuMP 0.21 and MOSEK 8.1.
## Usage

Assume that A is a tuple of matrices.
```Julia
julia> d = 2 # the relaxation order
julia> ub = SparseJSR(A, d, TS = "block") # ub is a upper bound of the joint spetral radius of A
```
By default, the initial lower bound and the initial upper bound for bisection are lb = 0 and ub = 2 respectively.  
The default tolerance is tol = 1e-5.  
One can set TS = "block" or TS = "MD" to use different chordal extentions.


## References
[1] [SparseJSR: A Fast Algorithm to Compute Joint Spectral Radius via Sparse SOS Decompositions](https://arxiv.org/abs/2008.11441), Jie Wang, Martina Maggio and Victor Magron, 2020.    

## Contact
[Jie Wang](https://wangjie212.github.io/jiewang/): jwang@laas.fr
