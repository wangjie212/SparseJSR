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

```Julia
julia> A = [[1 -1 0; -0.5 1 0; 1 1 0], [0.5 1 0; -1 1 0; -1 -0.5 0]]
julia> d = 2 # the relaxation order
julia> ub = SparseJSR(A, d, TS = "block") # computing an upper bound of the joint spetral radius of A via sparse SOS
julia> ub = JSR(A, d) # computing an upper bound of the joint spetral radius of A via dense SOS
```
By default, the initial lower bound and the initial upper bound for bisection are lb = 0 and ub = 2 respectively.  
The default tolerance is tol = 1e-5.  
The default sparse order is SparseOrder = 1.  
One can also set TS = "MD" to use approximately smallest chordal extentions (recommended when the relaxation order d = 1).

Gripenberg's algorithm can be used to produce a lower bound lb and an upper bound ub satisfying ub - lb < δ.
```Julia
julia> lb,ub = gripenberg(A, δ=0.2)
```

## References
[1] [SparseJSR: A Fast Algorithm to Compute Joint Spectral Radius via Sparse SOS Decompositions](https://arxiv.org/abs/2008.11441), Jie Wang, Martina Maggio and Victor Magron, ACC'21, 2021.    

## Contact
[Jie Wang](https://wangjie212.github.io/jiewang/): jwang@laas.fr
