# Poisson1D-PETSc
This code solves the one-dimensional Poisson equation using the finite element method written in C with the PETSc library. It uses quadratic shape functions and enforces homogeneous Dirichlet boundary conditions.
Compile using the given makefile

```
make Poisson1D
```
and run simply with 
```
 ./Poisson1D
```
By default, it uses 5 partitions. and provide $L_2$ and $H^1$ norm errors. A number of partitions can be changed using command line options. For example to find errors for 20 partitions can be found using 
```
./Poisson1D -fem_parts 20
```

reference: PETSc for Partial Differential Equations: Numerical Solutions in C and Python, by Ed Bueler.
