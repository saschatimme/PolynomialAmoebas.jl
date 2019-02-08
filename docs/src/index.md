# PolynomialAmoebas.jl

PolynomialAmoebas.jl is a package to compute and visualize the amoeba, coamoeba and imaginary projection
of bi- and trivariate polynomials as well as the contour and the spine of an two-dimensional amoeba.

## Getting started
To construct polynomials we export the macro `@polyvar` from the [DynamicPolynomials.jl](https://github.com/JuliaAlgebra/DynamicPolynomials.jl) package.

```julia
using PolynomialAmoebas
# Create variables
@polyvar x y
# construct a polynomial
f = x^2*y + y^2 + 3x^2*y^3 + y^4 + x^4*y^4
```

To compute the amoeba of `f` we then can simply do
```julia
A = amoeba(f)
```
To visualize the amoeba we use the plotting capabilities provided by [Plots.jl](http://docs.juliaplots.org/latest/).
Just do
```julia
using Plots

plot(A)
```
and you obtain
![amoeba text](./assets/amoeba.png)
