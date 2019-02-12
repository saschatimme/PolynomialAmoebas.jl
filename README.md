# PolynomialAmoebas.jl

| **Documentation** | **Build Status** |
|:-----------------:|:----------------:|
| [![][docs-stable-img]][docs-stable-url] | [![Build Status][build-img]][build-url] |
| [![][docs-latest-img]][docs-latest-url] | [![Codecov branch][codecov-img]][codecov-url] |

A package to compute [Amoebas](https://en.wikipedia.org/wiki/Amoeba_(mathematics)) in 2 and 3 variables.

## Installation
Just do
```julia
julia> using Pkg;
julia> pkg"add PolynomialAmoebas";
```
in a running Julia session.


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

There is also a Jupyter notebook [available](https://github.com/saschatimme/PolynomialAmoebas.jl/blob/master/docs/notebooks/Overview.ipynb) to get an overview over the capabilties of the package.
In order to set everything up for this follow the instructions in the [IJulia.jl](https://github.com/JuliaLang/IJulia.jl) repository.

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-stable-url]: https://saschatimme.github.io/PolynomialAmoebas.jl/stable
[docs-latest-url]: https://saschatimme.github.io/PolynomialAmoebas.jl/latest

[build-img]: https://travis-ci.org/saschatimme/PolynomialAmoebas.jl.svg?branch=master
[build-url]: https://travis-ci.org/saschatimme/PolynomialAmoebas.jl
[codecov-img]: https://codecov.io/gh/saschatimme/PolynomialAmoebas.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/saschatimme/PolynomialAmoebas.jl
