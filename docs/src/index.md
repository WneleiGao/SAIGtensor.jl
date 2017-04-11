# SAIGtensor Documentation

Tensor tools for seismic data processing

## Test for latex
inline math ``\frac{\partial^2 u}{\partial t^2}-\Delta^2 u = 0``

Here is an equation:
```math
\frac{n!}{k!(n-k)!} = \binom{n}{k}
```
## Include code to documentation
```@example
a = 1
b = 2
c = a+b
```

## Include figure
```@example
using PyPlot
a = rand(100)
plot(a)
savefig("plot.svg"); nothing # hide
```
![](plot.svg)

## Index
```@index
```
