# Amoeba

The amoeba $\mathcal{A}(𝑓)$ of a Laurent polynomial $f \in \mathbb{C}[z_1^\pm, \ldots, z_n^\pm]$ is the image
of the non-singular hypersurface $\mathcal{V}(f) \subset (\mathbb{C}^*)^n$
under the Log-absolute value map given by

$\operatorname{Log}|\cdot|:\; (\mathbb{C}^*)^n \rightarrow \mathbb{R}^n, \quad (z_1, \ldots, z_n) \mapsto (\log|z_1|, \ldots, \log|z_n|) \;.$

It was first introduced by Gelfand, Kapranov and Zelevinsky in their book
[Discriminants, Resultants, and Multidimensional Determinants](http://www.springer.com/de/book/9780817647704) in 1994.
```@docs
amoeba
```
