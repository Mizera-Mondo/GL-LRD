1. About alpha: it is suggested that alpha_real is constructed as alpha_ref/T, since the smoothness
measurement tr(DX'*L*DX) will expand as T grows, for what the normalization term ||L||_F^2
will eventually vanish when T is large, if alpha does not change accroding to T.

2. But wait, a too much small alpha will possibly make the update of X neglect the term alpha*Tr{DX'*L*DX}
since any change on X try to make this term small may greatly affect the term 1/2||D(Y - X)||_F^2

3. There may be some bug in update of P and Q snice the process does NOT actually 
lower the constraint X = P*Q'. This may because a too low choice of k, who knows.

# TODO:
GL_LRD_SVD.m, Solving 
```math
\begin{alignat}{1}
  &\mathrm{arg~min}~~~ &&\frac{1}{2}||D(\mathbf{Y}-\mathbf{X})||_F^2+\alpha \mathrm{Tr}\left(D(\mathbf{X})^\top*\mathbf{L}D(\mathbf{X})\right), \\
  &\mathrm{s.t.}~~~ &&\mathrm{rank}(\mathbf{X})\leq k
\end{alignat}
```
using 1-step SVD.
