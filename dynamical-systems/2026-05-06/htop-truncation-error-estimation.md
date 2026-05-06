# $h_{\rm top}$ approximation error due to Milnor-Thurston kneading determinant truncation

### 6 May 2026

Let $D(t) = \sum_{n \geq 0} a_nt^n$ denote the kneading determinant of a piecewise-continuous piecewise-monotone map $f: I \to I$, where $a_n \in \mathbb{Z}$. Then the topological entropy $h_{\rm top}$ satisfies $r = e^{-h_{\rm top}}$, where $r$ is the smallest zero of $D(t)$ (we have $r \in (0, 1]$, excluding the case of infinite topological entropy). For $N \in \mathbb{Z}_+$, let $D_N(t)$ denote the degree-$N$ truncated kneading determinant $\sum_{n = 0}^N a_nt^n$, and let $E_N(t) = D_N(t) - D(t)$ denote the truncation error.

Since $D(r) = 0$, we have $E_N(r) = D_N(r) - D(r) = -\sum_{n > N} a_nr^n = D_N(r)$. Assume we have positive topological entropy; then $r < 1$. It follows that $|D_N(r)| = \mathcal{O}(r^N)$.

Assume the root $r$ of $D(t)$ is simple (by The Real Teapot [1], for unimodal maps this is always the case); then $D'(r) \neq 0$. Let $r_N$ be the smallest root of $D_N(t)$; let $\delta = r_N - r$. Expanding $D_N(r_N)$ at $r$, we obtain
$$0 = D_N(r_N) = D_N(r + \delta) = D(r + \delta) + E_N(r + \delta)$$
$$= \cancel{D(r)} + \delta D'(r) + E_N(r) + \delta E_N'(r) + \mathcal{O}(\delta^2)$$
$$= \delta D'(r) + D_N(r) + \mathcal{O}(\delta^2)$$
so that
$$r_N - r = \delta \approx -\frac{D_N(r)}{D'(r)} = \mathcal{O}(r^N),$$
since $D'(r)$ is a nonzero constant. Letting $h_N = -\log r_N$, we expand to obtain
$$h_N = -\log(r + \delta) = -\log r + \sum_{k=1}^\infty \frac{(-\delta/r)^k}{k}$$
so that for $|\delta| \ll 1$ (e.g., for $N$ sufficiently large) we have the approximation
$$h_N - h_{\rm top} = -\sum_{k = 1}^\infty \frac{(-\delta/r)^k}{k} \approx -\frac{\delta}{r} = \mathcal{O}(r^N)$$
$$= \mathcal{O}(e^{-Nh_{\rm top}}).\tag{1}$$

What $(1)$ is ultimately saying is that the error in the truncation-derived approximation $h_N$ of $h_{\rm top}$ has the same asymptotic behavior as the reciprocal of the number of fixed points of $f^N$. This was first shown by Collet-Crutchfield-Eckmann [2].

#### *A posteriori* error estimation

Recall $D_N(t) - D(t) = -\sum_{n > N} a_nr^n$. We thus have
$$|D_N(r_N) - D(r_N)|  = \left|\sum_{n > N} a_nr_N^n\right| \leq \sum_{n > N} |a_n| r_N^n,$$
so, by root perturbation,
$$|h_N - h_{\rm top}| \lesssim \frac{\sum_{n > N} |a_n|r_N^n}{(1-r_N)|D_N'(r_N)|}.$$
The ansatz $|a_n| \sim 1$ (cf. [3], p. 332, which suggests $|a_n| = \mathcal{O}(1)$) results for large $N$ in the approximate upper bound
$$|h_N - h_{\rm top}| \lesssim \frac{r_N^{N}}{(1-r_N)|D_N'(r_N)|},\tag{2}$$
which is empirically a good estimate of the entropy error; see below the results for three maps:

![Logistic map entropy error estimate](https://hinsley.github.io/media/dynamical-systems/2026-05-06/logistic-map-entropy-error-estimate.png)

![Tent map entropy error estimate](https://hinsley.github.io/media/dynamical-systems/2026-05-06/tent-map-entropy-error-estimate.png)

![Cubic map entropy error estimate](https://hinsley.github.io/media/dynamical-systems/2026-05-06/cubic-map-entropy-error-estimate.png)

*Remark.* I would look at Rugh's 2015 paper [4] relating the kneading determinant to the regularized Fredholm determinant of a Ruelle transfer operator for additional clues about why this seems to work so well; questions of numerical conditioning should be clearer from that viewpoint. In particular, this may make clear why the ansatz $|a_n| \sim 1$ seems successful, or it may yield a more generally applicable ansatz.

If we assume that the initial estimate of the topological entropy has error equal to the maximal topological entropy possible given the number $\ell$ of laps in the domain of the map $f$, which is $\log(\ell)$, and we assume that the exponential error reduction per truncation degree adheres to the rate $r$ as suggested by $(1)$ (at least suggested in the limit for large $N$), then we have the stopping rule
$$N_{\rm stop} =\left\lceil\frac{\log\log\ell - \log\varepsilon}{h_N}\right\rceil,\tag{3}$$
where $\varepsilon$ is the desired error tolerance (e.g., $\varepsilon = 10^{-16}$ for double-precision floating point representation). This can be seen as the "running estimate" for each map in the following plot:

![Collet-Eckmann stopping rule plot](https://hinsley.github.io/media/dynamical-systems/2026-05-06/collet-eckmann-stopping-rule.png)

It should be noted that replacing $\log\ell$ in $(3)$ by $h_N$ results in slightly undershooting the true truncation degree required to meet the entropy-error tolerance $\varepsilon$.

On the other hand, letting $(2)$ yield a current estimate of the entropy error, we obtain a stopping rule
$$N_{\rm stop} = \left\lceil-\frac{\log(1-r_N) + \log|D_N'(r_N)| + \log\varepsilon}{h_N}\right\rceil\tag{4}$$
that matches the convergence rate of the above rule for the logistic and cubic maps, but converges by iterate $9$ for the tent map, whereas stopping rule $(3)$ only converges for the tent map by iterate $11$:

![Smart stopping rule plot](https://hinsley.github.io/media/dynamical-systems/2026-05-06/smart-stopping-rule.png)

#### References

[1] Alsedà Llu, Bobok J, Misiurewicz M, Snoha Ľubom. The Real Teapot. _Ergodic Theory and Dynamical Systems_. 2025;45(10):2945-2975. <https://doi.org/10.1017/etds.2025.15>

[2] Collet, P., Crutchfield, J.P. & Eckmann, J.P. Computing the topological entropy of maps. _Commun.Math. Phys._ **88**, 257-262 (1983). <https://doi.org/10.1007/BF01209479>

[3] Cvitanović, P., Artuso, R., Mainieri, R., Tanner, G., and Vattay, G. _Chaos: Classical and Quantum_, Chapter 20, Counting. <https://chaosbook.org/version17/chapters/count-2p.pdf>

[4] Rugh, H.H. The Milnor-Thurston Determinant and the Ruelle Transfer Operator. _Commun. Math. Phys._ **342**, 603-614 (2016). <https://doi.org/10.1007/s00220-015-2515-5>
