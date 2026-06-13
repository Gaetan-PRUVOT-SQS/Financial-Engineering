# Course 3 — Portfolio Optimization & Optimal Execution

Columbia "Financial Engineering and Risk Management" (Haugh & Iyengar) — *Optimization Methods in Asset Management*. Self-contained reference: optimization/probability review, Markowitz mean-variance optimization, the efficient frontier (with and without a risk-free asset), two-fund / one-fund theorems, the tangency (market) portfolio, Sharpe ratio, CAPM, practical MVO (estimation error, shrinkage, constraints, robustness, VaR/CVaR, leveraged ETFs), behavioral/statistical biases, and optimal trading/execution (market impact, mean-variance-liquidity, Almgren-Chriss style execution).

---

## Key formulas cheat-sheet

| Quantity | Formula |
|---|---|
| Portfolio expected return | $\mu_x = \sum_i \mu_i x_i = \mu^\top x$ |
| Portfolio variance | $\sigma_x^2 = \sum_{i,j}\sigma_{ij}x_i x_j = x^\top V x$ |
| Covariance | $\sigma_{ij}=\rho_{ij}\sigma_i\sigma_j$ |
| Budget constraint | $\sum_i x_i = 1$ (i.e. $\mathbf{1}^\top x=1$) |
| Min-risk MVO (target $r$) | $\min_x x^\top V x \;$ s.t. $\mu^\top x = r,\;\mathbf{1}^\top x=1$ |
| Risk-adjusted MVO | $\max_x \;\mu^\top x - \tau\, x^\top V x \;$ s.t. $\mathbf{1}^\top x=1$ |
| Lagrange stationarity (no $r_f$) | $2\sum_j\sigma_{ij}x_j - v\mu_i - u = 0\;\forall i$ |
| Frontier shape | $\sigma^2(r) = a r^2 + b r + c$ (parabola in $(\sigma^2,r)$, hyperbola in $(\sigma,r)$) |
| Two-fund theorem | $y^\ast(r) = (1-\beta)x^{(1)} + \beta x^{(2)},\;\beta=\tfrac{r-r_1}{r_2-r_1}$ |
| MVO with $r_f$ (excess) | $\hat\mu_i = \mu_i - r_f$; optimum $x(\tau)=\tfrac{1}{2\tau}V^{-1}\hat\mu$ |
| Tangency / Sharpe portfolio | $s^\ast = \dfrac{V^{-1}\hat\mu}{\mathbf{1}^\top V^{-1}\hat\mu}$ (one-fund theorem) |
| Capital Market Line | $\mu_x = r_f + \dfrac{\mu_m-r_f}{\sigma_m}\,\sigma_x$ |
| Sharpe ratio | $S=\dfrac{\mu_x-r_f}{\sigma_x}$; max achievable $=$ CML slope $=$ price of risk |
| CAPM | $\mu_j - r_f = \beta_j(\mu_m - r_f),\;\;\beta_j=\dfrac{\sigma_{jm}}{\sigma_m^2}$ |
| Risk decomposition | $\sigma_j^2 = \beta_j^2\sigma_m^2 + \mathrm{var}(\epsilon)$ (market + residual) |
| CAPM price of payoff $X$ | $P = \dfrac{E[X]}{1+r_f} - \dfrac{\mathrm{cov}(X,r_m)}{(1+r_f)\mathrm{var}(r_m)}(\mu_m-r_f)$ |
| Mean shrinkage | $\mu_i^{(sh)}=\alpha\mu_i^{(est)}+(1-\alpha)\big(\tfrac1d\sum_j\mu_j^{(est)}\big)$ |
| Covariance shrinkage | $V^{(sh)}=\alpha V^{(est)}+(1-\alpha)\big(\tfrac1d\sum_i\sigma_i^{2(est)}\big)I$ |
| VaR (Normal) | $\mathrm{VaR}_p(L)=\mu+\sigma\Phi^{-1}(p)$ |
| CVaR (Normal) | $\mathrm{CVaR}_p(L)=\mu+\sigma\,\tfrac{1}{1-p}\int_p^1\Phi^{-1}(\beta)\,d\beta$ |
| Execution slippage cost | $C(n)=\sum_k g(n_k)x_k + \sum_k h(n_k)n_k$ |
| Execution risk | $V(n)=\sigma^2\sum_k x_k^2$ |
| Optimal execution | $\min_n\; C(n) + \rho\, V(n)$ |

---

## 1. Optimization & linear algebra review

### 1.1 Linear programming and duality
General LP: $\max/\min_x\; c^\top x$ s.t. $A_{eq}x=b_{eq},\;A_{in}x\le b_{in}$.
**Hedging example (LP):** hold $\theta\in\mathbb{R}^d$ shares, cost $p^\top\theta$; payoff in state $i$ is $\sum_j S_{ij}\theta_j$; payoff $y=S\theta$ hedges obligation $X$ if $y\ge X$. Problem: $\min_\theta p^\top\theta$ s.t. $S\theta\ge X$.

**Primal/Dual pair:**
$$P=\min_x\{c^\top x:Ax\ge b\}\qquad D=\max_u\{b^\top u: A^\top u=c,\;u\ge 0\}$$
- **Weak duality:** $P\ge D$; any feasible $x,u$ bracket: $c^\top x\ge P\ge D\ge b^\top u$.
- **Strong duality:** if $P$ or $D$ finite, $P=D$. Dual of the dual is the primal.
- Constructed via *Lagrangian relaxation*: dualize and relax constraints.

### 1.2 Unconstrained nonlinear optimization
$\min_{x\in\mathbb{R}^n} f(x)$. Sufficient local-min conditions: gradient $\nabla f(x)=0$ (stationarity) **and** Hessian $\nabla^2 f(x)\succeq 0$ (PSD). For **convex** $f$, $\nabla f(x)=0$ alone is sufficient for a global min.
*Example* $f=x_1^2+3x_1x_2+x_2^3$: stationary points $x=0$ and $x=(-9/4,\,3/2)$; Hessian $\begin{bmatrix}2&3\\3&6x_2\end{bmatrix}$ is PSD only at $(-9/4,3/2)$ → that is the local min.

### 1.3 Lagrangian method
$\max_x f(x)$ s.t. $g(x)=0$ → form $L(x,v)=f(x)-v\,g(x)$, solve $\nabla_x L=0$ for $x(v)$, then fix $v$ from the constraint.
*Example* $\max 2\ln(1+x_1)+4\ln(1+x_2)$ s.t. $x_1+x_2=12$: $\nabla L=0\Rightarrow x_1=\tfrac2v-1,\,x_2=\tfrac4v-1$; constraint ⇒ $v=3/7$, $x=(11/3,\,25/3)$.

**Portfolio Lagrangian (preview):** $\max_x \mu^\top x-\lambda x^\top V x$ s.t. $\mathbf{1}^\top x=1$.
$L=\mu^\top x-\lambda x^\top V x - v(\mathbf{1}^\top x-1)$; $\nabla_x L=\mu-2\lambda Vx-v\mathbf{1}=0\Rightarrow x=\tfrac{1}{2\lambda}V^{-1}(\mu-v\mathbf{1})$; budget ⇒ $v=\dfrac{\mathbf{1}^\top V^{-1}\mu - 2\lambda}{\mathbf{1}^\top V^{-1}\mathbf{1}}$.

### 1.4 Matrices, vectors, norms
$A\in\mathbb{R}^{m\times d}$; transpose $(A^\top)_{ij}=A_{ji}$; product $C=AB$, $c_{ij}=\sum_k a_{ik}b_{kj}$. Inner product $v\cdot w=v^\top w=\sum_i v_iw_i$; $\ell_2$ norm $\|v\|_2=\sqrt{v\cdot v}$; $\cos\theta=\frac{v\cdot w}{\|v\|\|w\|}$. $\ell_1=\sum|x_i|$, $\ell_\infty=\max_i|x_i|$. Rank: row rank = column rank $\le\min(m,d)$; square $A$ with full rank $n$ is invertible. A function is linear iff $f(x)=Ax$.

### 1.5 Probability review
- **Discrete:** $E[X]=\sum_i x_i p(x_i)$; $\mathrm{Var}(X)=E[X^2]-E[X]^2$.
- **Binomial** $\mathrm{Bin}(n,p)$: $P(X=r)=\binom{n}{r}p^r(1-p)^{n-r}$, $E=np$, $\mathrm{Var}=np(1-p)$.
- **Poisson** $\lambda$: $P(X=r)=\frac{\lambda^r e^{-\lambda}}{r!}$, $E=\mathrm{Var}=\lambda$.
- **Bayes:** $P(A\mid B)=\frac{P(B\mid A)P(A)}{\sum_j P(B\mid A_j)P(A_j)}$.
- **Normal** $N(\mu,\sigma^2)$: $f(x)=\frac{1}{\sqrt{2\pi\sigma^2}}e^{-(x-\mu)^2/2\sigma^2}$. **Log-normal** $\log X\sim N(\mu,\sigma^2)$: $E[X]=e^{\mu+\sigma^2/2}$, $\mathrm{Var}(X)=e^{2\mu+\sigma^2}(e^{\sigma^2}-1)$.
- **Conditional identities:** $E[X]=E[E[X\mid Y]]$; $\mathrm{Var}(X)=\mathrm{Var}(E[X\mid Y])+E[\mathrm{Var}(X\mid Y)]$.
  - Random sum $W=\sum_{i=1}^N X_i$ ($X_i$ iid $\mu_x,\sigma_x^2$, $N$ random): $E[W]=\mu_x E[N]$, $\mathrm{Var}(W)=\mu_x^2\mathrm{Var}(N)+\sigma_x^2 E[N]$.
- **Multivariate normal** $X\sim MN_n(\mu,\Sigma)$: $f(x)=\frac{1}{(2\pi)^{n/2}|\Sigma|^{1/2}}e^{-\frac12(x-\mu)^\top\Sigma^{-1}(x-\mu)}$; MGF $e^{s^\top\mu+\frac12 s^\top\Sigma s}$. Marginals normal; conditional $X_2\mid X_1=x_1\sim MN(\mu_{2.1},\Sigma_{2.1})$ with $\mu_{2.1}=\mu_2+\Sigma_{21}\Sigma_{11}^{-1}(x_1-\mu_1)$, $\Sigma_{2.1}=\Sigma_{22}-\Sigma_{21}\Sigma_{11}^{-1}\Sigma_{12}$. Linear map: $\mathrm{Cov}(AX+a)=A\,\mathrm{Cov}(X)\,A^\top$.
- **Brownian motion** $B(\mu,\sigma)$: independent increments, $X_{t+s}-X_t\sim N(\mu s,\sigma^2 s)$, continuous paths; $X_t=x+\mu t+\sigma W_t$. SBM $W_t$: $W_0=0$, $E_0[W_{t+s}W_s]=s$.
- **Geometric BM** $X_t=e^{(\mu-\sigma^2/2)t+\sigma W_t}$: $\log(X_{t+s}/X_t)\sim N((\mu-\sigma^2/2)s,\sigma^2 s)$, $E_t[X_{t+s}]=e^{\mu s}X_t$ (positive prices → stock-price / Black-Scholes model).
- **Martingale** $\{X_n\}$ w.r.t. filtration $F_n$: $E[|X_n|]<\infty$ and $E[X_{n+m}\mid F_n]=X_n$ (submartingale $\ge$, supermartingale $\le$). E.g. $M_n=S_n-n\mu$ from a random walk; Polya's urn $M_n=X_n/(n+2)$.

---

## 2. Assets, portfolios and quantifying risk

- Asset price $P(t)$; gross return $R(t)=P(t+1)/P(t)$; net return $r(t)=R(t)-1$.
- Total wealth $W$, dollar amount $w_i$ in asset $i$ ($w_i>0$ long, $<0$ short). **Portfolio vector** $x=(x_1,\dots,x_d)$, $x_i=$ fraction in asset $i$, $\sum_i x_i=1$.
- Random portfolio net return $r_x=\sum_i r_i x_i$.
- Mean $\mu_i=E[r_i]$, variance $\sigma_i^2=\mathrm{var}(r_i)$, covariance $\sigma_{ij}=\rho_{ij}\sigma_i\sigma_j$ (assumed constant over time).
- Portfolio mean $\mu_x=\sum_i\mu_i x_i$; variance $\sigma_x^2=\sum_{i,j}\sigma_{ij}x_ix_j$.

**Diversification reduces uncertainty.** $d$ assets each $\mu_i\equiv\mu$, $\sigma_i\equiv\sigma$, $\rho_{ij}=0$. Concentrated $x=(1,0,\dots)$ vs equal-weight $y=\tfrac1d\mathbf{1}$: both have mean $\mu$, but $\sigma_x^2=\sigma^2$ while $\sigma_y^2=\sigma^2/d$.

**2-asset numerical example.** $r_1\sim N(1,0.1)$, $r_2\sim N(2,0.5)$, $\rho=-0.25$. Then $\sigma_{12}=-0.25\sqrt{0.05}=-0.0559$. Portfolio $(x,1-x)$: $\mu_x=x+2(1-x)$; $\sigma_x^2=0.1x^2+0.5(1-x)^2+2(0.0559)x(1-x)$ (sign per $\sigma_{12}<0$).

---

## 3. Markowitz mean-variance optimization

Markowitz (1954): "return" $\equiv\mu_x$, "risk" $\equiv\sigma_x$. Three equivalent formulations:

1. **Min risk, target return:** $\min_x x^\top V x$ s.t. $\mu^\top x\ge r,\;\mathbf{1}^\top x=1$.
2. **Max return, risk budget:** $\max_x \mu^\top x$ s.t. $x^\top V x\le\bar\sigma^2,\;\mathbf{1}^\top x=1$.
3. **Risk-adjusted return:** $\max_x \mu^\top x-\tau\, x^\top V x$ s.t. $\mathbf{1}^\top x=1$ ($\tau\ge0$ = risk-aversion).

The **efficient frontier** = max return for given risk (upper boundary of the feasible $(\sigma,\mu)$ region).

### 3.1 Two-asset frontier (closed form)
$\mu_x=\mu_1x+\mu_2(1-x)$, $\sigma_x^2=\sigma_1^2x^2+\sigma_2^2(1-x)^2+2\rho\sigma_1\sigma_2 x(1-x)$.
Target return ⇒ $x=\frac{r-\mu_2}{\mu_1-\mu_2}$. Substituting:
$$\sigma^2(r)=\sigma_1^2\Big(\tfrac{r-\mu_2}{\mu_1-\mu_2}\Big)^2+\sigma_2^2\Big(\tfrac{\mu_1-r}{\mu_1-\mu_2}\Big)^2+2\rho\sigma_1\sigma_2\Big(\tfrac{r-\mu_2}{\mu_1-\mu_2}\Big)\Big(\tfrac{\mu_1-r}{\mu_1-\mu_2}\Big)=ar^2+br+c.$$
A parabola in $(\sigma^2,r)$ / hyperbola in $(\sigma,r)$ with minimum-variance point $(\sigma_{\min},r_{\min})$. Only the **top half** ($r\ge r_{\min}$) is efficient; the bottom half is dominated.

### 3.2 General $d$-asset frontier — Lagrangian derivation
$$\sigma^2(r)=\min_x x^\top V x\quad\text{s.t. }\mu^\top x=r,\;\mathbf{1}^\top x=1.$$
Lagrangian $L=\sum_{i,j}\sigma_{ij}x_ix_j - v(\sum_i\mu_ix_i-r)-u(\sum_i x_i-1)$.
Stationarity $\partial L/\partial x_i=0$:
$$2\sum_{j}\sigma_{ij}x_j - v\mu_i - u = 0,\quad i=1,\dots,d.\qquad(\ast)$$
**Theorem.** $x$ is mean-variance optimal iff it is feasible and $\exists\,u,v$ satisfying $(\ast)$.
Solve the $d+2$ linear equations in $(x_1,\dots,x_d,v,u)$ in matrix form $A\,[x;v;u]=b$:
$$A=\begin{bmatrix}2V & -\mu & -\mathbf{1}\\ \mu^\top & 0 & 0\\ \mathbf{1}^\top & 0 & 0\end{bmatrix},\quad b=\begin{bmatrix}0\\ r\\ 1\end{bmatrix}\;\Rightarrow\;[x;v;u]=A^{-1}b.$$

### 3.3 Two-fund theorem
Fix $r_1\ne r_2$ with optimal portfolios $x^{(1)},x^{(2)}$ (multipliers $(v_1,u_1),(v_2,u_2)$). For any target $r$, set $\beta=\frac{r-r_1}{r_2-r_1}$ and $y=(1-\beta)x^{(1)}+\beta x^{(2)}$. Then $y$ is feasible ($\mathbf{1}^\top y=1$, $\mu^\top y=r$) and optimal (the linear combination of multipliers $v=(1-\beta)v_1+\beta v_2$, $u=(1-\beta)u_1+\beta u_2$ satisfies $(\ast)$).
**Theorem.** *Every efficient portfolio is a combination of any two distinct efficient portfolios.* Writing $y^\ast=rg+h$ with $g=\frac{x^{(2)}-x^{(1)}}{r_2-r_1}$, $h=\frac{r_2x^{(1)}-r_1x^{(2)}}{r_2-r_1}$, the $d$-asset frontier $\sigma^2(r)=r^2(g^\top Vg)+2r(g^\top Vh)+(h^\top Vh)$ has the same parabolic structure as the 2-asset case.

---

## 4. Mean-variance with a risk-free asset

Add a risk-free asset paying $r_f$ deterministically; $x_0$ = fraction in it, contributes no risk.
$$\max\;\Big(r_f x_0+\sum_i\mu_ix_i\Big)-\tau\sum_{i,j}\sigma_{ij}x_ix_j\quad\text{s.t. }x_0+\sum_i x_i=1.$$
Substitute $x_0=1-\sum_i x_i$ and define excess return $\hat\mu_i=\mu_i-r_f$:
$$\max\;r_f+\sum_i\hat\mu_i x_i-\tau\, x^\top V x.$$
First-order conditions $\hat\mu_i-2\tau\sum_j\sigma_{ij}x_j=0$ give the **family of frontier portfolios**
$$\boxed{\,x(\tau)=\tfrac{1}{2\tau}V^{-1}\hat\mu\,}$$
The set $\{(1-\sum_i x_i(\tau),\,x(\tau)):\tau\ge0\}$.

### 4.1 One-fund theorem & tangency portfolio
The risky weights $x(\tau)=\tfrac1{2\tau}V^{-1}\hat\mu$ do not sum to 1. Normalize:
$$\boxed{\,s^\ast=\frac{V^{-1}\hat\mu}{\mathbf{1}^\top V^{-1}\hat\mu}\,}$$
$s^\ast$ is **independent of $\tau$**. Since $\sum_i x_i=1-x_0$, $x=(1-x_0)s^\ast$.
**One-fund theorem.** *Every efficient portfolio with a risk-free asset is a combination of the riskless asset and the single risky portfolio $s^\ast$.* Family $=\{(x_0,(1-x_0)s^\ast):x_0\in\mathbb{R}\}$.

### 4.2 Efficient frontier = Capital Market Line
With $\mu_s^\ast=\sum_i\mu_i s_i^\ast$, $\sigma_s^\ast=\sqrt{(s^\ast)^\top V s^\ast}$:
$$\mu_x=x_0 r_f+(1-x_0)\mu_s^\ast,\qquad\sigma_x=(1-x_0)\sigma_s^\ast.$$
A straight line from $(0,r_f)$ with slope $m=\frac{\mu_s-r_f}{\sigma_s}$ — **tangent** to the risky-only frontier at $s^\ast$. $s^\ast$ maximizes $\tan\theta=\frac{\mu^\top x-r_f}{\sqrt{x^\top V x}}$ = expected excess return / volatility.

### 4.3 Sharpe ratio
$$S=\frac{\mu_x-r_f}{\sigma_x},\qquad s^\ast=\arg\max_{x:\,\mathbf{1}^\top x=1}\frac{\mu_x-r_f}{\sigma_x}.$$
$s^\ast$ is the **Sharpe-optimal (tangency) portfolio**. All investors diversify between $r_f$ and $s^\ast$ in fixed risky proportions → returns must be correlated → leads to CAPM.

---

## 5. CAPM

### 5.1 Market portfolio = tangency portfolio
Market caps $C_i$; **market portfolio** $x^{(m)}_i=C_i/\sum_j C_j$. If all investors are MV-optimizers they all hold $s^\ast$, so $C_i=\sum_k w^{(k)}(1-x_0^{(k)})s_i^\ast$ ⇒ $x^{(m)}=s^\ast$. The **Capital Market Line** is the efficient frontier through $(0,r_f)$ and $(\sigma_m,\mu_m)$:
$$m_{\mathrm{CML}}=\frac{\mu_m-r_f}{\sigma_m}=\text{max achievable Sharpe ratio}=\text{"price of risk".}$$

**Pipeline example.** Share price \$875, expected payoff \$1000 in 1yr, $\sigma=40\%$; $r_f=5\%$, $\mu_m=17\%$, $\sigma_m=12\%$. Actual return $r_{oil}=1000/875-1=14\%$. Required (CML) return $\bar r=r_f+\frac{\mu_m-r_f}{\sigma_m}\sigma=5\%+1\cdot 40\%=45\%\gg 14\%$ → **not worth considering.**

### 5.2 CAPM pricing relation
Treat asset $j$ as a portfolio and diversify with $x^{(m)}$: $\gamma x^{(j)}+(1-\gamma)x^{(m)}$. Equating the slope of this curve at $\gamma=0$ with $m_{\mathrm{CML}}$:
$$\boxed{\,\mu_j-r_f=\underbrace{\frac{\sigma_{jm}}{\sigma_m^2}}_{\beta_j}(\mu_m-r_f)\,}$$

### 5.3 CAPM via regression / Security Market Line
Regress excess returns: $(r_j-r_f)=\alpha+\beta(r_m-r_f)+\epsilon_j$, with $\beta_j=\frac{\mathrm{cov}(r_j-r_f,r_m-r_f)}{\mathrm{var}(r_m-r_f)}=\frac{\sigma_{jm}}{\sigma_m^2}$, $\alpha_j=(\mu_j-r_f)-\beta(\mu_m-r_f)$, $\mathrm{cor}(\epsilon_j,r_m-r_f)=0$. CAPM ⇒ $\alpha_j=0\;\forall j$.
**Risk decomposition:** $\sigma_j^2=\underbrace{\beta_j^2\sigma_m^2}_{\text{market risk}}+\underbrace{\mathrm{var}(\epsilon)}_{\text{residual}}$. Only market risk is compensated.
**Security Market Line:** plot returns vs $r_f+\beta(\mu_m-r_f)$; under CAPM all assets lie on it.
**Trading deviations:** Jensen's alpha $\alpha=(\hat\mu_j-r_f)-\beta_j(\hat\mu_m-r_f)$ (long if $>0$); stock Sharpe $s_j=\frac{\hat\mu_j-r_f}{\hat\sigma_j}$ (long if $>m_{\mathrm{CML}}$).
**Assumptions:** identical information; all MV-optimizers (or Normal returns); markets in equilibrium.

### 5.4 CAPM as a pricing formula
Payoff $X$ in 1yr, price $P$, $r_X=X/P-1$, $\beta_X=\frac1P\frac{\mathrm{cov}(X,r_m)}{\sigma_m^2}$. Putting $\mu_X$ on the SML:
$$\boxed{\,P=\frac{E[X]}{1+r_f}-\frac{\mathrm{cov}(X,r_m)}{(1+r_f)\,\mathrm{var}(r_m)}(\mu_m-r_f)\,}$$
(Risk-neutral discounted value minus a risk premium for market covariance.)

---

## 6. Practical MVO: estimation error, constraints, robustness

### 6.1 Estimation
True $(\mu,V)$ unknown. From $N$ historical returns:
$$\mu^{(est)}=\frac1N\sum_{t=1}^N r(t),\qquad V^{(est)}=\frac1{N-1}\sum_{t=1}^N (r(t)-\mu^{(est)})(r(t)-\mu^{(est)})^\top.$$
$V$ has $d(d+1)/2$ independent parameters → never enough data. With $N=60$ months of simulated data the estimated mean can lie far from the truth (true mean only inside the 95% confidence ellipse w.p. 0.95). The **estimated, realized and true efficient frontiers** diverge dramatically — the frontier is extremely unstable.

### 6.2 Why estimation error is so damaging
Two identical assets, mean $\mu$, var $\sigma^2$, $\rho=0$; true optimum $x^\ast=(0.5,0.5)$. Suppose $\mu_1^{(est)}=\mu+\epsilon$, $\mu_2^{(est)}=\mu-\epsilon$ (average error zero). MVO **overweights asset 1** — precisely the wrong asset — and overweighting worsens realized performance; the more shorting allowed, the worse it gets.
> R. Michaud: *"Mean-variance results in error-maximizing, investment-irrelevant portfolios."*

### 6.3 Constraints & robustness (each stabilizes the frontier)
- **No-short-sales** $x\ge 0$: the "corner" of $\{x\ge0\}$ can make solutions *more* sensitive at the boundary, but caps damage.
- **Leverage constraint** $\sum_i|x_i|\le 2$.
- **Robust target return:** with confidence region $S_m$, require $\min_{\mu\in S_m}\mu^\top x\ge r$ (worst-case over the estimation set).

### 6.4 Improving parameter estimates
- **Shrinkage** (James–Stein 1961; Ledoit–Wolf 2004):
  $\mu_i^{(sh)}=\alpha\mu_i^{(est)}+(1-\alpha)\big(\tfrac1d\sum_j\mu_j^{(est)}\big)$ (shrink to grand mean);
  $V^{(sh)}=\alpha V^{(est)}+(1-\alpha)\big(\tfrac1d\sum_i\sigma_i^{2(est)}\big)I$ (shrink to scaled identity).
- **Black–Litterman:** blend market-implied (equilibrium) returns with subjective views.
- **Non-parametric / nearest-neighbor:** with current return $r$, let $S=\{t:\|r-r_t\|\le\beta\}$ and predict next return from $\{r_{t+1}:t\in S\}$.
- (Implicitly, **resampling/averaging** over many estimated frontiers — the source of the "many independent samples" instability picture.)

### 6.5 Negative exposures & leveraged ETFs
Shorts can boost return but have **unlimited downside** (price can grow without bound). For negative exposure with limited liability use **inverse ETFs**.
- ETF return compounds daily: $R^{(ETF)}_T=\prod_t(1+r_t)$; $\beta$-leveraged: $R^{(\beta\text{-ETF})}_T=\prod_t(1+\beta r_t)$.
- **Daily-compounding decay (worked):** index $100\to105$ ($r_1=5\%$) $\to100$ ($r_2=-4.76\%$). A 2× bull ETF: $100\to110\to110(1-9.52\%)=99.52$ → *loses money*. Inverse $-1\times$ ETF: $100\to95\to99.52$ → also loses. Round-trip flat index ⇒ leveraged ETF ends below start.
- **Static leveraged position** gross return $=\frac{\beta S_T-(\beta-1)S_0(1+r_f T)}{S_0}$.
- **$\beta$-ETF** approx: $\big(\frac{S_T}{S_0}\big)^\beta e^{(1-\beta)r_fT-fT-\frac{\beta^2-\beta}{2}\sum_i\ln^2(S_{i+1}/S_i)}$; the term $-\frac{\beta^2-\beta}{2}\sum_i\ln^2(\cdot)$ is the **volatility drag (short-vol)**. → leveraged ETFs are short-term instruments; they erode over long, volatile horizons.

### 6.6 Beyond variance: VaR and CVaR
Variance is symmetric and only right for elliptic distributions. Tail measures of loss $L$:
- **VaR:** $\mathrm{VaR}_p(L)=$ $p$-th quantile $=F_L^{-1}(p)$. Increasing in $p$.
- **CVaR (Expected Shortfall):** $\mathrm{CVaR}_p(L)=E[L\mid L\ge\mathrm{VaR}_p(L)]\ge\mathrm{VaR}_p(L)$.
- **Normal $L\sim N(\mu,\sigma^2)$:** $\mathrm{VaR}_p=\mu+\sigma\Phi^{-1}(p)$; $\mathrm{CVaR}_p=\mu+\sigma\frac{1}{1-p}\int_p^1\Phi^{-1}(\beta)d\beta$.
- **From samples** (order statistics $L_{(1)}\le\dots\le L_{(N)}$, $K_p=\lceil pN\rceil$): $\mathrm{VaR}_p\approx L_{(K_p)}$; $\mathrm{CVaR}_p\approx\frac{1}{(1-p)N}\sum_{k=K_p}^N L_{(k)}$.
- **Distribution matters (worked, $N=10{,}000$ at 95%):** for the Sharpe-optimal portfolio's loss $-r_x$: Normal → VaR 1.67, CVaR 3.15; multivariate $t_{\nu=12}$ → VaR 2.36, CVaR 4.58; mixture $0.75\,N(\mu,0.76^2V)+0.25\,N(\mu,1.5^2V)$ → VaR 1.93, CVaR 5.60. Fat tails inflate tail risk.
- **VaR:** captures tail, robust to outliers; but ignores beyond the quantile, incentivizes "tail stuffing", and is **not sub-additive** (diversification can raise it). **CVaR:** sees beyond the quantile, **sub-additive**, mean-CVaR optimization solves efficiently; but is an expectation so sensitive to outliers.

---

## 7. Behavioral / statistical biases in asset allocation

### 7.1 Is manager "skill" (alpha) real? — three perspectives
Model: a no-skill manager outperforms each year w.p. $p=1/2$, independent across years; $X\sim\mathrm{Bin}(n,p)$, $P(X\ge r)=\sum_{i=r}^n\binom{n}{i}p^i(1-p)^{n-i}$.
1. **Single manager:** 9/10 years. $P(X\ge9)=0.0107$ → looks like skill.
2. **$M$ managers, best one:** $P(\text{best}\ge9/10)=1-[1-P(Z_1)]^M$ with $P(Z_1)=0.0107$; for $M=20$ this is **0.1942** — no longer impressive. → ask the *right* question given your prior; $P(V)\to1$ as $M\to\infty$.
3. **Survivorship thought-experiment:** if no manager has skill but losers are fired each year, the *survivors'* track records look great — but it does not reflect true skill. **Survivorship bias.**

### 7.2 How should average returns be computed?
+20% then −10%: naive time-average = **5%**. But with dollars invested (\$1m then \$10m): the **dollar-weighted** return $=\frac{1\cdot20\%-10\cdot10\%}{11}=-7.27\%$. From the investor's view the dollar-weighted figure is correct (return on dollars invested). Hedge funds prefer to report time-averages. Compounded by liquidity: expected returns fall as dollars invested rise (capacity limits).

### 7.3 Other averaging / sampling biases
- **Counting children** by sampling people and counting siblings+1: biased upward (large families over-sampled) — a length/size bias, same flavor as the dollar-weighting problem.
- **Waiting times** sampled one person per hour: biased (length-biased sampling).
- **Survivorship in backtests:** equi-weighting *today's* top-20 S&P stocks and back-testing 20 years over-states performance (winners were chosen ex-post).
- **Football-game scam:** correct "predictions" for 10 weeks just because you are the surviving recipient; you are the survivor.
- **Data snooping:** normalizing the *entire* dataset (full-sample mean/variance) before train/test split leaks future info into the training set → the strategy looks great in-sample, fails live. Keep test data strictly held out.
- Also: dubious claims ("market never fell over any 20-year period"), implausible "25-sigma moves several days running" (Normal model badly underestimates tails), and retail structured products that backtest well but hide volatility/credit risk, fees, and have no secondary market.
- **Monty Hall:** switching doors wins w.p. 2/3 (intuition trap).

---

## 8. Liquidity, trading costs and optimal execution

### 8.1 Liquidity & cost functions
A **liquid** security trades quickly, with low price impact (slippage), in large quantities (deep book). Measures: volume-based (trading volume, turnover = volume/shares outstanding); cost-based (% bid-ask spread, VWAP, Loeb, **Kissell–Glantz**).
**Kissell–Glantz** cost of trading $Q$ shares:
$$c(Q)=\Big(a_1\big(\tfrac{100|Q|}{V}\big)^\beta+a_2\sigma+a_3\Big)\,|PQ|,$$
$V$ = avg daily volume, $\sigma$ = volatility, $P$ = pre-trade price. Estimated by regressing $\frac{c(Q_t)}{|P_tQ_t|}=a_1(\tfrac{100|Q_t|}{V_t})^\beta+a_2\sigma_t+a_3+\epsilon$.

### 8.2 Mean-variance-liquidity optimization
Current position $y$, choose new $x$:
$$\max_x\;\mu^\top x-\tfrac{\lambda}{2}x^\top Vx-\eta\, C(x,y)\quad\text{s.t. }\mathbf{1}^\top x=\mathbf{1}^\top y,$$
$$C(x,y)=\sum_i\Big[a_1\big(\tfrac{100}{P_iV_i}\big)^\beta|x_i-y_i|^{1+\beta}+(a_2\sigma_i+a_3)|x_i-y_i|\Big].$$
(Solved in **trading.xlsx**, "mean-variance-liquidity" sheet.)

**Lo–Petrov–Wierzbicki** ("It's 11pm — do you know where your liquidity is?"): normalized liquidity $\ell_i=\frac1T\sum_t\frac{\tilde\ell_{it}-\min\tilde\ell}{\max\tilde\ell-\min\tilde\ell}\in[0,1]$. Three formulations (start all-cash):
- **Liquidity-filtered:** $\max\mu^\top x-\tfrac\lambda2 x^\top Vx$ s.t. $\mathbf{1}^\top x=1$, $x_i=0$ for illiquid $\ell_k\ge\bar L$.
- **MV-liquidity objective:** $\max\mu^\top x-\tfrac\lambda2 x^\top Vx+\eta\sum_i\ell_i|x_i|$.
- **Liquidity-constrained:** $\max\mu^\top x-\tfrac\lambda2 x^\top Vx$ s.t. $\mathbf{1}^\top x=1$, $\frac{\sum_i\ell_i|x_i|}{\sum_j|x_j|}\ge\bar L$ (ratio form prevents short illiquid positions cancelling long ones).

### 8.3 Optimal execution of a single stock (Almgren–Chriss style)
Sell $X$ shares over $\le T$ trades. $n_k$ = shares in trade $k$, $X=\sum_k n_k$; holdings $x_k=X-\sum_{j\le k}n_j$ ($x_0=X$). Two impacts:
- **Temporary** (own fill only): $\hat S_k=S_k-h(n_k)$.
- **Permanent** (moves future price): $S_{k+1}=S_k+\sigma z_k-g(n_k)$, $z_k$ iid $N(0,1)$.

Total revenue (capture):
$$\sum_k\hat S_kn_k=S_1\sum_k n_k+\sigma\sum_k\tilde z_k x_k-\sum_k g(n_k)x_k-\sum_k h(n_k)n_k.$$
Decompose into **expected cost (slippage)** and **risk (variance)**:
$$\boxed{\,C(n)=\sum_k g(n_k)x_k+\sum_k h(n_k)n_k\,}\qquad\boxed{\,V(n)=\sigma^2\sum_k x_k^2\,}$$
**Optimal-execution frontier** (trade off cost vs. timing risk; $T$ also optimized):
$$\min_n\; C(n)+\rho\, V(n).$$
Typical impact functions: linear permanent $g(v)=\gamma v$; nonlinear temporary (Kissell–Glantz) $h(n)=\big(a_1(\tfrac{100|n|}{V})^\beta+a_2\sigma_i+a_3\big)\mathrm{sign}(n)$. (Solved in **trading.xlsx**, "execution" sheet.)
- Larger $\rho$ (more risk-averse) → trade faster (front-loaded, less timing risk, more impact cost); smaller $\rho$ → trade slowly (lower cost, more risk).

### 8.4 Portfolio execution & microstructure
Move portfolio $x\to y$ over $T$ trades while staying balanced; complications: **cross-asset price impact** (a trade in one asset moves others — often *negative*, smaller when more speculators, larger under belief dispersion or shared specialist); **market vs. limit orders** (market = immediate but worse price; limit = better price but uncertain fill, choice driven by patience, volatility, liquidity-provider discount); **lit vs. dark pools** (dark pools hide size/identity to avoid impact but with uncertain fill volume; ~11–12% of volume in some names, raising price-discovery and *winner's-curse* concerns; SEC 2009 transparency measures). Similar issues in HFT: where to place orders in the book and how to manage inventory risk.

---

## 9. Excel workbooks (companion computations)

- **mvo.xlsx** — four sheets:
  - `data`: asset means $\mu$, covariance $V$, $r_f$.
  - `risky frontier`: efficient frontier with risky assets only — solves $\sigma^2(r)=\min x^\top Vx$ via the matrix system $A^{-1}b$ of §3.2; traces $(\sigma,r)$.
  - `riskfree frontier`: Capital Market Line — line from $(0,r_f)$ tangent at $s^\ast$.
  - `Sharpe portfolio`: computes the tangency/Sharpe-optimal $s^\ast=\frac{V^{-1}\hat\mu}{\mathbf{1}^\top V^{-1}\hat\mu}$ and its Sharpe ratio (also the input portfolio for the VaR/CVaR experiments).
- **trading.xlsx** — three sheets:
  - `mean-variance-liquidity`: solves the §8.2 MVL rebalancing problem (Kissell–Glantz costs with $1+\beta$ power term) from current position $y$.
  - `execution`: solves the §8.3 single-stock optimal-execution problem $\min_n C(n)+\rho V(n)$ with linear permanent + Kissell–Glantz temporary impact; produces the optimal trading trajectory $\{n_k\},\{x_k\}$.
  - `Sheet3`: scratch / auxiliary.

---

*Image-only companion PDFs (Trading_Friction_Optimization, Modélisation_Liquidité_et_Exécution, Finance_Beyond_Theory, Finance_Quant_Insights_and_Pitfalls) and TM*.pdf transcripts duplicate the deck content captured above; no unique quantitative material was lost.*
