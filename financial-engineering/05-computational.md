# Course 5 — Computational Methods: Transform Pricing & Model Calibration

Columbia FE&RM, Course 5: *Computational Methods in Pricing and Model Calibration* (Prof. Ali Hirsa / RAF). Three pillars:
1. **Transform (Fourier/FFT) option pricing** via the characteristic function — Carr–Madan damped call + FFT inversion (Module 2).
2. **Model calibration** — fitting parameters $\Theta$ to a market option surface by minimizing an objective (RMSE); brute-force, Nelder–Mead, BFGS (Module 3).
3. **Interest rates & rate-model calibration** — zero-coupon bonds, LIBOR/swap rates, empirical regression, Vasicek & CIR short-rate models calibrated to LIBOR+swap curves (Modules 4–5).

---

## Key formulas — cheat sheet

| Quantity | Formula |
|---|---|
| Risk-neutral price (density form) | $V_0 = e^{-rT}\int_0^\infty \text{Payoff}(S_T)\,f(S_T\mid S_0)\,dS_T$ |
| Characteristic function | $\phi(u)=E[e^{iuX}]$, $X=\ln S_T$ |
| Damped (modified) call | $c_T(k)=e^{\alpha k}\,C_T(k)$, $\alpha>0$, $k=\ln K$ |
| Fourier transform of damped call | $\psi_T(v)=\int_{-\infty}^{\infty} e^{ivk}\,c_T(k)\,dk = \dfrac{e^{-rT}\,\phi_T\!\big(v-(\alpha+1)i\big)}{\alpha^2+\alpha-v^2+i(2\alpha+1)v}$ |
| Call recovery (inverse FT) | $C_T(k)=\dfrac{e^{-\alpha k}}{\pi}\int_0^\infty e^{-ivk}\,\psi_T(v)\,dv$ |
| FFT grid constraint | $\lambda\cdot\eta = \dfrac{2\pi}{N}$, $N=2^n$ |
| GBM CF | $\phi_{\text{GBM}}(u)=\exp\!\big(iu\mu - \tfrac12\sigma^2u^2T\big)$, $\mu=\ln S_0+(r-q-\tfrac{\sigma^2}{2})T$ |
| Calibration objective (RMSE) | $\text{RMSE}(\Theta)=\sqrt{\frac1N\sum_{i=1}^N(\hat V_i^\Theta - V_i)^2}$ |
| Weighted least squares | $\min_{\Theta\in O}\sum_i w_i(\hat V_i^\Theta - V_i)^2$ |
| BFGS update | quasi-Newton: $\Theta_{n+1}=\Theta_n - H_n^{-1}\nabla f(\Theta_n)$ (Hessian estimated) |
| Gradient descent | $\Theta_{n+1}=\Theta_n - \gamma\,\nabla\ell(\Theta_n)$ |
| Simple spot (LIBOR) rate | $L(t,T)=\frac{1}{T-t}\big(\frac{1}{P(t,T)}-1\big)$ |
| Simple forward rate | $F(t;S,T)=\frac{1}{T-S}\big(\frac{P(t,S)}{P(t,T)}-1\big)$ |
| Continuous spot rate | $R(t,T)=-\frac{\ln P(t,T)}{T-t}$ |
| Instantaneous fwd rate | $f(t,T)=-\frac{\partial \ln P(t,T)}{\partial T}$, $P(t,T)=\exp(-\int_t^T f(t,u)du)$ |
| Swap rate | $s(t,T)=\dfrac{1-P(t,T)}{\Delta\sum_{j=1}^n P(t,T_j)}$ |
| Short-rate ZCB | $P(t,T)=E_t^Q[e^{-\int_t^T r_s\,ds}]=e^{A(t,T)-B(t,T)r_t}$ |
| Vasicek SDE | $dr_t=\kappa(\theta-r_t)dt+\sigma\,dW_t$ |
| CIR SDE | $dr_t=\kappa(\theta-r_t)dt+\sigma\sqrt{r_t}\,dW_t$ |

---

## 1. Transform (Fourier / FFT) Option Pricing

### 1.1 Motivation — why transforms

Direct risk-neutral pricing requires the density $f(S_T\mid S_0)$:
$$V_0 = e^{-rT}\int_0^\infty \text{Payoff}(S_T)\,f(S_T\mid S_0)\,dS_T.$$
**Bottleneck:** for advanced models (Heston, Variance Gamma) the density is *not available in closed form*, making numerical integration over the density slow/infeasible. However, the **characteristic function** (Fourier transform of the density) $\phi(u)=E[e^{iuX}]$ with $X=\ln S_T$ *is* known analytically for these Lévy / stochastic-vol processes. So: move to Fourier space (CF known) → invert to get the price.

### 1.2 Carr–Madan: damped call transform

The transform of the raw call $C_T(k)$ (as a function of log-strike $k=\ln K$) does **not** converge because the call payoff does not → 0 as $k\to-\infty$ (deep ITM call value → $S_0$). Carr–Madan introduce an **exponential damping** factor $\alpha>0$:
$$c_T(k)=e^{\alpha k}\,C_T(k)\quad(\text{modified / damped call}).$$
Its Fourier transform is
$$\psi_T(v)=\int_{-\infty}^{\infty} e^{ivk}\,c_T(k)\,dk
=\frac{e^{-rT}\,\phi_T\big(v-(\alpha+1)i\big)}{\alpha^2+\alpha-v^2+i(2\alpha+1)v},$$
where $\phi_T$ is the CF of the log-price $\ln S_T$. The denominator equals $(\alpha+iv)(\alpha+1+iv)$ — exactly the term in the code below. Recover the call by inverse transform:
$$C_T(k)=\frac{e^{-\alpha k}}{\pi}\int_0^\infty e^{-ivk}\,\psi_T(v)\,dv.$$

**Damping $\alpha$ choice:** must make $E[S_T^{\alpha+1}]<\infty$. For calls $\alpha>0$ (optimal $\approx 1.5$ in the course example); for **puts**, use $\alpha<0$ (the same `genericFFT` prices puts by choosing $\alpha\in\{-1.01,-1.25,-1.5,-1.75,-2,-5\}$). Same machinery, only $\alpha$'s sign differs.

### 1.3 FFT discretization (Carr–Madan)

Discretize the inversion integral on a grid $v_j=\eta j$, $j=0,\dots,N-1$, $N=2^n$, and evaluate the price at log-strikes $k_m=\beta+\lambda m$. Using Simpson/trapezoid weights $w_j$, the sum is recognized as a DFT computable by **FFT**:
$$C_T(k_m)\approx \frac{e^{-\alpha k_m}}{\pi}\,\text{Re}\!\left(\sum_{j=0}^{N-1} e^{-i\beta\eta j}\,e^{-i\frac{2\pi}{N}jm}\,e^{-rT}\,\psi(\nu_j)\,w_j\right),$$
where $\psi(\nu_j)=\dfrac{\phi_T(\nu_j-(\alpha+1)i)}{(\alpha+i\nu_j)(\alpha+1+i\nu_j)}$.

**Grid constraint (Nyquist):** the integration step $\eta$ (Fourier space) and log-strike step $\lambda$ (price space) are coupled:
$$\boxed{\lambda\cdot\eta=\frac{2\pi}{N}.}$$
You cannot choose them independently: a *fine* Fourier grid (small $\eta$) forces a *coarse* strike spacing (large $\lambda$) and vice versa. Trapezoid weight: $w_0=\eta/2$, $w_j=\eta$ otherwise.

**Complexity:** direct sequential integration is $O(N^2)$; FFT is $O(N\log N)$ — pricing *thousands* of strikes simultaneously. Measured speedup in the course: explicit loop $\approx0.5$ s vs FFT $\approx0.0005$ s ≈ **1000×** — essential for real-time calibration / HFT.

**Practical settings:** $N=2^{12}=4096$ (best precision/speed balance), $\eta\approx0.25$, $\alpha\approx1.5$ (calls). $\beta$ aligns the strike grid: `beta = log(S0) - N*lda/2` centers it at-the-money, or `beta = log(K)` puts the desired strike as the *first* array element (`cT_km[0]`). **Golden rule:** never write your own DFT loop — use NumPy / FFTW.

### 1.4 Pricing pipeline (6 steps)

1. **Input:** CF params $\Theta$, $\alpha$, $N$, $\eta$.
2. **Build vector $x$:** $\psi(\nu_j)$ with Simpson/trapezoid weights.
3. **FFT:** `y = np.fft.fft(x)`.
4. **Output:** complex array `y`.
5. **Post-process:** multiply by $e^{-\alpha k_m}/\pi$, take real part.
6. **Price:** $C(k_m)$ for all strikes.

### 1.5 Characteristic functions of the three models

These are the actual CFs implemented in the course code (`generic_CF`). $u$ is the Fourier variable.

**Black–Merton–Scholes / GBM** (constant vol $\sigma$ — benchmark, FFT matches analytic to $<10^{-4}$):
$$\mu=\ln S_0+(r-q-\tfrac{\sigma^2}{2})T,\quad a=\sigma\sqrt T,\quad \phi(u)=\exp\!\big(iu\mu-\tfrac12 a^2u^2\big).$$

**Heston stochastic volatility** — parameters $\Theta=\{\kappa,\theta,\sigma,\rho,v_0\}$ (mean-reversion speed, long-run var, vol-of-vol, correlation, initial var). Captures the **volatility smile** (BMS gives a flat line). With $\text{tmp}=\kappa-i\rho\sigma u$ and $g=\sqrt{\sigma^2(u^2+iu)+\text{tmp}^2}$:
$$\ln\phi = \underbrace{\frac{\kappa\theta T\,\text{tmp}}{\sigma^2}+iuTr+iu\ln S_0}_{\text{numer1}}
-\underbrace{\frac{2\kappa\theta}{\sigma^2}\ln\!\Big(\cosh\tfrac{gT}{2}+\tfrac{\text{tmp}}{g}\sinh\tfrac{gT}{2}\Big)}_{\text{log\_denum1}}
-\underbrace{\frac{(u^2+iu)v_0}{g\coth(gT/2)+\text{tmp}}}_{\text{tmp2}}.$$

**Variance Gamma (VG)** — pure-jump Lévy, params $\{\sigma,\nu,\theta\}$ (vol, kurtosis/variance-rate of the gamma clock, skew). Fat tails + sharp peak; as $\nu\to0,\theta\to0$ it → BMS. For $\nu\neq0$:
$$\mu=\ln S_0+\Big(r-q+\tfrac{1}{\nu}\ln(1-\theta\nu-\tfrac12\sigma^2\nu)\Big)T,\quad
\phi(u)=e^{iu\mu}\big(1-i\nu\theta u+\tfrac12\nu\sigma^2u^2\big)^{-T/\nu}.$$
For $\nu=0$: $\phi(u)=e^{iu\mu}\exp\big((i\theta u-\tfrac12\sigma^2u^2)T\big)$, $\mu=\ln S_0+(r-q-\theta-\tfrac12\sigma^2)T$.

### 1.6 Reference implementation — `genericFFT` (core code)

The actual course code (`modulesForCalibration.py`). This is the load-bearing FFT pricer.

```python
import numpy as np, cmath, math

def generic_CF(u, params, S0, r, q, T, model):
    if model == 'GBM':
        sig = params[0]
        mu  = np.log(S0) + (r - q - sig**2/2)*T
        a   = sig*np.sqrt(T)
        phi = np.exp(1j*mu*u - (a*u)**2/2)
    elif model == 'Heston':
        kappa, theta, sigma, rho, v0 = params
        tmp = (kappa - 1j*rho*sigma*u)
        g   = np.sqrt((sigma**2)*(u**2 + 1j*u) + tmp**2)
        pow1   = 2*kappa*theta/(sigma**2)
        numer1 = (kappa*theta*T*tmp)/(sigma**2) + 1j*u*T*r + 1j*u*math.log(S0)
        log_denum1 = pow1 * np.log(np.cosh(g*T/2) + (tmp/g)*np.sinh(g*T/2))
        tmp2 = ((u*u + 1j*u)*v0)/(g/np.tanh(g*T/2) + tmp)
        phi  = np.exp(numer1 - log_denum1 - tmp2)
    elif model == 'VG':
        sigma, nu, theta = params
        if nu == 0:
            mu  = math.log(S0) + (r-q-theta-0.5*sigma**2)*T
            phi = math.exp(1j*u*mu)*math.exp((1j*theta*u-0.5*sigma**2*u**2)*T)
        else:
            mu  = math.log(S0) + (r-q + math.log(1-theta*nu-0.5*sigma**2*nu)/nu)*T
            phi = cmath.exp(1j*u*mu)*((1-1j*nu*theta*u+0.5*nu*sigma**2*u**2)**(-T/nu))
    return phi

def genericFFT(params, S0, K, r, q, T, alpha, eta, n, model):
    N   = 2**n
    lda = (2*np.pi/N)/eta          # log-strike step  (lambda*eta = 2*pi/N)
    beta = np.log(K)               # so the wanted strike is km[0]
    km, xX = np.zeros(N), np.zeros(N)
    df  = math.exp(-r*T)           # discount factor
    nuJ = np.arange(N)*eta
    # damped, discounted CF / [(alpha+i nu)(alpha+1+i nu)]
    psi_nuJ = generic_CF(nuJ-(alpha+1)*1j, params, S0, r, q, T, model) \
              / ((alpha + 1j*nuJ)*(alpha + 1 + 1j*nuJ))
    for j in range(N):
        km[j] = beta + j*lda
        wJ = (eta/2) if j == 0 else eta        # trapezoid weights
        xX[j] = cmath.exp(-1j*beta*nuJ[j])*df*psi_nuJ[j]*wJ
    yY = np.fft.fft(xX)
    cT_km = np.zeros(N)
    for i in range(N):
        cT_km[i] = (math.exp(-alpha*km[i])/math.pi)*np.real(yY[i])
    return km, cT_km
```

Vectorized variant (from the notebook solutions) replaces the loop with:
```python
km = beta + lda*np.arange(N);  w = eta*np.ones(N); w[0] = eta/2
xX = np.exp(-1j*beta*nuJ)*df*psi_nuJ*w
cT_km = (np.exp(-alpha*km)/np.pi)*np.real(np.fft.fft(xX))
```

### 1.7 Numerical-integration alternative (lognormal density, for BMS puts)

When the density *is* known (BMS), price puts directly by quadrature (course quiz code). Step $\eta=K/N$, grid $S_j=j\eta$ (set $S_0\to10^{-8}$ to avoid div-by-zero), trapezoid weights $w_0=\eta/2$:
```python
def logNormal(S, r, q, sig, S0, T):  # lognormal density of S_T
    return np.exp(-0.5*((np.log(S/S0)-(r-q-sig**2/2)*T)/(sig*np.sqrt(T)))**2) \
           /(sig*S*np.sqrt(2*np.pi*T))
def numerical_integral_put(r,q,S0,K,sig,T,N):
    df=np.exp(-r*T); eta=K/N; S=np.arange(0,N)*eta; S[0]=1e-8
    w=np.ones(N)*eta; w[0]=eta/2
    priceP = df*np.sum((K-S)*logNormal(S,r,q,sig,S0,T)*w)
    return eta, priceP
```
Sample numbers (FFT puts, $S_0=100,K=80,r=.05,q=.01,T=1$): BMS($\sigma=.3$) $\approx2.7080$; Heston($\kappa=2,\theta=.05,\lambda=.3,\rho=-.7,v_0=.04$) $\approx1.3412$; VG($\sigma=.3,\nu=.5,\theta=-.4$) $\approx5.3137$.

---

## 2. Model Calibration

### 2.1 The option surface and the calibration problem

- **Model price:** $\hat V(S_0,K,r,q,T;\Theta)$, $\Theta=$ parameter set. Pick a model + $\Theta$ → for a fixed maturity it returns call prices across strikes; over $(K,T)$ this is the **option surface**.
- **Market prices** have bid/ask spread (buy at ask, sell at bid); tight for liquid/high-volume contracts, wide for illiquid. Use **mid** = (bid+ask)/2. Out-of-the-money options are the liquid reference.
- **Implied volatility (one price ↔ one parameter):** find the BMS $\sigma$ matching the market price for a given $(K,T)$. Market convention. The whole surface plotted as IV vs $(K,T)$ shows the smile/skew.
- **Calibration (definition):** adjust model parameters $\Theta$ so model prices are compatible with market prices. A perfect fit of the entire surface is generally impossible; the goal is a *decent* fit.

### 2.2 Objective / cost function

With market price $V_i$ and model price $\hat V_i^\Theta$ for instrument $i$:
$$\min_{\Theta\in O}\sum_{i=1}^n H(\hat V_i^\Theta - V_i),$$
common $H$ forms (weight $w_i$, power $p$): absolute $w_i|\hat V_i^\Theta-V_i|^p$, **relative** $w_i\big|\tfrac{\hat V_i^\Theta-V_i}{V_i}\big|^p$, log $w_i|\ln\hat V_i^\Theta-\ln V_i|^p$.

**Weighted least squares** (most common): $\min_\Theta\sum_i w_i(\hat V_i^\Theta-V_i)^2$.
**RMSE** (the course's choice — dollars, directly interpretable):
$$\text{RMSE}(\Theta)=\sqrt{\frac{1}{n}\sum_{i=1}^n(\hat V_i^\Theta-V_i)^2}.$$

**Parameter counting:** with $n_I$ instruments and $n_P$ parameters — $n_I=n_P$ (e.g. BMS implied vol), $n_I>n_P$ **under-parameterized** (typical; no unique $\Theta$), $n_I<n_P$ over-parameterized. Apple example: $N\approx300$ instruments vs 5 Heston params → under-parameterized, perfect fit impossible.

### 2.3 The calibration "recipe" (5 ingredients)

1. **Dataset** — all prices / only OTM calls+puts / user discretion.
2. **Model** — e.g. Heston $\Theta=\{\kappa,\theta,\sigma,\rho,v_0\}$ to fit (part of) the surface.
3. **Objective/cost function** — e.g. RMSE.
4. **Initial parameter set $\Theta_0$** — *critical*.
5. **Optimization routine**.

**Finding $\Theta_0$ / error-surface inspection:** interpolate between two arbitrary sets, $\Theta=\alpha\Theta_1+(1-\alpha)\Theta_2$, $\alpha\in(0,1)$, and plot the error along the line. A convex valley ⇒ promising start; an irregular/multi-modal landscape ⇒ risk of local-minimum traps. Course Apple start: $\Theta_0=(2.3,0.046,0.0825,-0.53,0.054)$, RMSE $=0.8659$.

### 2.4 Optimization routines (three families)

| Routine | Type | Pros | Cons |
|---|---|---|---|
| **Brute-force / grid** | exhaustive | global min regardless of shape | curse of dimensionality; $10^5$ evals for 5 params; infeasible high-dim |
| **Nelder–Mead simplex** | gradient-free | robust, no derivatives | slower; no global guarantee |
| **BFGS** (Broyden–Fletcher–Goldfarb–Shanno) | gradient-based, quasi-Newton | fast convergence on smooth/convex, precise | can stick at saddle points; needs smooth objective + good $\Theta_0$; no global guarantee |

**Brute-force grid algorithm** (Apple, Heston):
```
eMin = 1e12
for κ  in 1.8,…,2.8:
 for θ  in 0.036,…,0.056:
  for λ  in 0.0725,…,0.0925:
   for ρ  in -0.63,…,-0.43:
    for v0 in 0.044,…,0.064:
       e = RMSE(Θ); if e < eMin: eMin=e; Θ*=Θ
```
Result: $\Theta^*_{\text{Brute}}=(1.8,0.046,0.0925,-0.63,0.044)$, RMSE $=0.3670$.

**BFGS** — quasi-Newton, iterative; uses gradient + an estimated inverse Hessian, update $\Theta_{n+1}=\Theta_n-H_n^{-1}\nabla f$. `scipy.optimize.fmin_bfgs` in Python. Result: $\Theta^*_{\text{BFGS}}=(3.6941,0.0478,0.6059,-0.2186,0.0422)$, RMSE $=0.2661$.
**Nelder–Mead** result: $\Theta^*_{\text{NM}}=(1.9524,0.0469,0.1159,-0.7406,0.0397)$.

**Apple case-study scoreboard** (RMSE, lower=better): initial $0.86$ → Brute-force $0.36$ → Nelder–Mead $0.28$ → **BFGS $0.26$ (winner: fast + precise, conditioned on a good start)**. Stability check: interpolate $\Theta=\alpha\Theta^*_{\text{Brute}}+(1-\alpha)\Theta^*_{\text{BFGS}}$ over $\alpha\in(-0.5,1.5)$ — a single convex valley confirms BFGS refined the coarse grid solution without falling into an unstable local min.

### 2.5 Calibration objective implementation — `eValue`

Course objective function: loops over the $(T,K)$ grid, prices each cell via `genericFFT`, returns RMSE.
```python
def eValue(params, *args):
    marketPrices, maturities, strikes, r, q, S0, alpha, eta, n, model = args
    lenT, lenK = len(maturities), len(strikes)
    modelPrices = np.zeros((lenT, lenK)); count = 0; mae = 0
    for i in range(lenT):
        for j in range(lenK):
            count += 1; T = maturities[i]; K = strikes[j]
            km, cT_km = genericFFT(params, S0, K, r, q, T, alpha, eta, n, model)
            modelPrices[i,j] = cT_km[0]          # price at first log-strike = K
            mae += (marketPrices[i,j] - modelPrices[i,j])**2
    return math.sqrt(mae/count)                  # RMSE
```
**Parameter constraints via periodic linear extension** (`paramMapping`): keeps unconstrained-optimizer iterates inside admissible ranges by reflecting them back — e.g. Heston `kappa∈[0.1,20]`, `theta∈[0.001,0.4]`, `sigma∈[0.01,0.6]`, `rho∈[-1,1]`, `v0∈[0.005,0.25]`. This lets a box-free optimizer (BFGS/Nelder–Mead) respect bounds without explicit constraints.

### 2.6 Calibration driver — BFGS notebook (`exampleCalibration_BFGS.ipynb`)

```python
import readPlotOptionSurface, modulesForCalibration as mfc
from scipy.optimize import fmin_bfgs

alpha, eta, n, model = 1.5, 0.2, 12, 'Heston'
r, q, S0 = 0.0245, 0.005, 190.3                 # Apple, mkt data
maturities, strikes, callPrices = readPlotOptionSurface.readNPlot()
marketPrices = callPrices; maturities_years = maturities/365.0

params = [2.3, 0.046, 0.0825, -0.53, 0.054]      # Θ0
arg = (marketPrices, maturities_years, strikes, r, q, S0, alpha, eta, n, model)

[xopt, fopt, gopt, Bopt, fcalls, gcalls, flag] = fmin_bfgs(
    mfc.eValue, params, args=arg, fprime=None,   # numerical gradient
    callback=callbackF, maxiter=20, full_output=True, retall=False)
# xopt = Θ*  ;  fopt = optimal RMSE
```
- `fprime=None` ⇒ SciPy estimates the gradient by finite differences (no analytic Jacobian of the FFT pricer needed).
- The surface is read from `data_apple.xlsx`: Mid=(Bid+Ask)/2, strikes `arange(170,212.5,2.5)`, missing strikes linearly interpolated (`scipy.interpolate.interp1d`, `fill_value=0`) per maturity (`readPlotOptionSurface.py: readNPlot`).
- After calibration, model prices are recomputed on the grid and overlaid on market prices (Market vs. Model plot).

**Files:** `modulesForCalibration.py` (defines `myRange`, `paramMapping`, `eValue`, `generic_CF`, `genericFFT`); `readPlotOptionSurface.py` (`readNPlot` builds & plots the call surface); `exampleCalibration_BFGS.ipynb` (driver).

---

## 3. Interest Rates & Interest-Rate Instruments

### 3.1 Zero-coupon bonds and the building block $P(t,T)$

$P(t,T)$ = price at $t$ to receive \$1 at $T$; $P(t,t)=1$. Interest-rate derivatives depend on the *entire yield curve*, not a single underlying — models capturing curve evolution are **term-structure models**. $P(t,T)$ is bootstrapped/calibrated from liquid instruments: **LIBOR, swap rates, futures, caps/floors, swaptions**.

### 3.2 Forwards, FRAs and the rate zoo (all functions of $P(t,T)$)

A **forward contract** = obligation to buy/sell at a set forward price on a known date (linear payoff). A **Forward Rate Agreement (FRA)** locks a simple forward rate over $[S,T]$ (times $t<S<T$), payoff at $T$ is $1+(T-S)F(t;S,T)$ ≡ a forward \$1 investment at $S$. Replication (zero net cost at $t$): sell 1 $S$-bond, buy $\tfrac{P(t,S)}{P(t,T)}$ $T$-bonds ⇒ at $S$ pay \$1, at $T$ receive $\tfrac{P(t,S)}{P(t,T)}$. No-arbitrage ⇒

| Rate | Definition |
|---|---|
| Simple forward | $1+(T-S)F(t;S,T)=\dfrac{P(t,S)}{P(t,T)}\Rightarrow F(t;S,T)=\dfrac{1}{T-S}\!\big(\tfrac{P(t,S)}{P(t,T)}-1\big)$ |
| Simple spot (→ LIBOR) | $F(t,T)=\dfrac{1}{T-t}\big(\tfrac{1}{P(t,T)}-1\big)$ |
| Cont. compounded forward | $R(t;S,T)=-\dfrac{\ln P(t,T)-\ln P(t,S)}{T-S}$ |
| Cont. compounded spot | $R(t,T)=-\dfrac{\ln P(t,T)}{T-t}$ |
| Instantaneous forward | $f(t,T)=-\dfrac{\partial\ln P(t,T)}{\partial T}$, $P(t,T)=\exp(-\int_t^T f(t,u)du)$ |
| Instantaneous short rate | $r(t)=f(t,t)=\lim_{T\downarrow t}R(t,T)$ |

**Futures** vs forwards: futures are exchange-traded, standardized, cleared, and **marked-to-market daily** (re-settled each day); forwards are OTC/custom. Eurodollar futures reference LIBOR.

**LIBOR** = average interbank lending rate; 7 maturities (overnight–12M), 5 currencies, simply compounded, ACT/360, spot starting 2 business days forward:
$$L(t,T)=\frac{100}{T-t}\Big(\frac{1}{P(t,T)}-1\Big)\quad(\text{=100×simple spot rate}).$$

### 3.3 Swaps

An **interest-rate swap** exchanges fixed-leg vs floating-leg cashflows on a notional (principal not exchanged — "notional principal"; only interest differentials exchanged). Motivation: **comparative advantage** — A and B each borrow where they have an edge (fixed vs floating), then swap legs via a dealer (who earns a fee), both securing lower effective rates. A swap = portfolio of FRAs, each with fixed rate = swap rate $s$.

US vanilla conventions: terms 2/5/10/30y. **Fixed leg** = 2 semi-annual payments/yr; **floating leg** = 4 quarterly payments/yr at 3M LIBOR, resetting quarterly. Example (2y, \$1,000,000 notional, fixed 2.9894%): fixed payment $=10^6\times\tfrac{2.9894}{100}\times\tfrac12=\$14{,}947$; floating at 3M LIBOR 2.4012% $=10^6\times\tfrac{2.4012}{100}\times\tfrac14=\$6{,}003$.

**Swap rate as a function of $P(t,T)$** (no-arb: rate making the contract worth 0 at inception; $T_n=T$, $\Delta$ = accrual length):
$$s(t,T)=\frac{1-P(t,T)}{\Delta\sum_{j=1}^n P(t,T_j)}\quad(\times100\text{ in \%}).$$
$\Delta=\tfrac12$ for US semi-annual. The swap rate is a weighted average of forward rates derived entirely from the $P(t,T)$ curve.

### 3.4 Empirical / data-driven analysis (regression on rates)

**Cross-correlation findings** (LIBOR & swap time series, Jan-2014→Oct-2018): adjacent maturities are highly correlated (1M vs 2M LIBOR ≈ 0.77; swap 2y–30y block 0.64–0.99), but distant maturities **decorrelate** (1M vs 12M LIBOR scatter is pure noise) → you **cannot** linearly extrapolate the short end to the long end. Correlations are **not static** — they depend on the monetary-policy regime (stable-rates 2014–16 vs Fed-hiking 2016–18). QQ-plots: weak linearity in calm regimes, strong linearity under directional Fed moves.

**Linear regression of one rate on others** (synthetic reconstruction / arbitrage signal). Conjectured relations, e.g. $\text{ussw5}=a_0+a_1\text{ussw2}$, or multivariate $\text{ussw30}=a_0+a_1\text{ussw2}+a_2\text{ussw5}+a_3\text{ussw10}$, written $\hat y_t=f(x_t;\Theta)=a_0+\sum_k a_k x_{k,t}$. Minimize
$$\ell(\Theta)=\frac{1}{2T}\sum_{t=1}^T(\hat y_t-y_t)^2,\qquad R^2=\frac{\sum_t(\hat y_t-\bar y)^2}{\sum_t(y_t-\bar y)^2}.$$
Solvable by **OLS** (`sklearn.linear_model.LinearRegression`), **Nelder–Mead**, or **gradient descent**.
```python
from sklearn import linear_model
df = pandas.read_csv('swapLiborData.csv')
xX, yY = df.iloc[:,6:12], df.iloc[:,12:13]
regr = linear_model.LinearRegression(); regr.fit(xX, yY)
B, b0 = regr.coef_, regr.intercept_
```
Regime-dependence of fit quality: 5y-on-2y gives $R^2\approx4\%$ (2014–16, flat regime) but $R^2\approx97\%$ (2016–18) and $77\%$ over the full window; 30y-on-15y $R^2\approx99.5\%$; multivariate 30y-on-(2,5,10y) $R^2\approx99.4\%$. ⇒ Fit depends drastically on market regime (steepening vs flattening).

**Gradient descent** for $\ell(\Theta)$ — vanilla (batch) uses the whole dataset, mini-batch uses a subset $D$ per iteration:
$$\Theta_{n+1}=\Theta_n-\gamma\nabla\ell(\Theta_n),\quad \frac{\partial\ell}{\partial\theta_i}=\frac1T\sum_j(\hat y_t-y_t)\frac{\partial\hat y_t}{\partial\theta_i}.$$
**Learning rate $\gamma$ is decisive:** too small ⇒ premature stop / non-convergence (e.g. $\gamma=0.005$); optimal ($\gamma\approx0.1$) ⇒ smooth convergence to $\Theta^*=(1.036,0.628)$; too large ($\gamma\approx0.33$) ⇒ oscillation/divergence. Direction of steepest descent is perpendicular to the loss contours; steepness indicates proximity to the minimum (flatter near it ⇒ take smaller steps).

### 3.5 Short-rate models: Vasicek

Three IR-modeling frameworks: **short-rate models** (model instantaneous $r_t$), **HJM** (model instantaneous forward $f(t,T)$), **market models** (LMM models observable LIBOR; SMM models swap rates). Short-rate models remain popular for **parsimony** when the *level* (not the term-structure shape) matters most. In all short-rate models:
$$P(t,T)=E_t^Q\Big[e^{-\int_t^T r_s\,ds}\Big].$$

**Vasicek** — Ornstein–Uhlenbeck, Gaussian, mean-reverting (simplest, like BMS for rates):
$$dr_t=\kappa(\theta-r_t)dt+\sigma\,dW_t,$$
$\kappa$ = mean-reversion speed, $\theta$ = long-run mean, $\sigma$ = (constant) vol. Exact solution:
$$r_t=e^{-\kappa t}r_0+\theta(1-e^{-\kappa t})+\sigma e^{-\kappa t}\!\int_0^t e^{\kappa s}dW_s.$$
**Affine ZCB** (closed form): $P(t,T)=e^{A(t,T)-B(t,T)r_t}$ with
$$B(t,T)=\frac{1-e^{-\kappa(T-t)}}{\kappa},\quad
A(t,T)=\Big(\theta-\frac{\sigma^2}{2\kappa^2}\Big)\big[B(t,T)-(T-t)\big]-\frac{\sigma^2}{4\kappa}B^2(t,T).$$
Time-homogeneous (depends only on $T-t$). ZCB SDE: $dP(t,T)=r_tP\,dt-P\,B(t,T)\sigma\,dW_t$. **Drawback:** Gaussian ⇒ $r_t$ can go **negative**.

### 3.6 Short-rate models: CIR

**Cox–Ingersoll–Ross** adds a $\sqrt{r_t}$ (Feller / square-root) diffusion so vol vanishes as $r_t\to0$, forcing $r_t\ge0$:
$$dr_t=\kappa(\theta-r_t)dt+\sigma\sqrt{r_t}\,dW_t.$$
Non-Gaussian, so pricing is harder, but ZCB still closed-form: $P(t,T)=e^{A(t,T)-B(t,T)r_t}$ with $\gamma=\sqrt{\kappa^2+2\sigma^2}$ and
$$A(t,T)=\frac{2\kappa\theta}{\sigma^2}\ln\!\Bigg(\frac{e^{\kappa(T-t)/2}}{\cosh\!\tfrac{\gamma(T-t)}{2}+\tfrac{\kappa}{\gamma}\sinh\!\tfrac{\gamma(T-t)}{2}}\Bigg),\quad
B(t,T)=\frac{2}{\kappa+\gamma\coth\!\tfrac{\gamma(T-t)}{2}}.$$
Derivatives generally need Monte Carlo / PDE. Both Vasicek and CIR are 4-parameter, $\Theta=\{\kappa,\theta,\sigma,r_0\}$ — too parsimonious to fit the whole term structure exactly.

### 3.7 Calibrating Vasicek/CIR to LIBOR + swap rates

Four free params ⇒ no perfect fit; calibrate by minimizing the **sum of squared relative errors (SSRE)** over LIBOR ($I$) and swap ($J$) quotes:
$$\text{SSRE}=\sum_{i=1}^I\!\left(\frac{L^{(i)}_{\text{MODEL}}-L^{(i)}_{\text{MARKET}}}{L^{(i)}_{\text{MARKET}}}\right)^2
+\sum_{j=1}^J\!\left(\frac{S^{(j)}_{\text{MODEL}}-S^{(j)}_{\text{MARKET}}}{S^{(j)}_{\text{MARKET}}}\right)^2.$$
Model rates $L_{\text{MODEL}},S_{\text{MODEL}}$ come from the model $P(t,T)$ via the simple-spot/swap-rate formulas; minimize over $\Theta$ (e.g. Nelder–Mead / gradient methods, same optimization toolbox as §2.4). Course calibration dates: 14-Dec-2017 (low-rate) and 11-Oct-2018 (higher-rate), using LIBOR (1,2,3,6,12M) + swap (2,3,5,7,10,15,30y) quotes.

**Empirical conclusions:** Vasicek anchors well at the extremes (short end & 30y) but its single factor is too rigid to capture intermediate convexity. CIR's positivity constraint is theoretically vital (no negative rates) but, in a positive-rate environment, adds **marginal** pricing impact — Vasicek and CIR ZCB-price curves are visually indistinguishable on the 2018 data; the single-factor rigidity remains the dominant limitation for both. Overall arc: **Data (regression on residuals)** → **Optimization (gradient descent, the bridge between empirics and theory)** → **SDE (Vasicek/CIR)** → **calibrated curve / pricing**, trading off tractability (Vasicek) vs theoretical coherence (CIR).
