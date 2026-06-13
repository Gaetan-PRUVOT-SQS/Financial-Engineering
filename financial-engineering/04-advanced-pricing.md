# Course 4 — Advanced Derivative Pricing: Vol Surface, Structured Credit & Real Options

Columbia FE&RM "Advanced Topics in Derivative Pricing" (Haugh & Iyengar). This file is self-contained: it carries the actual definitions, formulas, algorithms and worked numbers from the lectures. Topics: Black–Scholes review and the Greeks, the implied volatility surface (skew/term structure, dynamics, no-arbitrage, what it tells us via Breeden–Litzenberger), structured credit (CDO tranches, one-factor Gaussian copula, tranche-loss pricing), and real options (Simplico gold-mine binomial valuation).

---

## Key formulas cheat-sheet

| Quantity | Formula |
|---|---|
| BS call | $C_0 = S_0 e^{-cT}N(d_1) - K e^{-rT}N(d_2)$ |
| $d_1,\,d_2$ | $d_1=\dfrac{\ln(S_0/K)+(r-c+\sigma^2/2)T}{\sigma\sqrt T}$, $\;d_2=d_1-\sigma\sqrt T$ |
| Put–call parity | $P_0 + S_0 e^{-cT} = C_0 + K e^{-rT}$ |
| GBM under $Q$ | $S_t=S_0 e^{(r-c-\sigma^2/2)t+\sigma W_t}$ |
| Delta (call) | $\partial C/\partial S = e^{-cT}N(d_1)$; put $= e^{-cT}(N(d_1)-1)$ |
| Gamma | $\partial^2C/\partial S^2 = \dfrac{e^{-cT}\varphi(d_1)}{\sigma S\sqrt T}$ |
| Vega | $\partial C/\partial\sigma = e^{-cT}S\sqrt T\,\varphi(d_1)$ |
| Theta (call) | $-e^{-cT}S\varphi(d_1)\dfrac{\sigma}{2\sqrt T} + c e^{-cT}S N(d_1) - rKe^{-rT}N(d_2)$ |
| Delta–gamma–vega P&L | $\text{P\&L}\approx \delta\,\Delta S + \tfrac12\Gamma(\Delta S)^2 + \text{vega}\,\Delta\sigma$ |
| Vol surface (implicit def) | $C_{\text{mkt}}(S,K,T)=C_{BS}(S,T,r,c,K,\sigma(K,T))$ |
| Risk-neutral density (Breeden–Litzenberger) | $f_T(K)=e^{rT}\,\partial^2 C/\partial K^2$ |
| Digital (cash-or-nothing) | $D(K,T)=-\partial C_{\text{mkt}}/\partial K = -\partial C_{BS}/\partial K - \text{vega}\times\text{skew}$ |
| Risk-neutral binomial prob | $q_u=(R-d)/(u-d)$, $q_d=1-q_u$ |
| Gaussian copula factor | $X_i = a_i M + \sqrt{1-a_i^2}\,Z_i$, $\;\text{Corr}(X_i,X_j)=a_i a_j$ |
| Default threshold | $\bar x_i(t)=\Phi^{-1}(q_i(t))$ |
| Conditional default prob | $q_i(t\,|\,M)=\Phi\!\big((\bar x_i(t)-a_iM)/\sqrt{1-a_i^2}\big)$ |
| Tranche loss function | $TL^{L,U}_t(l)=\max\{\min\{lA(1-R),U\}-L,\,0\}$ |
| Expected tranche loss | $E^Q_0[TL^{L,U}_t]=\sum_{l=0}^N TL^{L,U}_t(l)\,p(l,t)$ |
| Simplico lease recursion | $V_t(s)=\max\{s-C,0\}\,G + \dfrac{qV_{t+1}(us)+(1-q)V_{t+1}(ds)}{1+r}$ |

---

## 1. Binomial & Black–Scholes review

### 1.1 Risk-neutral pricing in the binomial lattice
For a European option in an $n$-period binomial tree with up/down factors $u,d$ and gross risk-free return $R$ per period:
$$q_u=\frac{R-d}{u-d},\qquad q_d=1-q_u,\qquad C_0=\frac{1}{R^{\,n}}\,E^Q_0\big[\max(S_T-K,0)\big].$$
Backward induction: at each node $C=\frac1R[q_u C_{\text{up}}+q_d C_{\text{down}}]$. Risk-neutral pricing avoids computing the price at every node — discount the terminal payoff under the binomial probabilities $\binom{n}{k}q^k(1-q)^{n-k}$.

**Worked example (deck):** $T=3$ periods, $K=100$, $R=1.01$, lattice $S_0=100$ with $u$-moves to $107,114.49,122.50$. Terminal call payoffs $(22.5,7,0,0)$. One step back: $15.48=\frac{1}{1.01}[q_u\cdot 22.5 + q_d\cdot 7]$. Root price $C_0\approx 6.57$.

### 1.2 Self-financing replication
A trading strategy $\theta_t=(x_t,y_t)$ ($x_t$ shares, $y_t$ cash units, $B_t=R^t$) with value process
$$V_t=\begin{cases}x_1S_0+y_1B_0 & t=0\\ x_tS_t+y_tB_t & t\ge1.\end{cases}$$
**Self-financing:** $V_t = x_{t+1}S_t+y_{t+1}B_t$ (value just before trading = value just after). Then $V_{t+1}-V_t=x_{t+1}(S_{t+1}-S_t)+y_{t+1}(B_{t+1}-B_t)$ — changes come only from capital gains/losses. A self-financing strategy that replicates the payoff = **dynamic replication**; its initial cost equals the option price (else arbitrage). At every node the option value equals the replicating-portfolio value.

### 1.3 Black–Scholes model and formula
Assumptions: continuously-compounded rate $r$; GBM $S_t=S_0 e^{(\mu-\sigma^2/2)t+\sigma W_t}$; dividend yield $c$; continuous frictionless trading, short-selling allowed. As $n\to\infty$ the binomial converges to GBM under $Q$: $S_t=S_0 e^{(r-c-\sigma^2/2)t+\sigma W_t}$, giving
$$C_0=S_0 e^{-cT}N(d_1)-Ke^{-rT}N(d_2),\quad d_1=\frac{\ln(S_0/K)+(r-c+\sigma^2/2)T}{\sigma\sqrt T},\ d_2=d_1-\sigma\sqrt T.$$
Note $\mu$ does **not** appear (just as $p$ never enters binomial pricing). Put via parity $P_0+S_0e^{-cT}=C_0+Ke^{-rT}$. Everybody knows GBM is **not** a good approximation to security prices.

---

## 2. The Greeks and hedging in practice

The "Greeks" are partial derivatives of the price w.r.t. model parameters ($\varphi$ = standard-normal PDF):

- **Delta** $=\partial C/\partial S = e^{-cT}N(d_1)$ (call); put delta $=e^{-cT}N(d_1)-e^{-cT}=e^{-cT}(N(d_1)-1)$. Sensitivity to underlying.
- **Gamma** $=\partial^2C/\partial S^2 = \dfrac{e^{-cT}\varphi(d_1)}{\sigma S\sqrt T}$. Sensitivity of delta to $S$. Always positive for European options (convexity). Same for call and put (from parity, the linear $S$, $K$ terms have zero second derivative).
- **Vega** $=\partial C/\partial\sigma = e^{-cT}S\sqrt T\,\varphi(d_1)$. Sensitivity to $\sigma$. (Conceptually inconsistent with BS, which assumes $\sigma$ constant — yet traders use it.) Vega is largest for ATM options and grows with time-to-maturity.
- **Theta** $=-\partial C/\partial T = -e^{-cT}S\varphi(d_1)\dfrac{\sigma}{2\sqrt T}+c e^{-cT}SN(d_1)-rKe^{-rT}N(d_2)$. Sensitivity to a negative change in time-to-maturity.

### 2.1 Delta–gamma–vega P&L approximation
Viewing $C=C(S,\sigma)$, Taylor's theorem gives
$$\text{P\&L}\approx \delta\,\Delta S + \tfrac12\Gamma(\Delta S)^2 + \text{vega}\,\Delta\sigma = \text{delta P\&L}+\text{gamma P\&L}+\text{vega P\&L}.$$
With $\Delta\sigma=0$ this is the **delta–gamma approximation** used in historical VaR. Rewritten in returns:
$$\text{P\&L}\approx \underbrace{\delta S}_{\text{ESP}}\Big(\tfrac{\Delta S}{S}\Big)+\tfrac{\Gamma S^2}{2}\Big(\tfrac{\Delta S}{S}\Big)^2+\text{vega}\,\Delta\sigma = \text{ESP}\times R + \$\Gamma\times R^2 + \text{vega}\,\Delta\sigma,$$
where ESP = equivalent stock position ("dollar delta"), $\$\Gamma$ = dollar gamma. Applies to whole portfolios; breaks down for large moves — then use **scenario analysis** (pivot table of P&L over stressed risk factors).

### 2.2 Delta-hedging (discrete, self-financing)
Continuous trading is infeasible; hedge at $0=t_0<\dots<t_n=T$, $\Delta t=t_{i+1}-t_i$. Self-financing rebalancing:
$$V_{i+1}=V_i+\delta_i\big(S_{i+1}+S_i c\,\Delta t-S_i\big)+(V_i-\delta_iS_i)\big(e^{r\Delta t}-1\big).$$
As $\Delta t\to0$ this replicates the payoff **if** the pricing $\sigma_0$ equals true $\sigma$. With the wrong $\sigma_0$, $V_0$ and all $\delta_i$ are wrong and replication can fail badly. Since true $\sigma$ is unknown and prices aren't even GBM, dynamic replication is only a theoretical concept — at best one hedges approximately.

---

## 3. The implied volatility surface

### 3.1 Definition and stylized facts
The vol surface $\sigma(K,T)$ is defined **implicitly** by inverting BS on market option prices:
$$C_{\text{mkt}}(S,K,T)=C_{BS}(S,T,r,c,K,\sigma(K,T)).$$
A unique solution always exists because $C_{BS}$ is strictly increasing in $\sigma$ (vega $>0$). If BS were correct the surface would be **flat and constant**: $\sigma(K,T)=\sigma$ for all $K,T$. In reality it is not flat and moves randomly. Stylized facts:

- **Skew / smile:** for a fixed $T$, lower strikes have **higher** implied vols (equity/index skew). Pronounced at short expirations; essentially absent before the 1987 crash.
- **Term structure:** for fixed $K$, $\sigma(K,T)$ can rise or fall with $T$; converges to a constant as $T\to\infty$. In stress, short-dated vols spike → **inverted** term structure.
- Single-stock options are usually American → calls and puts give different surfaces.
- BS is "the wrong number in the wrong formula to get the right price" — yet every equity/FX desk computes the BS implied surface and BS Greeks.

### 3.2 Why is there a skew?
1. **Risk aversion / supply-demand:** downside jumps are larger and more frequent; fear raises vol as markets fall; investors buy OTM puts for protection → extra demand at low strikes.
2. **Leverage effect (Merton structural view):** firm value $V=D+E$; equity = call on $V$ struck at $D$. From $V+\Delta V=(E+\Delta E)+(D+\Delta D)$, if equity absorbs losses and debt is near-riskless ($\Delta D\approx0$): $\sigma_V\approx \tfrac{E}{V}\sigma_E$, hence
$$\sigma_E \approx \frac{V}{E}\sigma_V = \Big(1+\frac{D}{E}\Big)\sigma_V.$$
With $\sigma_V$ constant, equity vol $\sigma_E$ **rises as $E$ falls** — producing the skew.

### 3.3 No-arbitrage constraints on the surface
1. $\sigma(K,T)\ge0$ for all $K,T$.
2. At any $T$, the **skew cannot be too steep** (else put-spread / butterfly arbitrage).
3. The **term structure cannot be too inverted** (else calendar-spread arbitrage).
These are hard to enforce when *stressing* the surface for risk management.

### 3.4 What the surface tells us — Breeden–Litzenberger
A butterfly centered at $K$ (long $K\pm\Delta K$ calls, short two $K$ calls) has value
$$B_0 = C(K-\Delta K,T)-2C(K,T)+C(K+\Delta K,T) \approx e^{-rT}f_T(K)(\Delta K)^2,$$
where $f_T$ is the risk-neutral PDF of $S_T$. Hence
$$f_T(K)\approx e^{rT}\frac{C(K-\Delta K,T)-2C(K,T)+C(K+\Delta K,T)}{(\Delta K)^2}\ \xrightarrow{\Delta K\to0}\ f_T(K)=e^{rT}\frac{\partial^2 C}{\partial K^2}.$$
So the surface gives the **marginal** risk-neutral density of $S_T$ for each fixed $T$. Therefore we can price **any** payoff $f(S_T)$ depending on a single fixed time: $P_0=E^Q_0[e^{-rT}f(S_T)]$. It tells us **nothing** about the *joint* distribution across multiple times $T_1,\dots,T_n$ — because it is built from European prices, which depend only on marginals. Hence path-dependent products (knockout/barrier options) **cannot** be priced from the surface alone.

### 3.5 Pricing with the surface — digital and range accrual
**Digital** (pays \$1 if $S_T>K$):
$$D(K,T)=\lim_{\Delta K\to0}\frac{C_{\text{mkt}}(K,T)-C_{\text{mkt}}(K+\Delta K,T)}{\Delta K}=-\frac{\partial C_{\text{mkt}}}{\partial K}.$$
Since $C_{\text{mkt}}(K,T)=C_{BS}(K,T,\sigma(K,T))$, the chain rule gives the practitioner formula
$$D(K,T)=-\frac{\partial C_{BS}}{\partial K}-\frac{\partial C_{BS}}{\partial\sigma}\frac{\partial\sigma}{\partial K}=-\frac{\partial C_{BS}}{\partial K}-\text{vega}\times\text{skew}.$$

**Worked example (Gatheral):** $r=c=0$, $T=1$, $S_0=K=100$ (ATM digital), skew $=2.5\%$ per $10\%$ change in strike, $\sigma_{\text{atm}}=25\%$.
$$D(100,1)=N\!\Big(-\tfrac{\sigma_{\text{atm}}}{2}\Big)-S_0\,\varphi\!\Big(\tfrac{\sigma_{\text{atm}}}{2}\Big)\times\frac{-.025}{.1\,S_0}=N(-.125)+.25\,\varphi(.125)\approx .45+.25\times.4=.55.$$
Correct price = **55%** of notional. Ignoring the skew (just the ATM BS digital) gives only **45%** — materially too low. The skew term is decisive.

**Range accrual** (pays $X\%$ of notional, $X$ = fraction of days the index stays in $[1500,1550]$): replicate as a portfolio of **pairs of digitals** for every observation date → priceable from the surface.

### 3.6 Beyond the surface — same marginals, different joints
Two non-dividend securities $A,B$ with
$$S_{T_1}=e^{(r-\sigma^2/2)T_1+\sigma\sqrt{T_1}Z_1},\quad S_{T_2}=e^{(r-\sigma^2/2)T_2+\sigma\sqrt{T_2}(\rho Z_1+\sqrt{1-\rho^2}Z_2)},$$
$Z_1,Z_2$ iid $N(0,1)$; $\rho>0$ = momentum, $\rho<0$ = mean-reversion. Because $\rho Z_1+\sqrt{1-\rho^2}Z_2\sim N(0,1)$, $A$ and $B$ have **identical marginals** at $T_1$ and $T_2$ → identical European prices → identical vol surfaces, regardless of $\rho_A,\rho_B$. But a knock-in put with payoff $\max(1-S_{T_2},0)\,\mathbf 1_{\{S_{T_1}\ge1.2\}}$ depends on the joint law and its price **differs** with $\rho$. Conclusion: the surface cannot price barrier/path-dependent options.

### 3.7 Pricing in practice (the honest summary)
Dynamic replication is elegant but not achievable in practice; **supply and demand** set prices of liquid options (volatility is itself an asset class). Models remain needed to (1) price illiquid/exotic securities and (2) risk-manage via Greeks / scenario analysis. Models are arbitrage-free by construction and calibrated to liquid prices, but are only approximations (witness frequent re-calibration) and ignore market endogeneity — an "almost fatal flaw."

---

## 4. Structured credit: CDOs and the Gaussian copula

### 4.1 Securitization and CDOs
Securitization builds new securities from cash-flows of a pool, enabling a broad range of risk profiles and broadening demand → lowers issuers' cost of capital. **CDOs** are built from a pool of fixed-income securities (first issued mid-1990s, motivated by **regulatory arbitrage**); **CLOs** from loans. Mechanics mirror those of credit default swaps (CDS).

### 4.2 The one-factor Gaussian copula model
$N$ credits, notional $A_i$, recovery $R_i$; loss on default of credit $i$ is $A_i(1-R_i)$. The risk-neutral marginal default-time law is assumed known → $q_i(t)=Q(\tau_i\le t)$ (from CDS spreads or bond prices). Define the normalized asset value
$$X_i = a_i M + \sqrt{1-a_i^2}\,Z_i,\qquad M,Z_1,\dots,Z_N\ \text{iid}\ N(0,1),\ a_i\in[0,1].$$
Each $X_i\sim N(0,1)$; $\text{Corr}(X_i,X_j)=a_i a_j$; correlation matrix $P_{ij}=a_ia_j$ ($i\ne j$). $M$ is the single common factor driving dependence. Credit $i$ defaults by $t_i$ if $X_i<\bar x_i(t_i)$ with threshold
$$\bar x_i(t_i)=\Phi^{-1}(q_i(t_i)).$$
Joint default-time distribution (the copula):
$$F(t_1,\dots,t_N)=\Phi_P\big(\Phi^{-1}(q_1(t_1)),\dots,\Phi^{-1}(q_N(t_N))\big).$$

### 4.3 Portfolio loss distribution
**Conditional on $M$, defaults are independent**, with conditional default probability
$$q_i(t\,|\,M)=\Phi\!\left(\frac{\bar x_i(t)-a_i M}{\sqrt{1-a_i^2}}\right).$$
Let $p_N(l,t)$ = prob. of exactly $l$ defaults before $t$. Then
$$p_N(l,t)=\int_{-\infty}^{\infty}p_N(l,t\,|\,M)\,\varphi(M)\,dM.$$
$p_N(l,t\,|\,M)$ is built by a simple **iterative (recursive) convolution** over credits; the outer integral over $M$ is done by numerical quadrature. If $A_i\equiv A$ and $R_i\equiv R$ are constant, each credit's loss is $0$ or $A(1-R)$, so the number-of-defaults distribution = total-loss distribution.

### 4.4 Worked example — 1-period CDO (Donnelly–Embrechts)
$T=1$ yr; $N=125$ bonds; each pays coupon 1 if no default; recovery $R=0$; identical default prob $q$ and identical pairwise correlation $\rho$:
$$X_i=\sqrt\rho\,M+\sqrt{1-\rho}\,Z_i,\qquad \bar x_1=\dots=\bar x_N=\Phi^{-1}(q).$$
Conditional default prob
$$q_M=\Phi\!\left(\frac{\Phi^{-1}(q)-\sqrt\rho\,M}{\sqrt{1-\rho}}\right),$$
and conditional on $M$ the number of defaults is **Binomial**$(N,q_M)$:
$$p_N(l\,|\,M)=\binom{N}{l}q_M^{\,l}(1-q_M)^{N-l},\qquad p(l)=\int_{-\infty}^\infty p_N(l\,|\,M)\,\varphi(M)\,dM.$$
(When correlations/default probs are not identical, the simple Binomial form is invalid; one must build $p_N(l\,|\,M)$ by the iterative convolution.)

Three equal tranches of notional 3 units each (equity = defaults 1–3, mezzanine = 4–6, senior = 7–9):
$$E^Q_0[\text{Equity loss}]=3\,P(\ge3\text{ defaults})+\sum_{k=1}^{2}k\,P(k\text{ defaults}),$$
$$E^Q_0[\text{Mezzanine loss}]=3\,P(\ge6)+\sum_{k=1}^{2}k\,P(k+3),\qquad E^Q_0[\text{Senior loss}]=3\,P(\ge9)+\sum_{k=1}^{2}k\,P(k+6).$$

**Observations on correlation $\rho$:**
- Always $E^Q_0[\text{Equity}]\ge E^Q_0[\text{Mezzanine}]\ge E^Q_0[\text{Senior}]$ (when tranches share equal notional).
- Equity-tranche expected loss is **decreasing in $\rho$**.
- Mezzanine is often **insensitive to $\rho$** (complicates calibration).
- Super-senior (upper attachment 100% / 125 units) expected loss is **increasing in $\rho$**.
- The **total** expected loss across all tranches (= loss on the index) is **independent of $\rho$** (not an accident — correlation reshuffles losses across tranches but does not change the mean total loss).

### 4.5 Synthetic CDO tranche mechanics
A tranche is defined by lower/upper attachment points $L,U$. With $l$ defaults the **tranche loss function** is
$$TL^{L,U}_t(l)=\max\{\min\{lA(1-R),U\}-L,\,0\}.$$
*Example:* $L=3\%$, $U=7\%$, portfolio loss $lA(1-R)=5\%$ → tranche loss $=2\%$ of portfolio notional $=50\%$ of the $4\%$ tranche notional. Selling protection = reimbursing realized tranche losses (default leg) in exchange for a periodic premium (premium leg), with possible upfront (common for equity tranches, $L=0$). Fair value sets premium-leg PV = default-leg PV (like a swap, initial value 0). Expected tranche loss at time $t$:
$$E^Q_0[TL^{L,U}_t]=\sum_{l=0}^{N}TL^{L,U}_t(l)\,p(l,t),\qquad p(l,t)=\int p_N(l,t\,|\,M)\varphi(M)\,dM.$$

### 4.6 Fair-value of a CDO tranche (premium/default legs)
Premium leg (paid on the *remaining* tranche notional; unlike a CDS it does not terminate on a single default):
$$PL^{L,U}_0 = S\sum_{t=1}^{n} d_t\,\Delta t\,\Big((U-L)-E^Q_0[TL^{L,U}_{t-1}]\Big),$$
with $n$ periods, discount factor $d_t$, annualized spread $S$, accrual factor $\Delta t$ (e.g. $1/4$ quarterly). Default leg (PV of incremental losses):
$$DL^{L,U}_0=\sum_{t=1}^{n} d_t\Big(E^Q_0[TL^{L,U}_t]-E^Q_0[TL^{L,U}_{t-1}]\Big).$$
Fair spread:
$$S^*=\frac{DL^{L,U}_0}{\sum_{t=1}^{n} d_t\,\Delta t\,\big((U-L)-E^Q_0[TL^{L,U}_{t-1}]\big)}.$$
The model's popularity is its **speed** — all $E^Q_0[TL]$ compute quickly. In practice $S^*$ was observed and (8) inverted for an **implied correlation** — but a different $\rho$ is implied per tranche (correlation skew), and sometimes no $\rho$ solves it at all.

### 4.7 Cash vs synthetic CDOs; CDO² and copula caveats
- **Cash CDOs:** reference portfolio physically exists (bonds on balance sheet); tranched to cut regulatory capital; issuer often keeps the equity tranche. Require waterfalls, credit enhancement, ratings, lengthy legal docs.
- **Synthetic CDOs:** reference portfolio is fictitious (CDS-based). Can sell a single tranche (e.g. 3%–7%) without placing the whole deal. Issuer no longer hedged by owning bonds → must **dynamically hedge** via the CDS markets — very hard ("too many moving parts").
- **CDO²:** pool the mezzanine tranches of ~100 CDOs into a new CDO's collateral; a mezz tranche plays the role of a bond. Incurs loss only when an underlying mezz tranche does. Same bonds re-appear across many underlying CDOs → enormous correlation/overlap complexity (~$100\times150=15{,}000$ pages of legal docs; thousands of lines of code). Extends to ABS-CDOs, CDO³, etc.
- **Copula caveat:** marginals + correlation do **not** determine the joint distribution; correlation captures only *linear* dependence. A Bivariate-Normal and a Meta-Gumbel copula can share standard-normal marginals and correlation $0.7$ yet differ dramatically in joint-tail behavior (Meta-Gumbel produces far more large joint moves) — tail dependence the Gaussian copula misses. Risk management challenges: scenario design, model-dependent Greeks (cf. Ford/GM downgrade May 2005), liquidity risk, market endogeneity, over-reliance on ratings/models.

---

## 5. Real options — valuing operational flexibility

### 5.1 The paradigm
Option theory applied to non-financial investment/decision problems: gold-mine lease, equipment-upgrade option, drug development, contracting out manufacturing capacity, power-plant tolling contracts, gas-storage leases. Projects face market, industry, technical, organizational and political uncertainty; **add flexibility via options** and manage that flexibility. Investments with negative standalone NPV can be positive once embedded options (e.g. option to develop follow-on technology) are valued.

### 5.2 Simplico gold-mine (Luenberger, Investment Science 12.7)
Value a **10-year lease** on a gold mine. Gold price binomial model: current $S_0=\$400$/oz; each year $\times u=1.2$ w.p. $p=0.75$ or $\times d=0.9$ w.p. $1-p=0.25$. Extraction cost $C=\$200$/oz; max rate $G=10{,}000$ oz/yr; interest $r=10\%$. Risk-neutral prob
$$q=\frac{(1+r)-d}{u-d}=\frac{1.1-0.9}{1.2-0.9}=\frac{0.2}{0.3}=\frac23,\qquad 1-q=\frac13.$$
(Implicitly assume gold can be bought/short-sold.) Convention: cash flows at year-end; revenue set by price at year's start.

**State** $s$ = gold price at the start of the year; $V_t(s)$ = lease value in state $s$ at start of year $t$. Lease worthless at expiry: $V_{10}(s)=0$. Backward recursion:
$$V_t(s)=\underbrace{\max\{s-C,0\}\,G}_{\text{year-}t\text{ revenue}}+\frac{q\,V_{t+1}(us)+(1-q)\,V_{t+1}(ds)}{1+r}.$$
Optimal to produce at full $G$ when $s>C$; when $s<\$200$ extract nothing (no production cost — free shutdown assumed).

*Gold price lattice (selected):* root $400$; up-spine $480,576,691.2,829.4,995.3,1194.4,1433.3,1719.9,2063.9,2476.7$; bottom-spine at $t=10$ down to $139.5$.

*Lease-value lattice (in \$millions):* terminal column all $0$. Top node at $t=9$: $16.94=10{,}000\,(2063.9-200)/1.1$. Mid node at $t=6$ (value $12.0$M):
$$12{,}000{,}000=\frac{10{,}000(503.9-200)}{1.1}+\frac{q\times 11{,}500{,}000+(1-q)\times 7{,}400{,}000}{1.1}.$$
**Lease value at $t=0$: $\approx\$24.1$ million.**

### 5.3 Equipment-upgrade (American) option
Option, exercisable any time over the lease, to buy new equipment for **\$4 M**: new rate $G_{\text{up}}=12{,}500$ oz/yr, new cost $C_{\text{up}}=\$240$/oz. Once installed it stays for all future years and reverts to the mine owner at lease end.

Step 1 — value the lease *assuming the enhancement is in place from $t=0$* (same recursion with $G_{\text{up}},C_{\text{up}}$): $V^{\text{up}}_t(s)$. Result $V^{\text{up}}_0\approx\$27$ M (before subtracting the \$4 M cost).

Step 2 — value the option as an American exercise comparing "exercise now" vs "continue":
$$U_t(s)=\max\Big\{\underbrace{V^{\text{up}}_t(s)-\$4\text{M}}_{\text{exercise}},\ \underbrace{\max\{s-C,0\}\,G+\tfrac{q\,U_{t+1}(us)+(1-q)\,U_{t+1}(ds)}{1+r}}_{\text{continue with old }C=200,\,G=10{,}000}\Big\}.$$
Working backward from $t=10$ (start with original equipment), at each node compare lease value $A$ (continue) to $B-\$4$M where $B=V^{\text{up}}$ at the matching node; install if $B-\$4\text M\ge A$, and store $\max(B-\$4\text M, A)$. Nodes where installation is optimal are marked `*` (an early-exercise/installation boundary; e.g. top-left region of the lattice).

**Lease value with option at $t=0$: $\approx\$24.6$ M** → **option value $\approx 24.6-24.1=\$0.5$M** (precisely \$558,595 per the workbook). All uncertainty here is *financial* (gold price); real cases add non-financial uncertainty (random extractable quantity, stochastic storage/extraction costs).

### 5.4 Related real options (binomial-lattice / DP framework)
**Gas-storage lease** (cavern capacity $C$, inventory $I_t$, pump rate $z_t$, loss $f(z,I)$, price $P_t$):
$$\max_{\{z_t\}}E^Q_0\!\Big[\sum_{t=1}^T e^{-rt}P_t\big(z_t-f(z_t,I_t)\big)\Big]\ \text{s.t.}\ I_{t+1}=I_t-z_t,\ I_t\in[0,C].$$
DP: $V_t(I,P)=\max_z\{P(z-f(z,I))+e^{-r}E^Q_t[V_{t+1}(I-z,P_{t+1})]\}$. Use a price lattice for $P$ but enumerate inventory levels $I$ (or approximate DP / simulation-based optimization).

**Tolling agreement on a gas-fired power plant** (states $s_t\in\{0,1\}$ off/on; actions $a\in\{0,1\}$; ramp-up cost $C_u$, ramp-down $C_d$, rent $K$):
$$c(0,0)=K,\ c(0,1)=C_u+K,\ c(1,0)=C_d+K,\ c(1,1)=\max\{P_tQ-G_tH-K,\ \dots\};$$
$$V_t(s_t,P_t,G_t)=\max_a\big\{c(s_t,a)+e^{-r}E^Q_t[V_t(u(s_t,a),P_{t+1},G_{t+1})]\big\}.$$
Solve on a joint binomial lattice for electricity price $P_t$ and gas price $G_t$, with two plant states per node.

---

## 6. Excel workbooks (purpose)

- **`EquityDerivsPractice_PSet3.xlsx`** — sheets `DemoSheet`, `StockPricePaths`. Simulates GBM stock-price paths and demonstrates discrete **delta-hedging** P&L: builds the self-financing replicating strategy of §2.2 and shows replication error when rebalancing at finite $\Delta t$ and/or hedging with a mis-specified $\sigma_0$.
- **`StructuredCredit_PSet4.xlsx`** — sheet `LossDistributionCalculations`. Implements the one-factor Gaussian copula: conditional binomial default probabilities $q_M$, numerical integration over the common factor $M$ to get the portfolio loss distribution $p(l)$, and the resulting **expected tranche losses** (equity/mezzanine/senior) for the 125-name, 1-year CDO of §4.4.
- **`simplico.xlsx`** — sheets `Lease`, `Equipment` plus named-parameter sheets (`cost`, `d`, `interest`, `q`, `rate`, `S0`, `u`). Builds the gold-price lattice, the lease-value lattice ($\approx\$24.1$M), the enhancement-in-place lattice ($\approx\$27$M), and the option-to-upgrade lattice ($\approx\$24.6$M) — yielding the exact upgrade-option value **\$558,595**.

---

## 7. Big-picture takeaways
1. Black–Scholes is known to be wrong (no jumps, thin tails, IID returns) but is the lingua franca: prices/Greeks are quoted via BS implied vols.
2. The implied vol surface encodes the **marginal** risk-neutral densities (Breeden–Litzenberger) and thus prices any single-date payoff — but says nothing about the joint/path distribution, so it cannot price barrier/path-dependent exotics.
3. Liquid derivative prices are set by **supply and demand**; models exist to interpolate to exotics and to risk-manage (Greeks + scenario analysis), always as imperfect, re-calibrated approximations.
4. The one-factor Gaussian copula made CDO-tranche pricing fast and standard, but a single correlation $\rho$ cannot match all tranches (correlation skew), correlation ignores tail dependence, and CDO²/ABS-CDO complexity made structured credit nearly impossible to risk-manage — central to the 2007–08 crisis.
5. Real options value operational flexibility with the same risk-neutral binomial machinery: where uncertainty is purely financial (Simplico gold price), the lease and its embedded American upgrade option are priced exactly by backward induction.
