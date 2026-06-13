# Course 2 — Term Structure & Credit Derivatives

Columbia "Financial Engineering and Risk Management" (Haugh & Iyengar), Course 2.
Covers: the binomial short-rate lattice, calibration (Ho-Lee, Black-Derman-Toy), fixed-income
derivatives (bond options/forwards/futures, caps/floors, swaps/swaptions, elementary/state prices),
defaultable bonds and the hazard-rate model, Credit Default Swaps (mechanics + par-spread valuation),
and mortgage mathematics / MBS / securitization (pass-throughs, PO/IO strips, CMOs).

Throughout we work **directly under the risk-neutral measure $Q$**: we specify $Q$-dynamics for the
short rate (no reference to true probabilities), price by risk-neutral expectation (guarantees
no-arbitrage **regardless of the probabilities used**), and then **calibrate** parameters to liquid
market prices. By convention we take $q_u=q_d=\tfrac12$ unless otherwise specified.

---

## Key formulas cheat-sheet

| Quantity | Formula |
|---|---|
| Risk-neutral pricing, no coupon | $Z_{i,j}=\dfrac{1}{1+r_{i,j}}\big[q_u Z_{i+1,j+1}+q_d Z_{i+1,j}\big]$ |
| With coupon $C$ | $Z_{i,j}=\dfrac{1}{1+r_{i,j}}\big[q_u(Z_{i+1,j+1}+C_{i+1,j+1})+q_d(Z_{i+1,j}+C_{i+1,j})\big]$ |
| Cash account | $B_{t}=\prod_{k=0}^{t-1}(1+r_k)$, $B_0=1$, $B_t/B_{t+1}=1/(1+r_t)$ |
| Martingale (deflated) pricing | $\dfrac{Z_t}{B_t}=E_t^Q\!\Big[\sum_{j=t+1}^{t+s}\dfrac{C_j}{B_j}+\dfrac{Z_{t+s}}{B_{t+s}}\Big]$ |
| Spot rate from ZCB | $Z_0^k(1+s_k)^k = 100$ (per-period compounding) |
| Forward price (delivery $n$) | $G_0=\dfrac{E_0^Q[S_n/B_n]}{E_0^Q[1/B_n]}$ |
| Futures price (expiry $n$) | $F_0=E_0^Q[S_n]$ |
| Caplet payoff (in arrears, mat $\tau$) | $(r_{\tau-1}-c)^+$ paid at $\tau$; floorlet $(c-r_{\tau-1})^+$ |
| Swap payment | $\pm(r_{i,j}-K)$ paid at $i{+}1$ |
| Elementary/state price (forward eqn) | $P^e_{k+1,s}=\dfrac{P^e_{k,s-1}}{2(1+r_{k,s-1})}+\dfrac{P^e_{k,s}}{2(1+r_{k,s})}$ |
| BDT short rate | $r_{i,j}=a_i e^{b_i j}$ (Ho-Lee: $r_{i,j}=a_i+b_i j$) |
| Survival prob (state-indep hazard) | $q(t)=\prod_{k=0}^{t-1}(1-h_k)$ |
| Defaultable ZCB, no recovery | $\bar Z^T_{i,j,0}=\dfrac{1-h_{ij}}{1+r_{ij}}\big[q_u\bar Z^T_{i+1,j+1,0}+q_d\bar Z^T_{i+1,j,0}\big]\approx e^{-(r_{ij}+h_{ij})}E[\cdot]$ |
| CDS par spread | $S_{par}=\dfrac{(1-R)\sum_k (q(t_{k-1})-q(t_k))\,d(0,t_k)}{\tfrac{\delta}{2}\sum_k (q(t_{k-1})+q(t_k))\,d(0,t_k)}\approx\dfrac{(1-R)h}{1-h/2}$ |
| Level mortgage payment | $B=\dfrac{c(1+c)^n M_0}{(1+c)^n-1}$ |
| Remaining principal | $M_k=M_0\dfrac{(1+c)^n-(1+c)^k}{(1+c)^n-1}$ |
| Scheduled principal in period $k$ | $P_k=(B-cM_0)(1+c)^{k-1}$; interest $I_k=cM_{k-1}$ |
| PSA benchmark | $\text{CPR}=6\%\cdot(t/30)$ for $t\le30$, else $6\%$ |
| CPR / SMM | $\text{SMM}=1-(1-\text{CPR})^{1/12}$, $\text{CPR}=1-(1-\text{SMM})^{12}$ |

---

## 1. The binomial short-rate lattice

Fixed income is harder than equities because we must model the **evolution of the entire term
structure**, not a single price. The state variable is the **short rate** $r_t$ — the one-period
risk-free rate applying between $t$ and $t+1$; it is random but **known by time $t$**.

**Lattice geometry.** Nodes $(i,j)$: date $i=0,\dots,n$, state $j=0,\dots,i$. From $(i,j)$ the rate
moves up to $(i+1,j+1)$ with prob $q_u$ or down to $(i+1,j)$ with prob $q_d$, $q_u+q_d=1$, both $>0$.
We specify the **short rate $r_{i,j}$ at each node** rather than all bond prices directly (specifying
bond prices directly is awkward to keep arbitrage-free).

**Risk-neutral pricing (no coupon):**
$$Z_{i,j}=\frac{1}{1+r_{i,j}}\big[q_u Z_{i+1,j+1}+q_d Z_{i+1,j}\big].$$
**With a coupon $C_{i+1,\cdot}$ paid at $i+1$** (here $Z_{i+1,\cdot}$ is the *ex-coupon* price):
$$Z_{i,j}=\frac{1}{1+r_{i,j}}\big[q_u(Z_{i+1,j+1}+C_{i+1,j+1})+q_d(Z_{i+1,j}+C_{i+1,j})\big].$$
Pricing this way is **arbitrage-free for any choice of probabilities** — this is the key reason we are
free to set $q_u=q_d=1/2$ and calibrate the $r_{i,j}$ instead.

**Cash account.** $B_t=(1+r_0)(1+r_1)\cdots(1+r_{t-1})$, $B_0=1$. It is *locally* risk-free
($B_{t+1}$ known at $t$) but not risk-free ($B_{t+s}$ unknown for $s>1$). Pricing is equivalent to
the **deflated-price martingale property**:
$$\frac{Z_t}{B_t}=E_t^Q\!\Big[\frac{Z_{t+s}}{B_{t+s}}\Big],\qquad
\frac{Z_t}{B_t}=E_t^Q\!\Big[\sum_{j=t+1}^{t+s}\frac{C_j}{B_j}+\frac{Z_{t+s}}{B_{t+s}}\Big].$$

### Sample short-rate lattice (used in all examples)
$r$ moves by $u=1.25$ or $d=0.9$ each period, starting at $r_{0,0}=6\%$:

| $t=0$ | $t=1$ | $t=2$ | $t=3$ | $t=4$ | $t=5$ |
|---|---|---|---|---|---|
| 6.00% | 7.50% | 9.38% | 11.72% | 14.65% | 18.31% |
| | 5.40% | 6.75% | 8.44% | 10.55% | 13.18% |
| | | 4.86% | 6.08% | 7.59% | 9.49% |
| | | | 4.37% | 5.47% | 6.83% |
| | | | | 3.94% | 4.92% |
| | | | | | 3.54% |
(Each row = state $j$, top row highest state.)

### Pricing a 4-period ZCB (face 100)
Backward induction. Terminal $Z^4_{4,\cdot}=100$. Example node:
$83.08=\frac{1}{1.0938}\big[\tfrac12\cdot89.51+\tfrac12\cdot92.22\big]$. Result $Z^4_0=77.22$.
Spot rate: $77.22(1+s_4)^4=100\Rightarrow s_4=6.68\%$. Pricing ZCBs of all maturities backs out the
**whole term structure**; re-pricing at $t=1$ gives a *new* term structure → the short-rate model
**is** a term-structure model.

---

## 2. Fixed-income derivatives on the lattice

All priced by the same backward-induction risk-neutral recursion.

**European bond option.** Call on the 4-yr ZCB, strike 84, expiry $t=2$:
payoff $\max(0,Z^4_{2,\cdot}-84)$, then discount back. e.g. $1.56=\frac{1}{1.075}[\tfrac12\cdot0+\tfrac12\cdot3.35]$. Price $=2.97$.

**American option.** Same recursion but take $\max(\text{intrinsic},\text{continuation})$ at every node.
Put strike 88, expiry $t=3$: $4.92=\max\big(88-83.08,\ \frac{1}{1.0938}[\tfrac12\cdot0+\tfrac12\cdot0]\big)$.

**Forward on a coupon bond.** Delivery of underlying at $t=n$, price known at 0 so from $0=E_0^Q[(Z^T_n-G_0)/B_n]$:
$$G_0=\frac{E_0^Q[S_n/B_n]}{E_0^Q[1/B_n]},\qquad E_0^Q[1/B_n]=Z_0^n \text{ (ZCB price)}.$$
Worked: forward at $t=4$ on a 2-yr 10% coupon bond. First find ex-coupon $Z^6_4$ at $t=4$ by backward
induction (e.g. $98.44=\frac{1}{1.1055}[\tfrac12\cdot107.19+\tfrac12\cdot110.46]$), then discount
$Z^6_4/B_4$ back to $0$ giving $E_0^Q[Z^6_4/B_4]=79.83$. With $Z_0^4=0.7722$:
$G_0=79.83/0.7722=103.38$.

**Futures.** Marked to market, zero entry value each period $\Rightarrow F_k=E_k^Q[F_{k+1}]$, so by
iterated expectations $F_0=E_0^Q[S_n]$ (**no discounting**), holds with or without coupons.
On the same bond: $F_0=103.22$ — close to but **not equal** to the forward $103.38$ (difference due to
the correlation of the underlying with rates / the $1/B_n$ weighting in the forward).

**Caplets / floorlets / caps / floors.** A caplet $\approx$ European call on the short rate. Maturity
$\tau$, strike $c$, settled **in arrears**: payoff $(r_{\tau-1}-c)^+$ at $\tau$ (floorlet $(c-r_{\tau-1})^+$).
A **cap/floor** is a sequence of caplets/floorlets sharing one strike. Trick: record the $t=\tau$
cash-flow at its $t=\tau-1$ predecessor discounted once: $(r_{\tau-1}-c)^+/(1+r_{\tau-1})$, then roll
back. e.g. caplet strike 2%, expiry $t=6$: $0.015=\frac{\max(0,.0354-.02)}{1.0354}$; rolling back
$0.021=\frac{1}{1.0394}[\tfrac12\cdot.028+\tfrac12\cdot.015]$.

**Swaps.** Fixed rate (strike) $K$; payment $\pm(r_{i,j}-K)$ at $t=i+1$ (in arrears). Record at
predecessor node discounted: $\pm(r_{5,5}-K)/(1+r_{5,5})$. Example $K=5\%$, payments $t=1..6$:
$0.1686=\frac{1}{1.0938}\big[(.0938-.05)+\tfrac12\cdot.1793+\tfrac12\cdot.1021\big]$.
(A payer swap pays fixed/receives floating; values shown are per \$1 notional.)

**Swaptions.** Option on a swap. Payoff at option expiry $T_o$ is $\max(0,S_{T_o})$ where $S_{T_o}$ is
the underlying swap value; below $T_o$ discount back **without** including the swap's own cash flows.
Example: swaption strike 0%, expiry $t=3$, underlying 5% swap expiring $t=6$:
$.0908=\frac{1}{1.075}[\tfrac12\cdot.1286+\tfrac12\cdot.0665]$.

### Elementary (state) prices and the forward equations
$P^e_{i,j}$ = time-0 price of a security paying \$1 only at node $(i,j)$, 0 elsewhere. They satisfy the
**forward equations** (iterate forward from $P^e_{0,0}=1$):
$$P^e_{k+1,s}=\frac{P^e_{k,s-1}}{2(1+r_{k,s-1})}+\frac{P^e_{k,s}}{2(1+r_{k,s})},\quad 0<s<k+1,$$
$$P^e_{k+1,0}=\frac{P^e_{k,0}}{2(1+r_{k,0})},\qquad P^e_{k+1,k+1}=\frac{P^e_{k,k}}{2(1+r_{k,k})}.$$
Once computed, many prices become a single dot product. e.g.
$Z^4_0=100\cdot(P^e_{4,0}+\cdots+P^e_{4,4})=100(.0449+.1868+.2901+.1992+.0511)=77.22$.

**Forward-start swap example (using elementary prices).** Swap starting $t=1$ ending $t=3$, notional
\$1M, fixed 7%, payments at $t=2,3$:
$$V_0=\tfrac{.07-.0938}{1.0938}.2194+\tfrac{.07-.0675}{1.0675}.4432+\tfrac{.07-.0486}{1.0486}.2238+\tfrac{.07-.075}{1.075}.4717+\tfrac{.07-.054}{1.054}.4717=\$5{,}800.$$

---

## 3. Calibration of the term-structure model

A model is useless unless its prices match liquid market prices (caps, floors, swaptions, ZCBs). With
$q\equiv1/2$ fixed, we calibrate a **parametric form** for $r_{i,j}$:

- **Ho-Lee:** $r_{i,j}=a_i+b_i j$ ($2n$ params; SD of one-period rate $=b_i/2$; can go negative).
- **Black-Derman-Toy (BDT):** $r_{i,j}=a_i e^{b_i j}$. $\log a_i$ = drift of $\log r$; $b_i$ =
  volatility of $\log r$. Rates stay positive.

**Calibrating BDT to the observed term structure** $(s_1,\dots,s_n)$ (per-period compounded). For each
maturity $i$ the model ZCB price must equal the market price, i.e. (with $b_i=b$ assumed known):
$$\frac{1}{(1+s_i)^i}=\sum_{j=0}^{i}P^e_{i,j}
=\frac{P^e_{i-1,0}}{2(1+a_{i-1})}+\sum_{j=1}^{i-1}\!\Big[\frac{P^e_{i-1,j}}{2(1+a_{i-1}e^{bj})}+\frac{P^e_{i-1,j-1}}{2(1+a_{i-1}e^{b(j-1)})}\Big]+\frac{P^e_{i-1,i-1}}{2(1+a_{i-1}e^{b(i-1)})}.$$
Solve iteratively for the $a_i$ (in Excel via **Solver** + the forward equations).

**Why volatility $b$ matters — payer-swaption example.** A 2×8 payer swaption (option to enter an
8-yr swap in 2 yrs; payments yrs 3–10; "payer" = pay fixed/receive floating), fixed 11.65%, 10-period
annual lattice. Calibrate $a_i$ to market ZCBs:
- $b=0.005 \Rightarrow$ swaption price **\$13,339**
- $b=0.01 \Rightarrow$ swaption price **\$19,497** (must **recalibrate** when $b$ changes)

Both models reproduce the *same* ZCB prices, yet swaption prices differ hugely: ZCB calibration does
not pin down volatility. **Lesson:** calibration securities must be *close* to the securities you want
to price.

**Calibration in practice.** Minimize weighted squared pricing error with a regularizer:
$$\min_\theta \sum_i \omega_i\big(P_i(\text{model})-P_i(\text{market})\big)^2+\lambda\,\|\theta-\theta_{\text{prev}}\|^2.$$
$\omega_i$ = confidence weight, $\theta_{\text{prev}}$ = previous calibration, $\lambda$ = stickiness.
The problem is **non-convex** (many local minima), requires re-calibration many times a day. Derivatives
pricing in practice = arbitrage-free **inter/extrapolation** from liquid prices to illiquid ones; true
probabilities and risk aversion enter only **through the calibration**.

---

## 4. Defaultable bonds and the hazard-rate model

A defaultable bond: coupon $c$, face $F$, **random recovery** $\tilde R$ = fraction of face recovered
on default. Model the term structure of default via the **one-step (conditional) default probability**
= **hazard rate** $h(t)=Q(\text{default in }[t,t+1)\mid \mathcal F_t)$.

**Lattice with defaults.** Split each node $(i,j)$ by a default indicator $\eta\in\{0,1\}$:
$(i,j,0)$ = survived ($\tau>i$), $(i,j,1)$ = defaulted ($\tau\le i$, **absorbing**). From a no-default
node $(i,j,0)$:
$$Q\big((i{+}1,s,\eta)\mid(i,j,0)\big)=\begin{cases}q_u h_{ij}&s=j{+}1,\eta=1\\ q_u(1-h_{ij})&s=j{+}1,\eta=0\\ q_d h_{ij}&s=j,\eta=1\\ q_d(1-h_{ij})&s=j,\eta=0.\end{cases}$$

**Default-free ZCB** (default events irrelevant): $Z^T_{i,j}=\frac{1}{1+r_{ij}}[q_u Z^T_{i+1,j+1}+q_d Z^T_{i+1,j}]$. Calibrate the short-rate lattice to these.

**Defaultable ZCB, no recovery** ($\bar Z^T_{i,j,1}\equiv0$):
$$\bar Z^T_{i,j,0}=\frac{1-h_{ij}}{1+r_{ij}}\big[q_u\bar Z^T_{i+1,j+1,0}+q_d\bar Z^T_{i+1,j,0}\big]\approx e^{-(r_{ij}+h_{ij})}E^{\bar Q}_i[\bar Z^T_{i+1,\cdot,\cdot}].$$
So a defaultable bond is priced by **discounting at $r_{ij}+h_{ij}$**: the hazard rate $h_{ij}$ **is the
one-period credit spread**.

**With recovery** $R=E[\tilde R]$ (independent of rate/default dynamics): add the recovered amount paid
on default:
$$\bar Z^T_{i,j,0}=\frac{1}{1+r_{ij}}\big[q_u(1-h_{ij})\bar Z^T_{i+1,j+1,0}+q_d(1-h_{ij})\bar Z^T_{i+1,j,0}\big]+\frac{1}{1+r_{ij}}\big[q_u h_{ij}R+q_d h_{ij}R\big].$$

### State-independent hazard rates
Assume $h_{ij}=h_i$ (default independent of rates). **Survival probability**
$$q(t)=Q(\text{survive to }t)=\prod_{k=0}^{t-1}(1-h_k),\qquad q(t+1)=(1-h_t)q(t).$$
Let $I(t)$ = survival indicator (1 if not defaulted by $t$). Then $E_0^Q[I(t)]=q(t)$, and the default
indicator for period $t$ is $I(t-1)-I(t)$.

**Pricing a defaultable coupon bond** (let $d(0,t)=Z_0^t=E_0^Q[1/B(t)]$, recovery $\tilde R(t_k)F$ paid
on default in period $k$):
$$\bar P(0)=\sum_{k=1}^{n}c\,q(t_k)d(0,t_k)+F\,q(t_n)d(0,t_n)+RF\sum_{k=1}^{n}\big(q(t_{k-1})-q(t_k)\big)d(0,t_k).$$
The three terms = surviving coupons + surviving face + recovery on the default in each period.
(Uses independence of $\tilde R$, default, and rates to split each expectation into $q(\cdot)\times d(0,\cdot)$.)

**Calibrating hazard rates.** With $r$ deterministic and known, the model price $P_i(h)$ of bond $i$
depends on $h=(h_0,\dots,h_{n-1})$. Given $m$ market prices $P^{mkt}_i$ (all driven by the same credit
event):
$$\min_{h\ge0}\ f(h)=\sum_{i=1}^{m}\big(P^{mkt}_i-P_i(h)\big)^2.$$
(Numerical example in **Bonds_and_cds_v2.xlsx**, sheets *Pricing*/*Calibration* — bonds e.g. 1yr c=5%
R=10%, 2yr c=8% R=25%, 3yr c=5% R=50%, 5yr c=10% R=20%.)

---

## 5. Credit Default Swaps (CDS)

The **seller** compensates the **buyer** on a credit event of a reference entity, in return for
periodic premium payments. Notional $N$; spread/coupon $S$; on default seller pays $(1-R)N$; buyer pays
**accrued interest** from last coupon to default.

**Mechanics.** Coupon dates $t_k=\delta k$, typically $\delta=\tfrac14$ (quarterly, Mar/Jun/Sep/Dec 20).
If no default at $t_k$, buyer pays premium $\delta S N$. If default $\tau\in(t_{k-1},t_k]$, contract
terminates at $t_k$: buyer pays accrued interest $(t_k-\tau)SN$, buyer **receives** $(1-R)N$.

**Numerical example.** 2-yr CDS, $N=\$1\text{M}$, $S=160\,$bps, quarterly. Default in month 16,
$R=45\%$:
- Buyer fixed premiums in months 3,6,9,12,15: $\tfrac{SN}{4}=\$4{,}000$ each.
- Buyer accrued interest, month 18: $\tfrac{SN}{12}=\$1{,}333.33$.
- Seller default-contingent payment, month 18: $(1-R)N=\$550{,}000$.

CDS standardized by ISDA in 1999 (amended 2003, 2009). Developed by Blythe Masters (JPMorgan, 1994) —
JPM hedged a \$4.8B Exxon credit line via a CDS with the EBRD. Market notional reached ~\$62T (end-2007);
DTCC gross notional ~\$25T (2012). A CDS is **not** insurance: you can buy protection without holding
the debt ("**naked CDS**"); used for hedging, speculation/leverage, and basis arbitrage between the
reference bond coupon and the CDS spread. Role in 2008 crisis: opacity (OTC), AIG sold ~\$500B of CDS
on its rating, contagion fears, naked CDS pushing spreads (and thus borrowing costs) higher.

### CDS valuation (par spread)
$$\text{Value to buyer}=\underbrace{\text{PV(protection)}}_{\text{protection leg}}-\underbrace{\text{PV(premiums+accrued)}}_{\text{premium leg}}.$$
With $d(0,t_k)=Z_0^{t_k}=E_0^Q[1/B(t_k)]$ and default uniform over each interval:

**Premium leg** (premiums + accrued, the accrued $\approx \delta/2$ of a coupon on the defaulting interval):
$$\delta S N\sum_{k=1}^{n}q(t_k)d(0,t_k)+\frac{\delta S N}{2}\sum_{k=1}^{n}\big(q(t_{k-1})-q(t_k)\big)d(0,t_k)=\frac{\delta S N}{2}\sum_{k=1}^{n}\big(q(t_{k-1})+q(t_k)\big)d(0,t_k).$$

**Protection leg** (contingent $(1-R)N$ on the default in each interval):
$$(1-R)N\sum_{k=1}^{n}\big(q(t_{k-1})-q(t_k)\big)d(0,t_k).$$

**Par spread** $S_{par}$ = spread making contract value zero:
$$\boxed{\,S_{par}=\frac{(1-R)\sum_{k=1}^{n}\big(q(t_{k-1})-q(t_k)\big)d(0,t_k)}{\tfrac{\delta}{2}\sum_{k=1}^{n}\big(q(t_{k-1})+q(t_k)\big)d(0,t_k)}\,}$$
If $q(t_k)\approx(1-h)q(t_{k-1})$ then $S_{par}\approx\dfrac{(1-R)h}{1-h/2}\approx(1-R)h$.
So **CDS spread $\approx$ (loss-given-default)$\times$(hazard rate)** — for fixed $R$ the spread is
directly proportional to $h$, the market's view of default risk. (Computed in **Bonds_and_cds_v2.xlsx**,
sheet *CDS pricing*: columns for PV of premium $N\!\times\!F\!\times\!B$, PV of accrued $N\!\times\!I\!\times\!B$,
PV of expected protection $(1-R)H$.)

---

## 6. Mortgage mathematics

**Level-payment, fully-amortizing mortgage.** Initial principal $M_0=M$, coupon $c$/period, equal
payments $B$, $n$ periods, $M_n=0$. Recursion: $M_k=(1+c)M_{k-1}-B$. Iterating:
$$M_k=(1+c)^kM_0-B\frac{(1+c)^k-1}{c}.$$
Imposing $M_n=0$:
$$\boxed{B=\frac{c(1+c)^nM_0}{(1+c)^n-1}},\qquad M_k=M_0\,\frac{(1+c)^n-(1+c)^k}{(1+c)^n-1}.$$

**Present value** (deterministic, no default/prepay, discount rate $r$/period):
$$F_0=\sum_{k=1}^{n}\frac{B}{(1+r)^k}=\frac{c(1+c)^nM_0}{(1+c)^n-1}\cdot\frac{(1+r)^n-1}{r(1+r)^n}.$$
If $r=c$ then $F_0=M_0$. In practice $r<c$ (default, prepayment, servicing fees, profit, payment
uncertainty).

**Principal/interest split of payment $k$:** interest $I_k=cM_{k-1}$; principal
$P_k=B-cM_{k-1}=(B-cM_0)(1+c)^{k-1}$. Early payments are mostly interest; the relationship reverses
later. (Modeled in **MBS_Structure.xlsx**, sheet *SingleMortgageCashFlows*: beginning balance, monthly
payment, interest, scheduled principal, ending balance.)

---

## 7. Mortgage-Backed Securities & securitization

MBS are asset-backed securities created by **pooling** mortgages (securitization). US issuers: GSAs
(Ginnie Mae, Freddie Mac, Fannie Mae — agency MBS are **guaranteed against default**) or non-agency
(commercial banks — not guaranteed). We assume agency / default-free MBS.

### Pass-through MBS
Pool mortgages; investors receive monthly the underlying interest+principal. The **pass-through coupon
rate < average mortgage coupon** (servicing fees).
- **WAC** (weighted average coupon): weighted avg of pool coupons, weights = outstanding amounts.
- **WAM** (weighted average maturity): weighted avg of remaining months, weights = outstanding amounts.

**Prepayment** = payments above scheduled (home sale, refinancing, default with insurer payoff,
destruction). Conventions:
- **CPR** = annual prepayment rate (% of current outstanding principal).
- **SMM** = single-month mortality = monthly version: $\text{SMM}=1-(1-\text{CPR})^{1/12}$,
  $\text{CPR}=1-(1-\text{SMM})^{12}$.
- **PSA benchmark** (30-yr): $\text{CPR}=6\%\cdot(t/30)$ for $t\le30$ months, $6\%$ thereafter; speeds
  quoted as multiples, e.g. 100 PSA, 300 PSA.

**Average life** (years), $P_k$ = principal (scheduled + projected prepay) at month $k$, $TP$ = total
principal, $T$ months:
$$\text{Average Life}=\sum_{k=1}^{T}\frac{k P_k}{12\cdot TP}.$$
Average life **decreases as PSA speed increases**. **Bond-equivalent yield** = annual rate, semiannual
compounding, that sets PV(expected cash flows at a given PSA) = market price. Yields are limited; the
market standard is the **option-adjusted spread (OAS)**.

**Prepayment risks** for a pass-through (on top of ordinary interest-rate risk):
- **Contraction risk**: rates fall → prepayments rise → prepaid principal reinvested at lower rates.
- **Extension risk**: rates rise → prepayments fall → less principal to reinvest at higher rates.

### Principal-Only (PO) and Interest-Only (IO) strips
Split the pass-through into principal and interest cash-flow streams.

**PO present value** ($P_k=(B-cM_0)(1+c)^{k-1}$):
$$V_0=(B-cM_0)\frac{(1+r)^n-(1+c)^n}{(r-c)(1+r)^n}=\frac{cM_0}{(1+c)^n-1}\cdot\frac{(1+r)^n-(1+c)^n}{(r-c)(1+r)^n}.$$
Limits: $\lim_{c\to r}V_0=\frac{n(B-rM_0)}{1+r}$; for $r=c$, $V_0=\frac{rnM_0}{(1+r)[(1+r)^n-1]}$;
$\lim_{n\to\infty}V_0=0$ (early payments are mostly interest).

**IO present value:** easiest via $W_0=F_0-V_0$:
$$W_0=\frac{cM_0}{[(1+c)^n-1](1+r)^n}\Big[\frac{(1+c)^n(1+r)^n-1}{r}-\frac{(1+r)^n-(1+c)^n}{r-c}\Big],$$
and for $r=c$, $W_0=M_0-\frac{rnM_0}{(1+r)[(1+r)^n-1]}$ (consistent with $F_0=M_0$).

**Durations** (annualized; PO has the **longer** duration):
$$D_P=\frac{1}{12V_0}\sum_{k=1}^{n}\frac{kP_k}{(1+r)^k},\qquad D_I=\frac{1}{12W_0}\sum_{k=1}^{n}\frac{kI_k}{(1+r)^k}=\frac{1}{12W_0}\sum_{k=1}^{n}\frac{kB}{(1+r)^k}-\frac{V_0}{W_0}D_P.$$

**With prepayments** (path-by-path): $I_k=cM_{k-1}$ and
$M_k=M_{k-1}-\text{ScheduledPrincipal}_k-\text{Prepayment}_k$.
**Risk profiles are opposite:** the **PO investor wants prepayments to increase** (gets principal sooner
at par); the **IO investor wants prepayments to decrease** (earns interest only on the remaining
balance). The IO is the rare fixed-income security whose price tends to **move with** the level of
rates (rates up → fewer prepayments → more expected interest, partly offsetting the lower discount factor).

### Collateralized Mortgage Obligations (CMOs)
CMOs **redirect** cash flows from other MBS (often pass-throughs, sometimes POs) to **mitigate
prepayment risk** and tailor securities to investors. Types: **sequential**, accrual (Z) bonds,
floater/inverse-floater tranches, **PAC** (planned amortization class).

**Sequential-pay CMO** — tranches A, B, C, D retired in order:
1. Periodic **coupon interest** paid to each tranche on its outstanding principal at the start of the period.
2. **All principal** (scheduled + prepay) goes to **A** until retired, then to **B**, then **C**, then
   **D**. So A has the shortest, D the longest average life — earlier tranches absorb prepayment/
   contraction risk, later tranches absorb extension risk.
(Modeled in **MBS_Structure.xlsx**, sheet *Sequential-Pay CMO*: par amounts and coupon rates per
tranche A–D, with the *Pass-Through* sheet feeding cash flows.)

### Pricing MBS in practice
The **prepayment model** is the most important and least-standardized part of an MBS pricing engine
(little public data → hard to calibrate). **Richard–Roll (1989)** model the CPR multiplicatively:
$$\text{CPR}_k=\text{RI}_k\times\text{AGE}_k\times\text{MM}_k\times\text{BM}_k,$$
- **Refinancing incentive** $\text{RI}_k=0.28+0.14\tan^{-1}\!\big(-8.57+430(\text{WAC}-r_k(10))\big)$,
  where $r_k(10)$ = prevailing 10-yr spot rate at $k$.
- **Seasoning** $\text{AGE}_k=\min(1,\,t/30)$.
- **Monthly multiplier** $\text{MM}_k$ (seasonal).
- **Burnout** $\text{BM}_k=0.3+0.7\,M_{k-1}/M_0$ (pools that have already prepaid heavily prepay less).

Also need a **term-structure model** (calibrated to market rates + liquid rate derivatives) to (i)
discount cash flows risk-neutrally and (ii) supply $r_k(10)$ for the refinancing incentive — it must
give the relevant rates **analytically**. Actual MBS pricing then requires **Monte-Carlo simulation**
(no analytic prices; computationally intensive).

**Financial crisis note.** Sub-prime mortgages (weak-credit borrowers, hidden quality, ARMs with
"teaser" rates) plus the MBS–ABS–CDO–"ABS-CDO" chain (too complex to model) contributed to 2008–09,
alongside moral hazard (brokers, rating agencies), weak regulation, and poor risk management/governance.

---

## Excel workbooks — what each computes

- **Bonds_and_cds_v2.xlsx** — sheets *Pricing*, *Calibration*, *CDS pricing*. Prices defaultable
  coupon bonds from hazard rates $h$, discount factors $d(0,t)$, survival probs $q(t)$ and recovery $R$;
  **calibrates** the hazard-rate vector by minimizing squared bond pricing error across several bonds
  (1–5 yr, various $c,R$); and prices a CDS by building the **premium leg** (PV of premium
  $N\!F\!B$ + accrued $N\!I\!B$) vs **protection leg** (PV of $(1-R)$ contingent payment), and solves
  for the par spread.
- **MBS_Structure.xlsx** — sheets *SingleMortgageCashFlows* (level-payment amortization: beginning
  balance, payment, interest, scheduled principal, ending balance), *Pass-Through* (PSA-driven pool
  cash flows: CPR/SMM, scheduled + prepaid principal, total cash flow, average life), and
  *Sequential-Pay CMO* (allocates the pool's principal/interest across tranches A–D with given par
  amounts and coupon rates).
- *(Course also references `bonds_cds.xlsx` in the lecture notes — same hazard-rate calibration example
  as Bonds_and_cds_v2.xlsx.)*
