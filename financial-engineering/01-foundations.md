# Course 1 — Foundations: No-Arbitrage, Binomial Model & Black-Scholes

Columbia "Financial Engineering and Risk Management" (Haugh & Iyengar). This file is a self-contained, formula-dense reference covering: no-arbitrage pricing, forwards/futures, fixed income & swaps, options & put-call parity, the 1-period and multi-period binomial model, risk-neutral pricing, American options & replicating strategies, dividends, and Black-Scholes.

---

## Key formulas — cheat sheet

| Concept | Formula |
|---|---|
| Compound value, $n$ periods | $A(1+r)^n$ |
| Annual rate $r$, $n$ comp./yr, $y$ yrs | $A\left(1+\tfrac{r}{n}\right)^{yn}$ |
| Continuous compounding | $Ae^{ry}$ |
| Present value | $PV(\mathbf{c};r)=\sum_{k=0}^{T}\frac{c_k}{(1+r)^k}$ |
| Discount factor | $d(0,t)=\frac{1}{(1+s_t)^t}$ |
| Perpetuity ($c_k=A$) | $p=\frac{A}{r}$ |
| Annuity ($n$ payments of $A$) | $p=\frac{A}{r}\left(1-\frac{1}{(1+r)^n}\right)$ |
| Forward price (no dividends) | $F=\frac{S_0}{d(0,T)}$ |
| Forward value at $t>0$ | $f_t=(F_t-F_0)\,d(t,T)$ |
| Bond price / YTM $\lambda$ | $P=\sum_{k=1}^{2T}\frac{c}{(1+\lambda/2)^k}+\frac{F}{(1+\lambda/2)^{2T}}$ |
| Floating-rate bond price | $P_f=F$ (par) |
| Swap fixed rate $X$ | $X=\frac{1-d(0,T)}{\sum_{t=1}^{T}d(0,t)}$ |
| Forward rate | $f_{uv}=\left(\frac{(1+s_v)^v}{(1+s_u)^u}\right)^{1/(v-u)}-1$ |
| Spot–forward link | $(1+s_t)^t=\prod_{k=0}^{t-1}(1+f_{k,k+1})$ |
| Call payoff | $\max(S_T-K,0)$ |
| Put payoff | $\max(K-S_T,0)$ |
| Put-call parity (European) | $p_E+S_t=c_E+K\,d(t,T)$ |
| Put-call parity w/ dividends | $p_E+S_t-D=c_E+K\,d(t,T)$ |
| No-arbitrage (binomial) | $d<R<u$ |
| Risk-neutral prob. | $q=\frac{R-d}{u-d}$ |
| 1-period derivative price | $C_0=\frac{1}{R}\left[qC_u+(1-q)C_d\right]=\frac{1}{R}\mathbb{E}^{Q}_0[C_1]$ |
| Multi-period price | $C_0=\frac{1}{R^n}\mathbb{E}^Q_0[\text{Payoff}]$ |
| Risk-neutral prob. w/ dividend $c$ | $q=\frac{R-d-c}{u-d}$ |
| Forward / futures (binomial) | $F_0=G_0=\mathbb{E}^Q_0[S_n]$ |
| BS call | $C_0=S_0e^{-cT}N(d_1)-Ke^{-rT}N(d_2)$ |
| BS $d_1$ | $d_1=\frac{\ln(S_0/K)+(r-c+\sigma^2/2)T}{\sigma\sqrt{T}}$ |
| BS $d_2$ | $d_2=d_1-\sigma\sqrt{T}$ |
| BS calibration | $u=e^{\sigma\sqrt{\Delta t}},\ d=1/u,\ R=e^{r\Delta t}$ |

---

## 1. No-arbitrage principle

**Setup.** A contract: pay price $p$ at $t=0$, receive cash flows $c_k$ at $t=k$, $k=1,\dots,T$ (the $c_k$ may be negative). No-arbitrage bounds $p$ ("no free lunch").

- **Weak no-arbitrage:** if $c_k\ge 0$ for all $k\ge1$ then $p\ge 0$.
- **Strong no-arbitrage:** if $c_k\ge0$ for all $k\ge1$ and $c_\ell>0$ for some $\ell$, then $p>0$.

*Rationale (weak):* if $p<0$, the buyer collects $-p>0$ now and never loses later → free lunch. Competition pushes $p$ up to $\ge0$.

**Implicit assumptions:** liquid markets (enough buyers/sellers), prices visible to all, competition corrects deviations.

### Pricing a simple bond (replication template)
Contract pays $A$ in 1 year; can borrow/lend at rate $r$.
- Portfolio A: buy contract + borrow $A/(1+r)$. Time-1 cash flow $=A-A=0$, price $z=p-\frac{A}{1+r}\ge0\Rightarrow p\ge\frac{A}{1+r}$.
- Portfolio B (reverse): sell contract + lend $A/(1+r)$. $\Rightarrow p\le\frac{A}{1+r}$.
- Together: $p=\frac{A}{1+r}$. This relies on borrowing **and** lending at the same $r$.

If borrow rate $r_B\ge$ lend rate $r_L$: price only bounded, $PV(\mathbf{c};r_B)\le p\le PV(\mathbf{c};r_L)$.

### Linear pricing theorem
If price of $\mathbf{c}_A$ is $p_A$ and price of $\mathbf{c}_B$ is $p_B$, then price of $\mathbf{c}_A+\mathbf{c}_B$ is $p_A+p_B$. (Proof: if cheaper, buy combined / sell parts → costless portfolio with $>0$ income now, zero future flows = free lunch.) Implies $PV(\mathbf{c};r)=\sum_t \frac{c_t}{(1+r)^t}$ by decomposing into single-date cash flows.

---

## 2. Interest rates & fixed income

**Simple interest:** $A(1+nr)$. **Compound:** $A(1+r)^n$. **$n$ comp/yr over $y$ yrs:** $A(1+r/n)^{yn}$. **Continuous:** $\lim_{n\to\infty}A(1+r/n)^{yn}=Ae^{ry}$.

**Present value:** $PV(\mathbf{c};r)=\sum_{k=0}^{N}\frac{c_k}{(1+r)^k}$ (proved via borrow/lend replication on both sides).

**Fixed-income risks** (not truly risk-free): default risk, inflation risk, market risk.

**Perpetuity** ($c_k=A,\ \forall k\ge1$): $p=\sum_{k=1}^{\infty}\frac{A}{(1+r)^k}=\frac{A}{r}$.

**Annuity** ($A$ for $k=1,\dots,n$) = perpetuity − perpetuity starting year $n+1$:
$$p=\frac{A}{r}-\frac{1}{(1+r)^n}\cdot\frac{A}{r}=\frac{A}{r}\left(1-\frac{1}{(1+r)^n}\right).$$

**Bonds.** Face value $F$ (usu. 100/1000), coupon rate $\alpha$ pays $c=\alpha F/2$ semiannually, maturity $T$, S&P quality ratings (AAA…CC). **Yield to maturity** $\lambda$ solves
$$P=\sum_{k=1}^{2T}\frac{c}{(1+\lambda/2)^k}+\frac{F}{(1+\lambda/2)^{2T}}.$$
Lower quality → lower price → higher YTM. YTM summarizes coupon/maturity/quality but is a crude single number.

### Term structure
**Spot rate** $s_t$: rate for a loan maturing in $t$ years; $PV=\frac{A}{(1+s_t)^t}$; discount $d(0,t)=(1+s_t)^{-t}$. Rates rise with term (liquidity preference, expectations, market segmentation).

**Forward rate** $f_{uv}$ (quoted today, lending years $u\to v$): from $(1+s_v)^v=(1+s_u)^u(1+f_{uv})^{v-u}$,
$$f_{uv}=\left(\frac{(1+s_v)^v}{(1+s_u)^u}\right)^{1/(v-u)}-1,\qquad (1+s_t)^t=\prod_{k=0}^{t-1}(1+f_{k,k+1}).$$

### Floating-rate bonds
Rate $r_k$ over $[k,k+1)$ is known only at time $k$. Floating bond pays coupon $r_{k-1}F$ at each $k$ and face $F$ at $n$. Price the random coupon $r_{k-1}F$ by a portfolio (borrow/relend to cancel the random term): set $\alpha=\frac{F}{(1+r_0)^{k-1}}$, giving deterministic $c_k=Fr_0$, so
$$p_k=\frac{Fr_0}{(1+r_0)^k},\qquad P_f=\frac{F}{(1+r_0)^n}+\sum_{k=1}^{n}\frac{Fr_0}{(1+r_0)^k}=F.$$
**A floating-rate bond is always worth its face value $F$ (trades at par).**

---

## 3. Forwards & futures

**Forward contract:** buyer has the right *and obligation* to buy a specified amount of an asset at time $T$ at a forward price $F$ set at $t=0$. Long payoff at $T$: $f_T=S_T-F$. $F$ is set so $f_0=0$.

**Short selling:** borrow shares from broker, sell, later "cover" (buy back). Profit if price drops; loss unbounded if price rises.

**No-arbitrage forward price (no dividends/storage):** buy contract, short the asset, lend $S_0$ to $T$. Net cost 0, deterministic time-$T$ flow $\frac{S_0}{d(0,T)}-F=0$:
$$\boxed{F=\frac{S_0}{d(0,T)}}.$$
$F>S_0$ = **cost of carry**. Example: $S_0=50$, 6-mo rate 4% semiannual, $d(0,.5)=\frac{1}{1+0.02}=0.9804$, $F=50/0.9804=51.0$.

**Forward value at $t>0$:** with $F_0,F_t$ the forward prices at $0$ and $t$ (same delivery $T$), short an $F_t$-contract, long the $F_0$-contract → deterministic flow $F_t-F_0$:
$$f_t=(F_t-F_0)\,d(t,T).$$

**Futures vs forwards (binomial result, see §7):** a futures contract always costs nothing; its "price" $F_t$ only sets the marked-to-market cash flow $\pm(F_t-F_{t-1})$ each period. $F_n=S_n$. In the binomial model with constant $R$, $F_0=G_0=\mathbb{E}^Q_0[S_n]$ (forward price = futures price); not equal in general.

---

## 4. Swaps

**Definition:** contracts that transform one kind of cash flow into another (plain-vanilla fixed↔floating, commodity, currency). Used to change cash-flow nature and exploit comparative advantage across markets.

**Comparative advantage example.** A: fixed 4.0%, floating LIBOR+0.3%; B: fixed 5.2%, floating LIBOR+1.0%. A borrows fixed at 4%, swaps with B (A pays LIBOR, receives 3.95%); B borrows floating at LIBOR+1%, swaps (B pays 3.95%, receives LIBOR). Effective: A pays LIBOR+0.05%, B pays 4.95% — both gain. With an intermediary (A receives 3.93%, B pays 3.97%) the intermediary keeps 0.04% (4 bps) for counterparty risk/service.

**Pricing an interest-rate swap.** Notional $N$, floating $r_{t-1}$. Long (A) receives $Nr_{t-1}$, pays $NX$ for $t=1,\dots,T$. Value to A = (floating-bond cash flow − face value):
$$V_A=N\big(1-d(0,T)\big)-NX\sum_{t=1}^{T}d(0,t).$$
Par swap rate ($V_A=0$):
$$\boxed{X=\frac{1-d(0,T)}{\sum_{t=1}^{T}d(0,t)}}.$$

---

## 5. Options

**Definitions.** European call/put: right (not obligation) to buy/sell 1 unit of underlying at strike $K$ at expiration $T$. American: same but at *any* time up to $T$.

**Payoffs at $T$:** call $\max(S_T-K,0)$; put $\max(K-S_T,0)$ — nonlinear, so a model for the underlying is needed.

**Intrinsic value** at $t\le T$: call $\max(S_t-K,0)$, put $\max(K-S_t,0)$.
- Call: in/at/out of money for $S_t>K$ / $=K$ / $<K$. Put: opposite.

**Notation:** $c_E,p_E$ European; $c_A,p_A$ American prices, written $\cdot(t;K,T)$.

### Put-call parity (European, non-dividend)
$$p_E(t;K,T)+S_t=c_E(t;K,T)+K\,d(t,T).$$
*Proof (zero-cost portfolio):* buy call, sell put, short underlying, lend $K\,d(t,T)$. Time-$T$ flow $\max(S_T-K,0)-\max(K-S_T,0)-S_T+K=0$; hence time-$t$ flow $0$.

With dividends (PV of dividends $=D$): $\ p_E+S_t-D=c_E+K\,d(t,T)$.

### Bounds
$c_A\ge c_E$, $p_A\ge p_E$.
- $c_E=\max\{S_t+p_E-K d(t,T),0\}\ge\max\{S_t-K d(t,T),0\}$; and $c_E\le S_t$.
- $p_E=\max\{K d(t,T)+c_E-S_t,0\}\ge\max\{K d(t,T)-S_t,0\}$; and $p_E\le K d(t,T)$.

**American call, non-dividend:** $c_A\ge c_E\ge\max\{S_t-Kd(t,T),0\}>\max\{S_t-K,0\}$ = exercise value. Since the price always exceeds the exercise value, **never optimal to early-exercise an American call on a non-dividend stock**: $c_A=c_E$. For American puts the bound $p_A\ge\max\{Kd(t,T)-S_t,0\}$ does *not* dominate exercise value $\max\{K-S_t,0\}$ — early exercise can be optimal (esp. deep in-the-money puts).

---

## 6. The 1-period binomial model

**Dynamics.** $S_0\to uS_0$ (prob. $p$) or $dS_0$ (prob. $1-p$). Cash account: \$1 → \$$R$ per period ($R$ = gross risk-free rate). Short-sales allowed. Canonical numbers: $S_0=100$, $u=1.07$ ($uS_0=107$), $d=0.9346$ ($dS_0=93.46$), $R=1.01$.

**Type A / Type B arbitrage** (needed once randomness enters):
- *Type A:* portfolio with $V_0<0$ and $V_1\ge0$ (immediate reward, no future loss).
- *Type B:* portfolio with $V_0\le0$, $V_1\ge0$, $V_1\ne0$ (positive payoff prob., zero loss prob.).

**No-arbitrage theorem:** no arbitrage $\iff d<R<u$.
- If $R<d<u$: borrow $S_0$, buy stock → type-B arbitrage.
- If $d<u<R$: short stock, invest proceeds in cash → type-B arbitrage.

### Replicating portfolio & pricing
Hold $x$ shares + \$$y$ cash. To replicate payoff $C_1(S_1)$ (values $C_u,C_d$):
$$uS_0 x+Ry=C_u,\qquad dS_0 x+Ry=C_d.$$
Then $C_0=xS_0+y$. Worked example: option payoff $\max(S_1-102,0)$, so $C_u=5$, $C_d=0$:
$$107x+1.01y=5,\quad 93.46x+1.01y=0\ \Rightarrow\ x=0.3693,\ y=-34.1708.$$
$C_0=0.3693\times100-34.1708=2.76$ (arbitrage-free price). $y<0$ = borrowing; $x<0$ = short stock. Price does **not** depend on buyer/seller utility.

### Risk-neutral pricing
Solving generally:
$$C_0=\frac{1}{R}\left[\frac{R-d}{u-d}C_u+\frac{u-R}{u-d}C_d\right]=\frac{1}{R}\big[qC_u+(1-q)C_d\big]=\frac{1}{R}\mathbb{E}^Q_0[C_1],$$
$$\boxed{q=\frac{R-d}{u-d}}.$$
No-arbitrage ($d<R<u$) $\iff q\in(0,1)$. $(q,1-q)$ are the **risk-neutral probabilities** ($Q$-measure). The real up-probability $p$ does **not** enter the price.

**The "paradox":** stock ABC ($p=0.99$) and XYZ ($p=0.01$) with the same $S_0=100$, $u$-node 110, $d$-node 90 have the *same* call price ($\approx\$4.80$ in the worked example, $K=100$). Reason: $S_0$ already embeds the market's expectations; pricing rests on hedging/replication, not speculation; only $q$ (from $u,d,R$) matters.

*Note:* naive expected-discounted-payoff pricing $\mathbb{E}^P_0[R^{-n}\,\text{payoff}]$ is wrong — illustrated by the **St. Petersburg paradox** (fair coin, payoff $\$2^n$ on first head at toss $n$): $\mathbb{E}[\text{payoff}]=\sum 2^n 2^{-n}=\infty$, yet nobody pays $\infty$. Bernoulli's log-utility fix gives $\sum \log(2^n)2^{-n}=\log2\sum n2^{-n}<\infty$, but the right answer is replication, not utility.

---

## 7. Multi-period binomial model

**Concept:** an $n$-period tree is a sequence of embedded 1-period models. $S_0$ moves by $u$ or $d$; $R$ constant; no-arbitrage $d<R<u$. Recombining lattice example ($S_0=100$, $u=1.07$, $d=0.9346$, $R=1.01$): $t=3$ terminal prices $122.50,\ 107,\ 93.46,\ 81.63$.

**Constant $q$:** since $u,d,R$ are constant, $q=\frac{R-d}{u-d}$ is identical at every node.

**Terminal distribution** (binomial): $P(k\text{ up-moves})=\binom{n}{k}q^k(1-q)^{n-k}$. For $n=3$: $S_{uuu}$ 1 path ($q^3$), $S_{uud}$ 3 paths, $S_{udd}$ 3 paths, $S_{ddd}$ 1 path.

### Pricing algorithm — backward induction (recurrence à rebours)
1. At maturity $T=n$: compute payoff at each node, e.g. $\max(S_T-K,0)$.
2. Work backward, at each node:
$$C_t=\frac{1}{R}\big[q\,C_u+(1-q)\,C_d\big].$$
3. Equivalent one-shot formula:
$$C_0=\frac{1}{R^n}\mathbb{E}^Q_0[\text{Payoff}]=\frac{1}{R^n}\sum_{j=0}^{n}\binom{n}{j}q^j(1-q)^{n-j}\,\text{Payoff}(S_0u^jd^{n-j}).$$

**First Fundamental Theorem of Asset Pricing:** No arbitrage $\iff$ a risk-neutral measure $Q$ exists ($q>0$). The binomial model is consistent because it admits such a $Q$.

### Forwards & futures in the binomial model
**Forward** expiring at $n$, price $G_0$ chosen so contract worth 0:
$$0=\mathbb{E}^Q_0\!\left[\frac{S_n-G_0}{R^n}\right]\Rightarrow G_0=\mathbb{E}^Q_0[S_n]\quad(8).$$
**Futures.** $F_n=S_n$. Futures always costs zero; characterized as a security worth zero that pays "dividend" $(F_t-F_{t-1})$ at each $t$. From $0=\mathbb{E}^Q_{t}[(F_{t+1}-F_t)/R]$: $F_t=\mathbb{E}^Q_t[F_{t+1}]$, and by iterated expectations $F_t=\mathbb{E}^Q_t[F_n]$ — the futures price is a **$Q$-martingale**. Hence
$$F_0=\mathbb{E}^Q_0[S_n]\quad(9),\qquad F_0=G_0\text{ in the binomial model (not in general)}.$$
Both (8) and (9) hold with or without dividends — dividends enter only through $Q$.

---

## 8. American options & optimal stopping

**Valuation:** price like a European option by backward induction, but at each node take the max of *continuation value* and *immediate exercise value*:
$$V_t=\max\Big(\underbrace{\tfrac{1}{R}[qV_u+(1-q)V_d]}_{\text{continuation}},\ \underbrace{\text{exercise value}}_{\text{e.g. }S_t-K\text{ (call) or }K-S_t\text{ (put)}}\Big).$$

**Worked American put** ($K=100$, $R=1.01$, expiry $t=3$, same lattice). Terminal put values: node 122.50→0, 107→0, 93.46→6.54, 81.63→18.37. Backward induction gives node values (top→bottom at each date): $t=2$: $0,\ 0,\ 12.66$ — where the lower node is $\max[12.66,\ \frac{1}{R}(q\cdot6.54+(1-q)\cdot18.37)]$ (early exercise binds, $87.34$ node: intrinsic $12.66$ > continuation); $t=1$: $1.26,\ 7.13$; $t=0$: **price $=3.82$**.

**Optimal-stopping analogy (dice game).** Throw a fair die up to 3 times; "stop" to collect the face value.
- 1 throw left: fair value $=3.5$.
- 2 throws left: keep $\{4,5,6\}$, reroll $\{1,2,3\}$: $\frac{1}{6}(4+5+6)+\frac{1}{2}\cdot3.5=4.25$.
- 3 throws: reroll if $<4.25$, i.e. keep $\{5,6\}$ on first throw, else compare to 4.25.
Same logic as American exercise: compare intrinsic (stop) to expected future (continue).

---

## 9. Replicating strategies (dynamic replication)

**Trading strategy** $\theta_t=(x_t,y_t)$: $x_t$ shares and $y_t$ units of cash account held between $t-1$ and $t$ (known at $t-1$, a random/adapted process). Cash account $B_t=R^t$ (with $B_0=1$).

**Value process:**
$$V_t=\begin{cases}x_1S_0+y_1B_0 & t=0\\ x_tS_t+y_tB_t & t\ge1.\end{cases}$$

**Self-financing:** changes in $V_t$ come only from trading gains/losses (no cash injected/withdrawn):
$$V_t=x_{t+1}S_t+y_{t+1}B_t,\quad t=1,\dots,n-1.$$
**Proposition:** if self-financing, $V_{t+1}-V_t=x_{t+1}(S_{t+1}-S_t)+y_{t+1}(B_{t+1}-B_t)$.

**Risk-neutral price ≡ price of replicating strategy.** A self-financing strategy can be built to replicate the option payoff (*dynamic replication*); its initial cost must equal the option price (else arbitrage), and equals the backward-induction / risk-neutral price. At every node, option value = replicating-portfolio value.

**Worked European call replication** (lattice above, payoffs at $t=3$: 22.5, 7, 0, 0; option prices $t=0$ down to nodes: $C_0=6.57$). Replicating $[x_t,y_t]$ at sample nodes: $t=0$: $[0.598,-53.25]$; $t=1$ up $[0.802,-74.84]$, down $[0.305,-26.11]$; $t=2$: $[1,-97.06]$, $[0.517,-46.89]$, $[0,0]$. Check: $0.802\times107+(-74.84)\times1.01=10.23$ = option value at the upper $t=1$ node.

---

## 10. Dividends in the binomial model

**Proportional dividend** $cS_0$ paid at $t=1$. Stock node values become $uS_0+cS_0$ and $dS_0+cS_0$. No-arbitrage: $d+c<R<u+c$. Replicate:
$$uS_0x+cS_0x+Ry=C_u,\qquad dS_0x+cS_0x+Ry=C_d,$$
$$C_0=\frac{1}{R}\left[\frac{R-d-c}{u-d}C_u+\frac{u+c-R}{u-d}C_d\right]=\frac{1}{R}\big[qC_u+(1-q)C_d\big],\quad \boxed{q=\frac{R-d-c}{u-d}}.$$
Dividends reduce $q$ → calls cheaper, puts more expensive. Multi-period: proportional dividend $cS_i$ each period gives identical $q$ in every embedded model.

**No-dividend identity:** $S_0=\mathbb{E}^Q_0[S_n/R^n]$ (risk-neutral pricing of a $K=0$ call). **With dividends** $D_i$ at time $i$ (ex-dividend price $S_n$):
$$S_0=\mathbb{E}^Q_0\!\left[\frac{S_n}{R^n}+\sum_{i=1}^{n}\frac{D_i}{R^i}\right].$$
*(Proof:* view each dividend $D_i$ as a separate security of value $P_i=\mathbb{E}^Q_0[D_i/R^i]$; the stock owner holds the portfolio $\sum_i P_i+\mathbb{E}^Q_0[S_n/R^n]=S_0$.)*

---

## 11. Black-Scholes model

**Assumptions:**
1. Continuously-compounded interest rate $r$.
2. Geometric Brownian motion: $S_t=S_0\,e^{(\mu-\sigma^2/2)t+\sigma W_t}$, $W_t$ standard Brownian motion.
3. Continuous dividend yield $c$.
4. Continuous trading, no transaction costs, short-selling allowed.

**Black-Scholes call formula** (strike $K$, maturity $T$):
$$C_0=S_0e^{-cT}N(d_1)-Ke^{-rT}N(d_2),$$
$$d_1=\frac{\ln(S_0/K)+(r-c+\sigma^2/2)T}{\sigma\sqrt{T}},\qquad d_2=d_1-\sigma\sqrt{T},$$
where $N(d)=P(\mathcal{N}(0,1)\le d)$ is the standard normal CDF. The drift $\mu$ does **not** appear (analogue of $p$ being irrelevant in the binomial model).

**European put** via put-call parity: $P_0+S_0e^{-cT}=C_0+Ke^{-rT}$.

**Risk-neutral form:** $C_0=\mathbb{E}^Q_0\!\left[e^{-rT}\max(S_T-K,0)\right]$, where under $Q$: $S_t=S_0e^{(r-c-\sigma^2/2)t+\sigma W_t}$. Derived by a replicating-strategy argument analogous to the binomial case.

### Binomial → Black-Scholes (calibration & convergence)
Specify the binomial model via BS parameters $r,\sigma$ over $T$ with $n$ periods, $\Delta t=T/n$:
$$R_n=e^{r\,T/n},\quad u_n=e^{\sigma\sqrt{T/n}},\quad d_n=1/u_n,\quad R_n-c_n\approx 1+r\tfrac{T}{n}-c\tfrac{T}{n},$$
$$q_n=\frac{e^{(r-c)T/n}-d_n}{u_n-d_n}.$$
European price in calibrated binomial model:
$$C_0=\frac{1}{R_n^n}\sum_{j=0}^{n}\binom{n}{j}q_n^j(1-q_n)^{n-j}\max(S_0u_n^jd_n^{n-j}-K,0),$$
$$=\frac{S_0}{R_n^n}\sum_{j=\eta}^{n}\binom{n}{j}q_n^j(1-q_n)^{n-j}u_n^jd_n^{n-j}-\frac{K}{R_n^n}\sum_{j=\eta}^{n}\binom{n}{j}q_n^j(1-q_n)^{n-j},$$
with $\eta=\min\{j:S_0u_n^jd_n^{n-j}\ge K\}$. As $n\to\infty$ ($\Delta t\to0$), $C_0\to$ Black-Scholes.

**History:** Bachelier (1900, Brownian motion for Paris Bourse, predated Einstein by 5 yrs); Samuelson (1965, GBM, Merton's advisor); Itô (1950s, stochastic calculus / Itô's Lemma; Doeblin 1940 independently); Black-Scholes-Merton (early 1970s); also Thorp, Cox-Ross, Harrison-Kreps.

---

## 12. Worked example — European put on a futures contract

Many of the most liquid options are options on futures (S&P 500, Eurostoxx 50, FTSE 100, Nikkei 225) where the underlying index isn't directly traded.

**Parameters:** $S_0=100$, $n=10$ periods, $r=2\%$, $c=1\%$, $\sigma=20\%$, futures expiry = option expiry = $T=0.5$ yr, strike $K=100$.

**Method:**
1. Build the futures-price lattice: set $F_n=S_n$ at maturity, then roll back via $F_t=\mathbb{E}^Q_t[F_{t+1}]$.
2. Apply backward induction on $F_t$ with the put payoff $\max(K-F_T,0)$.

**Result:** put value $\approx 5.21$.

*Practical note:* liquid options are priced by market supply/demand (which amounts to setting the **implied volatility** $\sigma$). Models are needed to **hedge** these and to price exotic/illiquid derivatives.

---

## 13. The `BinomialModel_Public.xlsx` workbook

Companion spreadsheet implementing this course's binomial pricing. It calibrates a binomial lattice from Black-Scholes parameters (inputs: spot price, strike, expiration, periods, volatility, $r$, dividend yield) using $u=e^{\sigma\sqrt{\Delta t}}$, $d=1/u$, $R=e^{r\Delta t}$, $q=\frac{e^{(r-c)\Delta t}-d}{u-d}$, and prices by backward induction so that prices converge to Black-Scholes as $n\to\infty$. Sheets:
- **EuropeanCall_EG**, **AmericanPut_EG**, **OptionsOnFuturesEG** — worked examples for each option type.
- **10PeriodBinomialModel** — the calibrated 10-period engine.
- **StockLattice / OptionLattice** — underlying-price and option-value trees (backward induction).
- **FuturesLattice / FuturesOptionLattice** — futures-price tree ($F_t=\mathbb{E}^Q_t[F_{t+1}]$, $F_n=S_n$) and the option-on-futures valuation (the §12 example, put $\approx5.21$).
