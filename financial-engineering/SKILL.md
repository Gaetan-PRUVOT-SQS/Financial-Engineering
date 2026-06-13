---
name: financial-engineering
description: Use when pricing or risk-managing derivatives and fixed income with classical (non-ML) quant methods — no-arbitrage and binomial pricing, Black-Scholes, term-structure / interest-rate lattices, CDS and credit derivatives, mortgages and MBS, mean-variance portfolio optimization, optimal execution, the implied volatility surface, CDO / structured credit, real options, or transform (FFT) pricing and model calibration. Self-contained reference with formulas and worked examples.
---

# Financial Engineering & Risk Management

## Overview

Self-contained quantitative reference distilled from Columbia's *Financial Engineering and Risk Management* specialization (Haugh & Iyengar). Covers the classical (pre-ML) toolkit: arbitrage-free pricing, the binomial method, Black-Scholes, fixed income and term-structure models, credit derivatives, mortgages/MBS, portfolio optimization, optimal execution, the volatility surface, structured credit, real options, and computational pricing/calibration.

**Core principle:** price by **no-arbitrage / risk-neutral expectation** — replicate the payoff, or take the expectation under the risk-neutral measure $Q$ and discount. Everything below is a special case of this idea.

**This skill is self-contained.** The formulas, algorithms, and worked numerical examples are written into the five reference files — there is no external repo to open. Distinct from and complementary to `quant-ml-finance` (which covers *machine-learning* methods in finance).

## When to Use

- Pricing/hedging an option, forward, swap, or fixed-income derivative analytically or on a lattice
- Building or calibrating a term-structure / short-rate model; pricing caps, floors, swaptions
- Valuing a CDS, CDO tranche, or mortgage / MBS
- Mean-variance portfolio optimization, efficient frontier, CAPM, practical allocation
- Designing an optimal-execution / market-impact trading schedule
- Interpreting or constraining an implied-volatility surface
- Pricing via Fourier/FFT transform methods or calibrating a model to market option prices

**When NOT to use:** ML-based finance (use `quant-ml-finance`); pure data/backtesting infrastructure.

## Topic Router

Load the matching reference file for full formulas, algorithms, and worked examples.

| Topic | Reference file |
|---|---|
| No-arbitrage, PV/bonds/YTM, forwards & futures, swaps, put-call parity | **01-foundations.md** |
| 1-period & multi-period binomial model, risk-neutral pricing, backward induction | **01-foundations.md** |
| American options & optimal stopping, dynamic replication, dividends | **01-foundations.md** |
| Black-Scholes formula, put-on-futures worked example | **01-foundations.md** |
| Binomial short-rate lattice, elementary/state prices, martingale pricing of rates | **02-credit-derivatives.md** |
| Fixed-income derivatives on the lattice; Ho-Lee / BDT calibration; swaptions | **02-credit-derivatives.md** |
| Defaultable bonds, hazard-rate model, **CDS mechanics & par-spread formula** | **02-credit-derivatives.md** |
| Mortgage math (level payment, PSA/CPR/SMM), pass-throughs, PO/IO, CMO tranching | **02-credit-derivatives.md** |
| Mean-variance optimization, efficient frontier, two-fund/one-fund, tangency, Sharpe | **03-portfolio-optimization.md** |
| CAPM, Capital Market Line; practical MVO (estimation error, shrinkage, constraints) | **03-portfolio-optimization.md** |
| Behavioral/statistical allocation biases; VaR/CVaR | **03-portfolio-optimization.md** |
| Optimal execution, market impact, Almgren-Chriss trading trajectory | **03-portfolio-optimization.md** |
| Greeks & delta/gamma/vega hedging in practice; discrete hedging P&L | **04-advanced-pricing.md** |
| Implied volatility surface (skew, term structure), Breeden-Litzenberger density, no-arb constraints | **04-advanced-pricing.md** |
| CDOs, tranches, **one-factor Gaussian copula**, expected tranche loss | **04-advanced-pricing.md** |
| Real options, operational flexibility, the Simplico gold-mine lattice | **04-advanced-pricing.md** |
| Transform / Fourier / **FFT pricing** (Carr-Madan); char. functions (GBM, Heston, VG) | **05-computational.md** |
| **Model calibration** (RMSE objective, BFGS/Nelder-Mead), with code | **05-computational.md** |
| Interest-rate instruments, $P(t,T)$/LIBOR/forward/swap, Vasicek & CIR calibration | **05-computational.md** |

## How to Use

1. Identify whether the task is **pricing**, **risk/portfolio**, or **calibration/computation**.
2. Open the matching reference file above for the precise formula, lattice algorithm, or code.
3. Default to the **risk-neutral pricing** recipe: identify the replicating portfolio or risk-neutral probabilities $q$, take $V_0 = \frac{1}{(1+r)^T}\,\mathbb{E}^Q[\text{payoff}]$ (or its lattice/continuous analogue).
4. For calibration tasks, set up the objective (sum of squared model-vs-market errors) and minimize over the parameter vector — see 05-computational.md for the worked BFGS loop.

## Common Pitfalls

- **Confusing real-world $p$ with risk-neutral $q$.** Price under $q$; $q$ comes from no-arbitrage ($q=\frac{(1+r)-d}{u-d}$ in the binomial model), not from estimated drift.
- **Forgetting early exercise.** American options need backward induction with an exercise check at every node (01-foundations.md §8).
- **Term-structure model must be calibrated to the observed curve** before pricing derivatives on it (02-credit-derivatives.md §3).
- **Naive mean-variance is estimation-error-amplifying.** Use shrinkage/constraints; tiny changes in expected returns swing weights wildly (03-portfolio-optimization.md §6).
- **Treating implied vol as constant.** The surface has skew and term structure; sticky-strike vs sticky-delta dynamics change your hedge (04-advanced-pricing.md §3).
- **FFT grid constraint.** Carr-Madan requires $\lambda\eta = 2\pi/N$; mismatching the grid corrupts prices (05-computational.md §1).

## Provenance

Distilled from the 5-course Columbia FE&RM specialization (≈99 lecture decks + Excel models). The original source folder is not retained; this skill is the durable artifact. Notation follows the course conventions.
