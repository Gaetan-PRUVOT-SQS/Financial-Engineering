# Financial Engineering & Risk Management

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT) [![Status: Active](https://img.shields.io/badge/Status-Active-success.svg)]()

Ce dépôt sert de référence technique et de boîte à outils pour la modélisation quantitative, la tarification de produits dérivés et la gestion des risques. L'objectif est de fournir une base mathématique et algorithmique solide pour l'ingénierie financière, en privilégiant une approche purement logique et dénuée de biais émotionnels.

## 🧠 Philosophie du Projet

La conception et l'exécution de stratégies financières nécessitent une analyse froide et algorithmique. Ce projet est piloté par les spécifications (Spec-Driven) et se concentre sur les méthodes quantitatives classiques (pré-Machine Learning) : tarification sans arbitrage, optimisation moyenne-variance, et calibration de modèles par transformée de Fourier rapide (FFT). 

Ces fondations mathématiques ont vocation à être implémentées et testées rigoureusement dans des environnements de développement modernes (Python pour le backend quantitatif, React/Next.js/Angular pour les interfaces de monitoring).

## 📚 Table des Matières et Contenu

Le dépôt est structuré en cinq modules de référence autonomes contenant les formules, algorithmes et exemples numériques :

* **[01-foundations.md](./01-foundations.md) : Fondations et Modèles Binomiaux**
    * Principe de non-arbitrage, contrats à terme (forwards/futures), swaps.
    * Parité put-call, modèle binomial à n-périodes, et réplication dynamique.
    * Modèle de Black-Scholes et évaluation neutre au risque.
* **[02-credit-derivatives.md](./02-credit-derivatives.md) : Structure par Terme et Dérivés de Crédit**
    * Modèles de taux à court terme (treillis binomiaux, Ho-Lee, BDT).
    * Credit Default Swaps (CDS) et tarification du spread au pair.
    * Mathématiques hypothécaires et titrisation (MBS, CMO, PO/IO).
* **[03-portfolio-optimization.md](./03-portfolio-optimization.md) : Optimisation de Portefeuille et Exécution**
    * Frontière efficiente de Markowitz, CAPM et ratio de Sharpe.
    * Gestion du risque en pratique : erreur d'estimation, contraintes, VaR/CVaR.
    * Biais comportementaux et statistiques en allocation d'actifs.
    * Exécution optimale et impact sur le marché (trajectoires d'Almgren-Chriss).
* **[04-advanced-pricing.md](./04-advanced-pricing.md) : Options Avancées et Crédit Structuré**
    * Les "Grecques", couverture (hedging) delta-gamma-vega en pratique.
    * La surface de volatilité implicite (skew, structure par terme, Breeden-Litzenberger).
    * Tranches de CDO et modèle de copule gaussienne à un facteur.
    * Options réelles (flexibilité opérationnelle via les arbres binomiaux).
* **[05-computational.md](./05-computational.md) : Méthodes Computationnelles et Calibration**
    * Tarification d'options par transformée de Fourier (FFT) et méthode de Carr-Madan.
    * Fonctions caractéristiques pour les modèles Heston, Variance Gamma et GBM.
    * Algorithmes de calibration de surface (optimisation RMSE via BFGS, Nelder-Mead).
    * Calibration des modèles de taux d'intérêt (Vasicek, CIR) sur les courbes LIBOR/Swap.

## 🛠 Cas d'Usage (Quand utiliser ce dépôt)

**Utilisez ces références pour :**
* Tarifer de manière analytique ou sur un treillis une option, un forward, ou un dérivé de taux.
* Calibrer un modèle de structure par terme ou évaluer un CDS / une tranche de CDO.
* Concevoir une allocation de portefeuille robuste et calculer une trajectoire d'exécution optimale.
* Implémenter des algorithmes de tarification via FFT dans un backend de trading systématique.

**À ne pas utiliser pour :**
* Le Machine Learning financier pur (qui nécessite une approche différente centrée sur les données).
* L'infrastructure pure de backtesting sans base de pricing.

## 📝 Licence

MIT

---

**Auteur :** Gaetan PRUVOT  
*Spec-Driven Quant Engineer | Systematic Trader | Software Tester | Certified IBM AI Developer*
