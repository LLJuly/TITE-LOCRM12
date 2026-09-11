# TITE-LOCRM12

**TITE-LOCRM12** is an R implementation of the Bayesian phase I/II design proposed in

> **TITE-LOCRM12: A Local Continual Reassessment Method for Drug Combination Optimization Based on Late-Onset Toxicity and Efficacy Outcomes**

**Authors**

- Li Liu
- Ruitao Lin
- Nolan A. Wages

Accepted for publication in **Pharmaceutical Statistics** (2026).

---

## Description

The primary goal of phase I/II combination trials is to identify the optimal biological
dose combination (OBDC), defined as the dose combination that achieves the highest
efficacy while maintaining an acceptable toxicity rate. Several designs have been
proposed for such trials. Among them, the local continual reassessment method for
phase I/II trials (LOCRM12) performs well when compared with other model-based
methods in most of the scenarios studied. However, LOCRM12 assumes that both
toxicity and efficacy outcomes are fully observed before enrolling the next cohort,
an assumption that may not hold in real-world settings, where late-onset outcomes
are increasingly common with the emergence of novel therapies such as molecularly
targeted agents and immunotherapies. To address this limitation, we propose a time-to-event 
extension of LOCRM12 (TITE-LOCRM12) that accommodates both late-onset toxicity and 
efficacy outcomes. The local modeling framework, which uses data
from neighboring dose combinations, is robust and efficient when data are sparse
in complex combination settings. Simulation studies show that the proposed design
can achieve a high percentage of correctly selecting the OBDC while maintaining a
feasible trial duration.

---

## Repository Contents

| File | Description |
|------|-------------|
| `TITE-LOCRM12.R` | Main implementation of the TITE-LOCRM12 design |
| `README.md` | Repository documentation and usage example |

---

## Requirements

The code was developed in R. The following packages are required:

```r
install.packages(c(
  "pocrm",
  "rjags"
))
```
---

## Function Arguments

The main function is

```r
tite.locrm12()
```

The function arguments are summarized below.

| Argument | Description |
|----------|-------------|
| `p.true.tox` | True toxicity probabilities of the dose combinations. |
| `p.true.eff` | True efficacy probabilities of the dose combinations. |
| `target.tox` | Target DLT probability. |
| `pdc.eff` | Efficacy threshold for defining promising dose combinations (PDCs). Dose combinations with efficacy probabilities greater than or equal to `pdc.eff` are considered promising. |
| `ntrial` | Number of simulated trials. |
| `nmax` | Maximum sample size. |
| `cohort.size` | Number of participants enrolled per cohort. |
| `rate` | Participant arrival rate, defined as the expected number of cohorts enrolled during one DLT assessment window. For example, `rate = 3` and `obswin = 60` means that 3 cohorts are expected to arrive within 60 time units. |
| `surv` | Distribution of time to toxicity (`"uniform"` or `"Weibull"`). |
| `surv.eff` | Distribution of time to efficacy (`"uniform"` or `"Weibull"`). |
| `weight.scheme` | Toxicity weighting scheme (`"linear"` or `"piecewise.linear"`). |
| `weight.scheme.eff` | Efficacy weighting scheme (`"linear"`, `"adaptive"`, or `"interval"`). |
| `first.prob` | Probability that a DLT occurs during the first half of the DLT assessment window, conditional on a DLT occurring (used when `surv = "Weibull"`). |
| `first.prob.eff` | Probability that an efficacy outcome occurs during the first half of the efficacy assessment window, conditional on an efficacy event occurring (used when `surv.eff = "Weibull"`). |
| `obswin` | Length of the DLT assessment window. |
| `obswin.eff` | Length of the efficacy assessment window. |
| `cutoff.tox` | Cutoff probability used in the safety stopping rule. |
| `t1` | First cut-off time point for the piecewise linear toxicity weighting scheme. |
| `t2` | Second cut-off time point for the piecewise linear toxicity weighting scheme. |
| `q1` | Weight assigned at follow-up time `t1` in the piecewise linear weighting scheme. |
| `q2` | Weight assigned at follow-up time `t2` in the piecewise linear weighting scheme. |
| `random.seed` | Random seed for the simulation. |
| `var.a` | Variance of the prior distribution for parameter `a`. |


---

## Output

The function returns a list containing the operating characteristics of the TITE-LOCRM12 design across all simulated trials.

| Output | Description |
|--------|-------------|
| `pct.sel` | Selection percentage of each dose combination. |
| `tox` | Average number of DLTs observed at each dose combination. |
| `ntox` | Average total number of DLTs per trial. |
| `eff` | Average number of efficacy responses observed at each dose combination. |
| `neff` | Average total number of efficacy responses per trial. |
| `pts` | Average number of participants treated at each dose combination. |
| `npts` | Average total number of enrolled participants per trial. |
| `nstop` | Number of trials stopped early for excessive toxicity. |
| `stop.pct` | Percentage of trials stopped early for excessive toxicity. |
| `duration` | Average trial duration. |
| `PCS.OBDC` | Percentage of correct selection (PCS) of the optimal biological dose combination (OBDC). |
| `pts.OBDC` | Average number of participants treated at the OBDC(s). |
| `pct.pts.OBDC` | Percentage of participants treated at the OBDC(s). |
| `PCS.PDC` | Percentage of selecting promising dose combinations (PDCs). |
| `pts.PDC` | Average number of participants treated at the PDC(s). |
| `pct.pts.PDC` | Percentage of participants treated at the PDC(s). |
| `PS.overdose` | Percentage of selecting an overly toxic dose combination. |
| `pts.overdose` | Average number of participants treated at overly toxic dose combinations. |
| `pct.pts.overdose` | Percentage of participants treated at overly toxic dose combinations. |


## Example

The following example illustrates how to run the `tite.locrm12()` function under a 3 × 5 dose combination scenario.

```r
# True toxicity probabilities
p.true.tox <- rbind(
  c(0.05, 0.15, 0.30, 0.45, 0.55),
  c(0.15, 0.30, 0.45, 0.55, 0.65),
  c(0.30, 0.45, 0.55, 0.65, 0.75)
)

# True efficacy probabilities
p.true.eff <- rbind(
  c(0.05, 0.25, 0.50, 0.55, 0.60),
  c(0.25, 0.50, 0.55, 0.60, 0.65),
  c(0.50, 0.55, 0.60, 0.65, 0.70)
)

# Run the simulation
output <- tite.locrm12(
  p.true.tox = p.true.tox,
  p.true.eff = p.true.eff,
  target.tox = 0.35,
  pdc.eff = 0.45,
  ntrial = 500,
  nmax = 51,
  cohort.size = 3,
  rate = 3,
  surv = "uniform",
  surv.eff = "uniform",
  weight.scheme = "linear",
  weight.scheme.eff = "adaptive",
  first.prob = 0.3,
  first.prob.eff = 0.3,
  obswin = 60,
  obswin.eff = 180,
  cutoff.tox = 0.85,
  t1 = 20,
  t2 = 40,
  q1 = 0.1,
  q2 = 0.3,
  random.seed = 23219,
  var.a = 1.34
)
print(output)
```

The function returns a list containing the operating characteristics of the design.

Some key outputs are

```r
output$PCS.OBDC
# [1] 65.2

output$PCS.PDC
# [1] 65.2

output$PS.overdose
# [1] 27.8

output$stop.pct
# [1] 0.6

output$duration
# [1] 643.23
```

The complete simulation results can be explored using

```r
str(output)
```

---

## Reference

Liu L, Lin R, Wages NA. 
*TITE-LOCRM12: A Local Continual Reassessment Method for Drug Combination Optimization Based on Late-Onset Toxicity and Efficacy Outcomes.* 
**Pharmaceutical Statistics**.
Accepted for publication (2026).

---
