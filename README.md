# GSuS-MIS: Generalized Subset Simulation with Multiple Importance Sampling

MATLAB implementation of **generalized subset simulation with multiple importance sampling (GSuS-MIS)** for estimating small failure probabilities of multiple stochastic responses. This code accompanies the paper:

> Xia, W., & Liao, Z. (2025). Enhanced generalized subset simulation with multiple importance sampling for reliability estimation. *Computers & Structures, 313*, Article 107741. [https://doi.org/10.1016/j.compstruc.2025.107741](https://doi.org/10.1016/j.compstruc.2025.107741)

Related implementation: [Python](https://github.com/ZihanLiao-ZJU/Generalized-Subset-Simulation-with-Multiple-Importance-Sampling-Python).

## Contents

- [Overview](#overview)
- [Common sequential sampling interface](#common-sequential-sampling-interface)
- [Unified MCMC interface](#unified-mcmc-interface)
- [Requirements](#requirements)
- [Repository structure](#repository-structure)
- [Quick start](#quick-start)
- [Built-in examples](#built-in-examples)
- [Using a custom performance function](#using-a-custom-performance-function)
- [Outputs](#outputs)
- [Computational cost](#computational-cost)
- [Code-to-paper mapping](#code-to-paper-mapping)
- [Citation](#citation)
- [License](#license)

## Overview

GSuS-MIS estimates the failure probabilities

$$p_{f,j}=P\{g_j(\mathbf X)<0\}=\int f_{\mathbf X}(\mathbf x)\,I_{\{g_j(\mathbf x)<0\}}\,\mathrm d\mathbf x$$

by constructing a sequence of intermediate union events:

1. Performance-function quantiles define adaptive thresholds, and responses that reach the final failure threshold are removed from the active union.
2. Component-wise Metropolis MCMC generates samples in each intermediate event. Uniform proposal increments are scaled by the sample standard deviation of the current seeds.
3. Samples from all levels are combined with the balance-heuristic multiple importance sampling estimator. Historical thresholds, active responses, and actual level sample counts are retained in the weights.

Intermediate membership is strictly `g < threshold`, and final failure is strictly `g < 0`. The conditional quantile probability is `p = 0.1`; MIS calculations use logarithmic weights for numerical stability.

The common calling sequence is **distribution → sampler → SeMCS → estimator**. `SeMCS` owns the iteration loop; distribution classes define evaluations, intermediate thresholds and convergence; drivers select examples, seed repetitions and summarize results. `MIS` is a separate estimator, not a parent class.

`SuS.m` and `SeMCS.m` are identical in both reliability toolboxes. The SuS drivers choose their numerical policy explicitly: GSuS uses physical coordinates, midpoint thresholds and an unclipped terminal threshold; SeMIS uses standard-normal coordinates, Hazen quantiles and clipping at zero. The final failure event is always strictly `g < 0`.

## Common sequential sampling interface

All five MATLAB toolboxes include the same `SeMCS.m`: choose a sampler, construct `simulation = SeMCS(sampler)`, set a constant `simulation.NumSam`, then call `raw = simulation.RunIte()`. The model in `sampler.BayStc` supplies initialization, level updates, seed selection, and the order of its original stopping checks. The framework stores `x`, `y`, and the evolving distribution states `g`; the estimator processes this history separately.

Fixed MCMC settings, including initial step sizes, acceptance targets, adaptation-group counts and ES options, are properties of the sampler. Values that evolve, such as the adapted `lambda` and HMC trajectory time, are stored in the model's `G`. Different toolboxes retain their own numerical settings. `G.active` carries GSuS's evolving active-response mask together with each historical distribution state.

## Unified MCMC interface

Each standalone toolbox includes `CW_aMH`, `aCS`, `HMC`, and `ES`. The original
outer algorithm and default sampler remain unchanged; choose the sampler in the
driver (Python CLI: `--sampler`). Example: `drive_GSS_MIS(1, 440, 1, 2026, "aCS")`.
All four expose `BayStc` and the same call:
`[u, y, Ncal, BayStc] = SamGen(u_seed, y_seed, Ns)` (Python returns a tuple).
`SamGen(Nchains, Nsteps)` initializes from the Gaussian reference. Arrays are
`[Ndim, Nchains, Nsteps]`; packed outputs are `[log_likelihood; log_reference; Y]`.
Dimension comes from the model, not from the two-dimensional examples.

Distribution state `G` carries model parameters and sampler adaptation. Legacy
numeric thresholds are stored in `G.parameters` (Python `G["parameters"]`);
existing named distribution parameters retain their names. aCS uses initial
`sigma0=1`, `lambda0=0.6`, five adaptation groups, and target acceptance 0.44.
The final `G.lambda` of one level initializes the next, rather than resetting.
ES uses the coordinatewise initial-angle `[0, pi]` version, with no stored seed
at the start of a conditional chain. HMC uses analytic Gaussian rotation in
Gaussian/Nataf reference coordinates, followed by the likelihood-factor MH test;
it is not a generic physical-coordinate gradient-HMC implementation. `G.L` or
`G.L_pi` is the Gaussian reference Cholesky factor, not the covariance matrix.
HMC's final trajectory time is likewise carried in `G.HMC_time`.

Custom nonnormal physical models must be mapped into Gaussian/Nataf coordinates
before these kernels are used. Arbitrary dimensional arrays are supported; the
supplied benchmark formulas retain their documented dimensions. Interface and
smoke checks do not establish statistical reproduction of every paper result.


GSuS now samples and stores Gaussian coordinates `u` in every kernel. `Direct.U2X`
maps these into physical independent marginals before calling the performance
function; transform returned samples with that method when physical coordinates
are required. MATLAB accepts its original probability-distribution objects;
Python `Direct(marginals)` accepts a list of SciPy frozen distributions (for
example `scipy.stats.norm` and `scipy.stats.lognorm`) as `PF.Rdis`. Set `PF.Ndim`
and `PF.Nfun` (or Python `nfun`), and provide `EvlLSF` (or Python `evaluate`).
This is the independent-input Nataf case; no correlated-marginal extension is
implied. Standard-normal benchmark mappings remain the identity. The shared
SuS `physical` option selects this `PF.Rdis` mapping; it does not run MCMC in
untransformed physical space.

## Requirements

- MATLAB R2020b or later.
- Statistics and Machine Learning Toolbox (`makedist`, `quantile`).
- No Parallel Computing Toolbox is required.

## Repository structure

```text
.
├── drive_GSS_MIS.m          % GSuS-MIS examples and independent repetitions
├── drive_SuS.m              % scalar SuS examples and independent repetitions
├── SeMCS.m                  % shared sequential sampling engine
├── GSS_MIS.m                % standalone intermediate union distribution
├── SuS.m                    % standalone scalar intermediate distribution
├── Limit-state-functions/   % numerical examples Ex1 and Ex2
├── MCMC-samplers/           % CW_aMH, aCS, ES, HMC and the common MCMC interface
├── Support-Fun/             % MIS estimator and numerical utilities
├── Operations on Logarithmic values/  % log-domain arithmetic
├── LICENSE
└── README.md
```

## Quick start

Run from the repository root; the drivers add the required subfolders. This runs Example 1 once with a constant `N = 440` per level and the default `CW_aMH` sampler:

```matlab
results = drive_GSS_MIS(1, 440, 1, 2026, 'CW_aMH');
```

The first argument selects Example 1 or 2, the second is the constant `N`, the third is the number of independent runs, the fourth is `seed`, and the last selects `CW_aMH`, `aCS`, `ES` or `HMC`. For Example 2, the default `N` is 410. The scalar SuS comparison uses the same argument order:

```matlab
comparison = drive_SuS(1, 240, 1, 2026, 'CW_aMH');
```

Without arguments, both drivers run Example 1 for 100 repetitions with seed 2026; `drive_GSS_MIS` uses `N = 440` and `drive_SuS` uses `N = 240`. Use multiple runs to estimate the empirical CoV.

## Built-in examples

The two numerical examples correspond to the first two examples in the paper. Both use independent standard-normal inputs and samples arranged in a `2 × n` array.

| Example | Description | Dimension | Performance functions | GSuS-MIS samples/level | SuS samples/level |
|---|---|---|---|---|---|
| `Ex1` | Two circular failure regions, radius 0.1, centered at `(0, 0)` and `(0.4, 0.4)` | 2 | 2 | 440 | 240 |
| `Ex2` | Two quadratic and two linear boundaries, in the paper's equation order | 2 | 4 | 410 | 180 |

Only these two simple numerical examples are included. Integer chain lengths can change the realized number of samples in a level slightly; the MIS estimator uses the realized counts.

Each repetition uses `seed + run - 1` with the `twister` generator. The SuS baseline adds `10000 * PF_index` to the seed for each independent performance function.

## Using a custom performance function

1. Copy one of `Limit-state-functions/Ex*.m` and rename its class, for example `MyProblem`.
2. Define `Ndim`, `Nfun`, and `Rdis`, and implement `g = EvlLSF(obj, X)`. Input samples have shape `[Ndim, Nsamples]`; the output has shape `[Nfun, Nsamples]`. Failure is `g < 0`.
3. Follow the driver to construct `model = GSS_MIS(MyProblem)`, `sampler = CW_aMH(model)`, and `simulation = SeMCS(sampler)`, then set `NumSam` before calling `RunIte`.

The supplied `Direct` distribution helper handles independent marginals. The built-in examples use independent standard-normal variables.

## Outputs

The raw sampling history `g` contains `Nite` sampled states followed by one final candidate state. `Nsze` records the chain and step counts at each sampled level.

All drivers return the same summary fields: `example`, `algorithm`, `N`, `repeats`, `seed`, `mean`, `standard_deviation`, `cov`, `mean_joint_calls`, `unfinished`, and `records`.

Each record contains `seed`, `probability`, `joint_calls`, `layer_joint_calls`, `actual_level_samples`, and `completed`. `unfinished` is the number of repetitions that did not reach the failure domain. `cov` is the empirical sample standard deviation divided by the mean; it is distinct from optional single-run uncertainty bounds. With one repetition, sample standard deviation and empirical CoV are returned as `NaN`.

The intermediate sampling output from `SeMCS.RunIte` retains `x`, `y`, `g`, `Nsam`, `Nsze`, `Ncal`, `Nite`, `flg_cvg`, and the final distribution object `intBay`. GSuS stores the active response history in each `g{level}.active`.

## Computational cost

**One evaluation of the complete performance-function vector at one input point counts as one joint model call.** Evaluated rejected candidates count; unchanged candidates whose values are reused do not. MIS post-processing uses already stored values.

The SuS baseline runs each performance function independently and sums the actual calls from those separate chains. Its total is not divided by the number of performance functions.

## Code-to-paper mapping

| Step | Code | Paper component |
|---|---|---|
| Intermediate union events and level sampling | `SeMCS.RunIte`, `GSS_MIS.UpdObj` | Generalized subset-simulation construction |
| Component-wise Metropolis sampling | `MCMC-samplers/CW_aMH.m` | Conditional sampling within intermediate events |
| Recursive normalizing constants | `MIS` | Intermediate importance-distribution normalization |
| Balance-heuristic MIS estimator | `MIS`, `WgtMis_Blc`, `MisInt` | Combined failure-probability estimator, Eq. (18) |
| Numerical benchmarks | `Ex1`, `Ex2` | First two numerical examples |

## Citation

If you use this toolbox in your research, please cite:

> Xia, W., & Liao, Z. (2025). Enhanced generalized subset simulation with multiple importance sampling for reliability estimation. *Computers & Structures, 313*, Article 107741. [https://doi.org/10.1016/j.compstruc.2025.107741](https://doi.org/10.1016/j.compstruc.2025.107741)

```bibtex
@article{xia2025enhanced,
  title   = {Enhanced generalized subset simulation with multiple importance sampling for reliability estimation},
  author  = {Xia, Weili and Liao, Zihan},
  journal = {Computers & Structures},
  volume  = {313},
  pages   = {107741},
  year    = {2025},
  doi     = {10.1016/j.compstruc.2025.107741}
}
```

Algorithm-specific references are included in the source: Li, Ma and Cao (2015) for generalized subset simulation and numerical benchmarks; Au and Beck (2001) for subset simulation; Veach and Guibas (1995) for the balance heuristic.

## License

Licensed under the [MIT License](LICENSE).

