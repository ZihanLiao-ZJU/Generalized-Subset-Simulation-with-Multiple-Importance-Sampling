# GSuS-MIS: Generalized Subset Simulation with Multiple Importance Sampling

MATLAB implementation of **generalized subset simulation with multiple importance sampling (GSuS-MIS)** for estimating small failure probabilities of multiple stochastic responses. This code accompanies the paper:

> Xia, W., & Liao, Z. (2025). Enhanced generalized subset simulation with multiple importance sampling for reliability estimation. *Computers & Structures, 313*, Article 107741. [https://doi.org/10.1016/j.compstruc.2025.107741](https://doi.org/10.1016/j.compstruc.2025.107741)

Related implementation: [Python](https://github.com/ZihanLiao-ZJU/Generalized-Subset-Simulation-with-Multiple-Importance-Sampling-Python).

## Contents

- [Overview](#overview)
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

## Requirements

- MATLAB R2020b or later.
- Statistics and Machine Learning Toolbox (`makedist`, `quantile`).
- No Parallel Computing Toolbox is required.

## Repository structure

```text
.
├── drive_GSS_MIS.m          % GSuS-MIS examples and independent repetitions
├── drive_SuS.m              % scalar SuS examples and independent repetitions
├── SeMCS.m                  % sequential sampling engine
├── GSS_MIS.m                % standalone intermediate union distribution
├── SuS.m                    % standalone scalar intermediate distribution
├── Limit-state-functions/   % numerical examples Ex1 and Ex2
├── MCMC-samplers/           % CW_aMH, aCS, ES and HMC samplers
├── Support-Fun/             % MIS estimator and numerical utilities
├── Operations on Logarithmic values/  % log-domain arithmetic
├── LICENSE
└── README.md
```

## Quick start

Edit `example`, `N`, `Nrun` and `SamplerClass` inside `drive_GSS_MIS.m` or `drive_SuS.m`. Select `SamplerClass = @CW_aMH` (the default), or choose another supplied sampler: `CW_aMH`, `aCS`, `HMC` or `ES`.

Both scripts default to Example 1 and 100 runs. `drive_GSS_MIS` uses 440 samples per level and `drive_SuS` uses 240. For Example 2, set the sample counts to 410 and 180, respectively. For a single run, set `Nrun = 1` inside the script.

Run the chosen script from the repository root:

```matlab
drive_GSS_MIS
drive_SuS
```

Use `drive_GSS_MIS` for GSuS-MIS or `drive_SuS` for the scalar SuS comparison.

Use multiple runs to estimate statistics across independent runs.

## Built-in examples

The two numerical examples correspond to the first two examples in the paper. Both use independent standard-normal inputs and samples arranged in a `2 × n` array.

| Example | Description | Dimension | Performance functions | GSuS-MIS samples/level | SuS samples/level |
|---|---|---|---|---|---|
| `Ex1` | Two circular failure regions, radius 0.1, centered at `(0, 0)` and `(0.4, 0.4)` | 2 | 2 | 440 | 240 |
| `Ex2` | Two quadratic and two linear boundaries, in the paper's equation order | 2 | 4 | 410 | 180 |

Only these two simple numerical examples are included. Integer chain lengths can change the realized number of samples in a level slightly; the MIS estimator uses the realized counts.

## Using a custom performance function

1. Copy one of `Limit-state-functions/Ex*.m` and rename its class, for example `MyProblem`.
2. Define `Ndim`, `Nfun`, and `Rdis`, and implement `g = EvlLSF(obj, X)`. Input samples have shape `[Ndim, Nsamples]`; the output has shape `[Nfun, Nsamples]`. Failure is `g < 0`.
3. Follow the driver to construct `model = GSS_MIS(MyProblem)`, `sampler = CW_aMH(model)`, and `simulation = SeMCS(sampler)`, then set `NumSam` before calling `RunIte`.

Set `Rdis = Direct(marginals)`, where `marginals` contains MATLAB probability-distribution objects for independent inputs. `EvlLSF` receives physical coordinates. Returned samples are stored in Gaussian coordinates; use `Rdis.U2X(u)` to obtain physical samples. The built-in standard-normal examples use the identity mapping.

## Outputs

The raw sampling history `g` contains `Nite` sampled states followed by one final candidate state. `Nsze` records the chain and step counts at each sampled level.

The scripts store summary fields in `results`: `example`, `algorithm`, `N`, `repeats`, `mean`, `standard_deviation`, `cov`, `mean_joint_calls`, `unfinished`, and `records`.

Each record contains `probability`, `joint_calls`, `layer_joint_calls`, `actual_level_samples`, and `completed`. `unfinished` is the number of repetitions that did not reach the failure domain. `cov` is the empirical sample standard deviation divided by the mean; it is distinct from optional single-run uncertainty bounds. With one repetition, sample standard deviation and empirical CoV are returned as `NaN`.

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

