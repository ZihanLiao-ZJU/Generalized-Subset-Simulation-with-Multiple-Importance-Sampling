# GSS-MIS: Enhanced Generalized Subset Simulation with Multiple Importance Sampling

MATLAB and Python implementations of **generalized subset simulation with multiple importance sampling (GSS-MIS)** for estimating the failure probabilities of multiple performance functions. This repository accompanies:

> Weili Xia and Zihan Liao (2025). **Enhanced generalized subset simulation with multiple importance sampling for reliability estimation.** *Computers & Structures*, 313, 107741. [doi:10.1016/j.compstruc.2025.107741](https://doi.org/10.1016/j.compstruc.2025.107741).

## Overview

The method constructs adaptive intermediate union events, samples each event with component-wise Metropolis MCMC, and combines samples from all levels using the balance-heuristic multiple importance sampling estimator. Both intermediate membership and final failure use the strict inequality **`g < threshold`** and **`g < 0`**, respectively. Historical thresholds, active performance functions, and realized level sample counts are retained for MIS weighting.

The MCMC proposals are uniform component-wise increments, scaled by the sample standard deviation of the current seeds. The conditional quantile probability is `p = 0.1`. The initial samples are independent standard normal; subsequent levels do not retain seeds as an extra initial chain step.

## Repository structure

The layout and naming follow the related [SeMIS repository](https://github.com/ZihanLiao-ZJU/Sequential-multiple-importance-sampling-for-Reliability-Estimation). Each language contains only the first two simple numerical examples.

```text
MATLAB/
  drive_GSS_MIS.m              Experiment driver
  GSS_MIS.m                   Generalized subset simulation
  GSS_MIS_model.m             Intermediate union events and seed selection
  MIS_Rel.m                   Reliability estimator
  SuS.m                       Independent subset-simulation comparison
  Limit-state-functions/      Ex1 and Ex2
  Support-Fun/                MCMC, MIS weights, and numerical helpers
Python/
  drive_GSS_MIS.py             Experiment driver and command-line interface
  GSS_MIS.py                   Generalized subset simulation
  GSS_MIS_model.py             Intermediate events and sample records
  MIS_Rel.py                  Balance-heuristic MIS estimator
  SuS.py                      Independent subset-simulation comparison
  Limit_state_functions/      Ex1 and Ex2
  Support_Fun/                MCMC, random draws, and numerical helpers
  tests/                      Numerical conventions and counting tests
```

Python package directory names use underscores so they can be imported normally. MATLAB numerical helpers that are shared with SeMIS are kept byte-for-byte identical; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Requirements

### MATLAB

- MATLAB R2020b or later.
- Statistics and Machine Learning Toolbox (`makedist`, `quantile`).
- No Parallel Computing Toolbox is required.

### Python

- Python 3.10 or later.
- NumPy 1.22 or later.

```bash
cd Python
python -m pip install -r requirements.txt
```

## Quick start

### MATLAB

From the `MATLAB` directory:

```matlab
result1 = drive_GSS_MIS(1, 100, 2026);
result2 = drive_GSS_MIS(2, 100, 2026);
comparison = SuS(2, 100, 2026);
```

### Python

From the `Python` directory:

```bash
python drive_GSS_MIS.py --example 1 --repeats 100
python drive_GSS_MIS.py --example 2 --repeats 100 --output results_ex2.json
python drive_GSS_MIS.py --example 2 --algorithm SuS --repeats 100
```

Or call the driver directly:

```python
from drive_GSS_MIS import drive_GSS_MIS
result = drive_GSS_MIS(example=1, N=440, repeats=100, seed=2026)
print(result["mean"], result["cov"], result["mean_joint_calls"])
```

Default per-level sample sizes are **440 / 410 for GSS-MIS** and **240 / 180 for the separate SuS comparisons**, for Examples 1 / 2. Integer chain lengths may yield a realized level sample count slightly different from the requested value. The estimator uses the realized counts. Set repetitions to at least two to calculate an empirical CoV.

The languages use their respective random-number generators: the same numerical seed is not intended to produce identical MATLAB and Python trajectories.

## Numerical examples

Both examples have two independent standard-normal inputs and accept samples in a `2 × n` array.

- **Example 1:** two circular failure regions with radius 0.1 and centers `(0, 0)` and `(0.4, 0.4)`.
- **Example 2:** four performance functions comprising two quadratic boundaries and two linear boundaries, ordered as in the paper.

Outputs contain one failure probability per performance function. No vibration examples, structural-response solvers, datasets, or pretrained models are included.

## Outputs and computational cost

The drivers retain all repetitions and summarize mean probability, sample standard deviation, empirical CoV, completion status, and model-evaluation counts. An unfinished run is reported explicitly rather than silently removed. Empirical CoV is the sample standard deviation across repetitions divided by their mean; it is not a single-run error bound.

**One evaluation of the complete performance-function vector at one sample point counts as one joint model call.** Rejected candidates that were evaluated count as calls; unchanged candidates whose stored values are reused do not. MIS post-processing reuses evaluated values. For SuS, the independently evaluated single-performance-function chains are counted separately and their actual calls are summed.

## Tests

Python:

```bash
cd Python
python -m unittest discover -s tests -v
```

The tests cover the numerical examples, strict event membership, rounding and quantile conventions, seed allocation, nonfinite log-domain behavior, and joint-call accounting. The implementations also have trajectory-level regression checks against the source algorithm. MATLAB preparation checks were executed using GNU Octave 9.4 with a compatibility harness; this is not a claim of native MATLAB execution.

## Citation

```bibtex
@article{XiaLiao2025GSSMIS,
  author  = {Xia, Weili and Liao, Zihan},
  title   = {Enhanced generalized subset simulation with multiple importance sampling for reliability estimation},
  journal = {Computers & Structures},
  volume  = {313},
  pages   = {107741},
  year    = {2025},
  doi     = {10.1016/j.compstruc.2025.107741}
}
```

Algorithm-specific references are included next to the corresponding implementation: generalized subset simulation (Li, Ma and Cao, 2015), subset simulation (Au and Beck, 2001), and the balance heuristic (Veach and Guibas, 1995).

## License

The authors' implementation is provided under the [MIT License](LICENSE). Reused SeMIS files retain its [MIT notice](LICENSE-SeMIS); Mo Chen's numerical utilities and derivatives retain their [BSD notice](LICENSE-Mo-Chen). See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for scope and attribution.
