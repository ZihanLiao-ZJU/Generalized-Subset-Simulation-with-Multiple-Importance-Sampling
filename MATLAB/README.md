# GSS-MIS — MATLAB

Enhanced generalized subset simulation with multiple importance sampling for reliability estimation. This implementation includes two numerical benchmarks: two circular failure domains and four performance functions in two independent standard-normal variables.

## Requirements
MATLAB R2020b or later, with Statistics and Machine Learning Toolbox (`makedist`, `quantile`). No parallel toolbox is required. The code uses double precision and samples arranged in columns.

## Quick start
From this directory:
```matlab
results = drive_GSS_MIS(1);           % Example 1, 100 runs
results = drive_GSS_MIS(2, 100, 2025); % Example 2, explicit seed
comparison = SuS(1, 100, 2025);        % Scalar SuS comparison
test_GSS_MIS
```
Example 1 uses 440 samples per level; Example 2 uses 410. The level probability is 0.1. Chain rounding may change the actual number generated in a level. Each run uses seed + run - 1 with the `twister` generator.

The scalar `SuS` comparator uses 240 samples per level for Example 1 and 180 for Example 2, running each performance function separately. Its `Ncal` array counts scalar evaluations, and `mean_total_evaluations` sums those costs. Scalar seeds are seed + run - 1 + 10000 * PF index.

## Structure
- `drive_GSS_MIS.m`: configurable experiment driver.
- `GSS_MIS.m`: level-by-level sample generation.
- `GSS_MIS_model.m`: union-indicator thresholds and active response sets.
- `MIS_Rel.m`: joint evaluations and multiple-importance-sampling estimator.
- `Limit-state-functions/Ex1.m`, `Ex2.m`: numerical benchmarks.
- `Support-Fun/`: sampler, independent-marginal distribution, MIS weights and shared utilities.

## Outputs
`Pf` contains one column per run. `mean_Pf`, `std_Pf`, and `CoV` summarize the independent runs. `Ncal` counts **joint vector evaluations**: evaluating all performance functions at one input counts once. `converged` reports whether all response thresholds were reached; `termination_flag` preserves the framework's separate termination flag. Unconverged runs are retained and flagged, not silently retried or excluded.

To define a new model, provide `Ndim`, `Nfun`, `Rdis`, and `EvlLSF(X)`, with physical-space inputs `[Ndim,Nsamples]` and output `[Nfun,Nsamples]`. Failure is strictly `g < 0`.

## Reference
Xia, W., and Liao, Z. (2025). *Enhanced generalized subset simulation with multiple importance sampling for reliability estimation*. Computers and Structures, 313, 107741. https://doi.org/10.1016/j.compstruc.2025.107741

The benchmark definitions trace to Li, H.-S., Ma, Y.-Z., and Cao, Z. (2015), Computers and Structures, 153, 239–251, https://doi.org/10.1016/j.compstruc.2014.10.014. Balance-heuristic weights follow Veach and Guibas (1995), https://doi.org/10.1145/218380.218498.

## Shared utilities and licensing
`DspRst.m`, `VabSrk.m`, `MisInt.m`, `logmean.m`, and `logsum.m` are identical to their counterparts in [Sequential multiple importance sampling for Reliability Estimation](https://github.com/ZihanLiao-ZJU/Sequential-multiple-importance-sampling-for-Reliability-Estimation), commit `dc227892d7c202ac7c3974b9e83135d69e102cc5`. Equation numbers in shared-file comments refer to that repository's paper. The remaining files implement the GSS-MIS model and interfaces.

See `LICENSE` and retained upstream notices in `Support-Fun/LICENSE-SeMIS` and `Support-Fun/LICENSE-Mo-Chen`.
