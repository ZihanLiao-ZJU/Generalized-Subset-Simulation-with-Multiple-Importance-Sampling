function results = drive_GSS_MIS(example, N, repeats, seed, sampler_name)
% DRIVE_GSS_MIS Run GSS-MIS: distribution -> sampler -> SeMCS -> MIS.
% results = drive_GSS_MIS(example, N, repeats, seed)
% Defaults: example 1, N=440 (Ex1) / 410 (Ex2), 100 runs, seed 2026.
% Ncal counts joint performance-function vector evaluations.
if nargin < 1, example = 1; end
if nargin < 5 || isempty(sampler_name), sampler_name = 'CW_aMH'; end
assert(any(strcmp(sampler_name,{'CW_aMH','aCS','HMC','ES'})), 'Unknown sampler');
if nargin < 2 || isempty(N)
    if example == 1, N = 440; else, N = 410; end
end
if nargin < 3, repeats = 100; end
if nargin < 4, seed = 2026; end
validateattributes(example, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(N, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(repeats, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(seed, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'Limit-state-functions'));
addpath(fullfile(root, 'Support-Fun'));
addpath(fullfile(root, 'MCMC-samplers'));
addpath(fullfile(root, 'Operations on Logarithmic values'));
records = cell(1, repeats);
assert(example <= 2, 'The toolbox provides examples 1 and 2.');
algorithm = 'GSS-MIS';
for run = 1:repeats
    rng(seed + run - 1, 'twister');
    if example == 1, PF = Ex1; else, PF = Ex2; end
    model = GSS_MIS(PF);
    sampler = feval(sampler_name,model);
    simulation = SeMCS(sampler);
    simulation.NumSam = N;
    simulation.MaxIte = 1000;
    simulation.Verbose = false;
    out = simulation.RunIte;
    mis = MIS(out);
    probability = mis.Pf_hat;
    records{run} = struct('seed', seed + run - 1, ...
        'probability', probability, 'joint_calls', sum(out.Ncal), ...
        'layer_joint_calls', out.Ncal, 'actual_level_samples', out.Nsam, ...
        'completed', out.intBay.FlgCvg);
end
probabilities = cat(2, records{:});
probabilities = cat(2, probabilities.probability);
calls = cellfun(@(record) record.joint_calls, records);
completed = cellfun(@(record) record.completed, records);
sd = std(probabilities, 0, 2);
if repeats < 2, sd(:) = NaN; end
results = struct('example', example, 'algorithm', algorithm, 'N', N, ...
    'repeats', repeats, 'seed', seed, 'mean', mean(probabilities, 2), ...
    'standard_deviation', sd, ...
    'cov', sd ./ mean(probabilities, 2), ...
    'mean_joint_calls', mean(calls), 'unfinished', sum(~completed), ...
    'records', {records});
fprintf('%s example %d: %d runs, mean joint calls %.2f\n', ...
    algorithm, example, repeats, results.mean_joint_calls);
end
