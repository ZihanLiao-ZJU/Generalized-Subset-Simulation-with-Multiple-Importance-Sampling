function results = drive_SuS(example, N, repeats, seed, sampler_name)
% DRIVE_SUS Run scalar SuS on each performance function of Ex1 or Ex2.
% results = drive_SuS(example, N, repeats, seed)
% Defaults: example 1, N=240 (Ex1) / 180 (Ex2), 100 runs, seed 2026.
% Explicit policy preserves midpoint thresholds and physical coordinates.
if nargin < 1, example = 1; end
if nargin < 5 || isempty(sampler_name), sampler_name = 'CW_aMH'; end
assert(any(strcmp(sampler_name,{'CW_aMH','aCS','HMC','ES'})), 'Unknown sampler');
if nargin < 2 || isempty(N)
    if example == 1, N = 240; else, N = 180; end
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
addpath(fullfile(root, 'Support-Fun', 'Operations on Logarithmic values'));
records = cell(1, repeats);
assert(example <= 2, 'The toolbox provides examples 1 and 2.');
algorithm = 'SuS';
if example == 1, nfun = 2; else, nfun = 4; end
for run = 1:repeats
    probability = zeros(nfun, 1);
    calls = zeros(nfun, 1);
    completed = false(nfun, 1);
    layer_calls = cell(nfun, 1);
    level_samples = cell(nfun, 1);
    for j = 1:nfun
        rng(seed + run - 1 + 10000*j, 'twister');
        PF = ScalarExample(example, j);
        model = SuS(PF, 'physical', 'midpoint', false);
        simulation = SeMCS(feval(sampler_name,model));
        simulation.NumSam = N;
        simulation.MaxIte = 1000;
        simulation.Verbose = false;
        out = simulation.RunIte;
        calls(j) = sum(out.Ncal);
        assert(calls(j) == PF.count, 'Scalar call count mismatch.');
        z = 1;
        for level = 1:out.Nite-1
            z = z*mean(out.y{level}(4,:) < out.g{level+1}.parameters);
        end
        probability(j) = z*mean(out.y{end}(4,:) < 0);
        completed(j) = out.intBay.FlgCvg;
        layer_calls{j} = out.Ncal;
        level_samples{j} = out.Nsam;
    end
    records{run} = struct('seed', seed + run - 1, ...
        'probability', probability, 'joint_calls', sum(calls), ...
        'layer_joint_calls', {layer_calls}, 'actual_level_samples', {level_samples}, ...
        'completed', all(completed));
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
