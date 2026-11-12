%% Problem and sampling settings
example = 1;
N = 240;                       % Constant samples per level
Nrun = 100;                   % Independent runs; set to 1 for a single run
SamplerClass = @CW_aMH;      % CW_aMH, aCS, HMC or ES
% Example 2 sample counts: GSS-MIS 410; scalar SuS 180.

%% Run simulations
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'Limit-state-functions'));
addpath(fullfile(root, 'Support-Fun'));
addpath(fullfile(root, 'MCMC-samplers'));
addpath(fullfile(root, 'Operations on Logarithmic values'));
records = cell(1, Nrun);
assert(example <= 2, 'The toolbox provides examples 1 and 2.');
algorithm = 'SuS';
if example == 1, nfun = 2; else, nfun = 4; end
for run = 1:Nrun
    probability = zeros(nfun, 1);
    calls = zeros(nfun, 1);
    completed = false(nfun, 1);
    layer_calls = cell(nfun, 1);
    level_samples = cell(nfun, 1);
    for j = 1:nfun
        PF = ScalarExample(example, j);
        model = SuS(PF, 'physical', 'midpoint', false);
        simulation = SeMCS(SamplerClass(model));
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
    records{run} = struct('probability', probability, 'joint_calls', sum(calls), ...
        'layer_joint_calls', {layer_calls}, 'actual_level_samples', {level_samples}, ...
        'completed', all(completed));
end
probabilities = cat(2, records{:});
probabilities = cat(2, probabilities.probability);
calls = cellfun(@(record) record.joint_calls, records);
completed = cellfun(@(record) record.completed, records);
sd = std(probabilities, 0, 2);
if Nrun < 2, sd(:) = NaN; end
results = struct('example', example, 'algorithm', algorithm, 'N', N, ...
    'repeats', Nrun, 'mean', mean(probabilities, 2), ...
    'standard_deviation', sd, ...
    'cov', sd ./ mean(probabilities, 2), ...
    'mean_joint_calls', mean(calls), 'unfinished', sum(~completed), ...
    'records', {records});
fprintf('%s example %d: %d runs, mean joint calls %.2f\n', ...
    algorithm, example, Nrun, results.mean_joint_calls);
