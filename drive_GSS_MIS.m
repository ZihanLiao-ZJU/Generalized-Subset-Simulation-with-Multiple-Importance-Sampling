%% Problem and sampling settings
PF = Ex1();
% PF = Ex2();
N = 440;                       % Constant samples per level
Nrun = 100;                   % Independent runs; set to 1 for a single run
SamplerClass = @CW_aMH;      % CW_aMH, aCS, HMC or ES
% Example 2 sample counts: GSS-MIS 410; scalar SuS 180.

%% Model and sampler
model = GSS_MIS(PF);
simulation = SeMCS(SamplerClass(model));
simulation.NumSam = N;

%% Run simulations
records = cell(1, Nrun);
algorithm = 'GSS-MIS';
example = sscanf(class(PF), 'Ex%d');
for run = 1:Nrun
    out = simulation.RunIte;
    mis = MIS(out);
    probability = mis.Pf_hat;
    records{run} = struct('probability', probability, 'joint_calls', sum(out.Ncal), ...
        'layer_joint_calls', out.Ncal, 'actual_level_samples', out.Nsam, ...
        'completed', out.intBay.FlgCvg);
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
fprintf('%s %s: %d runs, mean joint calls %.2f\n', ...
    algorithm, class(PF), Nrun, results.mean_joint_calls);
