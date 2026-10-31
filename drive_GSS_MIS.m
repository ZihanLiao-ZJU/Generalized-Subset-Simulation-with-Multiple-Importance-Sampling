function results = drive_GSS_MIS(example, repetitions, seed)
% DRIVE_GSS_MIS Run the two numerical GSS-MIS examples.
% results = drive_GSS_MIS(example, repetitions, seed)
% Defaults: example=1, repetitions=100, seed=2025.
% Ex1 uses 440 samples per level; Ex2 uses 410. Conditional probability=.1.
% Ncal is the number of joint performance-function vector evaluations.
% Xia, W.; Liao, Z. (2025), Computers and Structures 313, 107741.
% DOI: 10.1016/j.compstruc.2025.107741.
if nargin<1, example=1; end
if nargin<2, repetitions=100; end
if nargin<3, seed=2025; end
validateattributes(example,{'numeric'},{'scalar','integer','>=',1,'<=',2});
validateattributes(repetitions,{'numeric'},{'scalar','integer','positive'});
validateattributes(seed,{'numeric'},{'scalar','integer','nonnegative'});
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'Limit-state-functions'));
addpath(fullfile(root,'Support-Fun'));
addpath(fullfile(root,'Support-Fun','Operations on Logarithmic values'));
if example==1, Tfun=Ex1; N=440; else, Tfun=Ex2; N=410; end
model=GSS_MIS_model(Tfun); sampler=CW_aMH(model);
simulation=GSS_MIS(sampler); simulation.NumSam=N;
pf=zeros(Tfun.Nfun,repetitions); ncal=zeros(1,repetitions);
converged=false(1,repetitions); termination_flag=converged;
for run=1:repetitions
    rng(seed+run-1,'twister');
    evalc('out=simulation.RunIte;');
    [pf(:,run),ncal(run)]=model.EvlPf(out);
    converged(run)=out.intBay.FlgCvg;
    termination_flag(run)=out.flg_cvg;
end
results=struct('example',example,'samples_per_level',N,'seed',seed, ...
    'Pf',pf,'Ncal',ncal,'converged',converged, ...
    'termination_flag',termination_flag,'mean_Pf',mean(pf,2), ...
    'std_Pf',std(pf,0,2),'CoV',std(pf,0,2)./mean(pf,2), ...
    'mean_joint_evaluations',mean(ncal));
fprintf('Example %d: %d runs; mean joint evaluations %.2f\n',example,repetitions,mean(ncal));
fprintf('PF        Mean estimate          Empirical CoV\n');
for j=1:Tfun.Nfun, fprintf('%2d        %.8g             %.6g\n',j,results.mean_Pf(j),results.CoV(j)); end
if any(~converged), warning('GSS_MIS:Unconverged','%d runs did not reach all failure thresholds.',sum(~converged)); end
end
