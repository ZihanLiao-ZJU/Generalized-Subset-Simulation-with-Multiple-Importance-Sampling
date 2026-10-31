function results = SuS(example, repetitions, seed)
% SUS Scalar Subset Simulation comparison, one run per performance function.
% results = SuS(example,repetitions,seed); defaults: (1,100,2025).
% Example 1 uses 240 samples per level; Example 2 uses 180.
% Au, S.K.; Beck, J.L. (2001), Probabilistic Engineering Mechanics 16,
% 263-277. DOI: 10.1016/S0266-8920(01)00019-4.
if nargin<1,example=1;end
if nargin<2,repetitions=100;end
if nargin<3,seed=2025;end
validateattributes(example,{'numeric'},{'scalar','integer','>=',1,'<=',2});
validateattributes(repetitions,{'numeric'},{'scalar','integer','positive'});
validateattributes(seed,{'numeric'},{'scalar','integer','nonnegative'});
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'Support-Fun'));
addpath(fullfile(root,'Support-Fun','Operations on Logarithmic values'));
if example==1,N=240;nfun=2;else,N=180;nfun=4;end
pf=zeros(nfun,repetitions);ncal=pf;converged=false(nfun,repetitions);termination_flag=converged;
for run=1:repetitions
 for j=1:nfun
  rng(seed+run-1+10000*j,'twister');
  tf=ScalarExample(example,j);m=SuS_model(tf);s=GSS_MIS(CW_aMH(m));s.NumSam=N;
  evalc('out=s.RunIte;');ncal(j,run)=sum(out.Ncal);
  assert(ncal(j,run)==tf.count,'Joint counter does not equal scalar evaluations.');
  z=1;
  for i=1:out.Nite-1,z=z*mean(out.y{i}(4,:)<out.g{i+1});end
  pf(j,run)=z*mean(out.y{end}(4,:)<0);
  converged(j,run)=out.intBay.FlgCvg;termination_flag(j,run)=out.flg_cvg;
 end
end
results=struct('example',example,'samples_per_level',N,'seed',seed,'Pf',pf, ...
 'Ncal',ncal,'converged',converged,'termination_flag',termination_flag, ...
 'mean_Pf',mean(pf,2),'std_Pf',std(pf,0,2),'CoV',std(pf,0,2)./mean(pf,2), ...
 'mean_total_evaluations',mean(sum(ncal,1)));
fprintf('SuS example %d: %d runs; mean summed scalar evaluations %.2f\n',example,repetitions,results.mean_total_evaluations);
if any(~converged(:)),warning('SuS:Unconverged','%d scalar runs did not converge.',sum(~converged(:)));end
end
