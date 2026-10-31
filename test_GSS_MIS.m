function test_GSS_MIS
% TEST_GSS_MIS Deterministic benchmark and estimator checks.
root=fileparts(mfilename('fullpath'));addpath(genpath(root));
a=Ex1;b=Ex2;
assert(max(abs(a.EvlLSF([0;0])-[-.01;.31]))<1e-12);
g=b.EvlLSF([0;0]);assert(max(abs(g-[3;7/sqrt(2);3;7/sqrt(2)]))<1e-12);
x=[0,.4,1;0,.4,1];m=GSS_MIS_model(a);m=m.EvlY(x);
y=[m.EvlLKF;m.EvlPDF;m.Y];
o=struct('Nite',1,'Ncal',3,'y',{{y}},'g',{{inf}},'h',{{true}});
[pf,ncal]=m.EvlPf(o);
assert(max(abs(pf-[1;1]/3))<1e-12);assert(ncal==3);
assert(logsum([0,-Inf],'all')==0);assert(logsum([-Inf,-Inf],'all')==-Inf);
fprintf('Benchmark definitions, one-level MIS and joint counts passed.\n');
end
