% WGTMIS_BLC Balance-heuristic log weights for all sampled proposal levels.
% Inputs: model, stored thresholds G, active masks H, log normalizers C,
% and packed normalized sample outputs Y. Supports unequal level sizes.
% Returns one [Nchains,Nsteps] log-weight array for each level.
% Veach, E.; Guibas, L.J. (1995), SIGGRAPH 95, 419-428.
% DOI: 10.1145/218380.218498.
% Xia, W.; Liao, Z. (2025), Enhanced generalized subset simulation with
% multiple importance sampling for reliability estimation. Computers and
% Structures 313, 107741. DOI: 10.1016/j.compstruc.2025.107741.
function w_mis = WgtMis_Blc(intBay,g,h,c,y)

beta = 1;
Nite = length(y);
w_mis = cell(Nite,1);
Nsam = zeros(Nite,1);
Nsze = zeros(Nite,2);
for ite = 1:Nite
    Nsze(ite,:) = size(y{ite},2:3);
    Nsam(ite) = prod(Nsze(ite,:));
    w_mis{ite} = zeros(Nsze(ite,:));
end

for ite = 1:Nite
    pdf_crs = zeros(Nite,Nsam(ite));
    for ite1 = 1:Nite
        intBay.G = g{ite1};
        intBay.H = h{ite1};
        intBay.Y = y{ite}(3:end,:);
        pdf_crs(ite1,:) = intBay.EvlLKF+intBay.EvlPDF-c(ite1);
    end
    w_mis_tep = beta*log(Nsam(ite)) + beta*pdf_crs(ite,:) - logsum(beta*pdf_crs+beta*log(Nsam),1);
    w_mis{ite}(:) = w_mis_tep(:);
end
end
