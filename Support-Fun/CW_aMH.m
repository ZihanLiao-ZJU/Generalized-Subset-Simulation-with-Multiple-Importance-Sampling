% CW_AMH Component-wise Metropolis-Hastings sampling.
% BayStc supplies marginal log densities and the union-indicator likelihood.
% Proposals use uniform increments scaled by the seed sample standard
% deviation (held fixed within each level). Unchanged states reuse outputs.
% SamGen(Nchains,Nsteps) initializes from the reference density;
% SamGen(Xseed,Yseed,Nsteps) advances seeded chains.
% Outputs: X [Ndim,Nchains,Nsteps], packed Y, and joint call count.
% Au, S.K.; Beck, J.L. (2001), Estimation of small failure probabilities
% in high dimensions by subset simulation. Probabilistic Engineering
% Mechanics 16, 263-277. DOI: 10.1016/S0266-8920(01)00019-4.
% Xia, W.; Liao, Z. (2025), Enhanced generalized subset simulation with
% multiple importance sampling for reliability estimation. Computers and
% Structures 313, 107741. DOI: 10.1016/j.compstruc.2025.107741.
classdef CW_aMH
    properties(Constant)
        SedSav = false;
    end

    properties
        BayStc
        Ndim
        Nfun
    end

    methods
        function obj = CW_aMH(BayStc)
            obj.BayStc = BayStc;
            obj.Nfun = BayStc.Nfun;
            obj.Ndim = BayStc.Ndim;
        end

        function [x,y,Ncal] = SamGen(obj,varargin)

            baystc = obj.BayStc;
            sedsav = obj.SedSav || nargin<=3;
            NDim = obj.Ndim;
            NFun = obj.Nfun;
            if nargin<=3
                Nc = varargin{1};
                Ns = varargin{2};
                x_sed = baystc.SamGen(Nc);
                baystc = baystc.EvlY(x_sed);
                y_sed = [baystc.EvlLKF;baystc.EvlPDF;baystc.Y];
                Ncal = Nc;
            else
                x_sed = varargin{1};
                y_sed = varargin{2};
                Ns = varargin{3};
                Nc = size(x_sed,2);
                Ncal = 0;
            end
            x = zeros(NDim,Nc,Ns);
            y = zeros(NFun,Nc,Ns);

            du = 2*rand(NDim,Nc,Ns)-1;
            p = rand(NDim,Nc,Ns);

            sigma0 = max(std(x_sed,0,2),eps);
            for k=1:Ns
                if k == 1 && sedsav
                    x(:,:,1) = x_sed;
                    y(:,:,1) = y_sed;
                else
                    for idim=1:NDim
                        x(idim,:,k) = x_sed(idim,:) + sigma0(idim)*du(idim,:,k);
                        baystc.X = x_sed(idim,:);
                        pdf_sed = baystc.EvlPDF(idim);
                        baystc.X = x(idim,:,k);
                        pdf_cad = baystc.EvlPDF(idim);
                        ind_rej =  exp(pdf_cad-pdf_sed) < p(idim,:,k);
                        x(idim,ind_rej,k) = x_sed(idim,ind_rej);
                    end
                    ind_rep = all(x(:,:,k) == x_sed,1);
                    ind_evl = ~ind_rep;
                    y(:,ind_rep,k) = y_sed(:,ind_rep);
                    dcal = sum(ind_evl);
                    if dcal>0
                        ind_sam = 1:Nc;
                        ind_evl = ind_sam(ind_evl);
                        baystc = baystc.EvlY(x(:,ind_evl,k));
                        y(:,ind_evl,k) = [baystc.EvlLKF;baystc.EvlPDF;baystc.Y];
                        ind_rej =  ind_evl(exp(y(1,ind_evl,k)-y_sed(1,ind_evl)) < rand(1,dcal));
                        x(:,ind_rej,k) = x_sed(:,ind_rej);
                        y(:,ind_rej,k) = y_sed(:,ind_rej);
                        Ncal = Ncal + dcal;
                    end
                    x_sed = x(:,:,k);
                    y_sed = y(:,:,k);
                end
            end
        end

    end
end
