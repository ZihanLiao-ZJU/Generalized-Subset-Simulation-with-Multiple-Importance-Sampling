% MIS_REL Reliability evaluations and balance-heuristic MIS postprocessing.
% Tfun implements Ndim, Nfun, Rdis and EvlLSF(X) in physical space.
% EvlY stores [log reference density; performance functions].
% EvlPf combines all sampled levels using their stored thresholds and active
% sets. Ncal counts joint calls evaluating the performance-function vector.
% Reference for balance weights: Veach, E.; Guibas, L.J. (1995),
% Optimally combining sampling techniques for Monte Carlo rendering,
% SIGGRAPH 95, 419-428. DOI: 10.1145/218380.218498.
% Xia, W.; Liao, Z. (2025), Enhanced generalized subset simulation with
% multiple importance sampling for reliability estimation. Computers and
% Structures 313, 107741. DOI: 10.1016/j.compstruc.2025.107741.
classdef MIS_Rel
    properties
        Tfun
        Nfun
        Ndim
        X
        Y
    end

    methods
        function obj = MIS_Rel(Tfun)

            obj.Tfun = Tfun;
            obj.Nfun = Tfun.Nfun+3;
            obj.Ndim = Tfun.Ndim;
        end

        function obj = EvlY(obj,x)

            f = obj.Tfun.EvlLSF(x);
            p = obj.Tfun.Rdis.EvlPDF(x);
            obj.X = x;
            obj.Y = [p;f];
        end

        function y = UpdY(obj,y)

            y = y(3:end,:);
            obj.Y = y;
            pi = obj.EvlPDF;
            Li = obj.EvlLKF;
            y = [Li;pi;y];
        end

        function [Pf,Ncal] = EvlPf(obj,in)
            Nite = in.Nite;
            Ncal = sum(in.Ncal,"all");
            Nlsf = obj.Tfun.Nfun;
            y = in.y;
            g = in.g;
            h = in.h;

            c = zeros(Nite,1);
            for ite = 1:Nite-1
                obj.G = g{ite+1};
                obj.H = h{ite+1};
                yi = y{ite};
                obj.Y = yi(3:end,:);
                pi = obj.EvlPDF;
                Li = obj.EvlLKF;
                w_res = Li+pi-yi(1,:)-yi(2,:);
                w_max = max(w_res,[],"all");
                w_res = w_res-w_max;
                w_res(isnan(w_res)) = -inf;
                c(ite+1) = c(ite) + logmean(w_res,"all");
                y{ite+1}(1,:) = y{ite+1}(1,:)-c(ite+1);
            end

            w_mis = WgtMis_Blc(obj,g,h,c,y);

            Pf = zeros(Nlsf,1);
            for ilsf = 1:Nlsf
                w_is = cell(Nite,1);
                for ite = 1:Nite
                    Nsze = size(y{ite});
                    Nsze = [Nsze(2:end),1];
                    I = y{ite}(3+ilsf,:)>=0;
                    w_is_tep = y{ite}(3,:) - y{ite}(1,:) - y{ite}(2,:);
                    w_is_tep(I) = -inf;
                    w_is{ite} = zeros(Nsze);
                    w_is{ite}(:) = w_is_tep(:);
                    w_is{ite}(isnan(w_is{ite})) = -inf;
                end
                Pf(ilsf) = exp(logsum(MisInt(w_mis,w_is),"all"));
            end
        end
    end
end
