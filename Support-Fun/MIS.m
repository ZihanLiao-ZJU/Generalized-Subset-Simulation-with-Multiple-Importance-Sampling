classdef MIS
    % Balance-heuristic MIS postprocessing of a SeMCS output.
    properties
        Pf_hat
        Ncal
        log_z
    end
    methods
        function obj = MIS(in)
            Nite = in.Nite;
            Ncal = sum(in.Ncal,"all");
            model = in.intBay;
            Nlsf = model.TF.Nfun;
            y = in.y;
            g = in.g;
            h = in.h;

            c = zeros(Nite,1);
            for ite = 1:Nite-1
                model.G = g{ite+1};
                model.H = h{ite+1};
                yi = y{ite};
                model.Y = yi(3:end,:);
                pi = model.EvlPDF;
                Li = model.EvlLKF;
                w_res = Li+pi-yi(1,:)-yi(2,:);
                w_max = max(w_res,[],"all");
                w_res = w_res-w_max;
                w_res(isnan(w_res)) = -inf;
                c(ite+1) = c(ite) + logmean(w_res,"all");
                y{ite+1}(1,:) = y{ite+1}(1,:)-c(ite+1);
            end

            w_mis = WgtMis_Blc(model,g,h,c,y);

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
            obj.Pf_hat = Pf;
            obj.Ncal = Ncal;
            obj.log_z = c;
        end
    end
end
