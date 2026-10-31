% SUS_MODEL Scalar Subset Simulation threshold model.
% Strict threshold support; thresholds use the midpoint around round(N*p).
% Au, S.K.; Beck, J.L. (2001), Probabilistic Engineering Mechanics 16,
% 263-277. DOI: 10.1016/S0266-8920(01)00019-4.
classdef SuS_model < MIS_Rel
    properties
        G = inf;
        H = nan;
    end

    properties(Dependent)
        FlgCvg
    end

    methods
        function obj = SuS_model(Tfun)
            obj@MIS_Rel(Tfun);
        end
        function Li = EvlLKF(obj)
            y = obj.Y;
            f = y(2,:);
            g = obj.G;
            if isinf(g)
                Li = zeros(size(f));
            else
                Li = log(double(f<g));
            end
        end

        function Pi = EvlPDF(obj,varargin)
            if nargin>1
                idim = varargin{1};
                x = obj.X;
                Pi = obj.Tfun.Rdis.EvlPDF(x,idim);
            else
                y = obj.Y;
                if ~isempty(y)
                    q = y(1,:);
                    Pi = q;
                else
                    x = obj.X;
                    Pi = obj.Tfun.Rdis.EvlPDF(x);
                end
            end
        end

        function [g,h,obj] = UpdObj(obj,y,p)

            f = y(4,:);
            Nsam = size(f,2);
            l = sort(f,"ascend");
            Nsed = round(Nsam*p);
            g = (l(Nsed)+l(Nsed+1))/2;

            h = obj.H;

            obj.G = g;
            obj.H = h;
        end

        function [x_sed,y_sed] = SltSed(obj,x,y,Nsam)

            f = y(4,:);
            g = obj.G;
            Ns = round(Nsam./(1:Nsam));
            Nsed = round(Nsam./Ns);
            ind_sed = f<g;
            x_sed_can = x(:,ind_sed);
            y_sed_can = y(:,ind_sed);
            Nsed_can = sum(ind_sed,"all");
            Nsed = Nsed(find(Nsed_can-Nsed>=0,1,"last"));
            if Nsed_can>=Nsed
                randind = randperm(Nsed_can,Nsed);
            elseif Nsed_can>0
                randind = randi(Nsed_can,1,Nsed);
            else
                randind = [];
            end
            x_sed = x_sed_can(:,randind);
            y_sed = y_sed_can(:,randind);
            if ~isempty(y_sed)
                y_sed = obj.UpdY(y_sed);
            end
        end

        function FlgCvg = get.FlgCvg(obj)
            FlgCvg = obj.G<=0;
        end

        function x = SamGen(obj,Nsze)
            x = obj.Tfun.Rdis.SamGen(Nsze);
        end
    end
end
