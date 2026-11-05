% GSS_MIS Union-indicator model with independent input marginals.
% G stores performance thresholds and active performance functions.
% Failure and intermediate support use strict g < threshold.
% Quantiles are computed separately for each performance function.
% A response is removed from the active union once its threshold reaches zero.
% Constructor: model = GSS_MIS(Tfun).
% Xia, W.; Liao, Z. (2025), Enhanced generalized subset simulation with
% multiple importance sampling for reliability estimation. Computers and
% Structures 313, 107741. DOI: 10.1016/j.compstruc.2025.107741.
classdef GSS_MIS
    properties
        TF
        Ncal = 0;
        Nfun
        Ndim
        X
        Y
        Nsam = 1000;
        G = struct("parameters",inf,"p",0.1,"active",true);

    end

    properties(Dependent)
        FlgCvg
    end

    methods
        function obj = GSS_MIS(Tfun)
            obj.TF = Tfun;
            obj.Nfun = Tfun.Nfun+3;
            obj.Ndim = Tfun.Ndim;
        end
        function obj = EvlY(obj,x)

            f = obj.TF.EvlLSF(obj.TF.Rdis.U2X(x));
            p = -.5*(obj.Ndim*log(2*pi)+sum(x.^2,1));
            obj.X = x;
            obj.Ncal = size(x,2);
            obj.Y = [p;f];
        end

        function y = UpdY(obj,y)

            y = y(3:end,:);
            obj.Y = y;
            pi = obj.EvlPDF;
            Li = obj.EvlLKF;
            y = [Li;pi;y];
        end

        function Li = EvlLKF(obj)
            y = obj.Y;
            f = y(2:end,:);
            g = obj.G.parameters;
            h = obj.G.active;
            ind = f<g;
            ind = any(ind(h,:),1);
            Li = log(double(ind));
        end

        function Pi = EvlPDF(obj,varargin)
            if nargin>1
                idim = varargin{1};
                x = obj.X;
                Pi = -.5*(log(2*pi)+x(idim,:).^2);
            else
                y = obj.Y;
                if ~isempty(y)
                    q = y(1,:);
                    Pi = q;
                else
                    x = obj.X;
                    Pi = -.5*(obj.Ndim*log(2*pi)+sum(x.^2,1));
                end
            end
        end

        function [y,g,obj,N_cal] = UpdObj(obj,~,y,p)

            f = y(4:end,:);
            g = max(quantile(f,p,2),0);

            h = obj.G.active & g>0;

            obj.G.parameters = g; obj.G.active = h; g = obj.G;
            N_cal = 0;
        end

        function [x_sed,y_sed] = SltSed(obj,x,y,Nsam)

            f = y(4:end,:);
            g = obj.G.parameters;
            h = obj.G.active;
            Ns = round(Nsam./(1:Nsam));
            Nsed = round(Nsam./Ns);
            ind_sed = f<g;
            ind_sed = any(ind_sed(h,:),1);
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
            FlgCvg = ~any(obj.G.active);
        end

        function x = SamGen(obj,Nsze)
            x = randn(obj.Ndim,Nsze);
        end
    end
    methods
        function defaults=run_defaults(obj)
            defaults=struct('NumSam',1000,'ConPrb',0.1,'MaxIte',10000,'Verbose',true);
            if isprop(obj,'Nsam'),defaults.NumSam=obj.Nsam;end
            if isfield(obj.G,'p'),defaults.ConPrb=obj.G.p;end
        end
        function c=start_run(~,owner,sampler)
            active=owner.IntBay.G.active;
            c=SeMCS.start_context(owner,sampler,owner.MaxIte-1,'replace');
            if ~isfield(c.model.G,'active'),c.model.G.active=active;end
            c.g{1}=c.model.G;c.sampler.BayStc=c.model;
        end
        function [c,r]=sample_level(~,c,level)
            [c,r]=SeMCS.sample_context(c,level);
        end
        function c=advance_level(~,c,r,level)
            if level+1<=numel(c.owner.G)
                active=c.model.G.active;
                c.model.G=SeMCS.prescribed_state(c.model.G,c.owner.G{level+1},false,'replace');
                if ~isfield(c.model.G,'active'),c.model.G.active=active;end
                extra=0;
            else
                [c.y{level},~,c.model,extra]=c.model.UpdObj(r.x,r.y,c.owner.ConPrb);
            end
            c.Ncal(level)=c.Ncal(level)+extra;c.g{level+1,1}=c.model.G;
            c.sampler.BayStc=c.model;
        end
        function c=select_run_seeds(~,c,r,level)
            if c.targets(level+1)<=0,error('SeMCS:ScheduleExhausted','Sample budget exhausted before seed selection.');end
            target=c.targets(level+1);[c.xs,c.ys]=c.model.SltSed(r.x,c.y{level},target);
            c.nc=size(c.xs,2);c.steps=round(target/c.nc);
            if c.owner.Verbose
                fprintf('Level %d: N=%d, chains=%d, steps=%d, calls=%d\n',level,r.count,r.shape(1),r.shape(2),c.Ncal(level));
            end
        end
        function [yes,c]=run_complete(~,c,~,phase)
            yes=false;if ~strcmp(phase,'after_seeds'),return;end
            if c.nc<=0 || numel(c.x)==c.limit,yes=true;return;end
            c.complete=c.model.FlgCvg;yes=c.complete;
        end
        function out=finish_run(~,c)
            out=SeMCS.pack_output(c);
        end
    end
end
