% DIRECT Independent marginal probability distributions.
% Marg_X is a column array of MATLAB probability distribution objects.
% EvlPDF(X) returns the joint log density; EvlPDF(X,idim) the marginal
% log density. SamGen(N) returns independent samples in columns.
% Xia, W.; Liao, Z. (2025), Enhanced generalized subset simulation with
% multiple importance sampling for reliability estimation. Computers and
% Structures 313, 107741. DOI: 10.1016/j.compstruc.2025.107741.
classdef Direct
    
    properties
        Ndim
        Marg_X
    end

    methods
        function obj = Direct(Marg_X)

            obj.Ndim = numel(Marg_X);
            obj.Marg_X = Marg_X;
        end

        function X = U2X(obj,U)
            % Independent Nataf map: Gaussian coordinates to physical marginals.
            assert(size(U,1)==obj.Ndim,'Direct:Dimension','Input dimension mismatch');
            X=zeros(size(U));
            for j=1:obj.Ndim
                if iscell(obj.Marg_X),marginal=obj.Marg_X{j};else,marginal=obj.Marg_X(j);end
                if strcmpi(marginal.DistributionName,'Normal')
                    X(j,:)=marginal.mu+marginal.sigma*U(j,:);
                elseif strcmpi(marginal.DistributionName,'Lognormal')
                    X(j,:)=exp(marginal.mu+marginal.sigma*U(j,:));
                else
                    X(j,:)=marginal.icdf(.5*erfc(-U(j,:)/sqrt(2)));
                end
            end
        end
        function pdv = EvlPDF(obj,varargin)

            marg = obj.Marg_X;
            if nargin>2
                X = varargin{1};
                idim = varargin{2};
                if iscell(marg),m=marg{idim};else,m=marg(idim);end
                pdv = log(m.pdf(X));
            else
                X = varargin{1};
                ndim = obj.Ndim;
                nsam = size(X,2);
                pdv = zeros(ndim,nsam);
                for idim = 1:ndim
                    if iscell(marg),m=marg{idim};else,m=marg(idim);end
                    pdv(idim,:) = log(m.pdf(X(idim,:)));
                end
                pdv(isnan(pdv)) = 0;
                pdv = sum(pdv,1);
            end
        end

        function X = SamGen(obj,varargin)
            ndim = obj.Ndim;
            marg = obj.Marg_X;
            if varargin{1}(1,1) == floor(varargin{1}(1,1))
                Nsam = varargin{1};
            else
                x_sed = varargin{1};
                Nsam = size(x_sed,2);
            end
            X = zeros(ndim,Nsam);
            for idim = 1:ndim
                if iscell(marg),m=marg{idim};else,m=marg(idim);end
                X(idim,:) = m.random(1,Nsam);
            end
        end
    end
end
