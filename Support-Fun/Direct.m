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

            obj.Ndim = size(Marg_X,1);
            obj.Marg_X = Marg_X;
        end

        function pdv = EvlPDF(obj,varargin)

            marg = obj.Marg_X;
            if nargin>2
                X = varargin{1};
                idim = varargin{2};
                pdv = log(marg(idim).pdf(X));
            else
                X = varargin{1};
                ndim = obj.Ndim;
                nsam = size(X,2);
                pdv = zeros(ndim,nsam);
                for idim = 1:ndim
                    pdv(idim,:) = log(marg(idim).pdf(X(idim,:)));
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
                X(idim,:) = marg(idim).random(1,Nsam);
            end
        end
    end
end
