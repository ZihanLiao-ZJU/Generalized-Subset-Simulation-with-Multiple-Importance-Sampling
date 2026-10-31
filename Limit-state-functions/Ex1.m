% EX1 Numerical benchmark; paper Eq. (19).
% Li, H.-S.; Ma, Y.-Z.; Cao, Z. (2015), Computers and Structures 153,
% 239--251. DOI: 10.1016/j.compstruc.2014.10.014.
% Inputs: X [2,N], independent standard normal. Outputs: G [2,N].
%
% Reference: Xia, W.; Liao, Z. (2025). Enhanced generalized subset simulation with
% multiple importance sampling for reliability estimation. Computers and
% Structures 313, 107741. DOI: 10.1016/j.compstruc.2025.107741.
% Samples are columns; failure is strictly g < 0.
classdef Ex1
    properties
        Nfun = 2;
        Ndim = 2;
        Rdis = Direct([makedist("Normal","mu",0,"sigma",1);makedist("Normal","mu",0,"sigma",1)]);
        Dist = [makedist("Normal","mu",0,"sigma",1);makedist("Normal","mu",0,"sigma",1)];
    end

    methods
        function L = EvlLSF(~,theta)

            L1 = theta(1,:).^2+theta(2,:).^2-0.01;
            L2 = (theta(1,:)-0.4).^2+(theta(2,:)-0.4).^2-0.01;
            L = [L1;L2];
        end
    end
end
