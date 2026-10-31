% EX2 Numerical benchmark; paper Eq. (20).
% Li, H.-S.; Ma, Y.-Z.; Cao, Z. (2015), Computers and Structures 153,
% 239--251. DOI: 10.1016/j.compstruc.2014.10.014.
% Inputs: X [2,N], independent standard normal. Outputs: G [4,N].
%
% Reference: Xia, W.; Liao, Z. (2025). Enhanced generalized subset simulation with
% multiple importance sampling for reliability estimation. Computers and
% Structures 313, 107741. DOI: 10.1016/j.compstruc.2025.107741.
% Samples are columns; failure is strictly g < 0.
classdef Ex2
    properties
        Nfun = 4;
        Ndim = 2;
        Rdis = Direct([makedist("Normal","mu",0,"sigma",1);makedist("Normal","mu",0,"sigma",1)]);
        Dist = [makedist("Normal","mu",0,"sigma",1);makedist("Normal","mu",0,"sigma",1)];
    end

    methods
        function L = EvlLSF(~,theta)

            L1 = 3+(theta(1,:)-theta(2,:)).^2/10-(theta(1,:)+theta(2,:))/sqrt(2);
            L2 = theta(1,:)-theta(2,:)+7/sqrt(2);
            L3 = 3+(theta(1,:)-theta(2,:)).^2/10+(theta(1,:)+theta(2,:))/sqrt(2);
            L4 = theta(2,:)-theta(1,:)+7/sqrt(2);
            L = [L1;L2;L3;L4];
        end
    end
end
