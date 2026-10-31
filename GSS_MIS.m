% GSS_MIS Generalized subset simulation sample-generation framework.
% Authors: Zihan Liao, Xiao He and Weili Xia (original framework).
% Samples are columns; each level is [Ndim,Nchains,Nsteps].
% RunIte returns x/y/g/h and sample, chain and joint-evaluation counts.
% Packed y rows: log intermediate likelihood, log prior, log reference
% density, then the performance-function values.
% Multi-response union thresholds and active sets are provided by the model.
% Xia, W.; Liao, Z. (2025), Enhanced generalized subset simulation with
% multiple importance sampling for reliability estimation. Computers and
% Structures 313, 107741. DOI: 10.1016/j.compstruc.2025.107741.
classdef GSS_MIS

    properties(Constant)
        MaxIte = 1000;
    end

    properties
        IntBay
        Sampler
        ConPrb = 0.1;
        NumSam = 1000;
        H = [];
        G = [];
        X = [];
        Y = [];
    end

    methods
        function obj = GSS_MIS(Sampler)

            if nargin<1
                error('1.IntBay 2.Sampler are needed in the inputs of SuS()!!')
            else
                obj.Sampler = Sampler;
                obj.IntBay = Sampler.BayStc;
            end
        end

        function out = RunIte(obj)


            intBay = obj.IntBay;
            sampler = obj.Sampler;
            p = obj.ConPrb;
            maxi = obj.MaxIte;
            flg_Nsam = numel(obj.NumSam)>1;
            flg_IntD = ~isempty(obj.H) && length(obj.H)==length(obj.G);

            Nsam = zeros(maxi,1);
            Nsze = zeros(maxi,2);
            Ncal = zeros(maxi,1);
            x = cell(maxi,1);
            y = cell(maxi,1);
            g = cell(maxi,1);
            h = cell(maxi,1);

            if flg_Nsam
                Nite_sam = length(obj.NumSam);
                Nsam(1:Nite_sam) = obj.NumSam;
            else
                Nsam(:) = obj.NumSam;
            end
            Nsze(1,:) = [Nsam(1),1];

            if flg_IntD
                Nite_int = length(obj.G);
                g(1:Nite_int) = obj.G;
                h(1:Nite_int) = obj.H;
            else
                Nite_int = 0;
                g{1} = intBay.G;
                h{1} = intBay.H;
            end

            for ite = 1:maxi
                if all(Nsze(ite,:)>0) && ~isnan(Nsze(ite,2))
                    if ite==1
                        [x{ite},y{ite},Ncal(ite)] = sampler.SamGen(Nsze(ite,1),Nsze(ite,2));
                    else
                        [x{ite},y{ite},Ncal(ite)] = sampler.SamGen(x_sed,y_sed,Nsze(ite,2));
                    end
                end
                Nsam(ite) = prod(Nsze(ite,:));

                if ite+1<=Nite_int
                    intBay.G = g{ite+1};
                    intBay.H = h{ite+1};
                else
                    [g{ite+1},h{ite+1},intBay] = intBay.UpdObj(y{ite},p);
                end
                sampler.BayStc = intBay;

                [x_sed,y_sed] = intBay.SltSed(x{ite},y{ite},Nsam(ite+1));
                Nsze(ite+1,1) = size(x_sed,2);
                Nsze(ite+1,2) = round(Nsam(ite+1)/Nsze(ite+1,1));

                DspRst(ite,Nsam(ite),Nsze(ite,:),Ncal(ite))


                if Nsze(ite+1,1)<=0 || ite+1==maxi
                    flg_cvg = false;
                    break
                elseif intBay.FlgCvg
                    flg_cvg = true;
                    break
                end
            end
            Nite = ite;
            [x,y,Nsam,Nsze,Ncal] = VabSrk(Nite,x,y,Nsam,Nsze,Ncal);
            [g,h] = VabSrk(Nite+1,g,h);
            out.x = x;
            out.y = y;
            out.g = g;
            out.h = h;
            out.Nsam = Nsam;
            out.Nsze = Nsze;
            out.Ncal = Ncal;
            out.Nite = Nite;
            out.flg_cvg = flg_cvg;
            out.intBay = intBay;
        end
    end
end
