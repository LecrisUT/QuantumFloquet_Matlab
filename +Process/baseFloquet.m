classdef baseFloquet < handle
    properties
        calcObj
    end
    methods
        function obj=baseFloquet(calcObj)
            obj.calcObj = calcObj;
        end
        function Res=Variational(obj,Input,Args)
            arguments
                obj
                Input
                Args.Ex_Psi
                Args.Norm   logical =false
            end
            Res = Input;
            NSol = length(Input(:));
            for iN = 1:NSol
                for iStep = 1:length(Input(iN).steps)
                    if isempty(Input(iN).steps(iStep).Psi)
                        if isfield(Args,'Ex_Psi')
                            Res(iN).steps(iStep).S = nan(size(Args.Ex_Psi,2),1);
                            Res(iN).steps(iStep).S2 = nan(size(Args.Ex_Psi,2),1);
                        end
                        if Args.Norm
                            Res(iN).steps(iStep).norm = nan;
                            Res(iN).steps(iStep).norm0 = nan;
                        end
                    else
                        if isfield(Args,'Ex_Psi')
                            tS = nan(size(Args.Ex_Psi,2),1);
                            tS2 = nan(size(Args.Ex_Psi,2),1);
                            tPsi1 = Input(iN).steps(iStep).Psi;
                            tPsi2 = tPsi1 ./ norm(tPsi1);
                            for ik = 1:obj.calcObj.k_max2
                                tPsi = circshift(tPsi1,obj.calcObj.N*ik,1);
                                tS = max(tS,abs(Args.Ex_Psi' * tPsi));
                                tPsi = circshift(tPsi2,obj.calcObj.N*ik,1);
                                tS2 = max(tS2,abs(Args.Ex_Psi' * tPsi));
                            end
                            Res(iN).steps(iStep).S = tS;
                            Res(iN).steps(iStep).S2 = tS2;
                        end
                        if Args.Norm
                            tPsi = Input(iN).steps(iStep).Psi;
                            tPsi0 = obj.calcObj.Psi0(tPsi);
                            Res(iN).steps(iStep).norm = tPsi' * tPsi;
                            Res(iN).steps(iStep).norm0 = tPsi0' * tPsi0;
                        end
                    end
                end
            end
        end
    end
end