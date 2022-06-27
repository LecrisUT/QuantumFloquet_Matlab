classdef baseFloquet < handle
    properties
        calcObj
        exactObj
    end
    methods
        function obj=baseFloquet(calcObj,Args)
            arguments
                calcObj
                Args.exactObj   = calcObj
            end
            obj.calcObj = calcObj;
            obj.exactObj = Args.exactObj;
        end
        function Res=Variational(obj,Input,Args)
            arguments
                obj
                Input
                Args.Ex_Psi
                Args.varEps    logical =false
                Args.Ex_varEps    logical =false
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
                        if Args.varEps
                            Res(iN).steps(iStep).varEps = nan;
                        end
                        if Args.Ex_varEps
                            Res(iN).steps(iStep).Ex_varEps = nan;
                        end
                    else
                        if isfield(Args,'Ex_Psi')
                            tS = nan(size(Args.Ex_Psi,2),1);
                            tS2 = nan(size(Args.Ex_Psi,2),1);
                            tPsi1 = Input(iN).steps(iStep).Psi;
                            tPsi2 = tPsi1 ./ norm(tPsi1);
                            % Reshape to exact object's size
                            tPsi1 = obj.exactObj.matchSize(obj.calcObj,tPsi1);
                            tPsi2 = obj.exactObj.matchSize(obj.calcObj,tPsi2);
                            for ik = obj.calcObj.k_range
                                tPsi = circshift(tPsi1,obj.exactObj.N*ik,1);
                                tS = max(tS,abs(Args.Ex_Psi' * tPsi));
                                tPsi = circshift(tPsi2,obj.exactObj.N*ik,1);
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
                        if Args.varEps
                            tPsi = Input(iN).steps(iStep).Psi;
                            Res(iN).steps(iStep).varEps = obj.calcObj.varEps(tPsi,normalize=true);
                        end
                        if Args.Ex_varEps
                            tPsi = Input(iN).steps(iStep).Psi;
                            tPsi = obj.exactObj.matchSize(obj.calcObj,tPsi);
                            Res(iN).steps(iStep).Ex_varEps = obj.exactObj.varEps(tPsi,normalize=true);
                        end
                    end
                end
            end
        end
    end
end