classdef TwoLevel < Process.baseFloquet
    methods
        function obj=TwoLevel(calcObj)
            obj@Process.baseFloquet(calcObj);
        end
        function Res=Variational(obj,Input,Args,Args2)
            arguments
                obj
                Input
                Args.Ex_Psi
                Args.Norm   logical =false
                Args2.th    logical =false
            end
            Args = namedargs2cell(Args);
            Res = Variational@Process.baseFloquet(obj,Input,Args{:});
            NSol = length(Input(:));
            for iN = 1:NSol
                for iStep = 1:length(Input(iN).steps)
                    if isempty(Input(iN).steps(iStep).Psi)
                        if Args2.th
                            Res(iN).steps(iStep).th = nan;
                        end
                    else
                        if Args2.th
                            tPsi = Input(iN).steps(iStep).Psi;
                            tPsi = tPsi ./ norm(tPsi);
                            tPsi0 = obj.calcObj.Psi0(tPsi);
                            tPsi0 = tPsi0 ./ norm(tPsi0);
                            if tPsi0(1) < 0; tPsi0 = -tPsi0; end
                            Res(iN).steps(iStep).th = asin(tPsi0(2));
                        end
                    end
                end
            end
        end
    end
end