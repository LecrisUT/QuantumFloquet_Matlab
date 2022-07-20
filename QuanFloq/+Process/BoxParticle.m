classdef BoxParticle < Process.baseFloquet
    methods
        function obj=BoxParticle(calcObj,Args)
            arguments
                calcObj
                Args.exactObj
            end
            Args = namedargs2cell(Args);
            obj@Process.baseFloquet(calcObj,Args{:});
        end
    end
end