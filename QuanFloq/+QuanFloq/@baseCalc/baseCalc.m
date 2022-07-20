classdef baseCalc < handle & matlab.mixin.Heterogeneous & matlab.mixin.CustomDisplay
    % Calc.baseCalc Base calculation class
    % Base class from which other calculation inherit
    %
    % See also baseFloquet, baseSystemBath

    methods

        function json = jsonencode(obj,varargin)
            S.class = class(obj);
            json = jsonencode(S,varargin{:});
        end
    end
end