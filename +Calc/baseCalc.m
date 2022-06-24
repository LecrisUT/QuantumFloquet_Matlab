classdef baseCalc < handle
    methods
        function json=jsonencode(obj,varargin)
            S.class = class(obj);
            json = jsonencode(S,varargin{:});
        end
    end
end