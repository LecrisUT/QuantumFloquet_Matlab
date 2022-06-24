classdef Details
    properties (SetAccess=private)
        script
        model
        calculation
        objects     (1,:)   Calc.baseCalc
        variables
        notes
    end
    methods
        function obj=Details(Args)
            arguments
                Args.script
                Args.model
                Args.calculation
                Args.objects
                Args.variables
                Args.notes      = {}
            end
            obj.script = Args.script;
            obj.model = Args.model;
            obj.calculation = Args.calculation;
            obj.objects = Args.objects;
            obj.variables = Args.variables;
            obj.notes = Args.notes;
        end
        function json=jsonencode(obj,varargin)
            S.script = obj.script;
            S.calculation = obj.calculation;
            S.model = obj.model;
            S.variables = obj.variables;
            S.notes = obj.notes;
            json = jsonencode(S,varargin{:});
        end
    end
end