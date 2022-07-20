classdef Details
    % Details Metadata of a calculation
    properties (SetAccess=private)
        % script - Script used to perform the current calculation
        % Should be a plain/live script, not a function to keep track of inputs
        % See also Calc.Details.variables
        script
        % model - Model system being calculated
        % See also Calc.Details.objects
        model
        % calculation - Summary of the calculation
        % See also Calc.Details.notes
        calculation
        % objects - Calc objects used
        % See also baseCalc
        objects     (1,:)   Calc.baseCalc
        % variables - Basic variables of calculation
        % Backup some basic variable as a struct in case script changes
        % See also Calc.Details.script
        variables   (1,1)   struct
        % notes - Additional notes
        % Store additional notes about the calculation in a cell array
        notes       (1,:)   cell    {mustBeText}    = {}
    end
    methods
        function obj = Details(Args)
            arguments
%                 Args.?Calc.Details
                Args.script
                Args.model
                Args.calculation
                Args.objects
                Args.variables
                Args.notes      = {}
            end
            % Details Constructor
            %
            % Syntax:
            %   obj = Details(Name,Value)
            % 
            % Description:
            %   obj = Details(Name,Value) Properties are instantiated on construction
            %   using name-value pairs
            %
            % Inputs:
            %   Name-Value pairs
            %
            % Outputs:
            %   obj - Floquet object
            %
            % Name-value arguments:
            %   script - Script used
            %   model - Summary of model
            %   calculation - Summary of calculation
            %   objects - Calc objects used
            %   variables - Variables to backup
            %   notes - [{}] Additional notes
            % 
            % See also Details
            
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