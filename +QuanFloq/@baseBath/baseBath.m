classdef (Abstract) baseBath < QuanFloq.baseCalc
    % baseBath Base class describing bath properties
    %
    % See also baseCalc, baseSystemBath

    properties
        % systemBath - Associated system-bath coupling
        % See also baseSystemBath
        systemBath  QuanFloq.baseSystemBath {mustBeScalarOrEmpty}
        % beta - Inverse temperature
        % See also baseBath.NB
        beta        (1,1)   double
    end
    properties (SetAccess=private)
        % type - Type of bath
        % Either 'fermion' or 'boson'
        % See also baseBath.NB
        type
    end
    methods (Abstract)
        nu = nu(obj,omega);
    end
    methods
        function obj = baseBath(type,Args)
            arguments
                type   {mustBeMember(type,{'fermion','boson'})}
                Args.beta   (1,1)   double  = inf
            end
            % baseBath Constructor
            %
            % Syntax:
            %   obj = baseBath(type)
            %   [___] = baeBath(___,Name,Value)
            % 
            % Description:
            %   obj = baseBath(type) Constructs the bath object of specified type
            %   [___] = baseBath(___,Name,Value) specifies options using name-value
            %   arguments in addition to any of the input arguments in previous
            %   syntaxes.
            %
            % Inputs:
            %   type - Type of bath. Either 'fermion' or 'boson'
            %   Name-Value pairs
            %
            % Outputs:
            %   obj - Floquet object
            %
            % Name-value arguments:
            %   beta - [inf] Inverse temperature
            % 
            % See also baseBath

            obj.beta = Args.beta;
            obj.type = type;
        end
        gamma = gamma(obj,omega)
        NB = NB(obj,omega)
    end
end