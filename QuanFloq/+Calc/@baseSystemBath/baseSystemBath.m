classdef baseSystemBath < Calc.baseCalc
    % baseSystemBath Base class for system-bath interaction
    % Defines the details of the system bath interaction
    %
    % See also baseFloquet, baseBath

    properties
        % system - Floquet system being coupled
        % See also baseFloquet
        system  Calc.baseFloquet    {mustBeScalarOrEmpty}
        % bath - Bath being coupled
        % See also baseBath
        bath    Calc.baseBath       {mustBeScalarOrEmpty}
        % hk_max - Cut-off Fourier coefficient
        % See also Calc.baseFloquet.hk_max
        hk_max  (1,1)   double  {mustBeInteger,mustBeNonnegative}
    end
    properties (Abstract)
        % V - Interaction operator
        % Contains only system part. Bath components are included via bath spectral
        % density.
        % Assumed to be static.
        %
        % See also Calc.baseFloquet.N, Calc.baseBath.nu
        V       (:,:)   double
    end
    properties (Dependent)
        % hk_range - Helper property for [-hk_max:hk_max]
        % See also Calc.baseSystemBath.hk_max
        hk_range
        % hk_max2 - Helper property for 2*hk_max+1
        % See also Calc.baseSystemBath.k_max
        hk_max2
    end
    methods
        function obj = baseSystemBath(system,bath,Args)
            arguments
                system  (1,1)   Calc.baseFloquet
                bath    (1,1)   Calc.baseBath
                Args.hk_max (1,1)   double
            end
            % baseSystemBath Constructor
            %
            % Syntax:
            %   obj = baseSystemBath(system,bath)
            %   [___] = baseSystemBath(___,Name,Value)
            % 
            % Description:
            %   obj = baseSystemBath(system,bath) Constructs the base coupling object
            %   between system and bath
            %   [___] = baseSystemBath(___,Name,Value) specifies options using
            %   name-value arguments in addition to any of the input arguments in
            %   previous syntaxes.
            %
            % Inputs:
            %   system - System object
            %   bath - Bath object
            %   Name-Value pairs
            %
            % Outputs:
            %   obj - System-bath coupling object
            %
            % Name-value arguments:
            %   hk_max - [2*system.hk_max] Cut-off Fourier coefficient
            % 
            % See also baseFloquet

            obj.system = system;
            obj.bath = bath;
            obj.system.systemBath = obj;
            obj.bath.systemBath = obj;
            if ~isfield(Args,'hk_max')
                Args.hk_max = 2 * obj.system.hk_max;
            end
            obj.hk_max = Args.hk_max;
        end
    end
    methods
        function val = get.hk_max2(obj)
            val = 2 * obj.hk_max + 1;
        end
        function val = get.hk_range(obj)
            val = -obj.hk_max:obj.hk_max;
        end
    end
    methods
        VSB = VSB(obj,Psi)
        [S,om_mnk] = SpectraOverlap(obj,Psi1,Psi2,Args2)
        rho = SteadyState(obj,Args)
        L = Lindblad(obj,Args)
        Res = Gamma(obj,Args)
    end
end