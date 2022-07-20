classdef BosonBath2 < Calc.baseBath
    % bosonBath1 Boson bath with Ohmic spectral density
    %
    % See also BosonBath1

    properties
        % nu0 - Bath's spectral density constant
        % See also Calc.BosonBath2.nu
        nu0
    end
    methods
        function obj = BosonBath2(nu0,Args2)
            arguments
                nu0
                Args2.beta
            end
            % BosonBath2 Constructor
            %
            % Syntax:
            %   obj = BosonBath2(nu0)
            %   [___] = BosonBath2(___,Name,Value)
            % 
            % Description:
            %   obj = BosonBath2(nu0) Construct a boson bath with ohmic spectral
            %   density with constant nu0
            %   [___] = BosonBath2(___,Name,Value) specifies options using name-value
            %   arguments in addition to any of the input arguments in previous
            %   syntaxes.
            %
            % Inputs:
            %   nu0 - Spectral density constant
            %   Name-Value pairs
            %
            % Outputs:
            %   obj - Floquet object
            %
            % Name-value arguments:
            %   beta - Inverse temperature
            % 
            % See also BosonBath2, Calc.baseBath.baseBath

            Args2 = namedargs2cell(Args2);
            obj@Calc.baseBath(Args2{:});
            obj.nu0 = nu0;
        end
        function val = nu(obj,omega)
            % nu Ohmic spectra density
            % nu = nu0*omega
            % 
            % See also Calc.baseBath.nu
            val = obj.nu0 * abs(omega);
        end
    end
end