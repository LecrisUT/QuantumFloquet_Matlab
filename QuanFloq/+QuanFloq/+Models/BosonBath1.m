classdef BosonBath1 < QuanFloq.baseBath
    % BosonBath1 Boson bath with constant spectral density
    %
    % See also QuanFloq.Models.BosonBath2

    properties
        % nu0 - Bath's spectral density constant
        % See also QuanFloq.Models.BosonBath1.nu
        nu0
    end
    methods
        function obj = BosonBath1(nu0,Args2)
            arguments
                nu0
                Args2.beta
            end
            % BosonBath1 Constructor
            %
            % Syntax:
            %   obj = BosonBath1(nu0)
            %   [___] = BosonBath1(___,Name,Value)
            % 
            % Description:
            %   obj = BosonBath1(nu0) Construct a boson bath with constant spectral
            %   density nu0
            %   [___] = BosonBath1(___,Name,Value) specifies options using name-value
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
            % See also BosonBath1, QuanFloq.baseBath.baseBath

            Args2 = namedargs2cell(Args2);
            obj@QuanFloq.baseBath('boson',Args2{:});
            obj.nu0 = nu0;
        end
        function nu = nu(obj,~)
            % nu Constant spectral density
            % nu = nu0
            % 
            % See also QuanFloq.baseBath.nu
            nu = obj.nu0;
        end
    end
end