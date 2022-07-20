classdef TwoLevel_Boson < Calc.baseSystemBath
    % TwoLevel_Boson System bath coupling between two-level and boson bath
    % $$V_{SB}=\mqty[\sigma_z&\sigma_x-i\sigma_y\\\sigma_x+i\sigma_y&\sigma_z]$$
    %
    % See also Calc.baseFloquet

    properties
        % sigma - Pauli matrices of the system bath interaction
        % Represented as
        % $$\sigma=(\sigma_x,sigma_y,sigma_z)$$
        % See also Calc.TwoLevel_Boson.V
        sigma   (1,3)   double  {mustBeReal}
    end
    properties (Dependent)
        % V - System bath interaction
        % $$V_{SB}=\mqty[\sigma_z&\sigma_x-i\sigma_y\\\sigma_x+i\sigma_y&\sigma_z]$$
        % See also Calc.TwoLevel_Boson.sigma, Calc.baseSystemBath.V
        V
    end
    methods
        function obj = TwoLevel_Boson(system,bath,sigma)
            arguments
                system
                bath
                sigma
            end
            % TwoLevel_Boson Constructor
            %
            % Syntax:
            %   obj = TwoLevel_Boson(system,bath,sigma)
            % 
            % Description:
            %   obj = TwoLevel_Boson(system,bath,sigma) Constructs the system-bath
            %   coupling with the form sigma
            %
            % Inputs:
            %   system - System object
            %   bath - Bath object
            %   sigma - Pauli matrices of the system-bath interaction
            %
            % Outputs:
            %   obj - System-bath coupling object
            % 
            % See also TwoLevel_Boson, baseSystemBath

            Args2=namedargs2cell(system,bath);
            obj@Calc.baseSystemBath(Args2{:});
            obj.sigma = sigma;
        end
        function val = get.V(obj)
            val = [obj.sigma(3) obj.sigma(1)-1i*obj.sigma(2);
                obj.sigma(1)+1i*obj.sigma(2) -obj.sigma(3)];
        end
    end
end