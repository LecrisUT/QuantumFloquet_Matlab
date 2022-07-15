function L = Lindblad(obj,Args)
    arguments
        obj     Calc.baseSystemBath
        Args.Gamma  (:,:,:) double  {mustBeReal,mustBePositive}
    end
    % Lindblad Calculate the Lindbladian master equation
    % 
    % Syntax:
    %   L = Lindblad
    %   [___] = Lindblad(___,Name,Value)
    % 
    % Description:
    %   L = Lindblad Calculates the full Lindbladian master equation
    %   [___] = Lindblad(___,Name,Value) specifies options using name-value
    %   arguments in addition to any of the input arguments in previous
    %   syntaxes.
    % 
    % Inputs:
    %   Name-Value pairs
    %
    % Outputs:
    %   L - Lindbladian master equation
    %
    % Name-value arguments:
    %   Gamma - Dissipation strength. If not provided it is calculated
    %   automatically in the average energy basis
    %   
    % See also baseSystemBath.Gamma, baseSystemBath.SteadyState

    if ~isfield(Args,'Gamma')
        Args.Gamma = obj.Gamma(type='Lindblad').Gamma;
    end
    Args.Gamma = sum(Args.Gamma,3);
    M = size(Args.Gamma,1);
    L = zeros(M);
    for in = 1:M
        % \partial_t \rho_n(t) -= 2 * \Gamma_{nm} \rho_n(t)
        L(in,in) = -2 * sum(Args.Gamma(in,:));
        for im = 1:M
            % \partial_t\rho_n(t) += 2 * \Gamma_{mn} \rho_m(t)
            L(in,im) = L(in,im) + 2 * Args.Gamma(im,in);
        end
    end
end