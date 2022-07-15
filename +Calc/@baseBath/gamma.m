function gamma = gamma(obj,omega)
    arguments
        obj     Calc.baseBath
        omega   (:,1)   double
    end
    % gamma Calculate the bath coupling coefficient
    % 
    % Syntax:
    %   gamma = gamma(omega)
    % 
    % Description:
    %   gamma = gamma(omega) Calculate the bath coupling coefficient with
    %   frequency omega
    % 
    % Inputs:
    %   omega - Transition frequency
    %
    % Outputs:
    %   gamma - Bath coupling coefficient
    %   
    % See also Calc.baseBath.nu, baseBath.NB

    gamma = obj.nu(abs(omega)) .* obj.NB(-omega);
end