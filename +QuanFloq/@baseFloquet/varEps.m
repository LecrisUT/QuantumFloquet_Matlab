function varEps = varEps(obj,Psi,Args)
    arguments
        obj     QuanFloq.baseFloquet
        Psi     double
        Args.eps        (:,1)   double  {mustBeReal}
        Args.normalize  logical = false
    end
    % varEps Calculate the quasi-energy variance
    % 
    % Syntax:
    %   varEps = varEps(Psi)
    %   [___] = varEps(___,Name,Value)
    % 
    % Description:
    %   varEps = varEps(Psi) Calculate the quasi-energy variance of states Psi
    %   [___] = varEps(___,Name,Value) specifies options using name-value
    %   arguments in addition to any of the input arguments in previous
    %   syntaxes.
    % 
    % Inputs:
    %   Psi - Wave functions
    %   Name-Value pairs
    %
    % Outputs:
    %   varEps - Quasi-energy variance
    %
    % Name-value arguments:
    %   normalize - [false] Whether to normalize the wave functions
    %   eps - Pre-calculated quasi-energies
    %   
    % See also QuanFloq.baseFloquet.hf2, QuanFloq.baseFloquet.eps

    Psi = obj.Psi_Floquet(Psi);
    if Args.normalize
        Psi = Psi ./ vecnorm(Psi);
    end
    if ~isfield(Args,'eps')
        Args.eps = obj.eps(Psi,normalize=false);
    end
    varEps = diag(Psi' * obj.hf2 * Psi) - Args.eps .* Args.eps;
end