function varEps = varEps(obj,Psi,Args)
    arguments
        obj     QuanFloq.baseFloquetHF
        Psi     double
        Args.eps        (:,1)   double  {mustBeReal}
        Args.normalize  (1,1)   logical = false
        Args.orbital    (1,1)   logical = true
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
    %   orbital - [true] Whether to calculate orbital or many-body quasi-energy
    %   variance
    %   
    % See also baseFloquet.varEps

    if ~Args.orbital
        error('Not implemented')
    end
    %% Make sure the wave function is in Floquet representation
    Psi = obj.Psi_Floquet(Psi);
    %% Normalize the wave function if necessary
    if Args.normalize
        Psi = Psi ./ vecnorm(Psi);
    end
    if ~isfield(Args,'eps')
        Args.eps = obj.eps(Psi);
    end
    %% Calculate the quasi-energy variance
    varEps = diag(Psi' * obj.Ff * obj.Ff * Psi) - Args.eps .* Args.eps;
end