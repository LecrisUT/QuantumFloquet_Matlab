function eps = eps(obj,Psi,Args)
    arguments
        obj     QuanFloq.baseFloquet
        Psi     double
        Args.normalize  (1,1)   logical = false
    end
    % eps Calculate the quasi-energies
    % 
    % Syntax:
    %   eps = eps(Psi)
    %   eps(___,Name,Value)
    % 
    % Description:
    %   eps = eps(Psi) Calculate the quasi-energies of the wave functions Psi
    %   eps(___,Name,Value) Alter the calculation with name-value pairs
    % 
    % Inputs:
    %   Psi - Wave functions
    %   Name-Value pairs
    %
    % Outputs:
    %   eps - Quasi-energies
    %
    % Name-value arguments:
    %   normalize - [false] Whether to normalize the wave function
    %   
    % See also baseFloquet.eigs

    % Make sure the wave function is in Floquet representation
    Psi = obj.Psi_Floquet(Psi);
    % Normalize the wave function if necessary
    if Args.normalize
        Psi = Psi ./ vecnorm(Psi);
    end
    % Calculate the quasi-energies
    eps = diag(Psi' * obj.hf * Psi);
end