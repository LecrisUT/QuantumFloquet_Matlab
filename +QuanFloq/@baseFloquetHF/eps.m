function eps = eps(obj,Psi,Args)
    arguments
        obj     QuanFloq.baseFloquetHF
        Psi     double
        Args.normalize  (1,1)   logical = false
        Args.orbital    (1,1)   logical = true
    end
    % eps Calculate the quasi-energies
    % Switch between orbital or many-body quasi-energies using name-value flags
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
    %   orbital - [true] Whether to calculate orbital or many-body
    %   quasi-energies
    %   
    % See also baseFloquet.eps, baseFloquetHF.eigs

    %% Make sure the wave function is in Floquet representation
    Psi = obj.Psi_Floquet(Psi);
    %% Normalize the wave function if necessary
    if Args.normalize
        Psi = Psi ./ vecnorm(Psi);
    end
    %% Calculate the quasi-energies
    if Args.orbital
        eps = diag(Psi' * obj.Ff * Psi);
    else
        eps = diag(Psi' * obj.Hf * Psi);
    end
end