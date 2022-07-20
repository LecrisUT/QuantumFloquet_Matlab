function Ebar = Ebar(obj,Psi,Args)
    arguments
        obj     Calc.baseFloquet
        Psi     double
        Args.effective  (1,1)   logical = false
        Args.normalize  (1,1)   logical = false
    end
    % Ebar Calculate the average energies
    % Switch between effective or exact using name-value flag
    % 
    % Syntax:
    %   Ebar = Ebar(Psi)
    %   Ebar(___,Name,Value)
    % 
    % Description:
    %   Ebar = Ebar(Psi) Calculate the average energies
    %   Ebar(___,Name,Value) Alter the calculation with name-value pairs
    % 
    % Inputs:
    %   Psi - Wave function
    %   Name-Value pairs
    %
    % Outputs:
    %   Ebar - (Effective) Average energy
    %
    % Name-value arguments:
    %   effective - [false] Whether to calculate effective average energy or
    %   exact one
    %   normalize - [false] Whether to normalize the wave function
    %   
    % See also baseFloquet.HBar, Calc.baseFloquet.xi

    %% Make sure the wave function is in Floquet representation
    Psi = obj.Psi_Floquet(Psi);
    %% Normalize the wave function if necessary
    if Args.normalize
        Psi = Psi ./ vecnorm(Psi);
    end
    %% Calculate the average energies
    if Args.effective
        %% Calculate the effective average energies
        Ebar = diag(Psi' * obj.h * Psi);
    else
        %% Calculate the actual average energies
        % Requires the quasi-energies
        eps = obj.eps(Psi,normalize=true);
        % Calculate the average energy operator
        HBar = obj.HBar(Psi./vecnorm(Psi),eps);
        % Calculate the average energy expecation values
        Ebar = diag(Psi' * HBar * Psi);
    end
end