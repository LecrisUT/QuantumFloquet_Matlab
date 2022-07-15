function Psi0 = Psi0(obj,Psi,Args)
    arguments
        obj     Calc.baseFloquet
        Psi     double
        Args.normalize  (1,1)   logical = false
    end
    % Psi0 Calculate the wavefunction at $t=0$
    % 
    % Syntax:
    %   Psi0 = Psi0(Psi)
    %   Psi0(___,Name,Value)
    % 
    % Description:
    %   Psi0 = Psi0(Psi) Calculate $\ket{\Psi(0)}$ of Psi
    %   Psi0(___,Name,Value) Alter the calculation with name-value pairs
    % 
    % Inputs:
    %   Psi - Wave function to evaluate
    %   Name-Value pairs
    %
    % Outputs:
    %   Psi0 - Wave function at $t=0$
    %
    % Name-value arguments:
    %   normalize - [false] Whether to renormalize the wavefunction. Not
    %   recommended because it will represent a different wavefunction

    Psi0 = sum(obj.Psi_Fourier(Psi),3);
    if Args.normalize
        Psi0 = Psi0 ./ vecnorm(Psi0);
    end
end