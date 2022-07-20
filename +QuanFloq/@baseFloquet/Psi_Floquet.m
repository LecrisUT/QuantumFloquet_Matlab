function Psi = Psi_Floquet(obj,Psi)
    arguments
        obj     QuanFloq.baseFloquet
        Psi     double
    end
    % Psi_Floquet Transform a wave function in Floquet representation
    % 
    % Syntax:
    %   Psi = Psi_Floquet(Psi)
    % 
    % Description:
    %   Psi = Psi_Floquet(Psi) Transform Psi into Floquet representation
    % 
    % Inputs:
    %   Psi - Wave function to transform
    %
    % Outputs:
    %   Psi - Equivalent wave function in Floquet representation
    %
    % See also baseFloquet.Psi_Fourier

    Psi = obj.Psi_Fourier(Psi);
    Psi = reshape(permute(Psi,[1 3 2]),obj.N * obj.k_max2, []);
end