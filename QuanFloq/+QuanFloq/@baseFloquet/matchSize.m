function Psi = matchSize(obj,obj2,Psi2,Args)
    arguments
        obj     QuanFloq.baseFloquet
        obj2    QuanFloq.baseFloquet
        Psi2    double
        Args.normalize  (1,1)   logical = false
    end
    % matchSize Match size of wavefunctions
    % Makes sure that the wave functions of a partner Floquet object match the
    % current Floquet object's sizes. Truncate or fill zeros.
    % 
    % Syntax:
    %   Psi = matchSize(obj2,Psi2)
    %   matchSize(___,Name,Value)
    % 
    % Description:
    %   Psi = matchSize(obj2,Psi2) Shrink/Expand the wavefunction Psi2 of
    %   object obj2 to match the sizes of the current object
    %   matchSize(___,Name,Value) Alter the calculation with name-value pairs
    % 
    % Inputs:
    %   obj2 - Partner Floquet object associated with Psi2
    %   Psi2 - Wave function to shrink/expand
    %   Name-Value pairs
    %
    % Outputs:
    %   Psi - Equivalent wave function matching the size of current object
    %
    % Name-value arguments:
    %   normalize - [false] Whether to renormalize the wavefunction. Not
    %   recommended because it will represent a different wavefunction
    %   
    % See also QuanFloq.baseFloquet.N, QuanFloq.baseFloquet.k_max

    % Find the minimum sizes from which to take the values from
    minN = min(obj.N,obj2.N);
    mink = min(obj.k_max,obj2.k_max);
    % Make sure the wave function is in Fourier representation
    Psi2 = obj2.Psi_Fourier(Psi2);
    M = size(Psi2,2);
    % Initialize the wave functions with zeros
    Psi = zeros(obj.N,M,obj.k_max2);
    % Assign the relevant values
    Psi(1:minN,:,obj.k_max+1+(-mink:mink)) = Psi2(1:minN,:,obj2.k_max+1+(-mink:mink));
    % Transform back to Floquet representation
    Psi = obj.Psi_Floquet(Psi);
    % Renormalize if specified
    if Args.normalize
        Psi = Psi ./ vecnorm(Psi);
    end
end