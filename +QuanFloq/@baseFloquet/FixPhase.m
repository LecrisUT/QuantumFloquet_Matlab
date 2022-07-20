function Psi = FixPhase(obj,Psi)
    arguments
        obj     QuanFloq.baseFloquet
        Psi     double
    end
    % FixPhase Fix the wave function phase
    % Transform the wave function to have a pre-defined global phase convention
    % 
    % Syntax:
    %   Psi = FixPhase(Psi)
    % 
    % Description:
    %   Psi = FixPhase(Psi) Fix the phase of the wave functions Psi
    % 
    % Inputs:
    %   Psi - Wave functions to alter
    %
    % Outputs:
    %   Psi - Transformed wave functions with appropriate global phase

    Psi = obj.Psi_Floquet(Psi);
    Psi0 = obj.Psi0(Psi);
    for iM = 1:size(Psi,2)
        if Psi0(1,iM) < 0
            Psi(:,iM) = -Psi(:,iM);
        end
    end
end