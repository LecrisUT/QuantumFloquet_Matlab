function Ebar = Ebar(obj,Psi,Args)
    arguments
        obj     QuanFloq.baseFloquetHF
        Psi     double
        Args.normalize  (1,1)   logical = false
        Args.orbital    (1,1)   logical = false
    end
    % Ebar Calculate the average energies
    % Switch between effective and/or orbital using name-value flags
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
    %   Psi - Orbital wave functions
    %   Name-Value pairs
    %
    % Outputs:
    %   Ebar - Average energy
    %
    % Name-value arguments:
    %   effective - [false] Whether to calculate effective average energy or
    %   exact one
    %   normalize - [false] Whether to normalize the wave function
    %   orbital - [false] Whether to calculate orbital or many-body average
    %   energy
    %   
    % See also baseFloquet.Ebar, QuanFloq.baseFloquetHF.H, baseFloquetHF.FBar

    %% Make sure the wave function is in Floquet representation
    Psi = obj.Psi_Floquet(Psi);
    %% Check the wave function are in proper form
    if ~Args.orbital
        switch obj.mode
            case 'closed-shell'
                if size(Psi,2) ~= obj.Ne/2
                    error('Dimensions are not compatible')
                end
            otherwise
                error('Not implemented');
        end
        if Args.effective
            error('Many-body effective average energy is ill-defined in Floquet Hartree-Fock')
        end
        if ~Args.normalize
            warning('Many-body average energy is ill-defined without normalization. Re-normalizing the wave functions');
            Psi = Psi ./ vecnorm(Psi);
        end
        % Update the self-consistent operators
        obj.Psi = Psi;
    end
    %% Normalize the wave function if necessary
    if Args.normalize
        Psi = Psi ./ vecnorm(Psi);
    end
    %% Calculate the average energies
    if Args.effective
        %% Calculate the effective average energies
        Ebar = diag(Psi' * obj.F * Psi);
    else
        %% Calculate the actual average energies
        if ~Args.orbital
            % Requires the quasi-energies
            eps = obj.eps(Psi,normalize=true,orbital=true);
            % Calculate the average energy operator
            FBar = obj.FBar(Psi./vecnorm(Psi),eps);
            % Calculate the average energy expecation values
            Ebar = diag(Psi' * FBar * Psi);
        else
            Ebar = sum(diag(Psi' * obj.H * Psi));
        end
    end
end