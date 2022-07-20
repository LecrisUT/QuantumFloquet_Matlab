function [Psi_Orb,Psi,eps,Ebar] = OccupiedOrbitals(obj,Psi,eps,Ebar,Args)
    arguments
        obj     Calc.baseFloquetHF
        Psi     double
        eps     (:,1)   double  {mustBeReal}    = []
        Ebar    (:,1)   double  {mustBeReal}    = []
        Args.occupation {mustBeMember(Args.occupation,{'aufbau','overlap'})}
        Args.Psi0_prev  (:,:)   double
    end
    % OccupiedOrbitals Filter the occupied orbitals
    % Depends on Hartree-Fock mode and occupation algorithm
    % 
    % Syntax:
    %   Psi_Orb = OccupiedOrbitals(Psi)
    %   [Psi_Orb,Psi,eps,Ebar] = OccupiedOrbitals(Psi,eps,Ebar)
    %   [___] = OccupiedOrbitals(___,Name,Value)
    % 
    % Description:
    %   Psi_Orb = OccupiedOrbitals(Psi) Filter the occupied orbitals. Assumes
    %   the orebitals are pre-ordered
    %   [Psi_Orb,Psi,eps,Ebar] = OccupiedOrbitals(Psi,eps,Ebar) Output
    %   reordered eigenpairs as well
    %   [___] = OccupiedOrbitals(___,Name,Value) specifies options using
    %   name-value arguments in addition to any of the input arguments in
    %   previous syntaxes.
    % 
    % Inputs:
    %   Psi - (Pre-ordered) orbital wave functions
    %   eps - Orbital quasi-enegies
    %   Ebar - Orbital average energies
    %   Name-Value pairs
    %
    % Outputs:
    %   Psi_Orb - Filtered orbital wave functions
    %   Psi - Reordered orbital wave functions
    %   eps - Reordered quasi-energies
    %   Ebar - Reordered average energies
    %
    % Name-value arguments:
    %   occupation - Algorithm determining how to occupy the orbitals for
    %   converging calculation:
    %      - aufbau: From smallest average energy onwards
    %      - overlap: The orbitals with highest overlap to previous trial one
    %   (If not specified assuming pre-ordered)
    %   Psi0_prev - Previous orbitals to compare in overlap algorithm
    %   
    % See also Calc.baseFloquetHF.mode
    
    Psi = obj.Psi_Floquet(Psi);
    if ~isempty(Ebar) && length(Ebar) ~= size(Psi,2)
        error('Dimension missmatch between Ebar and Psi');
    end
    if ~isempty(eps) && length(eps) ~= size(Psi,2)
        error('Dimension missmatch between eps and Psi');
    end
    %% Reorder the eigenstates
    if isfield(Args,'occupation')
        switch Args.occupation
            case 'aufbau'
                if isempty(Ebar)
                    error('Average energy not provided for aufbau occupation');
                end
                [Ebar,ind] = sort(Ebar);
                Psi = Psi(:,ind);
                if ~isempty(eps)
                    eps = eps(ind);
                end
                % TODO: Handle average energy degenerate case
            case 'overlap'
                if ~isfield(Args.Psi0_prev)
                    error('Previous orbitals not provided for overlap method')
                end
                [Psi,ind] = obj.AdiabaticContinue(Psi,Psi0_prev);
                if ~isempty(eps)
                    eps = eps(ind);
                end
                if ~isempty(Ebar)
                    Ebar = Ebar(ind);
                end
        end
    end
    %% Filter the occupied orbitals
    switch obj.mode
        case 'closed-shell'
            Psi_Orb = Psi(:,1:obj.Ne/2);
        otherwise
            error('Not implemented')
    end
    if nargout > 2 && isempty(eps)
        error('Quasi-energy was not provided for reordering');
    end
    if nargout > 3 && isempty(Ebar)
        error('Average energy was not provided for reordering');
    end
end