function [FBar, Psi, eps] = FBar(obj,Psi,eps,Args)
    arguments
        obj     Calc.baseFloquetHF
        Psi     double  = get_Psi(obj)
        eps     (1,:)   double  {mustBeReal}    = obj.eps(Psi,normalize=true,orbital=true)
        Args.ignoreCheck    (1,1)   logical     = false
        Args.projectBack    (1,1)   logical     = false
    end
    % FBar Calculate the orbital average energy operator
    % 
    % Syntax:
    %   FBar = FBar
    %   FBar = FBar(Psi)
    %   FBar = FBar(Psi,eps)
    %   [FBar,Psi,eps] = FBar(___)
    %   [___] = FBar(___,Name,Value)
    % 
    % Description:
    %   FBar = FBar Calculate the full orbital average enery operator
    %   FBar = FBar(Psi) Calculate the orbital average enery operator around
    %   the states Psi
    %   FBar = FBar(Psi,eps) Include the quasi-energies
    %   [FBar,Psi,eps] = FBar(___) Output the ordered wave function and
    %   eigenstates
    %   [___] = FBar(___,Name,Value) specifies options using name-value
    %   arguments in addition to any of the input arguments in previous
    %   syntaxes.
    % 
    % Inputs:
    %   Psi - Orbital quasi-energy eigenstates
    %   eps - Orbital quasi-energy eigenvalues
    %   Name-Value pairs
    %
    % Outputs:
    %   FBar - Orbital average energy operator
    %   Psi - (Re)Ordered eigenstates
    %   eps - (Re)Ordered quasi-energies
    %
    % Name-value arguments:
    %   ignoreCheck - [false] Whether to ignore check of the quasi-energy
    %   projectBack - [false] Whether to project back to original basis
    %   
    % See also Calc.baseFloquet.HBar, Calc.baseFloquetHF.F

    %% Checks
    % Make sure the wave function is in Floquet representation
    Psi = obj.Psi_Floquet(Psi);
    M = size(Psi,2);
    % Check that quasi-energy and wave function size match
    if length(eps) ~= M
        error('Quasi-energy and wave function size are different');
    end
    % Make sure quasi-energies are ordered
    if ~issorted(eps)
        % Reorder the states
        [eps,ind] = sort(eps);
        Psi = Psi(:,ind);
    end
    % Check that quasi-energies are within a single BZ
    if ~Args.ignoreCheck && (eps(end)-eps(1)) > obj.w
        error('Quasi-energies exceed a single BZ');
    end
    % Check that quasi-energies are sufficiently far appart to define the
    % average energy operator
    if ~Args.ignoreCheck && abs(eps(end)-(eps(1)+obj.w)) < obj.xi
        error('Difference in quasi-energy end-points is too small. Try using shiftBZ.');
    end
    %% Calculate the average energy operator
    FBar = Psi' * obj.F * Psi;
    for iM = 1:M-1
        for iN = iM+1:M
            % Check for near resonance
            if eps(iN)-eps(iM) < obj.xi; continue; end
            % Set non-resonant components to 0
            FBar(iM,iN) = 0;
            FBar(iN,iM) = 0;
        end
    end
    if Args.projectBack
        % Include all quasi-energy replicas
        Psi_full = zeros(obj.N * obj.k_max2,M * obj.k_max2);
        for k = obj.k_range
            Psi_full(:,M * (obj.k_max+k) + (1:M)) = circshift(Psi,obj.N * k,1);
        end
        % Block diagonalize the average energy operator
        FBar_full = repmat(FBar,1,obj.k_max2);
        FBar_full = mat2cell(FBar_full,N,ones(1,obj.k_max2));
        FBar_full = blkdiag(FBar_full{:});
        % Project the operator back to the original basis
        FBar = Psi_full * FBar_full * Psi_full';
    else
        %% Permute the indeces back to original representation
        if nargout < 2 && ~issorted(eps)
            % Reorder back the states
            warning('Quasi-energies were not sorted');
            [~,ind2] = sort(ind);
            FBar = FBar(ind2,ind2);
        end
    end
end
function Psi=get_Psi(obj)
    [Psi,~] = obj.eigs;
end