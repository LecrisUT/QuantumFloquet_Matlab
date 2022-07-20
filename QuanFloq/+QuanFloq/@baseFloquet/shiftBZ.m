function varargout = shiftBZ(obj,eps,Psi)
    arguments
        obj     QuanFloq.baseFloquet
        eps     (:,1)   double
        Psi     double
    end
    % shiftBZ Shift the Brillouin zone center 
    % Shifts the Brillouin zone so that the ends of it are maximally seperated
    % 
    % Syntax:
    %   BZc = shiftBZ(eps)
    %   [eps,Psi] = shiftBZ(eps)
    %   [eps,Psi] = shiftBZ(eps,Psi)
    % 
    % Description:
    %   BZc = shiftBZ(eps) Calculate the Brillouin zone center
    %   [eps,Psi] = shiftBZ(eps) Reorder the quasi-energies around the
    %   Brillouin zone center. Output a dummy Psi=[] to separate from former.
    %   [eps,Psi] = shiftBZ(eps,Psi) Reorder the quasi-energies and eigenstates
    %   around the Brillouin zone center.
    % 
    % Inputs:
    %   eps - Quasi-energy eigenvalues
    %   Psi - Quasi-energy eigenstates
    %
    % Outputs:
    %   eps - Quasi-energy eigenvalues
    %   Psi - Reordered Quasi-energy eigenstates. If not provided outputs empty
    %   BZc - Brillouin zone center

    % TODO: Relax this condition and check that all eigenvalues are properly
    % shifted
    if length(eps) ~= obj.N * 3
        error('The full quasi-energy eigenvalues of 3 BZ need to be provided');
    end
    %% Order eigenstates by quasi-energy
    [eps,ind] = sort(eps);
    %% Calculate the difference to next quasi-energy
    deps = diff(eps);
    deps = deps(obj.N + (1:obj.N));
    %% Find the maximum difference in eigenstates
    [~,ind2] = max(deps);
    %% Set state with maximum difference at the end
    % mod(eps(1)-eps(end),w) is maximum
    eps = eps(ind2 + (1:obj.N));
    BZc = (eps(1) + eps(end)) / 2;
    %% Reorder the wavefunction as well
    if nargin > 2
        %% Make sure the wave function is in Floquet representation
        Psi = obj.Psi_Floquet(Psi);
        if size(Psi,2) ~= obj.N * 3
            error('Wrong size of Psi');
        end
        %% Reorder according to the quasi-energies
        Psi = Psi(:,ind);
        %% Reorder according to the BZ center
        Psi = Psi(:,ind2 + (1:obj.N));
    end
    %% Format the output
    switch nargout
        case 1
            varargout = {BZc};
        case 2
            if nargin > 2
                varargout = {eps,Psi};
            else
                varargout = {eps,[]};
            end
    end
end