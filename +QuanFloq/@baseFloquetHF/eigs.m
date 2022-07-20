function varargout = eigs(obj,Args,Args2,Args3)
    arguments
        obj             (1,1)   QuanFloq.baseFloquetHF
        Args.iterator   (1,1)   QuanFloq.baseCalcIterator
        Args.Print      (1,1)   logical = true
        Args2.eps0      (1,1)   double
        Args2.Psi_prev  (:,:)   double
        Args2.converge  (1,1)   logical = false
        Args2.occupation    {mustBeMember(Args2.occupation,{'aufbau','overlap'})}   = 'aufbau'
        Args2.alpha     (1,1)   double  {mustBeReal,mustBePositive} = 0.8
        Args2.tol       (1,1)   double  {mustBeReal,mustBePositive} = 1E-10
        Args3.eps_prev  (:,1)   double
    end
    % eigs Calculate the Floquet eigenstates/eigenvalues
    % 
    % Syntax:
    %   eps = eigs
    %   [Psi,eps] = eigs
    %   [Psi,eps,Ebar] = eigs
    %   [___] = eigs(___,Name,Value)
    % 
    % Description:
    %   eps = eigs Calculate the quasi-energy eigenvalues
    %   [Psi,eps] = eigs Output the quasi-energy eigenstates as well
    %   [Psi,eps,Ebar] = eigs Calculate the average energy eigenstates
    %   [___] = eigs(___,Name,Value) specifies options using name-value
    %   arguments in addition to any of the input arguments in previous
    %   syntaxes.
    % 
    % Inputs:
    %   Name-Value pairs
    %
    % Outputs:
    %   Psi - Quasi-energy/Average energy eigenstates
    %   eps - Quasi-energy eigenvalues
    %   Ebar - Average energy eigenvalues
    %
    % Name-value arguments:
    %   iterator - Calculate adiabatically continued eigenstates
    %   Print - [true] Print progress for adiabatic continuation
    %   eps0 - Center the (initial) quasi-energies around eps0
    %   Psi_prev - Adiabatically match the eigenstates with Psi_prev. Can be
    %   static wave functions.
    %   eps_prev - Adiabatically match the quasi-energies as well
    %   converge - [false] Whether to try to self-consistently converge
    %   occupation - ['aufbau'] See baseFloquetHF.OccupiedOrbitals
    %   alpha - [0.8] Mixing factor of the newly calculated eigenstates in the
    %   converging calculation
    %   tol - [1E-10] Overlap tolerance for stopping the self-consistent loop
    %   
    % See also baseFloquet.eigs, baseFloquetHF.FBar,
    % baseFloquetHF.OccupiedOrbitals
    
    %% Calculate the eigenstates
    if isfield(Args,'iterator')
        %% Calculate adiabatic eigenstates
        ind = 0;
        L = Args.iterator.size - Args.iterator.ind;
        eps = nan(obj.N,L);
        if nargout > 1
            Psi = nan(obj.N * obj.k_max2,obj.N,L);
        end
        if nargout > 2
            Ebar = nan(obj.N,L);
        end
        % If first itteration:
        % Set eps0: for BZ centered around eps0, otherwise call shiftCenter
        % Set Psi_prev/eps_prev: Adiabatically continue with given states
        Args2Cells = [namedargs2cell(Args2) namedargs2cell(Args3)];
        while ~Args.iterator.done
            Args.iterator.next;
            if (Args.iterator.object ~= obj)
                error('Iterator linked to different object');
            end
            ind = ind + 1;
            switch nargout
                case 1
                    eps(:,ind) = obj.eigs(Args2Cells{:});
                case 2
                    [Psi(:,:,ind), eps(:,ind)] = obj.eigs(Args2Cells{:});
                case 3
                    [Psi(:,:,ind), eps(:,ind), Ebar(:,ind)] = obj.eigs(Args2Cells{:});
            end
            Args2.Psi_prev = Psi(:,:,ind);
            Args2.eps_prev = eps(:,ind);
            % BZ center at first state's quasi-energy
            if ~isfield(Args2,'eps0')
                Args2.eps0 = eps(1,ind);
            end
            Args2Cells = namedargs2cell(Args2);
            if Args.Print
                fprintf('Done [%d/%d]\n',ind,L);
            end
        end
    else
        %% Calculate single eigenstate
        maxSteps = 1000;
        if Args.converge
            nSteps = maxSteps;
        else
            nSteps = 1;
        end
        if isfield(Args2,'Psi_prev')
            Psi0_prev = obj.Psi0(Args2.Psi_prev);
        else
            Psi0_prev = [];
        end
        if isfield(Args2,'eps0')
            eps0 = Args2.eps0;
        else
            eps0 = 0;
        end
        for iStep = 1:nSteps
            %% Iterate until self-consistent
            tFf = obj.Ff;
            %% Calculate the QE eigenstates
            % TODO: Allow for multiple eps0 values
            if nargout == 1
                % Just calculate the eigenvalues centered around eps0
                eps = eigs(tFf,obj.N,eps0);
            else
                % Calculate the eigenstates around eps0 including the nearest BZs
                try
                    [Psi,eps] = eigs(tFf,obj.N * 3,eps0);
                catch err
                    switch err.identifier
                        % eigs fails if singular so have to shift eps0
                        case 'MATLAB:eigs:AminusBSingular'
                            [Psi,eps] = eigs(tFf,obj.N * 3,eps0 + 1E-12);
                        otherwise
                            rethrow(err);
                    end
                end
                eps = diag(eps);
                % Filter the BZ so that the biggest QE difference is at the end
                [eps,Psi] = obj.shiftBZ(eps,Psi);
            end
            %% Calculate the FAE eigenstates
            if nargout > 2 || Args2.occupation == "aufbau"
                %% Calculate the FAE eigenstates
                Fbar = obj.FBar(Psi,eps);
                [Psi2,Ebar] = eig(Fbar,'vector');
                Psi = Psi * Psi2;
                eps = diag(Psi' * tFf * Psi);
                %% Sort eigenstates by average energy
                [Ebar,ind] = sort(Ebar);
                Psi = Psi(:,ind);
                eps = eps(ind);
            end
            %% Centralize the spectrum
            if nargout > 1 && ~isfield(Args2,'eps0')
                Psi = obj.shiftCenter(Psi);
                if size(Psi,2) ~= obj.N
                    error('Eigenstates contained replicas');
                end
                eps = diag(Psi' * tFf * Psi);
            end
            %% Reorder according to occupation algorithm
            % Dummy assingment in case of occupation is overlap
            if isempty(Psi0_prev)
                Psi0_prev = obj.Psi0(Psi);
            end
            if nargout > 2 || Args2.occupation == "aufbau"
                [Psi_Orb,Psi,eps,Ebar] = obj.OccupiedOrbitals(Psi,eps,Ebar,...
                    occupation=Args2.occupation,Psi0_prev=Psi0_prev);
            else
                [Psi_Orb,Psi,eps] = obj.OccupiedOrbitals(Psi,eps,...
                    occupation=Args2.occupation,Psi0_prev=Psi0_prev);
            end
            Psi0_Orb = obj.Psi0(Psi_Orb);
            %% Check convergence
            if iStep > 1
                minS0 = min(diag(abs(Psi0_prev' * Psi0_Orb)));
                if abs(1 - minS0) < Args2.tol
                    break;
                elseif iStep == maxSteps
                    warning('Floquet Hartree-Fock did not converge');
                end
            end
            %% Update wave functions
            if nSteps > 1
                %% Get current orbitals
                Psi_Curr = obj.FixPhase(obj.Psi);
                Psi_Orb = obj.FixPhase(Psi_Orb);
                %% Find the best matching BZ
                for iM = 1:size(Psi_Curr,2)
                    S_prev = -1;
                    for k = obj.k_range
                        ttPsi = circshift(Psi_Curr(:,iM),obj.N * k,1);
                        S = ttPsi' * Psi_Orb(:,iM);
                        if abs(S) > S_prev
                            tPsi = ttPsi;
                            S_prev = S;
                        end
                    end
                    Psi_Curr(:,iM) = tPsi;
                end
                %% Mix the wave functions
                Psi_Orb = Args2.alpha * Psi_Orb + (1-Args2.alpha) * Psi_Curr;
                Psi_Orb = Psi_Orb ./ vecnorm(Psi_Orb);
                obj.Psi = Psi_Orb;
            end
        end
        %% Adiabatically continue
        if nargout > 1 && isfield(Args2,'Psi_prev')
            Args3 = namedargs2cell(Args3);
            [Psi,ind] = obj.AdiabaticContinue(Psi,Args2.Psi_prev,Args3{:},eps=eps);
            eps = obj.eps(Psi,normalize=true);
            if nargout > 2
                Ebar = Ebar(ind);
            end
        end
    end
    %% Output the results
    switch nargout
        case 1
            varargout = {eps};
        case 2
            varargout = {Psi,eps};
        case 3
            varargout = {Psi,eps,Ebar};
    end
end