function varargout = eigs(obj,Args,Args2)
    arguments
        obj             (1,1)   Calc.baseFloquet
        Args.iterator   (1,1)   Calc.baseCalcIterator
        Args.Print      (1,1)   logical = true
        Args2.eps0      (1,1)   double
        Args2.Psi_prev  (:,:)   double
        Args2.eps_prev  (:,1)   double
    end
    % eigs Calculate the Floquet eigenstates/eigenvalues
    % 
    % Syntax:
    %   eps = eigs
    %   [Psi,eps] = eigs
    %   [Psi,eps,Ebar] = eigs
    %   eigs(___,Name,Value)
    % 
    % Description:
    %   eps = eigs Calculate the quasi-energy eigenvalues
    %   [Psi,eps] = eigs Output the quasi-energy eigenstates as well
    %   [Psi,eps,Ebar] = eigs Calculate the average energy eigenstates
    %   eigs(___,Name,Value) Alter the calculation with name-value pairs
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
    %   Print - Print progress for adiabatic continuation
    %   eps0 - Center the (initial) quasi-energies around eps0
    %   Psi_prev - Adiabatically match the eigenstates with Psi_prev. Can be
    %   static wave functions.
    %   eps_prev - Adiabatically match the quasi-energies as well
    %   
    % See also baseFloquet.HBar, Calc.baseFloquet.xi
    
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
        Args2Cells = namedargs2cell(Args2);
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
        if isfield(Args2,'eps0')
            eps0 = Args2.eps0;
        else
            eps0 = 0;
        end
        thf = obj.hf;
        %% Calculate the QE eigenstates
        % TODO: Allow for multiple eps0 values
        if nargout == 1
            % Just calculate the eigenvalues centered around eps0
            eps = eigs(thf,obj.N,eps0);
        else
            % Calculate the eigenstates around eps0 including the nearest
            % BZs
            try
                [Psi,eps] = eigs(thf,obj.N * 3,eps0);
            catch err
                switch err.identifier
                    % eigs fails if singular so have to shift eps0
                    case 'MATLAB:eigs:AminusBSingular'
                        [Psi,eps] = eigs(thf,obj.N * 3,eps0 + 1E-12);
                    otherwise
                        rethrow(err);
                end
            end
            eps = diag(eps);
            % Filter the BZ so that the biggest QE difference is at the end
            [eps,Psi] = obj.shiftBZ(eps,Psi);
        end
        %% Calculate the FAE eigenstates
        if nargout > 2
            %% Calculate the FAE eigenstates
            Hbar = obj.HBar(Psi,eps);
            [Psi2,Ebar] = eig(Hbar,'vector');
            Psi = Psi * Psi2;
            eps = diag(Psi' * thf * Psi);
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
            eps = diag(Psi' * thf * Psi);
        end
        %% Adiabatically continue
        if nargout > 1 && isfield(Args2,'Psi_prev')
            % Get previous values
            switch size(Args2.Psi_prev,1)
                case obj.N
                    % Special case if input static WF
                    Psi0_prev = Args2.Psi_prev;
                case obj.N * obj.k_max2
                    Psi0_prev = obj.Psi0(Args2.Psi_prev);
                otherwise
                    error('Wrong size for Psi_prev (1)');
            end
            if size(Args2.Psi_prev,2) ~= obj.N
                error('Wrong size for Psi_prev (2)');
            end
            Psi0 = obj.Psi0(Psi);
            % Overlap with previous
            S0 = Psi0_prev' * Psi0;
            % Get the closest indeces with respect to the previous WF (row)
            [mval,ind] = max(abs(S0),[],1);
            % Check if overlap is too small to adiabatically
            % continue
            if ~isempty(find(abs(1-mval)>1E-1, 1))
                warning('Overlap is too small');
            end
            % Check if missing values
            if length(ind) ~= length(unique(ind))
                error('Missing permutation index');
            end
            % Reorder to permute properly
            [~,ind] = sort(ind);
            Psi = Psi(:,ind);
            eps = eps(ind);
            if nargout > 2
                Ebar = Ebar(ind);
            end
            % Shift to adiabatically continue the QE
            if isfield(Args2,'eps_prev')
                dk_eps = round((eps - Args2.eps_prev - (mod(eps-Args2.eps_prev+obj.w/2,obj.w) - obj.w/2))/obj.w);
                for iN = 1:obj.N
                    Psi(:,iN)=circshift(Psi(:,iN), obj.N * dk_eps(iN), 1);
                    eps(iN) = Psi(:,iN)' * thf * Psi(:,iN);
                end
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