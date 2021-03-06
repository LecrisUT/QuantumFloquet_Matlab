clc; clear all;


function Res=NonSCF(obj,U_range,Args)
    arguments
        obj
        U_range
        Args.saveFile
        Args.loadFile
        Args.res_HF0
        Args.tol    = 1E-12
        Args.th_range
        Args.Psi0   = [1;0]
        Args.Print  = true
    end
    %% Get the static HF solutions if needed
    if ~isfield(Args,'res_HF0')
        if Args.Print
            fprintf('Calculating the static HF solutions\n');tic;
        end
        Hub0 = Calc.HubardDimerOrbital2('t',obj.t,'v0',obj.v0);
        subArgs = {'Psi0', Args.Psi0};
        if isfield(Args,'th_range')
            subArgs = [subArgs(:)',{'th_range'},{Args.th_range}];
        end
        Args.res_HF0 = Hub0.eigs(U_range,'tol',1E-12,...
            'Print',false,subArgs{:});
        if Args.Print
            fprintf('Finished calculating the static HF solutions (%fs)\n',toc);
        end
    end
    %% Calculate the non-SCF solutions
    nU = length(U_range);
    if isfield(Args,'loadFile')
        iU = load(Args.loadFile,'iU');
        Res = load(Args.loadFile,'Res');
    else
        Res(nU,obj.N) = struct('Psi',[],'eps',[],'th',[]);
        iU = 1;
    end
    for iU = iU:nU
        %% Project the time-periodic potential on this basis
        psi_HF = reshape([Args.res_HF0(iU,:).psi],obj.N,obj.N);
        e_HF = [Args.res_HF0(iU,:).eps];
        blk_diag = repmat({psi_HF},obj.k_max2,1);
        T_mat = blkdiag(blk_diag{:});
        vt = psi_HF' * [obj.v1/4 0;0 -obj.v1/4] * psi_HF;
        %% Construct the Floquet Hamiltonian
        vt_diag = spdiags(vt,-1:1);
        vt_diagp = spdiags(vt',-1:1);
        % Don't need to add the matrices because e_HF is diagonal
        h_diag = repmat([vt_diagp e_HF(:,iU) vt_diag],obj.k_max2,1);
        th = spdiags(h_diag,-3:3,obj.N*obj.k_max2,obj.N*obj.k_max2);
        thf = th + obj.pt;
        %% Calculate the eigenstates
        [Psi,eps] = eigs(thf,obj.N,0);
        % Project to original basis
        Psi = T_mat * Psi;
        Psi0 = obj.Psi0(Psi);
        % Fix the phase
        for iN=1:N
            if Psi0(1,iN)<0
                Psi(:,iN) = -Psi(:,iN);
                Psi0(:,iN) = -Psi0(:,iN);
            end
        end
        %% Centralize the spectrum at initial step, otherwise use adiabatic
        if iU == 1
            Psi = shiftCenter(Psi,N,tol);
            if size(Psi,2) ~= N
                error('Eigenstates contained replicas');
            end
            % Project back to the previous basis before calculating the QE
            eps = diag(Psi' * T_mat * thf * T_mat' * Psi);
            % Order initially by the average energy
            Ebar = diag(Psi' * T_mat * th * T_mat' * Psi);
            [~,ind] = sort(Ebar);
            Psi = Psi(:,ind);
            eps = eps(ind);
            Psi0 = Psi0(:,ind);
        end
        %% Reorder to adiabatically continued
        if iU > 1
            % Get previous values
            Psi_prev = reshape([Res(iU,:).Psi],[],obj.N);
            Psi0_prev = obj.Psi0(Psi_prev);
            eps_prev = [Res(iU-1,:).eps];
            % Overlap with previous
            S0 = Psi0_prev' * Psi0;
            % Get the closest indeces with respect to the previous WF (row)
            [~,ind] = max(abs(S0),[],1);
            Psi = Psi(:,ind);
            Psi0 = Psi0(:,ind);
            eps = eps(ind);
            eps = eps(:);
            % Shift to adiabatically continue the QE
            dk_eps = round((eps - eps_prev - (mod(eps-eps_prev+obj.w/2,obj.w) - obj.w/2))/obj.w);
            for iN =1:N
                Psi(:,iN) = circshift(Psi(:,iN), obj.N * dk_eps(iN), 1);
                eps(iN) = Psi(:,iN)' * T_mat * thf * T_mat' * Psi(:,iN);
            end
        end
        %% Calculate th
        th = asin(Psi0(2,:));
        %% Save the results
        for iN = 1:obj.N
            Res(iU,iN).Psi = Psi(:,iN);
            Res(iU,iN).eps = eps(iN);
            Res(iU,iN).th = th(iN);
        end
        if isfield(Args,'saveFile')
            save(Args.saveFile,'Res','iU','obj','-v7.3');
        end
        if Args.Print
            fprintf('Done [%d/%d]\n',iU,nU);
        end
    end
end