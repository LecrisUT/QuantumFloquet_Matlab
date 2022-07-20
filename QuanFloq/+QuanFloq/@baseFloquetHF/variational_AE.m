function Res = variational_AE(obj,Psi0,Args,Args2)
    arguments
        obj     QuanFloq.baseFloquetHF
    end
    arguments (Repeating)
        Psi0    (:,:)   double
    end
    arguments
        Args.SlaterCons     (1,1)   logical = true
        Args.filter_nconv   (1,1)   logical = true
        Args2.opts
        Args2.ms
        Args2.TypX
        Args2.tol
        Args2.TrackPsi      (1,1)   logical = false
    end
    % variational_AE Perform average energy minimization
    % See baseFloquet.variational_AE
    %
    % Additional name-value arguments:
    %   SlaterCons - [true] Whether to include Slater determinant normalization
    %   constrain
    %   
    % See also baseFloquet.variational_AE

    %% Share variables
    Args2.filter_nconv = Args.filter_nconv;
    NOrb = obj.NOrbitals;
    %% Alter X0 variable
    % Variation should not vary the ends of the Fourier space because the
    % operators are ill-represented there
    if isempty(Psi0)
        error('No initial Psi0 provided')
    end
    for iM = 1:length(Psi0)
        Psi0{iM} = get_Psi0(obj,Psi0{iM});
    end
    % flag the endpoints
    Nk_max2 = obj.N * obj.k_max2;
    flag_end = false(Nk_max2,NOrb);
    flag_end(1:obj.N*obj.hk_max,:) = true;
    flag_end(end+1-(1:obj.N*obj.hk_max),:) = true;
    %% Specify default options
    if ~isfield(Args2,'opts')
        Args2.opts=optimoptions('fmincon',Algorithm='sqp',...
            SubproblemAlgorithm='cg',...
            FunctionTolerance=obj.xi,OptimalityTolerance=obj.xi,...
            StepTolerance=Args2.tol,ConstraintTolerance=Args2.tol,...
            MaxIterations=1E4,MaxFunctionEvaluations=1E5,...
            FiniteDifferenceType='central');
    end
    if ~isfield(Args2,'ms')
        Args2.ms=MultiStart(UseParallel=false,...
            Display="iter",...
            FunctionTolerance=obj.xi);
    end
    %% Run minimization
    Args2 = namedargs2cell(Args2);
    Res = obj.variational_base(Psi0,@Ebar,@cons,Args2{:});
    %% Post-process results
    for iN = 1:length(Res)
        if ~isempty(Res(iN).Psi)
            Res(iN).eps = Res(iN).Psi(end,:);
            Res(iN).Psi = Res(iN).Psi(1:end-1,:);
        else
            Res(iN).eps = nan;
        end
        for iStep = 1:length(Res(iN).steps)
            if isempty(Res(iN).steps(iStep).Psi)
                Res(iN).steps(iStep).eps = nan;
            else
                Res(iN).steps(iStep).eps = Res(iN).steps(iStep).Psi(end,:);
                Res(iN).steps(iStep).Psi = Res(iN).steps(iStep).Psi(1:end-1,:);
            end
        end
        % The gradient is larger than tolerance:
        % should be invalidated
        % TODO: Currently minimizing the aggresiveness
        if Res(iN).optim > max(1E-3,obj.xi)
            Res(iN).conv = false;
        end
    end
    if filter_nconv
        Res=Res([Res(:).conv]);
    end
    %%
    function Ebar = Ebar(x)
        x = reshape(x,[],NOrb);
        tPsi = zeros(Nk_max2,NOrb);
        % Force endpoints to be 0
        tPsi(~flag_end,:) = x(1:end-1,:);
        obj.Psi = tPsi;
        Ebar = obj.Ebar(tPsi,orbital=false,normalize=true);
    end
    %%
    function [c,ceq] = cons(x)
        x = reshape(x,[],NOrb);
        tPsi = zeros(Nk_max2,NOrb);
        % Force endpoints to be 0
        tPsi(~flag_end,:) = x(1:end-1,:);
        obj.Psi = tPsi;
        tnorm = norm(tPsi);
        grad = obj.Ff * tPsi - x(end,:) * tPsi;
        grad = grad(:);
        c = [grad - obj.xi;
            -grad - obj.xi];
        if Args.SlaterCons
            tPsi_Slater = obj.OrbitalToSlater(tPsi,normalize=false,normalizeSlater=false);
            tnorm_Slater = norm(tPsi_Slater);
            ceq = [tnorm - 1;
                tnorm_Slater - 1];
        else
            ceq = tnorm - 1;
        end
    end
end
function Psi0 = get_Psi0(obj,Psi0)
    M = size(Psi0,2);
    switch size(Psi0,1)
        case obj.N
            % Convert to full Floquet representation. Also missing quasi-energies
            Psi0 = [zeros(obj.N * obj.k_max,M);
                Psi0;
                zeros(obj.N * obj.k_max,M)];
            teps = obj.eps(Psi0);
            Psi0 = [Psi0;teps'];
        case obj.N + 1
            % Convert to full Floquet representation
            teps = Psi0(end,:);
            Psi0 = Psi0(1:end-1,:);
            Psi0 = [zeros(obj.N * obj.k_max,M);
                Psi0;
                zeros(obj.N * obj.k_max,M)];
            Psi0 = [Psi0;teps'];
        case obj.N * (obj.k_max2-2*obj.hk_max)
            % Actual variation uses this size, but append zeros for consistent
            % manipulation. Also missing quasi-energies
            Psi0 = [zeros(obj.N * obj.hk_max,M);
                Psi0;
                zeros(obj.N * obj.hk_max,M)];
            teps = obj.eps(Psi0);
            Psi0 = [Psi0;teps'];
        case obj.N * (obj.k_max2-2*obj.hk_max) + 1
            % Actual variation uses this size, but append zeros for consistent
            % manipulation
            teps = Psi0(end,:);
            Psi0 = Psi0(1:end-1,:);
            Psi0 = [zeros(obj.N * obj.hk_max,M);
                Psi0;
                zeros(obj.N * obj.hk_max,M)];
            Psi0 = [Psi0;teps'];
        case obj.N * obj.k_max2
            % Missing quasi-energy value
            teps = obj.eps(Psi0);
            Psi0 = [Psi0;teps'];
        case obj.N * obj.k_max2 + 1
            % Already proper sized, do nothing
        otherwise
            error('Wrong size')
    end
end