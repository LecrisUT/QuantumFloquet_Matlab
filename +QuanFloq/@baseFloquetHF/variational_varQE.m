function Res = variational_varQE(obj,Psi0,Args,Args2)
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
        Args2.TrackPsi      (1,1)   logical =false
    end
    % variational_varQE Perform quasi-energy variance minimization
    % See baseFloquet.variational_varQE
    %
    % Additional name-value arguments:
    %   SlaterCons - [true] Whether to include Slater determinant normalization
    %   constrain
    %   
    % See also baseFloquet.variational_varQE

    tol = Args2.tol;
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
            FunctionTolerance=Args2.tol,OptimalityTolerance=Args2.tol,...
            StepTolerance=Args2.tol,ConstraintTolerance=Args2.tol,...
            MaxIterations=1E4,MaxFunctionEvaluations=1E5,...
            FiniteDifferenceType='central');
    end
    if ~isfield(Args2,'ms')
        Args2.ms=MultiStart(UseParallel=false,...
            Display="iter",...
            FunctionTolerance=Args2.tol);
    end
    %% Run minimization
    Args2 = namedargs2cell(Args2);
    Res = obj.variational_base(@varHf,@cons,Psi0{:},Args2{:});
    %% Post-process results
    for iN = 1:length(Res)
        if Res(iN).Fval > 10 * tol
            Res(iN).conv = false;
        end
        % The gradient is larger than tolerance:
        % should be invalidated
        % TODO: Currently minimizing the aggresiveness
        if Res(iN).optim > max(1E-3,tol)
            Res(iN).conv = false;
        end
    end
    if Args.filter_nconv
        Res = Res([Res(:).conv]);
    end
    %%
    function vHf = varHf(x)
        x = reshape(x,[],NOrb);
        tPsi = zeros(Nk_max2,NOrb);
        % Force endpoints to be 0
        tPsi(~flag_end,:) = x(1:end-1,:);
        obj.Psi = tPsi;
        vHf = obj.varEps(tPsi,normalize=true,orbital=true);
    end
    %%
    function [c,ceq] = cons(x)
        x = reshape(x,[],NOrb);
        tPsi = zeros(Nk_max2,NOrb);
        % Force endpoints to be 0
        tPsi(~flag_end,:) = x(1:end-1,:);
        obj.Psi = tPsi;
        tnorm = norm(tPsi);
        c = [];
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
            % Convert to full Floquet representation
            tPsi0 = zeros(obj.N * obj.k_max2, M);
            tPsi0(obj.k_max * obj.N + (1:obj.N),:) = Psi0;
            Psi0 = tPsi0;
        case obj.N * (obj.k_max2-2*obj.hk_max)
            % Actual variation uses this size, but append zeros for consistent
            % manipulation
            Psi0 = [zeros(obj.N * obj.hk_max,M);
                Psi0;
                zeros(obj.N * obj.hk_max,M)];
        case obj.N * obj.k_max2
            % Already in full Floquet representation, do nothing
        otherwise
            error('Wrong size')
    end
end