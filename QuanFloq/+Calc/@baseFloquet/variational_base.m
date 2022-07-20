function Res = variational_base(obj,fval,cons,Psi0,Args)
    arguments
        obj     Calc.baseFloquet
        fval
        cons
    end
    arguments (Repeating)
        Psi0    (:,:)   double
    end
    arguments
        Args.opts
        Args.ms
        Args.TypX           (:,1)   double  {mustBeReal,mustBePositive}
        Args.tol	        (1,1)	double  {mustBeReal,mustBePositive} = 1E-12
        Args.filter_nconv   (1,1)   logical = true
        Args.TrackPsi       (1,1)   logical = false
    end
    % variational_base Perform a minimization method
    % Defines default optimization options and formats the output for a Floquet
    % eigenstate calculation.
    % 
    % Syntax:
    %   Res = variational_base(Psi0,fval,cons)
    %   [___] = variational_base(___,Name,Value)
    % 
    % Description:
    %   Res = variational_base(Psi0,fval,cons) Perform the minimization with
    %   default arguments
    %   [___] = variational_base(___,Name,Value) specifies options using
    %   name-value arguments in addition to any of the input arguments in
    %   previous syntaxes
    % 
    % Inputs:
    %   Psi0 - (Repeating) Starting wave function guesses
    %   fval - Function to be minimized
    %   cons - Constraints
    %   Name-Value pairs
    %
    % Outputs:
    %   Res - Structured output vector array containing all (convergent)
    %   solutions. Post processed by calling function
    %
    % Name-value arguments:
    %   opts - Override gneerated optimoptions of 'fmincon'. See implementation
    %   for defaults
    %   ms - Override gneerated multistart options. See impelemntation for
    %   defaults
    %   TypX - Override generated TypX for opts
    %   tol - [1E-12] Convergence tolerance
    %   filter_nconv - [true] Whether to filter non-convergent solutions
    %   TrackPsi - [false] Whether to save the trial wave functions
    %   
    % See also MultiStart, fmincon

    %% Get basic meta data
    if isempty(Psi0)
        error('No initial Psi0 provided')
    end
    NPsi0 = length(Psi0);
    tPsi0 = Psi0{1};
    %% Detach variables
    TrackPsi = Args.TrackPsi;
    %% Fix vectors
    if ~isfield(Args,'TypX')
        Args.TypX = get_TypX(obj,tPsi0);
    end
    % flag the endpoints
    % Use size of Psi0 for consistency
    flag_end = false(size(tPsi0));
    flag_end(1:obj.N*obj.hk_max,:) = true;
    flag_end(obj.N*obj.k_max2+1-(1:obj.N*obj.hk_max),:) = true;
    % Trim the vectors
    for iN = 1:NPsi0
        Psi0{iN} = Psi0{iN}(~flag_end);
    end
    Args.TypX = Args.TypX(~flag_end);
    %% Initialize minimization options
    if ~isfield(Args,'opts')
        Args.opts=optimoptions('fmincon',Algorithm='sqp',...
            SubproblemAlgorithm='cg',...
            FunctionTolerance=Args.tol,OptimalityTolerance=Args.tol,...
            StepTolerance=Args.tol,ConstraintTolerance=Args.tol,...
            MaxIterations=1E4,MaxFunctionEvaluations=1E5,...
            FiniteDifferenceType='central');
    end
    Args.opts=optimoptions(Args.opts,...
        TypicalX = Args.TypX);
    if NPsi0 == 1
        Args.opts=optimoptions(Args.opts,...
            Display="iter");
    end
    if ~isfield(Args,'ms')
        Args.ms=MultiStart(UseParallel=false,...
            Display="iter",...
            FunctionTolerance=Args.tol);
    end
    pool = gcp('nocreate');
    if ~isempty(pool)
        Args.ms = MultiStart(Args.ms,UseParallel=true);
    end
    if TrackPsi
        iPsi=1;
        Args.opts=optimoptions(Args.opts,...
            OutputFcn=@getSteps);
        Args.ms=MultiStart(Args.ms,OutputFcn=@getMSStep, ...
            UseParallel=false);
    end
    %% Prepare output
    Res(NPsi0) = struct(Psi=[],steps=[],...
        conv=false,Fval=nan,optim=nan,NSteps=nan);
    for iN = 1:NPsi0
        Res(iN).Psi = [];
        Res(iN).conv = false;
        Res(iN).Fval = nan;
        Res(iN).optim = nan;
        Res(iN).NSteps = nan;
        if TrackPsi
            Res(iN).steps(Args.opts.MaxIterations).Psi=[];
        end
    end
    %% Prepare problem
    prob=createOptimProblem('fmincon',...
        objective=fval,nonlcon=cons,...
        options=Args.opts,x0=Psi0{1});
    %% Run minimization
    if NPsi0 > 1
        Psi0 = cell2mat(Psi0)';
        sPoints = CustomStartPointSet(Psi0);
        [~,~,flag,output,solutions] = run(Args.ms,prob,sPoints);
        for sol=solutions
	        for solX0 = sol.X0
                ind = ismembertol(Psi0,solX0{:}',Args.tol,...
			        ByRows=true);
                Res(ind).Psi = sol.X;
                Res(ind).Fval = sol.Fval;
                Res(ind).conv = true;
                Res(ind).optim = sol.Output.firstorderopt;
	        end
        end
        if Args.filter_nconv
            Res = Res([Res(:).conv]);
        end
    else
        [Res.Psi,Res.Fval,flag,output]=fmincon(prob);
        Res.optim = output.firstorderopt;
        if flag <= 0
            Res.conv = false;
        else
            Res.conv = true;
        end
    end
    %% Post-Process
    % Reshape Psi to original shape of Psi0 and pad with zeros
    for iN = 1:length(Res)
        if isempty(Res(iN).Psi); continue; end
        tPsi0 = zeros(size(tPsi0));
        tPsi0(~flag_end) = Res(iN).Psi;
        Res(iN).Psi = tPsi0;
    end
    % Get number of steps
    if TrackPsi
        for iN = 1:length(Res)
            for iStep = 1:Args.opts.MaxIterations
                if isempty(Res(iN).steps(iStep).Psi)
                    iStep = iStep - 1;
                    break;
                else
                    % Reshape Psi to original shape of Psi0 and pad with zeros
                    tPsi0 = zeros(size(tPsi0));
                    tPsi0(~flag_end) = Res(iN).steps(iStep).Psi;
                    Res(iN).steps(iStep).Psi = tPsi0;
                end
            end
            Res(iN).NSteps = iStep;
        end
        maxSteps = max([Res(:).NSteps]);
        for iN = 1:length(Res)
            Res(iN).steps = Res(iN).steps(1:maxSteps);
        end

    end
    %%
    function stop=getSteps(tPsi,oV,state)
        stop=false;
        % Skip either start or end steps
        if ~strcmp(state,'iter'); return; end
        if TrackPsi
            Res(iPsi).steps(oV.iteration+1).Psi = tPsi;
        end
    end
    %%
    function stop=getMSStep(oV,state)
        stop=false;
        if ~strcmp(state,'iter'); return; end
        Res(iPsi).Psi=oV.localsolution.X;
        iPsi=iPsi+1;
    end
end
function TypX = get_TypX(obj,Psi0)
    TypX = ones(size(Psi0));
    Lk = zeros(obj.k_max+1,1);
    for ik = 0:obj.k_max
        Lk(ik+1) = obj.Ladder(ik);
    end
    for k = 3+obj.hk_max:obj.k_max
        tX = Lk(k+1);
        tX = max(tX,1E-16);
        tX = min(tX,1);
        TypX(obj.N*(obj.k_max+k)+(1:obj.N),:) = tX;
        TypX(obj.N*(obj.k_max-k+1)+(-1:-1:-obj.N)+1,:) = tX;
    end
end