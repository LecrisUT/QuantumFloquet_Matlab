classdef baseFloquet < Calc.baseCalc
    properties
        N
        w
        k_max
        hk_max
        xi
    end
    properties (Dependent)
        k_max2
        hk_max2
        pt
        hf
    end
    properties (Dependent,Abstract)
        h
    end
    properties (SetAccess=protected)
        sMat
    end
    properties (SetAccess=private)
        sInd
        cacheMat
    end
    methods
%         function obj = baseFloquet(N,w,hk_max,k_max,Args)
%             arguments
%                 N
%                 w
%                 hk_max      = 1
%                 k_max       = 100
%                 Args.xi     = 1E-6
%                 Args.diag   = true
%                 Args.cacheMat = true
%             end
        function obj = baseFloquet(N,w,k_max,hk_max,Args)
            arguments
                N
                w
                k_max       = 100
                hk_max      = 1
                Args.xi     = 1E-6
                Args.diag   = true
                Args.cacheMat = true
            end
            obj.N = N;
            obj.w = w;
            obj.k_max = k_max;
            obj.hk_max = hk_max;
            obj.xi = Args.xi;
            obj.cacheMat = Args.cacheMat;
            if Args.diag
                % TODO: Only migrated the helper, not generalized
                % Construct sparse Index helper
                tcol = [];
                trow = [];
                for irow = 1:obj.N * obj.k_max2
                    % U has to do U_diag'; ans(:)
                    tcol = [tcol irow + obj.N * (-obj.hk_max:obj.hk_max)];
                    trow = [trow irow*ones(1,obj.hk_max2)];
                end
                obj.sInd.val_flag = tcol > 0 & tcol <= obj.N * obj.k_max2;
                obj.sInd.row = trow(obj.sInd.val_flag);
                obj.sInd.col = tcol(obj.sInd.val_flag);
                obj.sInd.nnz = sum(obj.sInd.val_flag);
                obj.sMat=spalloc(obj.N*obj.k_max2,obj.N*obj.k_max2,obj.sInd.nnz);
            else
                error('Not implemented')
            end
        end
        function k_max2 = get.k_max2(obj)
            k_max2 = 2 * obj.k_max + 1;
        end
        function hk_max2 = get.hk_max2(obj)
            hk_max2 = 2 * obj.hk_max + 1;
        end
        function pt = get.pt(obj)
            pt=[obj.k_max:-1:-obj.k_max;obj.k_max:-1:-obj.k_max] * obj.w;
            pt=pt(:);
            pt=spdiags(pt,0,obj.N * obj.k_max2,obj.N * obj.k_max2);
        end
        function hf = get.hf(obj)
            hf = obj.h + obj.pt;
        end
        function Ebar = Ebar(obj,psi,Args)
            arguments
                obj
                psi
                Args.normalize  logical = false
            end
            if Args.normalize
                psi = psi ./ sqrt(norm(psi));
            end
            Ebar = psi' * obj.h * psi;
        end
        function eps = eps(obj,psi,Args)
            arguments
                obj
                psi
                Args.normalize  logical = false
            end
            if Args.normalize
                psi = psi ./ sqrt(norm(psi));
            end
            eps = psi' * obj.hf * psi;
        end
        function varEps = varEps(obj,psi,Args)
            arguments
                obj
                psi
                Args.normalize  logical = false
            end
            if Args.normalize
                psi = psi ./ sqrt(norm(psi));
            end
            thf = obj.hf;
            eps = obj.eps(psi);
            varEps = psi' * thf * thf * psi - eps * eps;
        end
        function psi0 = Psi0(obj,psi)
            psi0 = sum(permute(reshape(psi,obj.N, obj.k_max2, []),[1 3 2]),3);
        end
        function [eps,Psi,Ebar] = eigs(obj,Args)
            arguments
                obj
                Args.eps0
                Args.Psi_prev
                Args.eps_prev
            end
            % EIGS Calculate the Floquet eigenstates
            %   eps=EIGS('eps0',) Returns the quasi-energies centered around eps0
            %   [eps,Psi]=EIGS  Returns the quasi-energies and eigenstates
            %   centeted around the spectra center
            %   [eps,Psi,Ebar]=EIGS Returns the average-energy eigenstates

            if isfield(Args,'eps0')
                eps0 = Args.eps0;
            else
                eps0 = 0;
            end
            thf = obj.hf;
            %% Calculate the QE eigenstates
            [Psi,eps] = eigs(thf,obj.N,eps0);
            eps = diag(eps);
            %% Centralize the spectrum
            if nargout > 1 && ~isfield(Args,'eps0')
                Psi = obj.shiftCenter(Psi);
                if size(Psi,2) ~= obj.N
                    error('Eigenstates contained replicas');
                end
                eps = diag(Psi' * thf * Psi);
            end
            %% Calculate the FAE eigenstates
            if nargout > 2
                %% Construct the FAE operator
                Hbar = Psi' * obj.h * Psi;
                for iN = 1:obj.N-1
                    for jN = iN+1:obj.N
                        d_eps = eps(iN) - eps(jN);
                        if abs(d_eps) < obj.xi; continue; end
                        Hbar(iN,jN)=0;
                        Hbar(jN,iN)=0;
                    end
                end
                %% Calculate the FAE eigenstates
                [Psi2,Ebar] = eig(Hbar,'vector');
                Psi = Psi * Psi2;
                eps = diag(Psi' * thf * Psi);
            end
            %% Adiabatically continue
            if isfield(Args,'Psi_prev')
                % Get previous values
                Psi0_prev = obj.Psi0(Args.Psi_prev);
                Psi0 = obj.Psi0(Psi);
                % Overlap with previous
                S0 = Psi0_prev' * Psi0;
                % Get the closest indeces with respect to the previous WF (row)
                [~,ind] = max(abs(S0),[],1);
                Psi = Psi(:,ind);
                eps = eps(ind);
                if nargout > 2
                    Ebar = Ebar(ind);
                end
                % Shift to adiabatically continue the QE
                if isfield(Args,'eps_prev')
                    dk_eps = round((eps - eps_prev - (mod(eps-eps_prev+obj.w/2,obj.w) - obj.w/2))/obj.w);
                    for iN = 1:obj.N
                        Psi(:,iN)=circshift(Psi(:,iN), obj.N * dk_eps(iN), 1);
                        eps(iN) = Psi(:,iN)' * thf * Psi(:,iN);
                    end
                end
            end
        end
        function psi = shiftCenter(obj,psi, Args)
            arguments
                obj
                psi
                Args.tol    = 1E-6
            end
            % SHIFTCENTER Center the wavefunction in Fourier space
        
            for iM = 1:size(psi,2)
                tpsi = reshape(psi(:,iM),obj.N,[]);
                tnorm = vecnorm(tpsi,2,1);
                [peak,tk] = findpeaks(smoothdata(tnorm),...
                    "NPeaks",1, "SortStr","descend");
                tpsi = circshift(tpsi, obj.k_max+1 - tk, 2);
                % Normalize to initial WF
                tpsi0 = sum(tpsi,2);
                if tpsi0(1)<0; tpsi=-tpsi; end
                psi(:,iM) = tpsi(:);
            end
            iM = 2;
            S = psi' * psi;
            while iM <= size(psi,2)
                % Flag if found duplicate
                flag = false;
                for iM2 = 1:iM-1
                    % Check if close to duplicate
                    if (1 - abs(S(iM2,iM))) < Args.tol
                        % Remove duplicate
                        flag = true;
                        psi(:,iM) = [];
                        S(:,iM) = [];
                        S(iM,:) = [];
                        break;
                    end
                end
                if ~flag; iM = iM + 1; end
            end
        end
        function json=jsonencode(obj,varargin)
            j = jsonencode@Calc.baseCalc(obj);
            S = jsondecode(j);
            S.matrix = struct('N',obj.N);
            S.floquet = struct('w',obj.w,'k_max',obj.k_max,'hk_max',obj.hk_max,'xi',obj.xi);
            json = jsonencode(S,varargin{:});
        end
        function Res=variational(obj,method,Args,Flags)
            arguments
                obj
                method
                Args.opts
                Args.Psi0
                Args.TypX
		        Args.tol	(1,1)	double
		        Args.filter_nconv	logical
                Args.functions
                Flags.TrackPsi      logical
            end
            switch method
                case 'varQE'
                    Args = namedargs2cell(Args);
                    Flags = namedargs2cell(Flags);
                    Res=obj.variational_varQE(Args{:},Flags{:});
                case 'AE'
                    Args = namedargs2cell(Args);
                    Flags = namedargs2cell(Flags);
                    Res=obj.variational_AE(Args{:},Flags{:});
                otherwise
                    error('Not implemented');
            end
        end
        function Res=variational_varQE(obj,Args,Flags)
	        arguments
                obj
                Args.opts
                Args.ms
                Args.Psi0
                Args.TypX
		        Args.tol	(1,1)	double	=1E-10
		        Args.filter_nconv	logical	=true
                Args.functions
                Flags.TrackPsi       logical
            end
            tol = Args.tol;
            filter_nconv = Args.filter_nconv;
            %% Alter X0 variable
            if ~isfield(Args,'Psi0')
                Args.Psi0 = zeros(obj.N * obj.k_max2, 1);
                Args.Psi0(obj.k_max * obj.N + 1) = 1;
            end
            switch size(Args.Psi0,1)
                case obj.N
                    NPsi = size(Args.Psi0,2);
                    tPsi0 = zeros(obj.N * obj.k_max2, NPsi);
                    tPsi0(obj.k_max * obj.N + (1:obj.N),:) = Args.Psi0;
                    Args.Psi0 = tPsi0;
                case obj.N * (obj.k_max2-2*obj.hk_max)
                    % Append zeros for consistent manipulation
                    NPsi = size(Args.Psi0,2);
                    Args.Psi0 = [zeros(obj.N * obj.hk_max,NPsi);
                        Args.Psi0;
                        zeros(obj.N * obj.hk_max,NPsi)];
                otherwise
                    error('Wrong size')
            end
            % flag the endpoints
            Nk_max2 = obj.N * obj.k_max2;
            flag_end=false(Nk_max2,1);
            flag_end(1:obj.N*obj.hk_max)=true;
            flag_end(end+1-(1:obj.N*obj.hk_max))=true;
            %% Specify default options
            if ~isfield(Args,'opts')
	            Args.opts=optimoptions('fmincon',Algorithm='sqp',...
                    SubproblemAlgorithm='cg',...
		            FunctionTolerance=Args.tol,OptimalityTolerance=Args.tol,...
		            StepTolerance=Args.tol,ConstraintTolerance=Args.tol,...
		            MaxIterations=30000,MaxFunctionEvaluations=300000,...
		            FiniteDifferenceType='central');
            end
            if ~isfield(Args,'ms')
	            Args.ms=MultiStart(UseParallel=false,...
                    Display="iter",...
		            FunctionTolerance=Args.tol);
            end
            %% Run minimization
            if obj.cacheMat
                thf = obj.hf;
                thf2 = thf * thf;
                Args = namedargs2cell(Args);
                Flags = namedargs2cell(Flags);
                Res = obj.variational_base(@varHf_cached,@cons,Args{:},Flags{:});
            else
                Args = namedargs2cell(Args);
                Flags = namedargs2cell(Flags);
                Res = obj.variational_base(@varHf,@cons,Args{:},Flags{:});
            end
            %% Post-process results
            for iN = 1:length(Res)
                ttPsi = zeros(Nk_max2,1);
                ttPsi(~flag_end) = Res(iN).Psi;
                Res(iN).Psi = ttPsi;
                for iStep = 1:length(Res(iN).steps)
                    if ~isempty(Res(iN).steps(iStep).Psi)
                        ttPsi = zeros(Nk_max2,1);
                        ttPsi(~flag_end) = Res(iN).steps(iStep).Psi;
                        Res(iN).steps(iStep).Psi = ttPsi;
                    end
                end
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
            if filter_nconv
                Res=Res([Res(:).conv]);
            end
	        %%
	        function vHf=varHf_cached(x)
                tPsi=zeros(Nk_max2,1);
                % Force endpoints to be 0
		        tPsi(~flag_end)=x;
		        tnorm2=tPsi' * tPsi;
	            eps=tPsi' * thf * tPsi;
	            vHf=tPsi' * thf2 * tPsi + eps * eps * (tnorm2-2);
		        vHf=vHf/tnorm2;
	        end
	        function vHf=varHf(x)
                tPsi=zeros(Nk_max2,1);
                % Force endpoints to be 0
		        tPsi(~flag_end)=x;
		        tnorm2=tPsi' * tPsi;
                tthf = obj.hf;
	            eps=tPsi' * tthf * tPsi;
	            vHf=tPsi' * tthf * tthf * tPsi + eps * eps * (tnorm2-2);
		        vHf=vHf/tnorm2;
	        end
	        %%
	        function [c,ceq]=cons(x)
                tPsi=zeros(Nk_max2,1);
                % Force endpoints to be 0
		        tPsi(~flag_end)=x;
                tnorm=norm(tPsi);
		        ceq=tnorm-1;
                c=[];
            end
        end
        function Res=variational_AE(obj,Args,Flags)
	        arguments
                obj
                Args.opts
                Args.ms
                Args.Psi0
                Args.TypX
		        Args.tol	(1,1)	double	=1E-10
		        Args.filter_nconv	logical	=true
                Args.functions
                Flags.TrackPsi      logical
            end
            %% Share variables
            txi = obj.xi;
            filter_nconv = Args.filter_nconv;
            %% Alter X0 variable
            if ~isfield(Args,'Psi0')
                tPsi0 = zeros(obj.N * obj.k_max2, 1);
                tPsi0(obj.k_max * obj.N + 1) = 1;
                teps = tPsi0' * obj.hf * tPsi0;
                Args.Psi0 = [tPsi0;teps];
            end
            switch size(Args.Psi0,1)
                case obj.N
                    NPsi = size(Args.Psi0,2);
                    teps = zeros(1,NPsi);
                    tPsi0 = zeros(obj.N * obj.k_max2, NPsi);
                    tPsi0(obj.k_max * obj.N + (1:obj.N),:) = Args.Psi0;
                    for iN = 1:NPsi
                        teps(iN) = tPsi0(:,iN)' * obj.hf * tPsi0(:,iN);
                    end
                    Args.Psi0 = [tPsi0;teps];
                case obj.N * obj.k_max2
                    NPsi = size(Args.Psi0,2);
                    tPsi0 = Args.Psi0;
                    for iN = 1:NPsi
                        teps(iN) = tPsi0(:,iN)' * obj.hf * tPsi0(:,iN);
                    end
                    Args.Psi0 = [tPsi0;teps];
                case obj.N * (obj.k_max2-2*obj.hk_max)
                    NPsi = size(Args.Psi0,2);
                    tPsi0 = [zeros(obj.N * obj.hk_max,NPsi);
                        Args.Psi0;
                        zeros(obj.N * obj.hk_max,NPsi)];
                    for iN = 1:NPsi
                        teps(iN) = tPsi0(:,iN)' * obj.hf * tPsi0(:,iN);
                    end
                    Args.Psi0 = [tPsi0;teps];
                otherwise
                    error('Wrong size')
            end
            % flag the endpoints
            Nk_max2 = obj.N * obj.k_max2;
            flag_end=false(Nk_max2,1);
            flag_end(1:obj.N*obj.hk_max)=true;
            flag_end(end+1-(1:obj.N*obj.hk_max))=true;
            %% Specify default options
            if ~isfield(Args,'opts')
	            Args.opts=optimoptions('fmincon',Algorithm='sqp',...
                    SubproblemAlgorithm='cg',...
		            FunctionTolerance=obj.xi,OptimalityTolerance=obj.xi,...
		            StepTolerance=Args.tol,ConstraintTolerance=Args.tol,...
		            MaxIterations=30000,MaxFunctionEvaluations=300000,...
		            FiniteDifferenceType='central');
            end
            if ~isfield(Args,'ms')
	            Args.ms=MultiStart(UseParallel=false,...
                    Display="iter",...
		            FunctionTolerance=obj.xi);
            end
            %% Run minimization
            if obj.cacheMat
                thf = obj.hf;
                th = obj.h;
                Args = namedargs2cell(Args);
                Flags = namedargs2cell(Flags);
                Res = obj.variational_base(@Ebar_cached,@cons_cached,Args{:},Flags{:});
            else
                Args = namedargs2cell(Args);
                Flags = namedargs2cell(Flags);
                Res = obj.variational_base(@Ebar,@cons,Args{:},Flags{:});
            end
            %% Post-process results
            for iN = 1:length(Res)
                Res(iN).eps = Res(iN).Psi(end);
                ttPsi = zeros(Nk_max2,1);
                ttPsi(~flag_end) = Res(iN).Psi(1:end-1);
                Res(iN).Psi = ttPsi;
                for iStep = 1:length(Res(iN).steps)
                    if isempty(Res(iN).steps(iStep).Psi)
                        Res(iN).steps(iStep).eps = nan;
                    else
                        Res(iN).steps(iStep).eps = Res(iN).steps(iStep).Psi(end);
                        ttPsi = zeros(Nk_max2,1);
                        ttPsi(~flag_end) = Res(iN).steps(iStep).Psi(1:end-1);
                        Res(iN).steps(iStep).Psi = ttPsi;
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
	        function Ebar=Ebar_cached(x)
                tPsi=zeros(Nk_max2,1);
                % Force endpoints to be 0
		        tPsi(~flag_end)=x(1:end-1);
		        Ebar=tPsi' * th * tPsi;
	        end
	        function Ebar=Ebar(x)
                tPsi=zeros(Nk_max2,1);
                % Force endpoints to be 0
		        tPsi(~flag_end)=x(1:end-1);
		        Ebar=tPsi' * obj.h * tPsi;
	        end
	        %%
	        function [c,ceq]=cons_cached(x)
                tPsi=zeros(Nk_max2,1);
                % Force endpoints to be 0
		        tPsi(~flag_end)=x(1:end-1);
		        tnorm=norm(tPsi);
		        grad=thf * tPsi - x(end) * tPsi;
		        c=[grad-txi;-grad-txi];
		        ceq=tnorm-1;
	        end
            function [c,ceq]=cons(x)
                tPsi=zeros(Nk_max2,1);
                % Force endpoints to be 0
		        tPsi(~flag_end)=x(1:end-1);
		        tnorm=norm(tPsi);
		        grad=obj.hf*tPsi-x(end)*tPsi;
		        c=[grad-txi;-grad-txi];
		        ceq=tnorm-1;
            end
        end
    end
    methods (Access=protected)
        function Res=variational_base(obj,fval,cons,Args,Flags)
	        arguments
                obj
                fval
                cons
                Args.opts
                Args.ms
                Args.Psi0
                Args.TypX
		        Args.tol	(1,1)	double
		        Args.filter_nconv	logical
                Args.functions
                Flags.TrackPsi       logical =false
            end
            %% Fix vectors
            if isfield(Args,'TypX')
                tTypX = Args.TypX;
            else
	            tTypX=ones(size(Args.Psi0,1),1);
	            tX=1;
	            for k=obj.hk_max+3:obj.k_max-obj.hk_max
		            tX=tX/k/obj.w;
		            tTypX(obj.N*(obj.k_max+k)+(1:obj.N))=tX;
		            tTypX(obj.N*(obj.k_max-k+1)+(-1:-1:-obj.N)+1)=tX;
	            end
            end
            % flag the endpoints
            % Use size of Args.Psi0 for consistency
            flag_end=false(size(Args.Psi0,1),1);
            flag_end(1:obj.N*obj.hk_max)=true;
            flag_end(obj.N*obj.k_max2+1-(1:obj.N*obj.hk_max))=true;
            % Trim the vectors
            Args.Psi0=Args.Psi0(~flag_end,:);
            tTypX=tTypX(~flag_end);
            %% Initialize minimization options
            if ~isfield(Args,'opts')
	            Args.opts=optimoptions('fmincon',Algorithm='sqp',...
                    SubproblemAlgorithm='cg',...
		            FunctionTolerance=obj.tol,OptimalityTolerance=obj.tol,...
		            StepTolerance=Args.tol,ConstraintTolerance=Args.tol,...
		            MaxIterations=30000,MaxFunctionEvaluations=300000,...
		            FiniteDifferenceType='central');
            end
            Args.opts=optimoptions(Args.opts,...
                TypicalX=tTypX);
            if ~isfield(Args,'ms')
	            Args.ms=MultiStart(UseParallel=false,...
                    Display="iter",...
		            FunctionTolerance=Args.tol);
            end
            pool = gcp('nocreate');
            if ~isempty(pool)
                Args.ms = MultiStart(Args.ms,UseParallel=true);
            end
            if Flags.TrackPsi
		        iPsi=1;
		        Args.opts=optimoptions(Args.opts,...
                    OutputFcn=@getSteps);
		        Args.ms=MultiStart(Args.ms,OutputFcn=@getMSStep, ...
                    UseParallel=false);
            end
            %% Prepare output
            NPsi0 = size(Args.Psi0,2);
            Res(NPsi0) = struct(Psi=[],steps=[],...
                conv=false,Fval=nan,optim=nan);
            for iN = 1:NPsi0
                Res(iN).Psi = [];
                Res(iN).conv = false;
                Res(iN).Fval = nan;
                Res(iN).optim = nan;
                if Flags.TrackPsi
                    Res(iN).steps(Args.opts.MaxIterations).Psi=[];
                end
            end
            %% Prepare problem
	        prob=createOptimProblem('fmincon',...
		        objective=fval,nonlcon=cons,...
		        options=Args.opts,x0=Args.Psi0(:,1));
            %% Run minimization
            if NPsi0 > 1
	            sPoints=CustomStartPointSet(Args.Psi0');
	            [~,~,flag,output,solutions]=run(Args.ms,prob,sPoints);
		        for sol=solutions
			        for solX0=sol.X0
                        ind = ismembertol(Args.Psi0',solX0{:}',Args.tol,...
					        ByRows=true);
                        Res(ind).Psi = sol.X;
                        Res(ind).Fval = sol.Fval;
                        Res(ind).conv = true;
                        Res(ind).optim = sol.Output.firstorderopt;
			        end
		        end
                if Args.filter_nconv
                    Res=Res([Res(:).conv]);
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
            if Flags.TrackPsi
                for iN = 1:length(Res)
                    for iStep = 1:Args.opts.MaxIterations
                        if isempty(Res(iN).steps(iStep).Psi)
                            iStep = iStep - 1;
                            break;
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
                if Flags.TrackPsi
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
    end
end