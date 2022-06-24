classdef HubardDimerOrbital2 < handle
    properties
        t
        v0
        U       = 0
        psi     = [1;0]
    end
    properties (Dependent)
        h
        hU
        h0
        E
    end
    properties (SetAccess=private)
        N       = 2
    end
    methods
        function obj=HubardDimerOrbital2(Args)
            arguments
                Args.t      = 1
                Args.v0     = 2
            end
            obj.t = Args.t;
            obj.v0 = Args.v0;
        end
        function json = jsonencode(obj,varargin)
            S.matrix = struct('N',2);
            S.hubbard_dimer2 = struct('t',obj.t,'v0',obj.v0,'v1',obj.v1,...
                'xi',obj.xi,'hk_max',obj.hk_max);
            json = jsonencode(S,varargin{:});
        end
        function h = get.h(obj)
            h = obj.h0 + obj.hU;
        end
        function hU = get.hU(obj)
            hU = obj.U * [abs(obj.psi(1)*obj.psi(1)) 0;
                0 abs(obj.psi(2)*obj.psi(2))];
        end
        function h0 = get.h0(obj)
            h0 = [obj.v0/2 -obj.t;
                -obj.t -obj.v0/2];
        end
        function E = get.E(obj)
            E = obj.psi' * obj.h * obj.psi + obj.psi' * obj.h0 * obj.psi;
        end
        function Res=eigs(obj,U_range,Args)
            arguments
                obj
                U_range
                Args.saveFile
                Args.loadFile
                Args.tol    = 1E-12
                Args.th_range
                Args.Psi0   = [1;0]
                Args.Print  = true
            end
            %% Initialize optimizaton options
	        opts=optimoptions('fmincon','Algorithm','interior-point',...
		        'FunctionTolerance',Args.tol,'OptimalityTolerance',Args.tol,...
		        'StepTolerance',1E-2 * Args.tol,'ConstraintTolerance',Args.tol,...
		        'MaxIterations',30000,'MaxFunctionEvaluations',300000);
	        prob=createOptimProblem('fmincon',...
		        'objective',@EHF,'nonlcon',@cons,...
		        'options',opts,'x0',Args.Psi0(:,1));
            if isfield(Args,'th_range')
                Args.Psi0 = [cos(th_range);sin(th_range)];
            end
            if size(Args.Psi0,2) > 1
	            ms=MultiStart('UseParallel',false,...
		            'FunctionTolerance',2 * Args.tol, ...
                    'Display','off', ...
                    'XTolerance', 2 * Args.tol);
	            sPoints=CustomStartPointSet(Args.Psi0');
                if ~isempty(gcp('nocreate'))
                    ms = MultiStart(ms,"UseParallel",true);
                end
            end
            %% Run optimization
            nU = length(U_range);
            if isfield(Args,'loadFile')
                iU = load(Args.loadFile,'iU');
                Res = load(Args.loadFile,'Res');
            else
                Res(nU,obj.N) = struct('Psi',[],'eps',[],'E',[]);
                iU = 1;
            end
            for iU = iU:nU
                obj.U = U_range(iU);
                if isempty(whos('ms'))
                    obj.Psi = fimincon(prob);
                else
                    obj.Psi = run(ms,prob,sPoints);
                end
                [tPsi,teps] = eig(obj.h,'vector');
                for iN = 1:obj.N
                    Res(iU,iN).Psi=tPsi(:,iN);
                    Res(iU,iN).eps=teps(iN);
                    obj.Psi = tPsi(:,iN);
                    Res(iU,iN).E=obj.E;
                end
                if isfield(Args,'saveFile')
                    save(Args.saveFile,'Res','iU','obj','-v7.3');
                end
                if Args.Print
	                fprintf('Done [%d/%d]\n',iU,nU);
                end
            end
            %% Optimization functions
            function [c,ceq] = cons(tPsi)
                c = [];
                ceq = norm(tPsi) - 1;
            end
            function E = EHF(tPsi)
                obj.Psi = tPsi;
                E = obj.E;
            end
        end
    end
end