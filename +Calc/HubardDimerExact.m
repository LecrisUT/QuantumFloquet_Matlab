classdef HubardDimerExact < Calc.baseFloquet
    properties
        t
        v0
        v1
        U       = 0
    end
    properties (Dependent)
        h
        hU
        h0
    end
    methods
        function obj=HubardDimerExact(Args)
            arguments
                Args.t      = 1
                Args.w      = 1.5
                Args.k_max  = 100
                Args.v0     = 2
                Args.v1     = 5
                Args.xi     = 1E-3
            end
            obj@Calc.baseFloquet(3,Args.w,...
                k_max=Args.k_max,hk_max=1,xi=Args.xi);
            obj.t = Args.t;
            obj.v0 = Args.v0;
            obj.v1 = Args.v1;
        end
        function json = jsonencode(obj,varargin)
            j = jsonencode@Calc.baseFloquet(obj);
            S = jsondecode(j);
            S.hubbard_dimer = struct('t',obj.t,'v0',obj.v0,'v1',obj.v1,...
                'xi',obj.xi);
            json = jsonencode(S,varargin{:});
        end
        function h = get_h(obj)
            h = obj.h0 + obj.hU;
        end
        function ht = get_ht(obj)
            error('Not implemented');
        end
        function h0 = get.h0(obj)
            h0 = [obj.v1/2  -sqrt(2)*obj.t   obj.v0   0 obj.v1/2;
                 0      -sqrt(2)*obj.t   0   -sqrt(2)*obj.t   0;
                -obj.v1/2    0 -obj.v0  -sqrt(2)*obj.t  -obj.v1/2];
            h0 = repmat(h0,obj.k_max2,1);
            h0 = spdiags(h0,[-3 -1:1 3],obj.N*obj.k_max2,obj.N*obj.k_max2);
        end
        function hU = get.hU(obj)
            hU = [obj.U;0;obj.U];
            hU = repmat(hU,obj.k_max2,1);
            hU = spdiags(hU, 0, obj.N*obj.k_max2, obj.N*obj.k_max2);
        end
        function Res=eigs(obj,U_range,Args)
            arguments
                obj
                U_range
                Args.saveFile
                Args.loadFile
                Args.xi     = 1E-4
                Args.tol    = 1E-12
                Args.th_range
                Args.Psi0   = [1;0]
                Args.Print  = true
            end

            %% Calculate the non-SCF solutions
            nU = length(U_range);
            if isfield(Args,'loadFile')
                iU = load(Args.loadFile,'iU');
                Res = load(Args.loadFile,'Res');
            else
                Res(nU,obj.N) = struct('Psi',[],'eps',[],'Ebar',[],'th',[],'S',[]);
                iU = 1;
            end
            th=zeros(2,3);
            S=zeros(2,3);
            for iU = iU:nU
                obj.U = U_range(iU);
                %% Calculate the AE eigenstates
                % Add adiabatic arguments if possible
                subArgs={};
                if iU > 1
                    Psi_prev = reshape([Res(iU-1,:).Psi],[],obj.N);
                    eps_prev = [Res(iU-1,:).eps];
                    subArgs = [subArgs(:)',{'Psi_prev'},{Psi_prev},{'eps_prev'},{eps_prev}];
                end
                [eps,Psi,Ebar] = eigs@Calc.baseFloquet('eps0',0,subArgs{:});
                Psi0 = obj.Psi0(Psi);
                % Fix the phase
                for iN=1:obj.N
                    if Psi0(1,iN)<0
                        Psi(:,iN) = -Psi(:,iN);
                        Psi0(:,iN) = -Psi0(:,iN);
                    end
                end
                %% Find the closest slater determinant
                for iN = 1:N
                    [th(:,iN),S(:,iN)] = obj.SlaterDecomp(Psi0);
                    if iU == 1
                        [S(:,iN),ind] = sort(S(:,iN),'descend');
                        th(:,iN) = th(ind,iN);
                    else
                        th_prev = Res(iU-1,iN).th;
                        dth = abs(mod(th(:,iN) - th_prev(1) + pi/2,pi) - pi/2);
                        [~,ind]=sort(dth,'ascend');
                        S(:,iN) = S(ind,iN);
                        th(:,iN) = th(ind,iN);
                    end
                end

                %% Save the results
                for iN = 1:obj.N
                    Res(iU,iN).Psi = Psi(:,iN);
                    Res(iU,iN).eps = eps(iN);
                    Res(iU,iN).Ebar = Ebar(iN);
                    Res(iU,iN).th = th(:,iN);
                    Res(iU,iN).S = S(:,iN);
                end
                if isfield(Args,'saveFile')
                    save(Args.saveFile,'Res','iU','obj','-v7.3');
                end
                if Args.Print
	                fprintf('Done [%d/%d]\n',iU,nU);
                end
            end
        end
        function [th,S]=SlaterDecomp(~,Psi0)
	        opts=optimoptions('fmincon','Algorithm','interior-point');
	        ms=MultiStart('UseParallel',false,'Display','off',...
                'XTolerance',1E-4);
	        prob=createOptimProblem('fmincon',...
		        'objective',@objS,...
		        'options',opts,'x0',0, ...
                'lb',-pi/2,'ub',pi/2);
            pool = gcp('nocreate');
            if ~isempty(pool)
                ms = MultiStart(ms,"UseParallel",true);
            end
	        [~,~,flag,output,solutions]=run(ms,prob,41);
            th = [solutions.X];
            S = abs(overlap(th));
            fltr = true(1,length(th));
            fltr(th < -pi/2 + 1E-4) = false;
            fltr(th > pi/2 - 1E-4) = false;
            th = th(fltr);
            S = S(fltr);
            if length(th) == 1
                th(1,2) = mod(pi+th,pi)-pi/2;
                S(1,2) = real(sqrt(1-S^2));
            end
            if length(th) ~= 2
                error('Too many solutions');
            end
            function Fval = objS(th)
                Fval = 1 - abs(overlap(th));
            end
            function S = overlap(th)
                S = Psi0' * OrbToSlat([cos(th);sin(th)]);
            end
        end
    end
end