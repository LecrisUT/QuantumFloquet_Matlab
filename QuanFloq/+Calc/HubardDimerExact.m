classdef HubardDimerExact < Calc.baseFloquet
    % HubardDimerExact Driven Hubard dimer exactly solvable system
    % Exaclty solvable driven Hamiltonian
    % $$H(t)=v(t)(n_{1\sigma}-n_{2\sigma})-t(c^\dagger_{1\sigma}c_{2\sigma}+\text{h.c.})+U(n_{1\uparrow}n_{1\downarrow}+n_{2\uparrow}n_{2\downarrow})$$
    % $$v(t)=v_0+v_1\cos(\omega t)$$
    % The system is solvable in the undriven basis
    % <latex>
    % \begin{align}
    % \ket{\varphi_1}&=\ket{1\uparrow1\downarrow}\\
    % \ket{\varphi_2}&=\frac{1}{\sqrt{2}}(\ket{1\uparrow2\downarrow}+\ket{2\uparrow1\downarrow})\\
    % \ket{\varphi_3}&=\ket{2\uparrow2\downarrow}\\
    % H(t)=\mqty[U+v(t)&-\sqrt(2)t&0\\-\sqrt(2)t&0&-\sqrt(2)t\\0&-\sqrt(2)t&U-v(t)]
    % \end{align}
    % </latex>
    %
    % See also HubardDimerOrbital

    properties
        % t - Hopping parameter
        t       (1,1)   double  {mustBeReal}
        % v0 - Static interaction strength
        v0      (1,1)   double  {mustBeReal}
        % v1 - Driving strength
        v1      (1,1)   double  {mustBeReal}
        % U - Cuolomb interaction
        U       (1,1)   double  {mustBeReal}
    end
    properties (Dependent,SetAccess=private)
        ht
    end
    methods
        function val = get.ht(obj)
            val = zeros(3,3,3);
            sqt = sqrt(2) * obj.t;
            val(:,:,2) = [obj.U+obj.v0 -sqt 0;
                -sqt 0 -sqt;
                0 -sqt obj.U-obj.v0];
            val(:,:,3) = [obj.v1/2 0 0;0 0 0;0 0 -obj.v1/2];
            val(:,:,1) = val(:,:,3)';
        end
    end
    methods
        function obj = HubardDimerExact(Args,Args2)
            arguments
                Args.U      = 0
                Args.t      = 1
                Args.v0     = 2
                Args.v1     = 5
                Args2.w     = 1.5
                Args2.k_max = 100
                Args2.xi    = 1E-3
            end
            % HubardDimerExact Constructor
            %
            % Syntax:
            %   obj = HubardDimerExact(Name,Value)
            % 
            % Description:
            %   obj = HubardDimerExact(Name,Value) Specify the parameters of the
            %   two-level system via name-value pairs
            %
            % Inputs:
            %   Name-Value pairs
            %
            % Outputs:
            %   obj - Floquet object
            %
            % Name-value arguments:
            %   U - [0] Cuolomb interaction
            %   t - [1] Hopping parameter
            %   v0 - [2] Static interaction
            %   v1 - [5] Driving strength
            %   w - [1.5] Driving frequency
            %   k_max - [100] Fourier cut-off of Fourier coefficients
            %   xi - [1E-3] Acceptable error
            % 
            % See also HubardDimerExact, baseFloquet

            Args2.hk_max = 1;
            Args2 = namedargs2cell(Args2);
            obj@Calc.baseFloquet(3,Args2{:});
            obj.t = Args.t;
            obj.v0 = Args.v0;
            obj.v1 = Args.v1;
        end
    end
    methods
        function set.t(obj,val)
            obj.t = val;
            obj.dirty_cache = true;
            obj.dirty_cacheU = true;
        end
        function set.v0(obj,val)
            obj.v0 = val;
            obj.dirty_cache = true;
            obj.dirty_cacheU = true;
        end
        function set.v1(obj,val)
            obj.v1 = val;
            obj.dirty_cache = true;
            obj.dirty_cacheU = true;
        end
        function set.U(obj,val)
            obj.U = val;
            obj.dirty_cacheU = true;
        end
    end
    methods
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
    methods
        function json = jsonencode(obj,varargin)
            j = jsonencode@Calc.baseFloquet(obj);
            S = jsondecode(j);
            S.hubbard_dimer = struct(t=obj.t,v0=obj.v0,v1=obj.v1,U=obj.U);
            json = jsonencode(S,varargin{:});
        end
    end
    methods (Access=protected)
        function groups = getPropertyGroups(obj)
            import matlab.mixin.util.PropertyGroup
            groups = getPropertyGroups@Calc.baseFloquet(obj);
            if isscalar(obj)
                Dimer = struct(t=obj.t,v0=obj.v0,v1=obj.v1,U=obj.U);
                groups = [PropertyGroup(Dimer,'Hubard dimer system properties:'),...
                    groups];
            end
        end
    end
end