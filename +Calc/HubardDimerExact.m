classdef HubardDimerExact < Calc.baseFloquet
    properties
        t
        v0
        v1
        U       = 0
    end
    properties (Dependent)
        hU
        h0
    end
    methods (Access=protected)
        function h = get_h(obj)
            h = obj.h0 + obj.hU;
        end
        function ht = get_ht(obj)
            error('Not implemented');
        end
    end
    methods
        function obj=HubardDimerExact(Args,Args2)
            arguments
                Args.t      = 1
                Args.v0     = 2
                Args.v1     = 5
                Args2.w      = 1.5
                Args2.k_max  = 100
                Args2.xi     = 1E-3
            end
            Args2.hk_max = 1;
            Args2 = namedargs2cell(Args2);
            obj@Calc.baseFloquet(3,Args2{:});
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