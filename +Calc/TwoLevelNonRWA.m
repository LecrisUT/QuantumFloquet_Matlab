classdef TwoLevelNonRWA < Calc.baseFloquet
    properties
        w0
        v
        V   = 0
    end
    methods
        function obj = TwoLevelNonRWA(Args)
            arguments
                Args.w0     = 1
                Args.w      = 1.5
                Args.k_max  = 100
                Args.v      = 1E-6
                Args.xi     = 1E-3
            end
            obj@Calc.baseFloquet(2,Args.w,...
                k_max=Args.k_max,hk_max=10,xi=Args.xi);
            obj.w0 = Args.w0;
            obj.v = Args.v;
        end
        function json = jsonencode(obj,varargin)
            j = jsonencode@Calc.baseFloquet(obj);
            S = jsondecode(j);
            S.two_level = struct('w0',obj.w0,'v',obj.v,'xi',obj.xi,'V',obj.V);
            json = jsonencode(S,varargin{:});
        end
        function h = get_h(obj)
            hdiag=[obj.V/2  0  0        obj.w0/2  obj.V/2   0  0;
                   0        0  obj.V/2 -obj.w0/2  0         0  obj.V/2];
	        hdiag = repmat(hdiag,obj.k_max2,1);
	        h = spdiags(hdiag,-3:3,2*obj.k_max2,2*obj.k_max2);
            vdiag = obj.v/2 * rand(2,21);
            vdiag = repmat(vdiag,obj.k_max2,1);
	        vt = spdiags(vdiag,-10:10,2*obj.k_max2,2*obj.k_max2);
            vt = vt + vt';
            h = h + vt;
        end
        function ht = get_ht(obj)
            error('Not implemented');
        end
        function Res=variational(obj,method,Args,Flags,Args2)
            arguments
                obj
                method
                Args.Psi0
                Args.TypX
		        Args.tol	(1,1)	double	=1E-10
		        Args.filter_nconv	logical
                Args.functions
                Flags.TrackPsi       logical
                Args2.th_range      (1,:)   double
            end
            if isfield(Args2,'th_range')
                Args.Psi0=[cos(Args2.th_range);sin(Args2.th_range)];
            end
            Args = namedargs2cell(Args);
            Flags = namedargs2cell(Flags);
            Res = variational@Calc.baseFloquet(obj,method,Args{:},Flags{:});
        end
        function Res = AEEigs(obj,V_range,Args)
            arguments
                obj
                V_range
                Args.Print  = true
            end
            nV = length(V_range);
            Res(nV,obj.N) = struct('Psi',[],'eps',[],'Ebar',[]);
            for iV=1:nV
	            obj.V=V_range(iV);
                hf = obj.hf;
                [tPsi,teps] = eigs(hf,obj.N,0);
                teps = diag(teps);
                Hbar = obj.HBar(tPsi,teps);
                [ttPsi,tEbar] = eig(Hbar,'vector');
                tPsi = tPsi * ttPsi;
                tPsi0 = obj.Psi0(tPsi);
                teps = diag(tPsi' * hf * tPsi);
                if iV == 1
                    [tEbar,ind] = sort(tEbar);
                    teps = teps(ind);
                    tPsi = tPsi(:,ind);
                    tPsi0 = tPsi0(:,ind);
                else
                    [mval,perm] = max(abs(tPsi0_prev'*tPsi0),[],1);
                    % Check if overlap is too small to adiabatically
                    % continue
                    if ~isempty(find(abs(1-mval)>1E-1, 1))
                        warning('Overlap is too small');
                    end
                    % Check if missing values
                    if length(perm) ~= length(unique(perm))
                        error('Missing permutation index');
                    end
                    % Reorder to permute properly
                    [~,perm] = sort(perm);
	                tPsi = tPsi(:,perm);
	                teps = teps(perm);
	                tEbar=tEbar(perm);
	                tPsi0 = tPsi0(:,perm);
	                dk = round((teps-teps_prev) / obj.w);
                    for iN = 1:obj.N
		                tdk = dk(iN);
                        tPsi(:,iN) = circshift(tPsi(:,iN),obj.N * tdk,2);
		                teps(iN) = teps(iN) - tdk * obj.w;
                    end
                end
                for iN = 1:obj.N
                    Res(iV,iN).Psi=tPsi(:,iN);
                    Res(iV,iN).eps=teps(iN);
                    Res(iV,iN).Ebar=tEbar(iN);
                end
                teps_prev=teps;
                tPsi0_prev=tPsi0;
                if Args.Print
	                fprintf('Done [%d/%d]\n',iV,length(V_range));
                end
            end
        end
        function Res = QEEigs(obj,V_range,Args)
            arguments
                obj
                V_range
                Args.Print  = true
            end
            nV = length(V_range);
            Res(nV,obj.N) = struct('Psi',[],'eps',[]);
            for iV=1:nV
	            obj.V=V_range(iV);
                [tPsi,teps] = eigs(obj.hf,obj.N,0);
                teps = diag(teps);
                tPsi0 = obj.Psi0(tPsi);
                if iV == 1
                    [teps,ind] = sort(teps);
                    tPsi = tPsi(:,ind);
                    tPsi0 = tPsi0(:,ind);
                else
                    [mval,perm] = max(abs(tPsi0_prev'*tPsi0),[],1);
                    % Check if overlap is too small to adiabatically
                    % continue
                    if ~isempty(find(abs(1-mval)>1E-1, 1))
                        warning('Overlap is too small');
                    end
                    % Check if missing values
                    if length(perm) ~= length(unique(perm))
                        error('Missing permutation index');
                    end
                    % Reorder to permute properly
                    [~,perm] = sort(perm);
	                tPsi = tPsi(:,perm);
	                teps = teps(perm);
	                tPsi0 = tPsi0(:,perm);
	                dk = round((teps-teps_prev) / obj.w);
                    for iN = 1:obj.N
		                tdk = dk(iN);
                        tPsi(:,iN) = circshift(tPsi(:,iN),obj.N * tdk,2);
		                teps(iN) = teps(iN) - tdk * obj.w;
                    end
                end
                for iN = 1:obj.N
                    Res(iV,iN).Psi=tPsi(:,iN);
                    Res(iV,iN).eps=teps(iN);
                end
                teps_prev=teps;
                tPsi0_prev=tPsi0;
                if Args.Print
	                fprintf('Done [%d/%d]\n',iV,length(V_range));
                end
            end
        end
    end
end