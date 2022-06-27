classdef BoxParticle < Calc.baseFloquet
    properties
        V   = 0
    end
    methods (Access=protected)
        function ht = get_ht(obj)
            ht = zeros(obj.N,obj.N,3);
            ht(:,:,2) = diag((1:obj.N).^2);
            th1 = ones(obj.N-1,1);
            th1(2:2:end) = -1;
            th1 = obj.V/4 * th1;
            ht(:,:,3) = diag(th1,1)+diag(th1,-1);
            ht(:,:,1) = ht(:,:,3)';
        end
        function h = get_h(obj)
            % Build Hamiltonian primitives
            % Column1: Non-interacting diagonal E_n=n^2*E_0
            % Column2: Driving of form [1 -1 1 -1 ...]
            h = [(1:obj.N).^2' ones(obj.N,1)];
            h(2:2:end,2) = -1;
            % Build the Hamiltonian diagonals
            h = [obj.V/4*[h(1:end-1,2);0],...
                obj.V/4*[0;h(1:end-1,2)],...
                h(:,1),...
                obj.V/4*[h(1:end-1,2);0],...
                obj.V/4*[0;h(1:end-1,2)]];
            % Extend the Hamiltonian diagonals
            h = repmat(h,obj.k_max2,1);
            h = spdiags(h,[-obj.N-1 -obj.N+1 0 obj.N-1 obj.N+1],...
                obj.N*obj.k_max2,obj.N*obj.k_max2);
        end
    end
    methods
        function obj = BoxParticle(N,Args)
            arguments
                N
                Args.w      = 8.3
                Args.k_max  = 200
                Args.xi     = 1E-2
            end
            obj@Calc.baseFloquet(N,Args.w,...
                k_max=Args.k_max,hk_max=1,xi=Args.xi);
        end
        function json = jsonencode(obj,varargin)
            j = jsonencode@Calc.baseFloquet(obj);
            S = jsondecode(j);
            S.box_particle = struct('xi',obj.xi,'V',obj.V);
            json = jsonencode(S,varargin{:});
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
                try
                    [tPsi,teps] = eigs(hf,obj.N,0);
                catch
                    [tPsi,teps] = eigs(hf,obj.N,1E-12);
                end
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
	                tEbar = tEbar(perm);
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
                try
                    [tPsi,teps] = eigs(obj.hf,obj.N,0);
                catch
                    [tPsi,teps] = eigs(obj.hf,obj.N,1E-12);
                end
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