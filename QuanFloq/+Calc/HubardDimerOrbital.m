classdef HubardDimerOrbital < Calc.baseFloquetHF
    % HubardDimerOrbital Driven Hubard dimer single Slater approximation
    % See Calc.HubardDimerExact for the Hamiltonian.
    % The Slater determinant is limited to the form
    % $$\Phi(r_1,r_2,t)=\frac{1}{\sqrt{2}}(\phi_{\uparrow}(r_1,t)\phi_{\downarrow}(r_2,t))-\phi_{\downarrow}(r_1,t)\phi_{\uparrow}(r_2,t))$$
    % $$\phi_\sigma(r,t)=\phi_1(t)\varphi_1(r)+\phi_2(t)\varphi_2(r)$$
    % 
    %
    % See also HubardDimerExact

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
        hUt
    end
    methods
        function val = get.ht(obj)
            val = zeros(2,2,3);
            val(:,:,2) = [obj.v0/2 -obj.t;-obj.t -obj.v0/2];
            val(:,:,3) = [obj.v1/4 0;0 -obj.v1/4];
            val(:,:,1) = val(:,:,3)';
        end
        function val = get.hUt(obj)
            Psi = obj.Psi_Fourier(obj.Psi);
            val = zeros(2,2,obj.hUk_max2);
            for k = obj.hUk_range
                if k > 0
                    tU = conj(Psi(:,1,1:end-k)) .* Psi(:,1,1+k:end);
                else
                    tU = conj(Psi(:,1,1-k:end)) .* Psi(:,1,1:end+k);
                end
                tU = obj.U * diag(sum(tU,3));
                val(:,:,obj.hUk_max+1+k) = tU;
            end
        end
    end
    methods
        function obj = HubardDimerOrbital(Args,Args2)
            arguments
                Args.U          = 0
                Args.t          = 1
                Args.v0         = 2
                Args.v1         = 5
                Args2.w         = 1.5
                Args2.k_max     = 100
                Args2.hUk_max   = 30
                Args2.xi        = 1E-3
            end
            % HubardDimerOrbital Constructor
            %
            % Syntax:
            %   obj = HubardDimerOrbital(Name,Value)
            % 
            % Description:
            %   obj = HubardDimerOrbital(Name,Value) Specify the parameters of the
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
            % See also HubardDimerOrbital, baseFloquetHF

            Args2.hk_max = 1;
            Args2 = namedargs2cell(Args2);
            obj@Calc.baseFloquetHF(2,Args2{:},Ne=2,mode='closed-shell');
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
%         function newobj = Calc.HubardDimerExact(obj)
% %             newobj = Calc.HubardDimerExact('v0');
%         end
        function Psi_Slat = OrbitalToSlater(obj,Psi_Orb,Args)
            arguments
                obj         Calc.HubardDimerOrbital
                Psi_Orb     double
                Args.objSlater          (1,1)   Calc.HubardDimerExact
                Args.normalize          (1,1)   logical = true
                Args.normalizeSlater    (1,1)   logical = false
            end
            % OrbitalToSlater Calculate the corresponding Slater determinant
            % Slater determinant basis corresponds to that in Calc.HubarDimerExact
            %
            % See usage of Calc.baseFloquetHF.OrbitalToSlater
            %
            % See also baseFloquetHF.OrbitalToSlater, HubarDimerExact

            if Args.normalize
                Psi_Orb = obj.Psi_Floquet(Psi_Orb);
                Psi_Orb = Psi_Orb ./ vecnorm(Psi_Orb);
            end
            Psi_Orb = obj.Psi_Fourier(Psi_Orb);
            M = size(Psi_Orb,2);
            % Normally should check if orbital dimensions are acceptable, but the
            % single orbital is a special case
            if isfield(Args,'objSlater')
                Psi_Slat = zeros(Args.objSlater.N,M,Args.objSlater.k_max2);
                Slat_k_max = Args.objSlater.k_max;
            else
                Psi_Slat = zeros(3,M,obj.k_max2);
                Slat_k_max = obj.k_max;
            end
            for k = -obj.k_max:obj.k_max
                if k>0
                    % psi(k+l)*psi(l) with l in [-k_max:k_max]
                    tPsi = [Psi_Orb(1,:,1+k:obj.k_max2) .* Psi_Orb(1,:,obj.k_max2:-1:1+k);
                        Psi_Orb(1,:,1+k:obj.k_max2) .* Psi_Orb(2,:,obj.k_max2:-1:1+k);
                        Psi_Orb(2,:,1+k:obj.k_max2) .* Psi_Orb(2,:,obj.k_max2:-1:1+k)];
                else
                    tPsi = [Psi_Orb(1,:,1:obj.k_max2+k) .* Psi_Orb(1,:,obj.k_max2+k:-1:1);
                        Psi_Orb(1,:,1:obj.k_max2+k) .* Psi_Orb(2,:,obj.k_max2+k:-1:1);
                        Psi_Orb(2,:,1:obj.k_max2+k) .* Psi_Orb(2,:,obj.k_max2+k:-1:1)];
                end
                tPsi = sum(tPsi,3);
                tPsi(2,:) = tPsi(2,:) * sqrt(2);
                Psi_Slat(:,:,Slat_k_max+1+k) = tPsi;
            end
            % Transform to Floquet representation
            Psi_Slat = reshape(permute(Psi_Slat,[1 3 2]),[],M);
            if Args.normalizeSlater
                Psi_Slat = Psi_Slat ./ vecnorm(Psi_Slat);
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
            groups = getPropertyGroups@Calc.baseFloquetHF(obj);
            if isscalar(obj)
                Dimer = struct(t=obj.t,v0=obj.v0,v1=obj.v1,U=obj.U);
                groups = [PropertyGroup(Dimer,'Hubard dimer system properties:'),...
                    groups];
            end
        end
    end
end