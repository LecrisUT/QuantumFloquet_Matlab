classdef TwoLevel < Calc.baseFloquet
    % TwoLevel Driven two level system
    % Type of Hamiltonian is controlled by RWA and type of purturbation by
    % rand_v
    % Hamiltonian with rotating wave approximation (RWA=true)
    % $$H(t)=\mqty[\frac{\omega_0}{2}&\frac{V}{2}e^{-i\omega t}\\\frac{V}{2}e^{i\omega t}&&-\frac{\omega_0}{2}$$
    % Hamiltonian without rotating wave approximation (RWA=false)
    % $$H(t)=\mqty[\frac{\omega_0}{2}&V\cos(\omega t)\\V\cos(\omega t)&&-\frac{\omega_0}{2}$$
    % The perturbation with rand_v=false
    % $$v(t)=\mqty[0&v\\v&0]$$
    % The perturbation with rand_v=true is a random time-periodic Hermitian
    % operator
    %
    % See also baseFloquet

    properties
        % w0 - Undriven transition frequency $\omega_0$
        w0      (1,1)   double  {mustBeReal,mustBeNonnegative}
        % v - Perturbation strength
        % See also Calc.TwoLevel.rand_v
        v       (1,1)   double  {mustBeReal}
        % V - Driving strength
        % See also Calc.TwoLevel.RWA
        V       (1,1)   double  {mustBeReal,mustBeNonnegative}
        % rand_v - Flag to use random perturbation
        % See also Calc.TwoLevel.v
        rand_v  (1,1)   logical
        % RWA - Flag to use rotating wave approximation
        % Equivalent to using circularly polarized driving
        % See also Calc.TwoLevel.rand_v
        RWA     (1,1)   logical
    end
    properties (Dependent,SetAccess=private)
        ht
    end
    methods
        function val = get.ht(obj)
            if obj.rand_v
                val = zeros(2,2,obj.hk_max2);
                val(:,:,obj.hk_max+1) = [obj.w0/2 0;0 -obj.w0/2];
                if obj.RWA
                    val(:,:,obj.hk_max+1+1) = [0 obj.V/2;0 0];
                else
                    val(:,:,obj.hk_max+1+1) = [0 obj.V/2;obj.V/2 0];
                    val(:,:,obj.hk_max+1-1) = val(:,:,obj.hk_max+1+1)';
                end
                vt = obj.v/2 * rand(2,2,obj.hk_max2);
                vt = vt + flip(pagectranspose(vt),3);
                val = val + vt;
            else
                val = zeros(2,2,3);
                val(:,:,2) = [obj.w0/2 obj.v;obj.v -obj.w0/2];
                if obj.RWA
                    val(:,:,3) = [0 obj.V/2; 0 0];
                    val(:,:,1) = val(:,:,3)';
                else
                    val(:,:,3) = [0 obj.V/2; obj.V/2 0];
                    val(:,:,1) = val(:,:,3)';
                end
            end
        end
    end
    methods
        function obj = TwoLevel(Args,Args2)
            arguments
                Args.RWA    = true
                Args.rand_v = true
                Args.w0     = 1
                Args.V      = 0
                Args.v      = 1E-4
                Args2.w     = 1.5
                Args2.k_max = 100
                Args2.xi    = 1E-2
            end
            % TwoLevel Constructor
            %
            % Syntax:
            %   obj = TwoLevel(Name,Value)
            % 
            % Description:
            %   obj = TwoLevel(Name,Value) Specify the parameters of the two-level
            %   system via name-value pairs
            %
            % Inputs:
            %   Name-Value pairs
            %
            % Outputs:
            %   obj - Floquet object
            %
            % Name-value arguments:
            %   RWA - [true] Whether to use Rotating Wave Approximation
            %   rand_v - [true] Whether to use random perturbation
            %   w0 - [1] Undriven transition frequency
            %   V - [0] Driving strength
            %   v - [0] Pertiurbation strength
            %   w - [1.5] Driving frequency
            %   k_max - [200] Fourier cut-off of Fourier coefficients
            %   xi - [1E-2] Acceptable error
            % 
            % See also TwoLevel, baseFloquet
            
            if Args.v == 0
                Args.rand_v = false;
            end
            if Args.rand_v
                Args2.hk_max = 10;
            else
                Args2.hk_max = 1;
            end
            Args2 = namedargs2cell(Args2);
            obj@Calc.baseFloquet(2,Args2{:});
            obj.w0 = Args.w0;
            obj.v = Args.v;
            obj.rand_v = Args.rand_v;
            obj.RWA = Args.RWA;
        end
    end
    methods
        function set.w0(obj,val)
            obj.w0 = val;
            obj.dirty_cache = true;
        end
        function set.v(obj,val)
            obj.v = val;
            obj.dirty_cache = true;
        end
        function set.V(obj,val)
            obj.V = val;
            obj.dirty_cache = true;
        end
        function set.rand_v(obj,val)
            obj.rand_v = val;
            obj.dirty_cache = true;
        end
        function set.RWA(obj,val)
            obj.RWA = val;
            obj.dirty_cache = true;
        end
    end
    methods
        function json = jsonencode(obj,varargin)
            j = jsonencode@Calc.baseFloquet(obj);
            S = jsondecode(j);
            S.two_level = struct(RWA=obj.RWA,rand_v=obj.rand_v,w0=obj.w0,V=obj.V,v=obj.v);
            json = jsonencode(S,varargin{:});
        end
    end
    methods (Access=protected)
        function groups = getPropertyGroups(obj)
            import matlab.mixin.util.PropertyGroup
            groups = getPropertyGroups@Calc.baseFloquet(obj);
            if isscalar(obj)

                TwoLevel = struct(RWA=obj.RWA,rand_v=obj.rand_v,w0=obj.w0,V=obj.V,v=obj.v);
                groups = [PropertyGroup(TwoLevel,'Two-level system properties:'),...
                    groups];
            end
        end
    end
end