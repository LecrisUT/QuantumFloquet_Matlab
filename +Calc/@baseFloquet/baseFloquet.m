classdef (Abstract) baseFloquet < Calc.baseCalc
    % baseFloquet Base class for Floquet calculations
    %
    % See also TwoLevel

    properties (SetAccess=private)
        % N - Size of the base Hamiltonian Hilbert space
        N           (1,1)   double  {mustBeInteger}
        % hk_max - Fourier cutoff of the time-periodic Hamiltonian
        % See also Calc.baseFloquet.hk_max2, Calc.baseFloquet.hk_range
        hk_max      (1,1)   double  {mustBeInteger,mustBeNonnegative}
    end
    properties
        % k_max - Fourier space cutoff
        % See also Calc.baseFloquet.k_max2, Calc.baseFloquet.k_range
        k_max       (1,1)   double  {mustBeInteger,mustBeNonnegative}
        % w - Driving frequency $\omega$
        % See also Calc.baseFloquet.pt, Calc.baseFloquet.hf
        w           (1,1)   double  {mustBeReal,mustBeNonnegative}
        % xi - Acceptable error
        % See also baseFloquet.eigs, baseFloquet.Ebar
        xi          (1,1)   double  {mustBeReal,mustBeNonnegative}
        % systemBath - System bath coupling
        % See also baseSystemBath
        systemBath  Calc.baseSystemBath {mustBeScalarOrEmpty}
    end
    properties (Dependent)
        % k_range - Helper property for [-k_max:k_max]
        % See also Calc.baseFloquet.k_max
        k_range
        % hk_range - Helper property for [-hk_max:hk_max]
        % See also Calc.baseFloquet.hk_max
        hk_range
        % k_max2 - Helper property for 2*k_max+1
        % See also Calc.baseFloquet.k_max
        k_max2
        % hk_max2 - Helper property for 2*hk_max+1
        % See also Calc.baseFloquet.hk_max
        hk_max2
        % ht - Time-periodic Hamiltonian in Fourier representation
        % $$\hat{H}(t)=\sum_ke^{-ik\omega t}\hat{H}^{(k)}$$
        % See also baseFloquet.get_ht, Calc.baseFloquet.cache_ht
        ht
        % h - Time-periodic Hamiltonian in Floquet representation
        % $$\hat{\hat{H}}$$
        % See also Calc.baseFloquet.ht, baseFloquet.get_h,
        % Calc.baseFloquet.cache_h
        h
        % hf - Floquet Hamiltonian
        % $$\hat{H}^F(t)=\hat{H}(t)-i\hat{\partial}_t$$
        % See also Calc.baseFloquet.h, Calc.baseFloquet.pt,
        % Calc.baseFloquet.cache_hf, Calc.baseFloquet.hf2
        hf
    end
    properties (Hidden,Dependent,GetAccess=protected)
        % pt - Time derivative operator in Floquet representation
        % $$\hat{\partial}_t=-ik\omega\delta_{kl}\ket{\Phi^{(k)}}\bra{\Phi^{(l)}}$$
        % See also Calc.baseFloquet.hf, Calc.baseFloquet.w
        pt
        % hf2 - Helper property for hf*hf
        % See also Calc.baseFloquet.hf
        hf2
    end
    properties (Transient,NonCopyable,Access=private)
        % cache_ht - Cached matrix of ht
        % See also baseFloquet.get_ht, Calc.baseFloquet.ht,
        % Calc.baseFloquet.dirty_cache, Calc.baseFloquet.cacheMat
        cache_ht    (:,:)   double
        % cache_h - Cached matrix of h
        % See also baseFloquet.get_h, Calc.baseFloquet.h,
        % Calc.baseFloquet.dirty_cache, Calc.baseFloquet.cacheMat
        cache_h     (:,:)   double
        % cache_hf - Cached matrix of hf
        % See also Calc.baseFloquet.hf,
        % Calc.baseFloquet.dirty_cache, Calc.baseFloquet.cacheMat
        cache_hf    (:,:)   double
        % cache_hf2 - Cached matrix of hf2
        % See also Calc.baseFloquet.hf2,
        % Calc.baseFloquet.dirty_cache, Calc.baseFloquet.cacheMat
        cache_hf2   (:,:)   double
    end
    properties (Hidden,Transient,NonCopyable,SetAccess=protected)
        % dirty_cache - State of the cache
        %   Internal flag controlling when the cached matrices should be
        %   recalculated
        % See also Calc.baseFloquet.cacheMat
        dirty_cache (1,1)   logical = true
    end
    properties (Hidden,SetAccess=private)
        % cacheMat - Flag for caching the matrices
        % See also Calc.baseFloquet.ht, Calc.baseFloquet.h,
        % Calc.baseFloquet.hf, Calc.baseFloquet.hf2,
        % Calc.baseFloquet.dirty_cache
        cacheMat    (1,1)   logical
    end
    %% Constructor
    methods
        function obj = baseFloquet(N,Args)
            arguments
                N
                Args.w          = 0
                Args.k_max      = 100
                Args.hk_max     = 1
                Args.xi         = 1E-6
%                 Args.diag       (1,1)   logical = true
                Args.cacheMat   = true
            end
            % baseFloquet Constructor
            %
            % Syntax:
            %   obj = baseFloquet(N)
            %   [___] = baseFloquet(___,Name,Value)
            % 
            % Description:
            %   obj = baseFloquet(N) Constructs the base floquet object of size N
            %   [___] = baseFloquet(___,Name,Value) specifies options using name-value
            %   arguments in addition to any of the input arguments in previous
            %   syntaxes.
            %
            % Inputs:
            %   N - Size of the originial Hilbert space
            %   Name-Value pairs
            %
            % Outputs:
            %   obj - Floquet object
            %
            % Name-value arguments:
            %   w - [0] Driving frequency
            %   k_max - [100] Fourier cut-off of Fourier coefficients
            %   hk_max - [1] Fourier dimension of the time-periodic Hamiltonian
            %   xi - [1E-6] Acceptable error
            %   cacheMat - [true] Whether to cache Floquet matrices or not
            % 
            % See also baseFloquet

            obj.N = N;
            obj.w = Args.w;
            obj.k_max = Args.k_max;
            obj.hk_max = Args.hk_max;
            obj.xi = Args.xi;
            obj.cacheMat = Args.cacheMat;
%             if Args.diag
%                 % TODO: Only migrated the helper, not generalized
%                 % Construct sparse Index helper
%                 tcol = [];
%                 trow = [];
%                 for irow = 1:obj.N * obj.k_max2
%                     tcol = [tcol irow + obj.N * (-obj.hk_max:obj.hk_max)];
%                     trow = [trow irow*ones(1,obj.hk_max2)];
%                 end
%                 obj.sparseInd.val_flag = tcol > 0 & tcol <= obj.N * obj.k_max2;
%                 obj.sparseInd.row = trow(obj.sInd.val_flag);
%                 obj.sparseInd.col = tcol(obj.sInd.val_flag);
%                 obj.sparseInd.nnz = sum(obj.sInd.val_flag);
%                 obj.sparseMat=spalloc(obj.N*obj.k_max2,obj.N*obj.k_max2,obj.sInd.nnz);
%             else
%                 error('Not implemented')
%             end
        end
    end
    %% Basic overrides
    methods
        json = jsonencode(obj,varargin)
    end
    methods (Access=protected)
        groups = getPropertyGroups(obj)
    end
    %% Get/Setters
    methods
        function val = get.k_range(obj)
            val = -obj.k_max:obj.k_max;
        end
        function val = get.hk_range(obj)
            val = -obj.hk_max:obj.hk_max;
        end
        function val = get.k_max2(obj)
            val = 2 * obj.k_max + 1;
        end
        function val = get.hk_max2(obj)
            val = 2 * obj.hk_max + 1;
        end
        function val = get.pt(obj)
            val = -obj.k_range * obj.w;
            val = repmat(val,obj.N,1);
            val = val(:);
            val = spdiags(val,0,obj.N * obj.k_max2,obj.N * obj.k_max2);
        end
        function val = get.ht(obj)
            if obj.cacheMat
                if obj.dirty_cache; obj.CacheAll; end
                val = obj.cache_ht;
            else
                val = obj.get_ht;
            end
        end
        function val = get.h(obj)
            if obj.cacheMat
                if obj.dirty_cache; obj.CacheAll; end
                val = obj.cache_h;
            else
                val = obj.get_h;
            end
        end
        function val = get.hf(obj)
            if obj.cacheMat
                if obj.dirty_cache; obj.CacheAll; end
                val = obj.cache_hf;
            else
                val = obj.h + obj.pt;
            end
        end
        function val = get.hf2(obj)
            if obj.cacheMat
                if obj.dirty_cache; obj.CacheAll; end
                val = obj.cache_hf2;
            else
                thf = obj.hf;
                val = thf * thf;
            end
        end
        function set.w(obj,val)
            obj.w = val;
            obj.dirty_cache = true; %#ok<MCSUP> 
        end
        function set.k_max(obj,val)
            obj.k_max = val;
            obj.dirty_cache = true; %#ok<MCSUP> 
        end
    end
    %% Main methods
    % Basic definition methods
    methods (Hidden,Access=protected)
        h = get_h(obj)
        CacheAll(obj)
    end
    methods (Hidden,Abstract,Access=protected)
        get_ht(obj)
    end
    % Operation methods
    methods
        [S,om_mnk] = SpectraOverlap(obj,Psi1,Psi2,Args)
        Lk = Ladder(obj,k,Args)
        deps = epsDistance(obj,eps1,eps2)
        varargout = shiftBZ(obj,eps,Psi)
        [HBar,Psi] = HBar(obj,Psi,eps)
        Ebar = Ebar(obj,Psi,Args)
        eps = eps(obj,Psi,Args)
        varEps = varEps(obj,Psi,Args)
        Psi0 = Psi0(obj,Psi)
        Psi = Psi_Fourier(obj,Psi)
        Psi = Psi_Floquet(obj,Psi)
        Psit = Psit(obj,Psi,t)
        pk = Spectra(obj,Psi,Args)
        varargout = eigs(obj,Args)
        Psi = shiftCenter(obj,Psi,Args)
        Res = variational(obj,method,Psi0,Args)
        Psi = matchSize(obj,obj2,Psi2)
    end
    methods (Hidden)
        Res = variational_varQE(obj,Psi0,Args)
        Res = variational_AE(obj,Args)
    end
    methods (Hidden,Access=protected)
        Res = variational_base(obj,Psi0,fval,cons,Args)
    end
end