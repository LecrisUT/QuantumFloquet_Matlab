classdef (Abstract) baseFloquetHF < Calc.baseFloquet
    % baseFloquetHF Base class for Floquet Hartree-Fock calculations
    %
    % See also HubardDimerOrbital

    properties
        % Ne - Number of electrons
        Ne          (1,1)   double  {mustBeInteger,mustBePositive}
        % hUk_max - Fourier cutoff of the self-consistent Interaction
        % See also Calc.baseFloquetHF.hUk_max2, Calc.baseFloquetHF.hUk_range
        hUk_max     (1,1)   double  {mustBeInteger,mustBeNonnegative}
        % Psi - Orbital wave functions
        Psi         (:,:)   double
        % cacheMatU - Flag for caching the U matrices
        cacheMatU   (1,1)   logical
    end
    properties (SetAccess=private)
        % mode - Orbital occupation method
        % Currently only 'closed-shell` is implemented
        mode    {mustBeMember(mode,{'closed-shell','restricted open-shell','unrestricted open-shell'})}
    end
    properties (Abstract,Dependent,SetAccess=private)
        % hUt - Time-periodic self-consistent interaction in Fourier representation
        hUt     (:,:,:)     double
    end
    properties (Dependent)
        % NOrbitals - Number of electron orbitals
        % See also Calc.baseFloquetHF.Ne, Calc.baseFloquetHF.mode
        NOrbitals
        % hUk_range - Helper property for [-hUk_max:hUk_max]
        % See also Calc.baseFloquetHF.hUk_max
        hUk_range
        % hUk_max2 - Helper property for 2*hUk_max+1
        % See also Calc.baseFloquetHF.hUk_max
        hUk_max2
        % Ht - Many-body Hamiltonian in Fourier representation
        Ht
        % hU - Time-periodic self-consistent interaction in Floquet representation
        hU
        % H - Many-body Hamiltonian in Floquet representation
        H
        % Hf - Many-body Floquet Hamiltonian
        Hf
        % Ft - Fock operator in Fourier representation
        Ft
        % F - Fock operator in Floquet representation
        F
        % Ff - Floquet Fock operator
        Ff
    end
    properties (Transient,NonCopyable,Access=private)
        % cache_hU - Cached matrix of hU
        cache_hU    (:,:)   double
        % cache_H - Cached matrix of H
        cache_H     (:,:)   double
        % cache_Hf - Cached matrix of Hf
        cache_Hf    (:,:)   double
        % cache_F - Cached matrix of F
        cache_F     (:,:)   double
        % cache_Ff - Cached matrix of Ff
        cache_Ff    (:,:)   double
    end
    properties (Hidden,Transient,NonCopyable,SetAccess=protected)
        % dirty_cacheU - State of the U cache
        % Flag signaling if the orbital wave functions has changed
        dirty_cacheU    (1,1)   logical = true
    end
    %% Constructor
    methods
        function obj = baseFloquetHF(N,Args,Args2)
            arguments
                N
                Args.Ne         = 2
                Args.hUk_max    = 10
                Args.mode       = 'closed-shell'
                Args.cacheMatU  = true
                Args2.w
                Args2.k_max
                Args2.hk_max
                Args2.xi
                Args2.cacheMat
            end
            % baseFloquetHF Constructor
            % See Calc.baseFloquet.baseFloquet for usage
            %
            % Additional name-value arguments:
            %   Ne - [2] Number of electrons
            %   hUk_max - [10] Fourier cut-off of the self-consistent
            %   interaction
            %   mode - ['closed-shell'] Hartree-Fock orbital model
            %   cacheMatU - [true] Whether to cache the self-consistent
            %   interaction matrices
            % 
            % See also baseFloquetHF, Calc.baseFloquet.baseFloquet

            Args2 = namedargs2cell(Args2);
            obj@Calc.baseFloquet(N,Args2{:});
            obj.Ne = Args.Ne;
            obj.hUk_max = Args.hUk_max;
            obj.mode = Args.mode;
            obj.cacheMatU = Args.cacheMatU;
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
        function set.Psi(obj,val)
            if isempty(val)
                % Allow for initial initialization
                obj.Psi = val;
            else
                tPsi = obj.Psi_Floquet(val);
                switch obj.mode %#ok<MCSUP> 
                    case 'closed-shell'
                        if size(tPsi,2) ~= obj.Ne/2 %#ok<MCSUP> 
                            error('Dimensions are not compatible')
                        end
                    otherwise
                        error('Not implemented');
                end
                % If update is not within numerical error mark to recalculate U
                if ~all(abs(obj.Psi(:)-tPsi(:)) < 1E-15)
                    obj.dirty_cacheU = true; %#ok<MCSUP> 
                end
                obj.Psi = tPsi;
            end
        end
        function val = get.NOrbitals(obj)
            switch obj.mode
                case 'closed-shell'
                    val = obj.Ne / 2;
                case {'restricted open-shell','unrestricted open-shell'}
                    val = obj.Ne;
            end
        end
        function val = get.hUk_range(obj)
            val = -obj.hUk_max:obj.hUk_max;
        end
        function val = get.hUk_max2(obj)
            val = 2 * obj.hUk_max + 1;
        end
        function val = get.Ht(obj)
            max_hk = max(obj.hk_max,obj.hUk_max);
            val = zeros(obj.N,obj.N,2 * max_hk + 1);
            val(:,:,max_hk+1+obj.hk_range) = obj.ht;
            val(:,:,max_hk+1+obj.hUk_range) = val(:,:,max_hk+1+obj.hUk_range) + ...
                obj.hUt;
        end
        function val = get.hU(obj)
            if obj.cacheMatU
                if obj.dirty_cacheU; obj.CacheAllU; end
                val = obj.cache_hU;
            else
                val = obj.get_hU;
            end
        end
        function val = get.H(obj)
            if obj.cacheMatU
                if obj.dirty_cacheU; obj.CacheAllU; end
                val = obj.cache_H;
            else
                val = obj.get_H;
            end
        end
        function val = get.F(obj)
            if obj.cacheMatU
                if obj.dirty_cacheU; obj.CacheAllU; end
                val = obj.cache_F;
            else
                val = obj.get_F;
            end
        end
        function val = get.Hf(obj)
            if obj.cacheMatU
                if obj.dirty_cacheU; obj.CacheAllU; end
                val = obj.cache_Hf;
            else
                val = obj.H + obj.pt;
            end
        end
        function val = get.Ff(obj)
            if obj.cacheMatU
                if obj.dirty_cacheU; obj.CacheAllU; end
                val = obj.cache_Ff;
            else
                val = obj.F + obj.pt;
            end
        end
    end
    %% Main methods
    % Basic definition methods
    methods (Hidden,Access=protected)
        hU = get_hU(obj)
        H = get_H(obj)
        F = get_F(obj)
        CacheAllU(obj)
    end
    % Overloads
    methods
        Ebar = Ebar(obj,Psi,Args)
        eps = eps(obj,Psi,Args)
        varEps = varEps(obj,Psi,Args)
        function HBar(varargin)
            % HBar [Ill-defined]
            % Exact many-body average energy operator is ill-defined in Floquet
            % Hartree-Fock. Use H or Ebar instead
            %
            % See also Calc.baseFloquetHF.H, baseFloquetHF.Ebar

            error('Exact many-body average energy operator is ill-defined in Floquet Hartree-Fock.\nUse H or Ebar instead.')
        end
        varargout = eigs(obj,Args)
        Res = variational_varQE(obj,Psi0,Args)
        Res = variational_AE(obj,Args)
    end
    % Operation methods
    methods
        [Psi_Orb,Psi,eps,Ebar] = OccupiedOrbitals(obj,Psi,eps,Ebar,Args)
        [HBar,Psi] = FBar(obj,Psi,eps)
        function Psi_Slat = OrbitalToSlater(obj,Psi_Orb,Args)
            arguments
                obj         Calc.baseFloquetHF
                Psi_Orb     double
                Args.objSlater          (1,1)   Calc.baseFloquet
                Args.normalize          (1,1)   logical = true
                Args.normalizeSlater    (1,1)   logical = false
            end
            % OrbitalToSlater Calculate the corresponding Slater determinant
            % No common Slater representation is provided. This method has to be
            % individually implemented for each system
            % 
            % Syntax:
            %   Psi_Slat = OrbitalToSlater(Psi_Orb)
            %   [___] = OrbitalToSlater(___,Name,Value)
            % 
            % Description:
            %   Psi_Slat = OrbitalToSlater(Psi_Orb) Transform the Psi_Orb orbital into
            %   a corresponding Slater determinant
            %   [___] = OrbitalToSlater(___,Name,Value) specifies options using
            %   name-value arguments in addition to any of the input arguments in
            %   previous syntaxes.
            % 
            % Inputs:
            %   Psi_Orb - Orbital wave function
            %   Name-Value pairs
            %
            % Outputs:
            %   Psi_Slat - Slater determinant in Floquet representation with
            %   pre-defined basis
            %
            % Name-value arguments:
            %   objSlater - Slater Floquet object to match size with
            %   normalize - [true] Whether to normalize the orbital wave functions
            %   normalizeSlater - [false] Whether to normalize the Slater determinant.
            %   (Not recommended)

            if Args.normalize
                Psi_Orb = Psi_Orb ./ vecnorm(Psi_Orb);
            end
            error('Not implemented');
            if Args.normalizeSlater
                Psi_Slat = Psi_Slat ./ vecnorm(Psi_Slat);
            end
        end
    end
end