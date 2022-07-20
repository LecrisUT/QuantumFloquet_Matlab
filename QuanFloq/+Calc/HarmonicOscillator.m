classdef HarmonicOscillator < Calc.baseFloquet
    % HarmonicOscillator 1D Driven harmonic oscillator
    % $$H(x,t)=-\frac{1}{2}\partial_x^2+\frac{1}{2}\omega_0x^2-Vx\cos(\omega t)$$
    % The Hamiltonian and wave functions are projected onto the static
    % eigensates. See HarmonicOscillator.Phinx.
    %
    % See also HarmonicOscillatorB, HarmonicOscillatorC,
    % Calc.HarmonicOscillator.Phinx
    
    properties
        % w0 - Undriven harmonic frequency $\oemga_0$
        w0  (1,1)   double  {mustBeReal,mustBePositive}
        % V - Driving strength
        V   (1,1)   double  {mustBeReal}
    end
    properties (Dependent,SetAccess=private)
        ht
    end
    methods
        function val = get.ht(obj)
            % h0 = w_0 * (a^\dagger a + 0.5)
            % h(t) = h0 - V cos(wt) / \sqrt{2*w_0} * (a^\dagger + a)
            val = zeros(obj.N,obj.N,3);
            val(:,:,2) = obj.w0 * diag((0:obj.N-1) + 0.5);
            th1 = sqrt(1:(obj.N-1));
            th1 = -obj.V / 2 / sqrt(2 * obj.w0) * th1;
            val(:,:,3) = diag(th1,1) + diag(th1,-1);
            val(:,:,1) = val(:,:,3)';
        end
    end
    methods
        function obj = HarmonicOscillator(N,Args,Args2)
            arguments
                N
                Args.w0     = 1
                Args.V      = 0
                Args2.w     = 0
                Args2.k_max = 200
                Args2.xi    = 1E-10
            end
            % HarmonicOscillator Constructor
            %
            % Syntax:
            %   obj = HarmonicOscillator(N)
            %   [___] = HarmonicOscillator(___,Name,Value)
            % 
            % Description:
            %   obj = HarmonicOscillator(N) Construct the driven harmonic oscillator
            %   with Hilbert space cut off at N
            %   [___] = HarmonicOscillator(___,Name,Value) specifies options using
            %   name-value arguments in addition to any of the input arguments in
            %   previous syntaxes.
            %
            % Inputs:
            %   N - Hilbert space cut-off of a Harmonic oscillator
            %   Name-Value pairs
            %
            % Outputs:
            %   obj - Floquet object
            %
            % Name-value arguments:
            %   w0 - [1] Undriven harmonic frequency
            %   V - [0] Driving strength
            %   w - [0] Driving frequency
            %   k_max - [200] Fourier cut-off of Fourier coefficients
            %   xi - [1E-10] Acceptable error
            % 
            % See also HarmonicOscillator, baseFloquet

            Args2 = namedargs2cell(Args2);
            obj@Calc.baseFloquet(N,Args2{:});
            obj.w0 = Args.w0;
        end
        function set.V(obj,val)
            obj.V = val;
            obj.dirty_cache = true;
        end
        function json = jsonencode(obj,varargin)
            j = jsonencode@Calc.baseFloquet(obj);
            S = jsondecode(j);
            S.box_particle = struct('xi',obj.xi,'V',obj.V);
            json = jsonencode(S,varargin{:});
        end
        function Hx = Hx(obj,x,Args)
            arguments
                obj     Calc.HarmonicOscillator
                x       (:,1)   double
                Args.t  (:,1)   double
            end
            % Hx Time-periodic Hamiltonian in position basis
            % 
            % Syntax:
            %   Hx = Hx(x)
            %   [___] = Hx(___,Name,Value)
            % 
            % Description:
            %   Hx = Hx(x) Project the time-periodic Hamiltonian on the position basis
            %   x
            %   [___] = Hx(___,Name,Value) specifies options using name-value arguments
            %   in addition to any of the input arguments in previous syntaxes.
            % 
            % Inputs:
            %   x - Position basis
            %   Name-Value pairs
            %
            % Outputs:
            %   Hx - Classical solution
            %
            % Name-value arguments:
            %   t - Time parameter. Calculate the wave function at times t instead of
            %   Fourier representation

            nx = length(x);
            dx = diff(x);
            if length(uniquetol(dx)) > 1
                error('Different steps not implemented');
            end
            dx = dx(1);
            d2x = diag(ones(nx-1,1),1) - eye(nx);
            d2x = (d2x + d2x');
%             d2x(1,2) = 2;
%             d2x(end,end-1) = 2;
            d2x = d2x / dx^2;
            x2 = diag(x.^2);
            x = diag(x);
            if isfield(Args,'t')
                cwt = cos(obj.w*Args.t);
                cwt = reshape(cwt,1,1,[]);
                Hx = -0.5 * d2x + 0.5 * obj.w0^2 * x2 - obj.V * x .* cwt;
            else
                Hx = zeros(nx,nx,3);
                Hx(:,:,2) = -0.5 * d2x + 0.5 * obj.w0^2 * x2;
                Hx(:,:,3) = -obj.V/2 * x;
                Hx(:,:,1) = Hx(:,:,3)';
            end
        end
        function x = xc(obj,Args)
            arguments
                obj     Calc.HarmonicOscillator
                Args.t  (:,1)   double
            end
            % xc Calculate the classical solution
            % $$x_c=-\frac{V}{\omega^2-\omega_0^2}\cos(\omega t)$$
            % 
            % Syntax:
            %   x = xc
            %   [___] = xc(___,Name,Value)
            % 
            % Description:
            %   x = xc Calculate the classical solution
            %   [___] = xc(___,Name,Value) specifies options using name-value arguments
            %   in addition to any of the input arguments in previous syntaxes.
            % 
            % Inputs:
            %   Name-Value pairs
            %
            % Outputs:
            %   x - Classical solution
            %
            % Name-value arguments:
            %   t - Time parameter. Calculate the wave function at times t instead of
            %   Fourier representation
            %   
            % See also HarmonicOscillator.dxc

            if isfield(Args,'t')
                x = -obj.V / (obj.w^2 - obj.w0^2) * cos(obj.w * Args.t);
            else
                x = zeros(1,1,3);
                x(3) = -0.5 * obj.V / (obj.w^2 - obj.w0^2);
                x(1) = x(3);
            end
        end
        function dx = dxc(obj,Args)
            arguments
                obj     Calc.HarmonicOscillator
                Args.t  (:,1)   double
            end
            % dxc Calculate the classical solution time derivative $\dot{x}_c$
            % 
            % Syntax:
            %   x = dxc
            %   [___] = dxc(___,Name,Value)
            % 
            % Description:
            %   x = dxc Calculate the classical solution time derivative
            %   [___] = dxc(___,Name,Value) specifies options using name-value
            %   arguments in addition to any of the input arguments in previous
            %   syntaxes.
            % 
            % Inputs:
            %   Name-Value pairs
            %
            % Outputs:
            %   dx - Time-derivative
            %
            % Name-value arguments:
            %   t - Time parameter. Calculate the wave function at times t instead of
            %   Fourier representation
            %   
            % See also HarmonicOscillator.dxc
            % \dot{xc} = V*w/(w^2-w_0^2) * sin(w*t)
            if isfield(Args,'t')
                dx = obj.V * obj.w / (obj.w^2 - obj.w0^2) * sin(obj.w * Args.t);
            else
                dx = zeros(1,1,3);
                dx(3) = 0.5i * obj.V * obj.w / (obj.w^2 - obj.w0^2);
                dx(1) = dx(3)';
            end
        end
        function L = L(obj,Args)
            arguments
                obj     Calc.HarmonicOscillator
                Args.t  (:,1)   double
            end
            % L Classical Lagrangian
            % $$L(t,x_c,\dot{x}_c)=\frac{1}{2}\dot{x}_c^2-\frac{1}{2}\omega_0^2x_c^2+x_cV\cos(\omega t)$$
            % 
            % Syntax:
            %   L = L
            %   [___] = L(___,Name,Value)
            % 
            % Description:
            %   L = L Calculate the classical Lagrangian in Fourier representation
            %   [___] = L(___,Name,Value) specifies options using name-value arguments
            %   in addition to any of the input arguments in previous syntaxes.
            % 
            % Inputs:
            %   Name-Value pairs
            %
            % Outputs:
            %   L - Lagrangian
            %
            % Name-value arguments:
            %   t - Time parameter. Calculate the wave function at times t instead of
            %   Fourier representation
            %   
            % See also HarmonicOscillator.xc, HarmonicOscillator.dxc,
            % HarmonicOscillator.intL

            % 
            A = 0.5 * obj.V^2 / (obj.w^2 - obj.w0^2)^2;
            B = A * (obj.w0^2 - 3 * obj.w^2) / 2;
            C = A * obj.w^2 + B;
            % L(t) = C + B*cos(2wt)
            if isfield(Args,'t')
                L = C + B * cos(2 * obj.w * Args.t);
            else
                L = [B/2;0;C;0;B/2];
                L = reshape(L,1,1,5);
            end
        end
        function L = intL(obj,t)
            arguments
                obj     Calc.HarmonicOscillator
                t       (:,1)   double
            end
            % intL Integrated classical Lagrangian
            % $$\int_0^tL(\tau,x_c,\dot{x}_c)\dd{\tau}$$
            % 
            % Syntax:
            %   L = intL(t)
            % 
            % Description:
            %   L = intL(t) Calculate the intergrated classical Lagrangian from 0 to t
            % 
            % Inputs:
            %   t - Time parameter
            %
            % Outputs:
            %   L - Integrated Lagrangian
            %   
            % See also HarmonicOscillator.L

            A = 0.5 * obj.V^2 / (obj.w^2 - obj.w0^2)^2;
            B = A * (obj.w0^2 - 3 * obj.w^2) / 2;
            C = A * obj.w^2 + B;
            % L(t) = C*t + B/(2w)*sin(2wt)
            L = C*t + B / (2*obj.w) * sin(2 * obj.w * t);
        end
        function tPsix = Psix(obj,Psi,x,Args,Args2)
            arguments
                obj     Calc.HarmonicOscillator
                Psi     double
                x       (:,1)   double
                Args.t      (:,1)   double
                Args2.eps
            end
            % Psix Project the wave function onto position basis
            % $$\Psi(x,t) = \Psi_n(t)*\Phi_n(x)$$
            % 
            % Syntax:
            %   Psix = Psix(Psi,x)
            %   [___] = Psix(___,Name,Value)
            % 
            % Description:
            %   Psix = Psix(Psi,x) Project the wave function Psi onto the position
            %   basis set at points x
            %   [___] = Psix(___,Name,Value) specifies options using name-value
            %   arguments in addition to any of the input arguments in previous
            %   syntaxes.
            % 
            % Inputs:
            %   Psi - Wave function
            %   x - Position points
            %   Name-Value pairs
            %
            % Outputs:
            %   Psix - Projected wave function
            %
            % Name-value arguments:
            %   t - Time parameter. Calculate the wave function at times t instead of
            %   Fourier representation
            %   eps - Quasi-energies
            %   
            % See also baseFloquet.Psit, HarmonicOscillator.Phinx
            
            nx = length(x);
            if isfield(Args,'t')
                Args2 = namedargs2cell(Args2);
                Psi = obj.Psit(Psi,Args.t,Args2{:});
            else
                Psi = obj.Psi_Fourier(Psi);
            end
            Psi = permute(Psi,[1 2 4 3]);
            tPsix = zeros(obj.N,1,nx);
            for iN = 1:obj.N
                tPsix(iN,1,:) = obj.Phinx(iN,x);
            end
            tPsix = Psi .* tPsix;
            tPsix = sum(tPsix,1);
            tPsix = permute(tPsix,[3 2 4 1]);
        end
        function Phix = Phinx(obj,N,x)
            arguments
                obj     Calc.HarmonicOscillator
                N       (1,1)   double  {mustBeInteger}
                x       (:,1)   double
            end
            % Phinx Project the static eigenstates onto position basis
            % $$\varphi_n(x)=(\omega_0/\pi)^{\frac{1}{4}}\frac{1}{\sqrt{2^nn!}}H_n(\sqrt{\omega_0})e^{-\frac{\omega_0}{2}x^2}$$
            % 
            % Syntax:
            %   Phix = Phinx(N,x)
            % 
            % Description:
            %   Phix = Phinx(N,x) Project the static eigenstates $\varphi_n$ at
            %   positions x
            % 
            % Inputs:
            %   N - eigenstate number
            %   x - Position points
            %
            % Outputs:
            %   Phix - Projected wave function
            %   
            % See also HarmonicOscillator.Psix

            C = (obj.w0/pi)^0.25 / sqrt(2^N * factorial(N));
            Phix = C * hermiteH(N,sqrt(obj.w0) * x) .* exp(-obj.w0/2 * x.^2);
        end
    end
end