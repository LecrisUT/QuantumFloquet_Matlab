classdef BoxParticle < Calc.baseFloquet
    % BoxParticle 1D Driven particle in a box
    % <latex>
    % \begin{align}
    % H(x,t)&=-\frac{1}{2}\partial_x^2+V_0(x)+V\sin(\frac{\pi x}{2a})\cos(\omega t)\\
    % V_0(x)&=\begin{cases}
    % 0&\quad\text{if}\;\abs{x}\leq a\\
    % \inf&\quad\text{if}\;\abs{x}>a
    % \end{cases}
    % \end{align}
    % </latex>
    % The system is normalized so that _a_=1.
    % The Hamiltonian and wave functions are projected onto the static
    % eigensates. See BoxParticle.Psix.
    %
    % See also Calc.BoxParticle.Psix

    properties
        % V - Driving strength
        % See also BoxParticle
        V   (1,1)   double  {mustBeReal}
    end
    methods (Hidden,Access=protected)
        function ht = get_ht(obj)
            ht = zeros(obj.N,obj.N,3);
            ht(:,:,2) = diag((1:obj.N).^2);
            th1 = ones(obj.N-1,1);
            th1(2:2:end) = -1;
            th1 = obj.V/4 * th1;
            ht(:,:,3) = diag(th1,1)+diag(th1,-1);
            ht(:,:,1) = ht(:,:,3)';
        end
    end
    methods
        function obj = BoxParticle(N,Args,Args2)
            arguments
                N
                Args.V      = 0
                Args2.w     = 8.3
                Args2.k_max = 200
                Args2.xi    = 1E-2
            end
            % BoxParticle Constructor
            %
            % Syntax:
            %   obj = BoxParticle(N)
            %   [___] = BoxParticle(___,Name,Value)
            % 
            % Description:
            %   obj = BoxParticle(N) Constructs a driven particle in a box object with
            %   Hilbert space cut off at N
            %   [___] = BoxParticle(___,Name,Value) specifies options using name-value
            %   arguments in addition to any of the input arguments in previous
            %   syntaxes.
            %
            % Inputs:
            %   N - Hilbert space cut-off of a particle in a box Hamiltonian
            %   Name-Value pairs
            %
            % Outputs:
            %   obj - Floquet object
            %
            % Name-value arguments:
            %   V - [0] Driving strength
            %   w - [8.3] Driving frequency
            %   k_max - [200] Fourier cut-off of Fourier coefficients
            %   xi - [1E-2] Acceptable error
            % 
            % See also BoxParticle, baseFloquet

            Args2 = namedargs2cell(Args2);
            obj@Calc.baseFloquet(N,Args2{:});
            obj.V = Args.V;
        end
        function json = jsonencode(obj,varargin)
            j = jsonencode@Calc.baseFloquet(obj);
            S = jsondecode(j);
            S.box_particle = struct('xi',obj.xi,'V',obj.V);
            json = jsonencode(S,varargin{:});
        end
        function set.V(obj,val)
            obj.V = val;
            obj.dirty_cache = true;
        end
        function Psix = Psix(obj,Psi,x,Args,Args2)
            arguments
                obj     Calc.BoxParticle
                Psi     double
                x       (1,1,:) double
                Args.t
                Args.a      (1,1)   double  =1
                Args2.eps
            end
            % Psix Project the wave function onto position basis
            % <latex>
            % \begin{gather}
            % \varphi_n(t)=\begin{cases}
            % 1/\sqrt{a}\cos(\frac{n\pi x}{2a})&\quad\text{if}\;n\,\text{odd}\\
            % 1/\sqrt{a}\sin(\frac{n\pi x}{2a})&\quad\text{if}\;n\,\text{even}\\
            % \end{cases}
            % \end{gather}
            % </latex>
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
            %   a - [1] Position renormalization constant
            %   t - Time parameter. Calculate the wave function at times t instead of
            %   Fourier representation
            %   eps - Quasi-energies
            %   
            % See also baseFloquet.Psit

            nx = length(x);
            if isfield(Args,'t')
                Args2 = namedargs2cell(Args2);
                Psi = obj.Psit(Psi,t,Args2{:});
            else
                Psi = obj.Psi_Fourier(Psi);
            end
            Psi = permute(Psi,[1 2 4 3]);
            Psix = zeros(obj.N,1,nx);
            pi2a = pi / 2 / Args.a;
            for iN = 1:obj.N
                if mod(iN,2)
                    Psix(iN,1,:) = sqrt(1/Args.a) * cos(iN * pi2a * x);
                else
                    Psix(iN,1,:) = sqrt(1/Args.a) * sin(iN * pi2a * x);
                end
            end
            Psix = Psi .* Psix;
            Psix = sum(Psix,1);
            Psix = permute(Psix,[4 2 3 1]);
        end
    end
end