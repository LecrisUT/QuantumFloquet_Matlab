classdef HarmonicOscillatorB < Calc.HarmonicOscillator
    % HarmonicOscillatorB 1D Driven harmonic oscillator
    % The position parameter is shifted so that $y=x-x_c(t)$
    % $$H(y,t)=i\dot{x}_c\partial_y-\frac{1}{2}\partial_y^2+\frac{1}{2}\omega_0^2(y+x_c)^2-V(y+x_c)cos(\omega t)$$
    %
    % See also HarmonicOscillatorA, HarmonicOscillatorC

    methods (Access=protected)
        function ht = get_ht(obj)
            % h0 = w_0 * (a^\dagger a + 0.5)
            % h0(t) = h0 + 0.5*w_0^2*xc^2 - V*xc*cos(w*t)
            % h1(t) = \sqrt{w_0/2}*xc + 1i*\sqrt{w_0/2}*\dot{xc}-
            % V/\sqrt{2*w_0}*cos(w*t)
            % h(t) = h0(t) + (h1(t) * a + h.c.)
            A = 0.5 * obj.V^2 * (obj.w^2 - 0.5 * obj.w0^2) / (obj.w^2 - obj.w0^2)^2;
            B = obj.V * obj.w / sqrt(2 * obj.w0) / (obj.w^2 - obj.w0^2)^2;
            C = -obj.w/2 - obj.w0/2;
            Cp = -obj.w/2 + obj.w0/2;
            ht = zeros(obj.N,obj.N,5);
            ht(:,:,3) = obj.w0 * diag((0:obj.N-1) + 0.5 + A/2);
            th1 = sqrt(1:(obj.N-1));
            ht(:,:,5) = A/2 * eye(obj.N);
            ht(:,:,4) = B*C*diag(th1,1) + B*Cp*diag(th1,-1);
            ht(:,:,1) = ht(:,:,5)';
            ht(:,:,2) = ht(:,:,4)';
        end
    end
    methods
        function obj = HarmonicOscillatorB(N,Args2)
            arguments
                N
                Args2.w0
                Args2.w
                Args2.k_max
                Args2.xi
            end
            % HarmonicOscillatorB Constructor
            % 
            % See also Calc.HarmonicOscillator.HarmonicOscillator

            Args2 = namedargs2cell(Args2);
            obj@Calc.HarmonicOscillator(N,Args2{:});
            obj.hk_max = 2;
        end
        function Hy = Hy(obj,y,Args)
            arguments
                obj     Calc.HarmonicOscillatorB
                y       (:,1)   double
                Args.t  (:,1)   double
            end
            % Hy Time-periodic Hamiltonian in shifted position basis
            % See Calc.HarmonicOscillator.Hx
            % 
            % See also See Calc.HarmonicOscillator.Hx

            ny = length(y);
            dy = diff(y);
            if length(uniquetol(dy)) > 1
                error('Different steps not implemented');
            end
            dy = dy(1);
            d1y = diag(ones(ny-1,1),1) - eye(ny);
            d2y = (d1y + d1y');
%             d1y(end,end-1) = 1;
%             d2y(1,2) = 2;
%             d2y(end,end-1) = 2;
            d1y = d1y / dy;
            d2y = d2y / dy^2;
            y = diag(y);
            if isfield(Args,'t')
                xc = obj.xc(t=Args.t);
                xc = reshape(xc,1,1,[]);
                dxc = obj.dxc(t=Args.t);
                dxc = reshape(dxc,1,1,[]);
                cwt = cos(obj.w*Args.t);
                cwt = reshape(cwt,1,1,[]);
                Hy = 1i * dxc * d1y - 0.5 * d2y + 0.5 * obj.w0^2 * (y + diag(xc(:))).^2 - obj.V * (y + diag(xc(:))) .* cwt;
            else
                xc = obj.xc;
                dxc = obj.dxc;
                Hy = zeros(ny,ny,5);
                Hy(:,:,3) = - 0.5 * d2y + 0.5 * obj.w0^2 * y.^2 ...
                    + obj.w0^2 * xc(3)^2 * eye(ny) ...
                    - obj.V * xc(3) * eye(ny);
                Hy(:,:,4) = 1i * dxc(3) * d1y + obj.w0^2 * y * xc(3) - obj.V/2 * y;
                Hy(:,:,5) = 0.5 * obj.w0^2 * xc(3)^2 * eye(ny) - obj.V/2 * xc(3) * eye(ny);
                Hy(:,:,2) = Hy(:,:,4)';
                Hy(:,:,1) = Hy(:,:,5)';
            end
        end
        function tPsiy = Psiy(obj,Psi,y,Args2)
            arguments
                obj
                Psi
                y
                Args2.t
                Args2.eps
            end
            % Psiy Project the wave function onto the shifted position basis
            % Equivalent to Psix(___,centerx=true)
            %   
            % See also Calc.HarmonicOscillatorB.Psix
            
            Args2 = namedargs2cell(Args2);
            tPsiy = obj.Psix(Psi,y,Args2{:},centerx=true);
        end
        function tPsix = Psix(obj,Psi,x,Args,Args2)
            arguments
                obj     Calc.HarmonicOscillatorB
                Psi     double
                x       (:,1)   double
                Args.t          (:,1)   double
                Args.centerx    (1,1)   logical = false
                Args2.eps
            end
            % Psix Project the wave function onto position basis
            % See Calc.HarmonicOscillator.Psix.
            %
            % Additional name-value arguments:
            %   centerx - [false] Whether or not to shift the positions with respect to
            %   x_c
            %   
            % See also Calc.HarmonicOscillator.Psix

            if ~Args.centerx && isfield(Args,'t')
                if ~isfield(Args,'t')
                    error('Missing argument t');
                end
                Args2 = namedargs2cell(Args2);
                nx = length(x);
                nt = length(Args.t);
                M = size(obj.Psi_Floquet(Psi),2);
                tPsix = zeros(nx,M,nt);
                xc = obj.xc(t=Args.t);
                for it = 1:length(Args.t)
                    tt = Args.t(it);
                    tPsix(:,:,it) = Psix@Calc.HarmonicOscillator(obj,...
                        Psi,x-xc(it),Args2{:},t=tt);
                end
            else
                if isfield(Args,'t')
                    Args2.t = Args.t;
                end
                Args2 = namedargs2cell(Args2);
                tPsix = Psix@Calc.HarmonicOscillator(obj,Psi,x,Args2{:});
            end
        end
    end
end