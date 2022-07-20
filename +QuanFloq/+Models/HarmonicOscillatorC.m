classdef HarmonicOscillatorC < QuanFloq.Models.HarmonicOscillatorB
    % HarmonicOscillatorC 1D Driven harmonic oscillator
    % Projected onto an efficient basis. See QuanFloq.Models.HarmonicOscillatorC.Psix
    % $$H(y,t)=-\frac{1}{2}\partial_y^2+\frac{1}{2}\omega_0^2y^2-L(t,x_c,\dot{x}_c}$$
    %
    % See also QuanFloq.Models.HarmonicOscillator,
    % QuanFloq.Models.HarmonicOscillatorB,
    % QuanFloq.Models.HarmonicOscillatorC.Psix

    methods (Access=protected)
        function val = get_ht(obj)
            % h0 = w_0 * (a^\dagger a + 0.5)
            % h(t) = h0 - L(t)
            L=obj.L;
            val=zeros(obj.N,obj.N,5);
            val(:,:,3) = obj.w0 * diag((0:obj.N-1) + 0.5 - L(3));
            val(:,:,1) = -L(1) * eye(obj.N);
            val(:,:,5) = -L(5) * eye(obj.N);
        end
    end
    methods
        function obj = HarmonicOscillatorC(N,Args2)
            arguments
                N
                Args2.w0
                Args2.w
                Args2.k_max
                Args2.xi
            end
            % HarmonicOscillatorC Constructor
            % 
            % See also QuanFloq.HarmonicOscillator.HarmonicOscillatorB

            Args2 = namedargs2cell(Args2);
            obj@QuanFloq.Models.HarmonicOscillatorB(N,Args2{:});
            obj.hk_max = 2;
        end
        function tPsix = Psix(obj,Psi,x,Args,Args2)
            arguments
                obj
                Psi
                x               (:,1) double
                Args.t
                Args.centerx    (1,1)   logical = false
                Args2.eps
            end
            % Psix Project the wave function onto position basis
            % The basis set is transformed into
            % $$\Phi_n(y,t)=e^{i\dot{x}_cy}\varphi_n(y)\ket{\Phi_n(t)}$$
            % See QuanFloq.Models.HarmonicOscillatorB.Psix.
            %   
            % See also QuanFloq.Models.HarmonicOscillatorB.Psix

            if isfield(Args,'t'); Args2.t = Args.t; end
            Args2.centerx = Args.centerx;
            Args2 = namedargs2cell(Args2);
            tPsix = Psix@QuanFloq.HarmonicOscillatorB(obj,...
                Psi,x,Args2{:});
            if ~isfield(Args,'t')
                % TODO: transform exp(1i*sin(w*t)) to Fourier components
                error('Not implemented');
            else
                if Args.centerx
                    dxc = obj.dxc(t=Args.t);
                    dxc = reshape(dxc,1,1,[]);
                    tPsix = exp(1i * dxc .* x) .* tPsix;
                else
                    dxc = obj.dxc(t=Args.t);
                    dxc = reshape(dxc,1,1,[]);
                    xc = obj.xc(t=Args.t);
                    xc = reshape(xc,1,1,[]);
                    tPsix = exp(1i * dxc .* (x - xc)) .* tPsix;
                end
            end
        end
    end
end