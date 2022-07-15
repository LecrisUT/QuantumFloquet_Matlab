function Res = Gamma(obj,Args)
    arguments
        obj     Calc.baseSystemBath
        Args.Psi    double
        Args.eps    (:,1)   double
        Args.type   {mustBeMember(Args.type,{'Redfield','Lindblad','Secular'})} = 'Lindblad'
    end
    % Gamma Calculate the dissipation strength
    % 
    % Syntax:
    %   Res = Gamma
    %   [___] = Gamma(___,Name,Value)
    % 
    % Description:
    %   Res = Gamma Calculate the full average enery operator
    %   [___] = Gamma(___,Name,Value) specifies options using name-value
    %   arguments in addition to any of the input arguments in previous
    %   syntaxes.
    % 
    % Inputs:
    %   Name-Value pairs
    %
    % Outputs:
    %   Res - Structured output. Depending on type chosen may contain fields:
    %     - Gamma - Dissipation strengths
    %     - om_mnk - Transition frequencies
    %
    % Name-value arguments:
    %   Psi - Quasi-energy eigenstates. If not provided calculates the full
    %   average energy eigenstates
    %   eps - Quasi-energy eigenvalues. If not provided calculates from Psi
    %   type - ['Lindblad'] Type of master equation. Can be 'Redfield',
    %   'Lindblad', or 'Secular'.
    %   
    % See also baseSystemBath.SteadyState, baseSystemBath.VSB, basebath.gamma

    % Calculate the eigensates if not provided
    if ~isfield(Args,'Psi')
        [Args.Psi,Args.eps,~] = obj.system.eigs;
    end
    % Calculate the quasi-energies if not provided
    if ~isfield(Args,'eps')
        Args.eps = obj.system.eps(Args.Psi,normalize=true);
    end
    Args.Psi = obj.system.Psi_Fourier(Args.Psi);
    M = size(Args.Psi,2);
    VSB = obj.VSB(Args.Psi);
    Gamma = zeros(M,M,M,M,obj.hk_max2,obj.hk_max2);
    om_mnk = zeros(M,M,M,M,obj.hk_max2,obj.hk_max2);
    for jm = 1:M
        for jn = 1:M
            for jk = obj.hk_range
                jom_mnk = Args.eps(jn) - Args.eps(jm) + jk * obj.system.w;
                tgam = obj.bath.gamma(jom_mnk);
                for im = 1:M
                    for in = 1:M
                        for ik = obj.hk_range
                            iom_mnk = Args.eps(in) - Args.eps(im) + ik * obj.system.w;
                            om_mnk(im,in,jm,jn,...
                                obj.hk_max+1+ik,obj.hk_max+1+jk) = iom_mnk + jom_mnk;
                            Gamma(im,in,jm,jn,...
                                obj.hk_max+1+ik,obj.hk_max+1+jk) = VSB(im,in,obj.hk_max+1+ik) * ...
                                VSB(jm,jn,obj.hk_max+1+jk) * tgam;
                        end
                    end
                end
            end
        end
    end
    switch Args.type
        case 'Lindblad'
            Res.Gamma = zeros(M,M,obj.hk_max2);
            for im = 1:M
                for in = 1:M
                    for ik = obj.hk_range
                        Res.Gamma(im,in,obj.hk_max+1+ik) = Gamma(im,in,in,im,obj.hk_max+1+ik,obj.hk_max+1-ik);
                    end
                end
            end
        case 'Secular'
            error('Not implemented')
        case 'Redfield'
            Res.Gamma = Gamma;
            Res.omega_mnk = om_mnk;
        otherwise
            error('Unsupported type')
    end
end