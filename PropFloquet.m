function Psit = PropFloquet(t,w,Psi,Args)
    arguments
        t
        w
        Psi         (:,:,:)
        Args.eps    (1,1,:) double
        Args.Psi0   (:,1)   double
    end
    % PROPFLOQUET Propagate the Floquet wavefunction

    k_max = (size(Psi,3) - 1) / 2;
    % Fourier factor
    kwt(1,1,:) = exp(-1i*(-k_max:k_max)*w*t);
    if isfield(Args,'eps')
        % quasi-energy factor
        epst(1,:) = exp(-1i*eps_ex*t);
        % Time dependent wavefunction
        Psit = sum(epst .* kwt .* Psi,3);
    else
        % Time dependent wavefunction
        Psit = sum(kwt .* Psi,3);
    end
    if isfield(Args,'Psi0')
        Psi0 = sum(Psi,3);
        % Projection of initial wavefunction
        tPsi0 = Psi0' * Args.Psi0;
        % Project the wavefunction onto the time-dependent one
        Psit = Psit * tPsi0;
    end
end