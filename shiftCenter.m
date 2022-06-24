function psi = shiftCenter(psi, N, Args)
    arguments
        psi
        N
        Args.tol    = 1E-6
    end
    % SHIFTCENTER Center the wavefunction in Fourier space

    k_max = (size(psi,1)/N - 1) / 2;
    for iM = 1:size(psi,2)
        tpsi = reshape(psi(:,iM),N,[]);
        tnorm = vecnorm(tpsi,2,1);
        [peak,tk] = findpeaks(smoothdata(tnorm),...
            "NPeaks",1, "SortStr","descend");
        tpsi = circshift(tpsi, k_max+1 - tk, 2);
        % Normalize to initial WF
        tpsi0 = sum(tpsi,2);
        if tpsi0(1)<0; tpsi=-tpsi; end
        psi(:,iM) = tpsi(:);
    end
    iM = 2;
    S = psi' * psi;
    while iM <= size(psi,2)
        % Flag if found duplicate
        flag = false;
        for iM2 = 1:iM-1
            % Check if close to duplicate
            if (1 - abs(S(iM2,iM))) < Args.tol
                % Remove duplicate
                flag = true;
                psi(:,iM) = [];
                S(:,iM) = [];
                S(iM,:) = [];
                break;
            end
        end
        if ~flag; iM = iM + 1; end
    end
end