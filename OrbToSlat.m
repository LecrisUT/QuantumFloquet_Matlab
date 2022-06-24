function Psi=OrbToSlat(psi)
    % ORBTOSLAT Generate the corresponding Slater orbital
    %
    % Currently this is assuming HubbardDimer wavefunctions

    switch ndims(psi)
        case 2
            psi = reshape(psi,2,1,[]);
        case 3
        otherwise
            error('Not expected')
    end
    N = size(psi,1);
    M = size(psi,2);
    k_max = (size(psi,3) - 1) / 2;
    k_max2 = 2 * k_max+1;
    Psi=zeros(3,M,k_max2);
    for k = -k_max:k_max
        if k>0
            % psi(k+l)*psi(l) with l in [-k_max:k_max]
            tPsi = [psi(1,:,1+k:k_max2).*psi(1,:,k_max2:-1:1+k);
                psi(1,:,1+k:k_max2).*psi(2,:,k_max2:-1:1+k);
                psi(2,:,1+k:k_max2).*psi(2,:,k_max2:-1:1+k)];
        else
            tPsi = [psi(1,:,1:k_max2+k).*psi(1,:,k_max2+k:-1:1);
                psi(1,:,1:k_max2+k).*psi(2,:,k_max2+k:-1:1);
                psi(2,:,1:k_max2+k).*psi(2,:,k_max2+k:-1:1)];
        end
        tPsi=sum(tPsi,3);
        tPsi(2,:) = tPsi(2,:)*sqrt(2);
        Psi(:,:,k_max+1+k) = tPsi;
    end
%     % If Psi(0) is normalized, this shouldn't be necessary (affects the
%     % time-dependence)
%     nrm = norm(Psi(:));
%     Psi=Psi./nrm;
end