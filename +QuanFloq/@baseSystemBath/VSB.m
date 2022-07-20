function VSB = VSB(obj,Psi)
    arguments
        obj     QuanFloq.baseSystemBath
        Psi     double
    end
    % VSB Project the system-bath coupling onto a given basis
    % 
    % Syntax:
    %   VSB = VSB(Psi)
    % 
    % Description:
    %   VSB = VSB(Psi) Project the system-bath coupling onto the Psi basis
    % 
    % Inputs:
    %   Psi - Floquet basis
    %
    % Outputs:
    %   VSB - System-bath coupling in the provided basis
    %   
    % See also QuanFloq.baseSystemBath.V

    Psi = obj.system.Psi_Fourier(Psi);
    M = size(Psi,2);
    VSB = zeros(M,M,obj.hk_max2);
    for ik = obj.hk_range
        if ik < 0
            % Calculate V * Psi^{k-ik}
            ttV = pagemtimes(obj.V,Psi(:,:,1:obj.system.k_max2+ik));
            % Calculate Psi^{k}' * V * Psi^{k-ik}
            ttV = pagemtimes(Psi(:,:,1-ik:obj.system.k_max2),'ctranspose',...
                ttV,'none');
            % Sum over k
            VSB(:,:,obj.hk_max+1+ik) = sum(ttV,3);
        else
            % Calculate V * Psi^{k+ik}
            ttV = pagemtimes(obj.V,Psi(:,:,1+ik:obj.system.k_max2));
            % Calculate Psi^{k}' * V * Psi^{k+ik}
            ttV = pagemtimes(Psi(:,:,1:obj.system.k_max2-ik),'ctranspose',...
                ttV,'none');
            % Sum over k
            VSB(:,:,obj.hk_max+1+ik) = sum(ttV,3);
        end
    end
end