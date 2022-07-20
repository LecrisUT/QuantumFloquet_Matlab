function h = get_h(obj)
    arguments
        obj     Calc.baseFloquet
    end
    % get_h Basic getter for h
    % Allows to be subclassed for custom optimizations
    % 
    % Syntax:
    %   h = get_h
    % 
    % Description:
    %   h = get_h Calculates the corresponding time-periodic Hamiltonian in
    %   Floquet representation
    % 
    % Inputs:
    %   [none]
    %
    % Outputs:
    %   h - Time-periodic Hamiltonian in Floquet representation
    % 
    % See also Calc.baseFloquet.h, Calc.baseFloquet.cache_h

    tht = reshape(permute(obj.ht,[1 3 2]),obj.N*obj.hk_max2,obj.N);
    [h_diag,ind] = spdiags(tht);
    ind = ind + obj.N * obj.hk_max;
    h_diag = repmat(h_diag,obj.k_max2,1);
    h = spdiags(h_diag,ind,obj.N * obj.k_max2, obj.N * obj.k_max2);
end