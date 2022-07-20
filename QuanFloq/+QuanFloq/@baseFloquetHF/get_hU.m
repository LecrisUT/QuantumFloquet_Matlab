function hU = get_hU(obj)
    arguments
        obj     QuanFloq.baseFloquetHF
    end
    % get_hU Basic getter for hU
    % [To be documented]
    % 
    % Syntax:
    %   hU = get_hU
    % 
    % Description:
    %   hU = get_hU
    % 
    % Inputs:
    %   [none]
    %
    % Outputs:
    %   hU
    % 
    % See also QuanFloq.baseFloquet.get_h

    thUt = reshape(permute(obj.hUt,[1 3 2]),obj.N*obj.hUk_max2,obj.N);
    [hU_diag,ind] = spdiags(thUt);
    ind = ind + obj.N * obj.hUk_max;
    hU_diag = repmat(hU_diag,obj.k_max2,1);
    hU = spdiags(hU_diag,ind,obj.N * obj.k_max2, obj.N * obj.k_max2);
end