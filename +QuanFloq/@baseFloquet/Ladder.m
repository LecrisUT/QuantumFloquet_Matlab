function Lk = Ladder(obj,k,Args)
    arguments
        obj     QuanFloq.baseFloquet
        k       (1,1)   double  {mustBeInteger}
        Args.dagger (1,1)   logical = false
        Args.eps    (1,1)   double  {mustBeReal}
        Args.tol    (1,1)   double  {mustBeReal,mustBePositive} = 1E-10
    end
    % Ladder Calculate or estimate the ladder operator
    % (Currently only estimation of single frequency Hamiltoninan is implemented)
    % Ladder operator $L^{(\dagger)}_k$ lowers (raises) the Fourier index
    % of the Floquet eigenstate
    % $$L_k\ket{\Phi^{(k))}=\ket{\Phi^{(k-1)}}$$
    % $$L^{\dagger}_k\ket{\Phi^{(k))}=\ket{\Phi^{(k+1)}}$$
    % $$[H^{(0)}-(\epsilon+k\omega)]+\sum_l^{hk_{max}}H^{(l)}\prod_m=1^lL_{k+m}+\sum_l=1^{hk_{max}}H^{(-l)}\prod_m=1^lL^{\dagger}_{k-m}=0$$
    % 
    % Syntax:
    %   Lk = Ladder(k)
    %   [___] = Ladder(___,Name,Value)
    % 
    % Description:
    %   Lk = Ladder(k) Calculate the ladder operators from the _k_th index to
    %   convergence
    %   [___] = Ladder(___,Name,Value) specifies options using name-value
    %   arguments in addition to any of the input arguments in previous
    %   syntaxes.
    % 
    % Inputs:
    %   k - Fourier index from which to start calculating
    %   Name-Value pairs
    %
    % Outputs:
    %   Lk - All of the ladder operators from _k_th to convergence
    %
    % Name-value arguments:
    %   dagger - [false] Calculate the raising ladder operator instead of
    %   lowering
    %   eps - Quasi-energy eigenvalue (estimate). If not provided, calculating
    %   an estimate
    %   tol - Cut-off tolerance to converge the ladder operators

    if isfield(Args,'eps')
        if obj.hk_max > 1
            error('Not implemented for hk_max > 1');
        else
        end
        error('Not Implemented');
    else
        if obj.hk_max > 1
            warning('Not implemented for hk_max > 1, using only H^(1)')
        end
        tht = obj.ht;
        th = tht(:,:,obj.hk_max+1 + 1) + tht(:,:,obj.hk_max+1 - 1);
        max_E = eig(th,'vector');
        max_E = max(abs(max_E));
        Lk = abs(max_E / k / obj.w);
    end
end