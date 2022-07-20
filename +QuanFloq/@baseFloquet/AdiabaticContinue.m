function [Psi,ind] = AdiabaticContinue(obj,Psi,Psi_prev,Args)
    arguments
        obj     QuanFloq.baseFloquet
        Psi         double
        Psi_prev    (:,:)   double
        Args.eps        (:,1)   double
        Args.eps_prev   (:,1)   double
        Args.tol        (1,1)   double  {mustBeReal,mustBePositive} = 1E-1
    end
    % AdiabaticContinue Reorder the Floquet states to adiabatically continue
    % 
    % Syntax:
    %   Psi = AdiabaticContinue(Psi,Psi_prev)
    %   [Psi,ind] = AdiabaticContinue(Psi,Psi_prev)
    %   [___] = AdiabaticContinue(___,Name,Value)
    % 
    % Description:
    %   Psi = AdiabaticContinue(Psi,Psi_prev) Reorder the states Psi to
    %   adiabatically continue to Psi_prev
    %   [Psi,ind] = AdiabaticContinue(Psi,Psi_prev) Output the permutation
    %   indices used in the reordering
    %   [___] = AdiabaticContinue(___,Name,Value) specifies options using
    %   name-value arguments in addition to any of the input arguments in
    %   previous syntaxes.
    % 
    % Inputs:
    %   Psi - Wave functions to be reordered
    %   Psi_prev - Reference wave functions to adiabatically continue. Either
    %   static of Floquet representation
    %   Name-Value pairs
    %
    % Outputs:
    %   Psi - Reordered wave functions
    %   ind - Indices of the reordering permutation
    %
    % Name-value arguments:
    %   eps - Quasi-energies of the states to be reordered. If not provided can
    %   be calculated automatically
    %   eps_prev - Quasi-energies of the reference states. If provided will
    %   adiabatically continue the quasi-energies as well.

    Psi = obj.Psi_Floquet(Psi);
    switch size(Psi_prev,1)
        case obj.N
            % Special case if input static WF
            Psi0_prev = Psi_prev;
        case obj.N * obj.k_max2
            Psi0_prev = obj.Psi0(Args2.Psi_prev);
        otherwise
            error('Wrong size for Psi_prev (1)');
    end
    % TODO: Fill nan if incompatible sizes and use tol as required limit
    if size(Psi_prev,2) > size(Psi,2)
        error('Size missmatch between Psi and Psi_prev');
    end
    Psi0 = obj.Psi0(Psi);
    % Overlap with previous
    S0 = Psi0_prev' * Psi0;
    % Get the closest indeces with respect to the previous WF (row)
    [mval,ind] = max(abs(S0),[],1);
    % Check if overlap is too small to adiabatically continue
    if ~isempty(find(abs(1 - mval) > Args.tol, 1))
        warning('Overlap is too small');
    end
    % Check if missing values
    if length(ind) ~= length(unique(ind))
        error('Missing permutation index');
    end
    % Reorder to permute properly
    [~,ind] = sort(ind);
    % Append extra wave functions
    ind = ind(:)';
    ind = [ind setdiff(1:size(Psi,2),ind)];
    Psi = Psi(:,ind);
    % Shift to adiabatically continue the QE
    if isfield(Args,'eps_prev')
        if isfield(Args,'eps')
            Args.eps = Args.eps(ind);
        else
            Args.eps = obj.eps(Psi,normalize=true);
        end
        % Calculate the difference in quasi-energies
        dk_eps = round((Args.eps - Args.eps_prev - (mod(Args.eps-Args2.eps_prev+obj.w/2,obj.w) - obj.w/2)) / obj.w);
        for iN = 1:obj.N
            Psi(:,iN)=circshift(Psi(:,iN), obj.N * dk_eps(iN), 1);
        end
    end
end