function rho = SteadyState(obj,Args)
    arguments
        obj     QuanFloq.baseSystemBath
        Args.type   {mustBeMember(Args.type,{'Redfield','Lindblad','Secular'})} = 'Lindblad'
        Args.MEq    (:,:)   double
        Args.tol    (1,1)   double  {mustBeReal,mustBePositive} = 1E-12
    end
    % SteadyState Calculate the steady state density matrix
    % 
    % Syntax:
    %   rho = SteadyState
    %   [___] = SteadyState(___,Name,Value)
    % 
    % Description:
    %   rho = SteadyState Calculate the steady state
    %   [___] = SteadyState(___,Name,Value) specifies options using name-value
    %   arguments in addition to any of the input arguments in previous
    %   syntaxes.
    % 
    % Inputs:
    %   Name-Value pairs
    %
    % Outputs:
    %   rho - Steady state density matrix
    %
    % Name-value arguments:
    %   MEq - Master equation. If not provided it is automatically generated
    %   based on type
    %   type - ['Lindblad'] Type of master equation to generate
    %   tol - [1E-12] Tolerance where the staedy state 
    %   
    % See also baseSystemBath.Lindblad

    if ~isfield(Args,'MEq')
        switch Args.type
            case 'Lindblad'
                Args.MEq = obj.Lindblad;
            otherwise
                error('Not implemented')
        end
    end
    rho = null(Args.MEq,Args.tol);
    for iN = 1:size(rho,2)
        rho(:,iN) = rho(:,iN) / sum(rho(:,iN));
    end
    rho(abs(rho)<Args.tol) = 0;
end