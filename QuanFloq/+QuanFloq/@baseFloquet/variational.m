function Res = variational(obj,method,Psi0,Args2)
    arguments
        obj     QuanFloq.baseFloquet
        method  {mustBeMember(method,{'AE','varQE'})}
    end
    arguments (Repeating)
        Psi0
    end
    arguments
        Args2.opts
        Args2.ms
        Args2.Psi0
        Args2.TypX
        Args2.tol
        Args2.filter_nconv
        Args2.TrackPsi
    end
    % variational Calculate Floquet eigenstates using a variational method
    % 
    % Syntax:
    %   Res = variational(method)
    %   [___] = variational(method,Psi0)
    %   [___] = variational(___,Name,Value)
    % 
    % Description:
    %   Res = variational(method) Perform the minimization specified in method
    %   with default generated starting wave functions
    %   [___] = variational(method,Psi0) Specify the starting wave functions
    %   Psi0
    %   [___] = variational(___,Name,Value) specifies options using name-value
    %   arguments in addition to any of the input arguments in previous
    %   syntaxes. See specific implementations for default values.
    % 
    % Inputs:
    %   method - minimization method (either 'AE' for average energy
    %   minimization, or 'varQE' for quasi-energy variance)
    %   Psi0 - (Repeating) Starting wave function guesses
    %   Name-Value pairs
    %
    % Outputs:
    %   Res - Structured output vector array containing all (convergent)
    %   solutions (see specific implementations)
    %
    % Name-value arguments:
    %   opts - Override gneerated optimoptions of 'fmincon'
    %   ms - Override gneerated multistart options
    %   TypX - Override generated TypX for opts
    %   tol - Convergence tolerance. Not to be confused with acceptable error
    %   xi
    %   filter_nconv - Whether to filter non-convergent solutions
    %   TrackPsi - Whether to save the trial wave functions
    %   
    % See also QuanFloq.baseFloquet.xi, baseFloquet.variational_AE,
    % baseFloquet.variational_varQE

    if isempty(Psi0)
        tPsi0 = zeros(obj.N, 1);
        tPsi0(1) = 1;
        Psi0 = {tPsi0};
    end
    switch method
        case 'varQE'
            Args2 = namedargs2cell(Args2);
            Res=obj.variational_varQE(Psi0{:},Args2{:});
        case 'AE'
            Args2 = namedargs2cell(Args2);
            Res=obj.variational_AE(Psi0{:},Args2{:});
        otherwise
            error('Not implemented');
    end
end