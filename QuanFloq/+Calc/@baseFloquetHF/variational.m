function Res = variational(obj,method,Psi0,Args2)
    arguments
        obj     Calc.baseFloquetHF
        method  {mustBeMember(method,{'AE','varQE'})}
    end
    arguments (Repeating)
        Psi0
    end
    arguments
        Args2.SlaterCons
        Args2.opts
        Args2.ms
        Args2.Psi0
        Args2.TypX
        Args2.tol
        Args2.filter_nconv
        Args2.TrackPsi
    end
    % variational Calculate Floquet eigenstates using a variational method
    % See baseFloquet.variational
    %
    % Additional name-value arguments:
    %   SlaterCons - Whether to include Slater determinant normalization
    %   constrain
    %   
    % See also baseFloquet.variational

    if isempty(Psi0)
        error('No Psi0 argument passed, cannot automatically generate');
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