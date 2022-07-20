function [S,om_mnk] = SpectraOverlap(obj,Psi1,Psi2,Args2)
    arguments
        obj     Calc.baseSystemBath
        Psi1    double
        Psi2    double
        Args2.sk_max    (1,1)   double  {mustBeInteger,mustBePositive}  = obj.hk_max
        Args2.eps1
        Args2.eps2
    end
    % SpectraOverlap Calculate the spectra overlap/transition measure
    % Includes the effects of the bath
    %   
    % See also baseFloquet.SpectraOverlap, baseBath.gamma

    Args2 = namedargs2cell(Args2);
    [S,om_mnk] = obj.system.SpectraOverlap(Psi1,Psi2,Args2{:});
    for ind = 1:length(S(:))
        S(ind) = S(ind)^2 * obj.bath.gamma(-om_mnk(ind));
    end
    if nargout < 2
        S = max(S,[],3);
    end
end