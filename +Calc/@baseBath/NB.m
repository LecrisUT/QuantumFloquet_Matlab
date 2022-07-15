function NB = NB(obj,omega)
    arguments
        obj     Calc.baseBath
        omega   (:,1)   double
    end
    % NB Calculate bath occupation number
    % 
    % Syntax:
    %   NB = NB(omega)
    % 
    % Description:
    %   NB = NB(omega) Calculate bath occupation number with frequenct omega
    % 
    % Inputs:
    %   omega - Bath frequency
    %
    % Outputs:
    %   NB - Occupation number
    %   
    % See also Calc.baseBath.type, Calc.baseBath.beta

    switch obj.type
        case 'fermion'
            error('Not implemented')
        case 'boson'
            if isinf(obj.beta)
                if omega < 0
                    NB = 1;
                else
                    NB = 0;
                end
            else
                if omega == 0
                    NB = 0;
                elseif omega < 0
                    omega = -omega;
                    NB = 1 / (1 - exp(-obj.beta * omega));
                else
                    NB = 1 / (exp(obj.beta * omega) - 1);
                end
            end
        otherwise
            error('Unknown bath type');
    end
end