function F = get_F(obj)
    arguments
        obj     Calc.baseFloquetHF
    end
    % get_F Basic getter for F
    % [To be documented]
    % 
    % Syntax:
    %   F = get_F
    % 
    % Description:
    %   F = get_F
    % 
    % Inputs:
    %   [none]
    %
    % Outputs:
    %   F

    switch obj.mode
        case 'closed-shell'
            F = obj.h + obj.hU;
        otherwise
            error('Not implemented')
    end
end