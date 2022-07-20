function H = get_H(obj)
    arguments
        obj     QuanFloq.baseFloquetHF
    end
    % get_H Basic getter for H
    % [To be documented]
    % 
    % Syntax:
    %   H = get_H
    % 
    % Description:
    %   H = get_H
    % 
    % Inputs:
    %   [none]
    %
    % Outputs:
    %   H

    switch obj.mode
        case 'closed-shell'
            H = obj.h + 0.5 * obj.hU;
        otherwise
            error('Not implemented')
    end
end