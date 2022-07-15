function Psi = Psi_Fourier(obj,Psi)
    arguments
        obj     Calc.baseFloquet
        Psi     double
    end
    % Psi_Fourier Transform a wave function in Fourier representation
    % 
    % Syntax:
    %   Psi = Psi_Fourier(Psi)
    % 
    % Description:
    %   Psi = Psi_Fourier(Psi) Transform Psi into Fourier representation
    % 
    % Inputs:
    %   Psi - Wave function to transform
    %
    % Outputs:
    %   Psi - Equivalent wave function in Fourier representation
    %
    % See also baseFloquet.Psi_Floquer

    % Check the current size of wave function
    switch ndims(Psi)
        case 2
            switch size(Psi,1)
                case obj.N
                    if size(Psi,2) ~= obj.k_max2
                        error('Wrong size psi:2');
                    end
                    Psi = reshape(Psi,obj.N,1,obj.k_max2);
                case obj.N * obj.k_max2
                    Psi = permute(reshape(Psi,obj.N,obj.k_max2,[]),[1 3 2]);
                otherwise
                    error('Unsupported size');
            end
        case 3
            if size(Psi,1) ~= obj.N
                error('Wrong size psi:1');
            end
            if size(Psi,3) ~= obj.k_max2
                if size(Psi,2) ~= obj.k_max2
                    error('Wrong size psi:2');
                end
                Psi = permute(Psi,[1 3 2]);
            end
        otherwise
            error('Unsupported dimensions')
    end
end