function Psit = Psit(obj,Psi,t,Args)
    arguments
        obj     QuanFloq.baseCalc
        Psi     double
        t       (1,1,1,:)   double
        Args.eps    (1,:)   double
    end
    % Psit Calculate the wave function at time $t$
    % 
    % Syntax:
    %   Psit = Psit(Psi,t)
    %   Psit(___,Name,Value)
    % 
    % Description:
    %   Psit = Psit(Psi,t) Calculate $\ket{\Psi(t)}$ of wave function Psi at
    %   time t
    %   Psit(___,Name,Value) Alter the calculation with name-value pairs
    % 
    % Inputs:
    %   Psi - Floquet wave function to evaluate
    %   t - Time parameter
    %   Name-Value pairs
    %
    % Outputs:
    %   Psit - Wave function at time t
    %
    % Name-value arguments:
    %   eps - Quasi-energies to calculate physical wave functions instead of
    %   Floquet ones
    %   
    % See also baseFloquet.Psi_Fourier

    % Make sure the wave function is in Fourier representation
    Psi = obj.Psi_Fourier(Psi);
    nt = length(t);
    M = size(Psi,2);
    % Calculate the Fourier exponents
    expt = exp(-1i * obj.w * reshape(obj.k_range,1,1,[]) .* t);
    if isfield(Args,'eps')
        if length(Args.eps) ~= M
            error('Wrong size of eps');
        end
        % Include the quasi-energy exponents
        expt = exp(-1i * Args.eps .* t) .* expt;
    end
    % Include the wave function
    Psit = Psi .* expt;
    % Sum the Fourier components and reshape
    Psit = reshape(sum(Psit,3),obj.N,M,nt);
end