function pk = Spectra(obj,Psi,Args)
    arguments
        obj     QuanFloq.baseFloquet
        Psi     double
        Args.E_range    (1,:)   double  {mustBeReal}
        Args.t_max      (1,1)   double  {mustBeReal,mustBePositive}
        Args.eps        (1,:)   double  {mustBeReal}
    end
    % Spectra Calculate the energy spectra
    % (Currently only Fourier spectra is implemented)
    %
    % Syntax:
    %   pk = Spectra(psi)
    %   [___] = Spectra(___,Name,Value)
    % 
    % Description:
    %   pk = Spectra(Psi) Calculate the energy spectra of states Psi
    %   [___] = Spectra(___,Name,Value) specifies options using name-value
    %   arguments in addition to any of the input arguments in previous
    %   syntaxes.
    % Inputs:
    %    Psi - Wave function
    %
    % Outputs:
    %    pk - Energy spectra
    %
    % Name-value arguments:
    %   E_range - Calculate the continuous energy spectra instead of Fourier
    %   one
    %   t_max - Calculate the finite-time energy spectra
    %   eps - Quasi-energies to shift the continuous energy spectra
    %   
    % See also baseFloquet.eigs, QuanFloq.baseFloquet.xi

    if isfield(Args,'E_range')
        % TODO: Implement calculating energy spectra for specific
        % energy values
        error('Not implemented');
        if isfield(Args,'t_max')
            % TODO: Implement calulating finite energy spectra
            error('Not implemented');
        end
    else
        % Make sure the wave function is in Fourier representation
        Psi = obj.Psi_Fourier(Psi);
        pk = zeros(obj.k_max2,size(Psi,2));
        for ik = 1:obj.k_max2
            pk(ik,:) = vecnorm(Psi(:,:,ik));
        end
    end
end