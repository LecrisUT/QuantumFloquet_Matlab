function [S,om_mnk] = SpectraOverlap(obj,Psi1,Psi2,Args)
    arguments
        obj     QuanFloq.baseFloquet
        Psi1    double
        Psi2    double
        Args.sk_max (1,1)   double  {mustBeInteger,mustBePositive}  = 2*obj.hk_max
        Args.eps1   (:,1)   double  {mustBeReal}
        Args.eps2   (:,1)   double  {mustBeReal}
    end
    % SpectraOverlap Calculate the spectra overlap/transition measure
    % 
    % Syntax:
    %   S = SpectraOverlap(Psi1,Psi2)
    %   [S,om_mnk] = SpectraOverlap(___)
    %   [___] = SpectraOverlap(___,Name,Value)
    % 
    % Description:
    %   S = SpectraOverlap(Psi1,Psi2) Calculate the spectra overlap measure of
    %   wave functions Psi1 and Psi2
    %   [S,om_mnk] = SpectraOverlap(___) Calculate the spectra transition with
    %   corresponding transition frequencies
    %   [___] = SpectraOverlap(___,Name,Value) specifies options using
    %   name-value arguments in addition to any of the input arguments in
    %   previous syntaxes.
    % 
    % Inputs:
    %   Psi1 - First set of wave functions
    %   Psi2 - Second set of wave functions
    %   Name-Value pairs
    %
    % Outputs:
    %   S - Spectra overlap/transition measure
    %   om_mnk - Transition frequency
    %
    % Name-value arguments:
    %   sk_max - [2*hk_max] Fourier cut off of the spectra overlap 
    %   eps1 - Provide the quasi-energies of the states Psi1, otherwise
    %   calculated automatically
    %   eps2 - Provide the quasi-energies of the states Psi2, otherwise
    %   calculated automatically
    %   
    % See also baseFloquet.Spectra, baseSystemBath.SpectraOverlap

    % Limit the Fourier cut off
    if Args.sk_max > obj.k_max
        Args.sk_max = obj.k_max;
    end
    % Calculate the quasi-energies if not provided
    if ~isfield(Args,'eps1')
        Args.eps1 = obj.eps(Psi1);
    end
    if ~isfield(Args,'eps2')
        Args.eps2 = obj.eps(Psi2);
    end
    % Calculate the Fourier space energy spectra
    P1 = obj.Spectra(Psi1);
    P2 = obj.Spectra(Psi2);
    % Initialize the outputs
    M1 = size(P1,2);
    M2 = size(P2,2);
    sk_max2 = 2 * Args.sk_max + 1;
    S = zeros(M1,M2,sk_max2);
    om_mnk = zeros(size(S));
    % Calculate the Spectra transition measures and their frequencies
    for iM = 1:M1
        for jM = 1:M2
            for k = -Args.sk_max:Args.sk_max
                % Transition frequencies
                om_mnk(iM,jM,Args.sk_max+1+k) = Args.eps2(jM) - Args.eps1(iM) + k * obj.w;
                % Spectra transition measure
                if k < 0
                    S(iM,jM,Args.sk_max+1+k) = sum(P1(1-k:obj.k_max2,iM) .* P2(1:obj.k_max2+k,jM));
                else
                    S(iM,jM,Args.sk_max+1+k) = sum(P1(1:obj.k_max2-k,iM) .* P2(1+k:obj.k_max2,jM));
                end
            end
        end
    end
    % If ignoring the transition frequency output the spectra overlap measure
    if nargout < 2
        S(om_mnk>=0) = nan;
        S = max(S,[],3);
    end
end