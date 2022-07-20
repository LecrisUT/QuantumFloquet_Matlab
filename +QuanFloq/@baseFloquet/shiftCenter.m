function Psi = shiftCenter(obj,Psi,Args)
    arguments
        obj     QuanFloq.baseFloquet
        Psi     double
        Args.filter (1,1)   logical = true
        Args.tol    (1,1)   double  = 1E-6
    end
    % shiftCenter Center the wave function in Fourier space
    % 
    % Syntax:
    %   Psi = shiftCenter(Psi)
    %   [___] = shiftCenter(___,Name,Value)
    % 
    % Description:
    %   Psi = shiftCenter(Psi) Shift the wave function to center the Fourier
    %   spectra
    %   [___] = shiftCenter(___,Name,Value) specifies options using name-value
    %   arguments in addition to any of the input arguments in previous
    %   syntaxes.
    % 
    % Inputs:
    %   Psi - Wave function to center
    %   Name-Value pairs
    %
    % Outputs:
    %   Psi - Centered wave function
    %
    % Name-value arguments:
    %   filter - [true] Whether to filter equivalent wave functions
    %   tol - [1E-6] Tolerance for considering wave functions equivalent
    %   
    % See also baseFloquet.Spectra

    %% Shift the wave functions
    % Make sure Psi is in Floquet representation
    Psi = obj.Psi_Floquet(Psi);
    % Compute the energy Fourier space energy spectra
    pk = obj.Spectra(Psi);
    for iM = 1:size(Psi,2)
        [peak,tk] = findpeaks(smoothdata(pk(:,iM)),...
            NPeaks=1, SortStr="descend");
        Psi(:,iM) = circshift(Psi(:,iM), obj.k_max+1 - tk, 2);
    end
    %% Filter equivalent wave functions
    if Args.filter
        iM = 2;
        S = Psi' * Psi;
        while iM <= size(Psi,2)
            % Flag if found duplicate
            flag = false;
            for iM2 = 1:iM-1
                % Check if close to duplicate
                if (1 - abs(S(iM2,iM))) < Args.tol
                    % Remove duplicate
                    flag = true;
                    Psi(:,iM) = [];
                    S(:,iM) = [];
                    S(iM,:) = [];
                    break;
                end
            end
            if ~flag; iM = iM + 1; end
        end
    end
end