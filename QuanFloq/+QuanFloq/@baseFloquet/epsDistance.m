function deps = epsDistance(obj,eps1,eps2)
    arguments
        obj     QuanFloq.baseFloquet
        eps1    (:,1)   double  {mustBeReal}
        eps2    (:,1)   double  {mustBeReal}
    end
    % epsDistance Calculate the quasi-energy distance within a BZ
    % 
    % Syntax:
    %   deps = epsDistance(eps1,eps2)
    % 
    % Description:
    %   deps = epsDistance(eps1,eps2) Calculate the closest quasi-energy
    %   distance from eps1 to eps2
    % 
    % Inputs:
    %   eps1 - First quasi-energy
    %   eps2 - Second quasi-energy
    %
    % Outputs:
    %   deps - Quasi-energy distance
    %   
    % See also QuanFloq.baseFloquet.w

    eps1 = mod(eps1,obj.w);
    eps2 = mod(eps2,obj.w);
    deps = eps1 - eps2;
    % Wrap around the top of the BZ
    deps(deps>obj.w/2) = deps(deps>obj.w/2) - obj.w;
    % Wrap around the bottom of the BZ
    deps(deps<-obj.w/2) = deps(deps<-obj.w/2) + obj.w;
end