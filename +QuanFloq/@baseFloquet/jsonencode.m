function json = jsonencode(obj,varargin)
    j = jsonencode@Calc.baseCalc(obj);
    S = jsondecode(j);
    S.matrix = struct(N=obj.N);
    S.floquet = struct(w=obj.w,k_max=obj.k_max,hk_max=obj.hk_max,xi=obj.xi);
    json = jsonencode(S,varargin{:});
end