function json = jsonencode(obj,varargin)
    j = jsonencode@Calc.baseFloquet(obj);
    S = jsondecode(j);
    S.floquet_hf = struct(Ne=obj.Ne,hUk_max=obj.hUk_max,mode=obj.mode);
    json = jsonencode(S,varargin{:});
end