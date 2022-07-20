function CacheAllU(obj)
    % CacheAllU Cache all self-consistent matrix objects
    %
    % See also Calc.baseFloquetHF.dirty_cacheU
    
    obj.dirty_cacheU = false;
    obj.cache_hU = obj.get_hU;
    obj.cache_H = obj.get_H;
    obj.cache_F = obj.get_F;
    obj.cache_Hf = obj.H + obj.pt;
    obj.cache_Ff = obj.F + obj.pt;
end