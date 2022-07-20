function CacheAll(obj)
    % CacheAll Cache all fundamental matrix objects
    %
    % See also QuanFloq.baseFloquet.dirty_cache
    
    obj.dirty_cache = false;
    obj.cache_h = obj.get_h;
    obj.cache_hf = obj.h + obj.pt;
    obj.cache_hf2 = obj.hf * obj.hf;
end