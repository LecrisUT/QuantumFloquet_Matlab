function groups = getPropertyGroups(obj)
    import matlab.mixin.util.PropertyGroup
    if ~isscalar(obj)
        groups = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
    else
        BaseMatrix = struct(N=obj.N,k_max=obj.k_max);
        BaseFloquet = struct(w=obj.w,hk_max=obj.hk_max,xi=obj.xi);
        BaseCoupling = struct(systemBath=obj.systemBath);
        groups = [PropertyGroup(BaseMatrix,'Floquet matrix properties:'),...
            PropertyGroup(BaseFloquet,'Floquet physical properties:'),...
            PropertyGroup(BaseCoupling,'Couplings:')];
    end
end