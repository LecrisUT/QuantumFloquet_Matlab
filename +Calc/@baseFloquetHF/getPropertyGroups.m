function groups = getPropertyGroups(obj)
    import matlab.mixin.util.PropertyGroup
    groups = getPropertyGroups@baseFloquet(obj);
    if isscalar(obj)
        BaseFloquetHF = struct(Ne=obj.Ne,mode=obj.mode,hUk_max=obj.hUk_max);
        groups = [ groups,...
            PropertyGroup(BaseFloquetHF,'Floquet Hartree-Fock properties:')];
    end
end