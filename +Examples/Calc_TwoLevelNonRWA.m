clc; clear all;
w0=1; w=1.5;
Unpert = Calc.TwoLevel(w0=w0,w=w,k_max=100,v=0,xi=1E-12,...
    RWA=false,rand_v=true);
Pert = Calc.TwoLevel(w0=w0,w=w,k_max=100,v=1E-6,xi=1E-3,...
    RWA=false,rand_v=true);
fine_tune(1) = struct(Vc=1.7035769013796129, delV=4E-4, nV=1001,...
    V_range=[]);
fine_tune(2) = struct(Vc=4.0843758839750000, delV=4E-4, nV=1001,...
    V_range=[]);
nV=501;
V_range=linspace(0,5,nV);
for iF = 1:length(fine_tune)
    ft = fine_tune(iF);
    fine_tune(iF).V_range = linspace(-ft.delV,ft.delV,ft.nV) + ft.Vc;
    ft = fine_tune(iF);
    V_range=unique(sort([V_range ft.V_range]));
end
nV = length(V_range);
%% Save calculation metadata
Details = Calc.Details(script='Calc_TwoLevelNonRWA.m',...
    model='Two-level system noRWA [cos(w*t)]',...
    calculation='Adiabatic Floquet eigenstate calculation',...
    objects=[Pert,Unpert],...
    variables=struct(V_range=V_range,fine_tune=fine_tune));
%% Calculate eigenstates
% Adiabatically continue from static eigenstates
Psi0 = [0 1;1 0];
eps0 = [-w0/2;w0/2];
iter=Calc.GenericCalcIterator(Unpert,data=V_range,updatefcn=@AdiabaticV);
iter.reset;
[tPsi,teps,tEbar] = Unpert.eigs(iterator=iter,...
    Psi_prev=Psi0,eps_prev=eps0);
Res_Unpert_AE = PackRes(Unpert,tPsi,teps,Ebar=tEbar);
iter.reset;
[tPsi,teps] = Unpert.eigs(iterator=iter,...
    Psi_prev=Psi0,eps_prev=eps0);
Res_Unpert_QE = PackRes(Unpert,tPsi,teps);
iter=Calc.GenericCalcIterator(Pert,data=V_range,updatefcn=@AdiabaticV);
iter.reset;
[tPsi,teps,tEbar] = Pert.eigs(iterator=iter,...
    Psi_prev=Psi0,eps_prev=eps0);
Res_Pert_AE = PackRes(Pert,tPsi,teps,Ebar=tEbar);
iter.reset;
[tPsi,teps] = Pert.eigs(iterator=iter,...
    Psi_prev=Psi0,eps_prev=eps0);
Res_Pert_QE = PackRes(Pert,tPsi,teps);
%% Save results
% save('Calc_TwoLevelNonRWA.mat',...
%     'Res_Unpert_QE','Res_Unpert_AE',...
%     'Res_Pert_QE','Res_Pert_AE',...
%     'Details',...
%     '-v7.3');
%% Helper functions
function AdiabaticV(obj)
    obj.object.V = obj.data(obj.ind);
end
function Res = PackRes(obj,Psi,eps,Args)
    arguments
        obj
        Psi
        eps
        Args.Ebar
    end
    M = size(Psi,2);
    L = size(Psi,3);
    Res(L,M) = struct(Psi=[],eps=[],P=[]);
    for iL = 1:L
        tPsi = obj.Psi_Fourier(Psi(:,:,iL));
        for iM = 1:M
            %% Pack calculated results
            Res(iL,iM).Psi = Psi(:,iM,iL);
            Res(iL,iM).eps = eps(iM,iL);
            if isfield(Args,'Ebar')
                Res(iL,iM).Ebar = Args.Ebar(iM,iL);
            end
            %% Other PostProcess
            Res(iL,iM).P = zeros(obj.k_max2,1);
            for ik = 1:obj.k_max2
                Res(iL,iM).P(ik) = tPsi(:,iM,ik)' * tPsi(:,iM,ik);
            end
        end
    end
end